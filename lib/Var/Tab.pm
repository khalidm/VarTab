package Var::Tab;
use strict;
#use warnings;
use Carp;
use Data::Dumper;
use Vcf;

use Misc::PrintFile qw(print_to_file print_to_var_type_single print_to_var_type);
use Misc::Calc qw( add multiply percent_to_count_threshold maf_filter );

use Exporter qw(import);
our @EXPORT_OK = qw(get_gene_counts get_ann bed_annotate_cpg bed_annotate_tfbs 
    bed_annotate_polyphen bed_annotate_cadd bed_annotate_gwascatalog 
    bed_annotate_hapmap get_effect_type_region get_variant_contingency_table);



#
# gene burden
#
sub get_gene_counts {
    my ( $opts ) = @_;
    # my $opts = parse_params();
    my $getseq = "F";
    my $sample_count = 0;

    # output file names
    my $output = $$opts{output};
    my $outputburden = $output.".burden";
    # OPEN OUTPUT FILE PREFIX.tab
    open(OUTBUR, ">$outputburden") || die "Can't open the file: $outputburden: $!\n";    

    my $freq_threshold;
    my $annotation_bed = "False";
    
    if ( exists($$opts{frequency}) )
    {
        $freq_threshold = $$opts{frequency} * 1.0;
    }
    else 
    {
        # set to default
        $freq_threshold = 5.0;
    }

    # check print flanking sequence : requires genome fasta file
    if (exists($$opts{sequence}))
    {
        $getseq = "T";
    }

    my $input = $$opts{input};
    #my $vcf = Vcf->new(fh=>\*STDIN);
    my $vcf = Vcf->new(file=>$input);
    $vcf->parse_header();

    my $header_printed=0;
    my $total = 0;
    my $snpeff_gene;
    my $snpeff_type;
    my %count;

    while (my $x=$vcf->next_data_hash())
    {
        $total++;
        # print Dumper($x);
        if ( !$header_printed ) 
        {
            #print "dbSNP_ID\tChr\tPOS\tREF\tALT\tGENE\tTYPE\tNETWORK\tTF_binding_peak\tDNASE";
            #print "\tCLINICAL\tAA_CHANGE\tConservation\tRMSK\tCpG\tGWAS\tSampleFreq.HOM_REF\tSampleFreq.HET\tSampleFreq.HOM_ALT\tFS";
            if($getseq eq "T") { print "\tFASTA"; }
            for my $col (sort keys %{$$x{gtypes}})
            {
                #print "\t$col";
                $sample_count++;
            }
            #print "\n";
            $header_printed = 1;
        }

        if( $$opts{nondbsnp} )
        {
            if($$x{ID} =~ m/^\./)
            {         
                #do_count($x, $freq_threshold, $getseq, $vcf);
                #print $x."\n";
                $snpeff_gene = get_ann($x, "SNPEFF_GENE_NAME");

                ($snpeff_type, $snpeff_gene) = get_snpeff_ann($x);
                
                #$vargenecounts{$snpeff_gene}++;
                if($snpeff_type eq "NON_SYNONYMOUS_CODING"){
                    $count{$snpeff_gene}++;
                }
            }
        }
        else
        {            
            #do_count($x, $freq_threshold, $getseq, $vcf);
            #print $x."\n";
            #$snpeff_gene = get_ann($x, "SNPEFF_GENE_NAME");
            #my ($snpeff_type, $snpeff_gene, $snpeff_aa_change) = get_snpEffannotations($x);
            
            #$snpeff_gene = get_snpeff_ann($x);
            ($snpeff_type, $snpeff_gene) = get_snpeff_ann($x);
            
            #$vargenecounts{$snpeff_gene}++;
            if($snpeff_type eq "NON_SYNONYMOUS_CODING"){
                $count{$snpeff_gene}++;
            }
        }
    }
    $vcf->close();

    my $gene_var_freq = 0.00;

    foreach my $str (sort keys %count) 
    {
        #printf OUTBUR "%-31s %s\t%s\n", $str, $count{$str}, scalar(keys %{$$x{gtypes}});
        $gene_var_freq = $count{$str}/$sample_count;
        #printf OUTBUR "%-31s %s\t%s\t%s\n", $str, $gene_var_freq, $count{$str}, $sample_count;
        printf OUTBUR "%s\t%0.5f\t%s\t%s\n", $str, $gene_var_freq, $count{$str}, $sample_count;
        #printf OUTBUR "%s\t%s\n", $$x{CHROM},$$x{POS};
    }

    #print OUTBUR "\n";
    #foreach my $word (reverse sort { $count{$a} <=> $count{$b} } keys %count) 
    #{
    #    printf OUTBUR "%-31s %s\n", $word, $count{$word};
    #}
    close OUTBUR;
}

#
# variant contingency table
#
sub get_variant_contingency_table {
    my ($opts) = @_;
    my $getseq = "F";

    my $freq_threshold;
    my $annotation_bed = "False";
    
    if ( exists($$opts{frequency}) )
    {
        $freq_threshold = $$opts{frequency} * 1.0;
    }
    else 
    {
        # set to default
        $freq_threshold = 5.0;
    }

    # check print flanking sequence : requires genome fasta file
    if (exists($$opts{sequence}))
    {
        $getseq = "T";
    }
   
    my $iupac;
    if ( $$opts{iupac} ) { $iupac=$$opts{iupac}; }

    my $input = $$opts{input};
    my $output = $$opts{output};
    
    # output file names
    my $outputtab = $output.".matrix";
    
    # output array
    # my @frame_events = (((1) x 10), ((1) x 10));
    my @arr;    
    # my @gt_array = (0);
    my $gt_counts = 0;
        
    # OPEN OUTPUT FILE PREFIX.tab
    # open(OUTTAB, ">$outputtab") || die "Can't open the file: $outputtab: $!\n";
    open my $OUTTAB, '>', $outputtab or die "...$!\n";
   
    #my $vcf = Vcf->new(fh=>\*STDIN);
    my $vcf = Vcf->new(file=>$input);
    $vcf->parse_header();
    # print length($vcf);

    my $header_printed=0;
    my $total = 0;
    my $print_string = "";
    my @maf_output = ();
    my $header = "";    

    while (my $x=$vcf->next_data_hash())
    {
        my $r = $$x{REF};
        my $a = join("", @{$$x{ALT}});
        #print $r."\t".$a."\n";
        $total++;

        # print Dumper($x);
        if ( !$header_printed ) 
        {
            my @gt_array = (0);
            my $gt_counts = 0;
            for my $col (sort keys %{$$x{gtypes}}) {
                # print $col."\t";
                $gt_counts++;               
            }    
            @gt_array = (0) * ($gt_counts + 1);
            # $header = "chr:pos\tREF\tALT\tGENE\tAA\tTYPE\tREGION\tCADD\tFUNSEQ\tPOLY_SIFT\tCONS\tTFBS\tDNASE\tCpG\tGWAS\tRMSK\tncRNA_p-gene\tCLINICAL\tMAF\tSampleFreq.HOM_REF\tSampleFreq.HET\tSampleFreq.HOM_ALT";
            $header = "chr:pos:id";
            $gt_array[0] = $header;
            
            # ADD FUNCTION TO PRINT TO FILE
            print_to_file($OUTTAB, $header); 
            
            if($getseq eq "T") { print OUTTAB "\tFASTA"; }
            my $index = 1;
            for my $col (sort keys %{$$x{gtypes}})
            {
                $gt_array[$index] = $col;
                print_to_file($OUTTAB, "\t".$col);
            }
            print_to_file($OUTTAB, "\n");
            
            $header_printed = 1;
            push @arr, [@gt_array];
        }

        # non dbsnp variants only
        if( $$opts{nondbsnp} )
        {
            if($$x{ID} =~ m/^\./)
            {         
                # print_info($x, $freq_threshold, $getseq, $vcf, $output);
                if( $$opts{maf1kg} ) {
                    my $gt_string = "";
                    # get the ALT
                    my $alt;
                    for my $alt (@{$$x{ALT}}) {
                        if ( $alt eq '.' ) { $alt=$$x{REF}; }                        
                        for my $col (sort keys %{$$x{gtypes}}){                            
                            # my ($al1,$sep,$al2) = exists($$x{gtypes}{$col}{GT}) ? $vcf->parse_alleles($x,$col) : ('.','/','.');
                            # my $gt = $al1.'/'.$al2;
                            # my ($current_gt, $gt_index) = get_gt_type($$x{gtypes}{$col}{GT});
                            # $gt_string = "$gt_string\t$current_gt";
                            my $current_gt = $$x{gtypes}{$col}{GT};
                            if($current_gt ne "."){
                                $gt_string = "$gt_string\t1";
                            } else {
                                $gt_string = "$gt_string\t.";
                            }
                        }
                    }
                    @maf_output = maf_filter($$x{CHROM},$$x{POS},$$x{POS}+1,7,$$opts{maf1kg});
                    my $hapmap = bed_annotate_hapmap($$x{CHROM},$$x{POS},$$x{POS}+1, 13, $$x{REF}, $alt, $$opts{maf1kg});
                    # ONLY PRINT VARIANT TO OUTPUT IF 1KG MAF IS BELOW THRESHOLD OR NOT REPORTED IN THE BED FILE
                    # ALSO ONLY PRINT THOSE VARIANTS WITH HAPMAP MAF > 0.05
                    if ( $maf_output[0] == 1 && $hapmap == 0 ) {
                        # $print_string = print_info($x, $freq_threshold, $getseq, $vcf, $output);
                        $print_string = "$$x{CHROM}:$$x{POS}:$$x{ID}\t$gt_string\n";
                        print_to_var_type_single($r, $a, $OUTTAB, $print_string);
                        
                    } elsif ( $maf_output[0] == 2 && $hapmap == 0 ) {
                        # $print_string = print_info($x, $freq_threshold, $getseq, $vcf, $output);
                        $print_string = "$$x{CHROM}:$$x{POS}:$$x{ID}\t$gt_string\n";
                        print_to_var_type_single($r, $a, $OUTTAB, $print_string);
                    }
                } else {
                    my $gt_string = "";
                    # get the ALT
                    my $alt;
                    for my $alt (@{$$x{ALT}}) {
                        if ( $alt eq '.' ) { $alt=$$x{REF}; }                        
                        for my $col (sort keys %{$$x{gtypes}}){
                            my $current_gt = $$x{gtypes}{$col}{GT};
                            if($current_gt ne "."){
                                $gt_string = "$gt_string\t1";
                            } else {
                                $gt_string = "$gt_string\t.";
                            }
                        }
                    }
                    $print_string = "$$x{CHROM}:$$x{POS}:$$x{ID}\t$gt_string\n";
                    print_to_var_type_single($r, $a, $OUTTAB, $print_string);
                }
            }
        } else {

            # if(@{$$x{ALT}} == 1 && length(@{$$x{ALT}}[0]) == 1) {
            # if(@{$$x{ALT}} == 1) {
            my @gt_array = (0);
            $gt_counts = 0;
            for my $col (sort keys %{$$x{gtypes}}) {
                # print $col."\t";
                $gt_counts++;               
            }
            # reset array
            @gt_array = (0) * ($gt_counts + 1);            

            if( $$opts{maf1kg} ) {
                # print_info($x, $freq_threshold, $getseq, $vcf, $output);
                my $gt_string = "";
                # get the ALT
                my $alt;
                my $index = 1;  
                $gt_array[0] = "**$$x{CHROM}:$$x{POS}:$$x{ID}";
                for my $alt (@{$$x{ALT}}) {
                    if ( $alt eq '.' ) { $alt=$$x{REF}; }
                    for my $col (sort keys %{$$x{gtypes}}) {
                        my $current_gt = $$x{gtypes}{$col}{GT};
                        if($current_gt ne "."){
                            $gt_string = "$gt_string\t1";
                            # $frame_events[0][$index] = 1;
                            $gt_array[$index] = 1;
                        } else {
                            $gt_string = "$gt_string\t.";
                            # $frame_events[0][0] = 0;
                            $gt_array[$index] = 0;
                        }
                        $index++;
                    }
                }
                @maf_output = maf_filter($$x{CHROM},$$x{POS},$$x{POS}+1,7,$$opts{maf1kg});
                my $hapmap = bed_annotate_hapmap($$x{CHROM},$$x{POS},$$x{POS}+1, 13, $$x{REF}, $alt, $$opts{maf1kg});
                my $sum = 0;
                for(my $i = 1; $i < @gt_array; $i++) {
                    $sum = $sum + $gt_array[$i];
                }
                $sum = ($sum/$gt_counts) * 100.0;
                #my $sum = eval join '+', splice(@gt_array,1);                    
                # ONLY PRINT VARIANT TO OUTPUT IF 1KG MAF IS BELOW THRESHOLD OR NOT REPORTED IN THE BED FILE
                # ALSO ONLY PRINT THOSE VARIANTS WITH HAPMAP MAF > 0.05                    
                if ( $maf_output[0] == 1 || $maf_output[0] == 2 || $hapmap == 1 ) {                        
                    # $print_string = "$$x{CHROM}:$$x{POS}:$$x{ID}\t$gt_string\n";
                    # $frame_events[0][0] = "$$x{CHROM}:$$x{POS}:$$x{ID}";
                    # DEBUG $gt_array[0] = "$$x{CHROM}:$$x{POS}:$$x{ID}:$sum:$gt_counts:$freq_threshold";
                    # $gt_array[0] = "**$$x{CHROM}:$$x{POS}:$$x{ID}";
                    # Only print if
                    if( length($r) == 1 && length($a) == 1 && $sum < $freq_threshold) {
                        push @arr, [@gt_array];
                        # print_to_var_type_single($r, $a, $OUTTAB, $print_string);
                        # ORIGINAL print $OUTTAB join(",", @gt_array);
                        # ORIGINAL print $OUTTAB "\n";
                    }
                }
                # for (my $i = 1; $i < ($gt_counts+1); $i++) {
                #     for (my $j = 1; $j < ($gt_counts+1); $j++) {
                #         print $OUTTAB "$arr[$i]_[$j] ";
                #     }
                #     print $OUTTAB "\n";
                # }

                #elsif ( $maf_output[0] == 2 && $hapmap == 0 ) {
                #     # $print_string = print_info($x, $freq_threshold, $getseq, $vcf, $output);
                #     $print_string = "-$$x{CHROM}:$$x{POS}:$$x{ID}\t$gt_string\n";
                #     print_to_var_type_single($r, $a, $OUTTAB, $print_string);
                # }
            } else {
                my $gt_string = "";
                # get the ALT
                my $alt;
                for my $alt (@{$$x{ALT}}) {
                    if ( $alt eq '.' ) { $alt=$$x{REF}; }
                    for my $col (sort keys %{$$x{gtypes}}){
                        # my ($al1,$sep,$al2) = exists($$x{gtypes}{$col}{GT}) ? $vcf->parse_alleles($x,$col) : ('.','/','.');
                        # my $gt = $al1.'/'.$al2;
                        # my ($current_gt, $gt_index) = get_gt_type($$x{gtypes}{$col}{GT});
                        # $gt_string = "$gt_string\t$current_gt";
                        my $current_gt = $$x{gtypes}{$col}{GT};
                        if($current_gt ne "."){
                            $gt_string = "$gt_string\t1";
                        } else {
                            $gt_string = "$gt_string\t.";
                        }
                    }
                }
                $print_string = "-$$x{CHROM}:$$x{POS}:$$x{ID}\t$gt_string\n";
                # print_to_var_type_single($r, $a, $OUTTAB, $print_string);
            }
            # }
            # last;
        }        

    }
    print $OUTTAB "\ngt_counts=".$gt_counts."\n";
    print $OUTTAB "total=".$total."\n\n";
    # for (my $i = 0; $i < ($gt_counts+1); $i++) {
    #     for (my $j = 0; $j < ($total+1); $j++) {
    #         print $OUTTAB "$arr[$i] ";
    #     }
    #     print $OUTTAB "\n";
    # }

    # TRANSPOSE
    my @rows = ();
    my @transposed = ();
    for my $row (@arr) {
        for my $column (0 .. $#{$row}) {
            push(@{$transposed[$column]}, $row->[$column]);
        }
    }
    # PRINT    
    for my $new_row (@transposed) {
        for my $new_col (@{$new_row}) {
            print $OUTTAB $new_col, " ";
        }
        print $OUTTAB "\n";
    }
    # for (my $i = 0; $i < @transposed; $i++) {
    #     for (my $j = 0; $j < @transposed[$i]; $j++) {
    #         print "$transposed[$i][$j] ";
    #     }
    #     print "\n";
    #     last;
    # }
    
    close $OUTTAB;    
    $vcf->close();
}


sub get_ann
{
    my $x = shift;
    my $annotation_name = shift;

    my $str = "";
    #$str = $$x{"$annotation_name"};
    $str = $$x{EFF}{$annotation_name};
    #$str = $$x{"MQRankSum"};

    #print Dumper($x);
    #print $annotation_name."\t*".$str."*\n";
    #exit;
    if($str eq "")
    {
        $str = ".";
    }
    return $str;
}

sub get_snpeff_ann
{
    my $x = shift;

    my $str = "";
    $str = $$x{INFO}{"EFF"};

    #print Dumper($x);
    #print $annotation_name."\t*".$str."*\n";
    #exit;
    #NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gAg/gGg|E791G|925|CEP104||CODING|NM_014704.3|19|1)
    
    my @eff_values = split('\|', $str);
    my $eff_type = $eff_values[0];
    $eff_type =~ s/\(.*//g;

    my $eff_gene = $eff_values[5];
    my $eff_aa_change = $eff_values[3];
    #my $eff_gene = $eff_values[6];

    if($eff_gene eq "")
    {
        $eff_gene = ".";
    }
    if($eff_type eq "")
    {
        $eff_type = ".";
    }
    if($eff_aa_change eq "")
    {
        $eff_aa_change = ".";
    }
    #return ($eff_type, $eff_gene, $eff_aa_change);
    return ($eff_type, $eff_gene);
}


sub bed_annotate_cpg
{
    # my $annotation_tabix = shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $annotation_id = shift;
    #my @tabix_current = shift;
    my @var = ();
    #my @var_array = ();
    #my $out_str = "";
	
    #@snp_tabix = split /\n/, `tabix $tgp_snp -B $input | awk '{FS="\t";OFS="\t"} \$4 >= $maf_cut_off'`;	
    
    #@var = split('\n', $tabix_current[$annotation_id]->read($tabix_current[$annotation_id]->query( $chr, $start, $end)));
    @var = split('\n', $main::tabix[$annotation_id]->read($main::tabix[$annotation_id]->query( $chr, $start, $end)));
    
    if (0+@var > 0){
        return "T";
    } else {
        return ".";
    }
}

sub bed_annotate_tfbs
{
    # my $annotation_tabix = shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $annotation_id = shift;
    #my @tabix_current = shift;
    my @var = ();
    my @var_array = ();
    my @tfbs_array = ();
    #my $out_str = "";

    #my $ann_bed = $main::annotation_beds[$annotation_id];
    my $ann_bed = $annotation_id;
    my $input = $chr.":".$start."-".$end;

    #print "ANN BED = ".$ann_bed."\n";
    #print "ANN INPUT = ".$input."\n";
	
    #@var = split /\n/, `tabix $ann_bed -B $input | awk ' { print \$4 }awk '{FS="\t";OFS="\t"} \$4 >= $maf_cut_off'`;	
    @var = split /\n/, `tabix $ann_bed -B $input | awk ' { print \$4"\t"\$5 } '`;

    #print "ANN SIZE = ".0+@var."\n";
    
    my $tfbs = "";
    my $current_cell_count = 0;
    if (0+@var <= 0){
        return ".";
    } else {
        # split tabix returned string - return string with highest cell count for the transcription factor (ENCODE - wgEncodeRegTfbsClusteredV2)
        foreach my $snp (@var) {
            @var_array = split('\t', $snp);
            if (0+@var_array > 0) {
                #$current_cell_count = $var_array[1];
                if($var_array[1]>$current_cell_count){
                    $current_cell_count = $var_array[1];
                    if($var_array[0] ne ""){
                        #return ($var_array[3]);
                        #$tfbs = $tfbs.";".$var_array[0];
                        $tfbs = $tfbs.";".$var_array[0]."(".$var_array[1].")";
                    }
                    else{
                        return (".");
                    }
                }
            } else {
                return (".");
            }
        }
    }
    #print "ANN = ".$tfbs."\n";
    return $tfbs;
}

# get polyphen scores
sub bed_annotate_polyphen {
    # my $annotation_tabix = shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $annotation_id = shift;
    my $ref = shift;
    my $alt = shift;
    #my @tabix_current = shift;
    my @var = ();
    my @var_array = ();
    my @annotation_array = ();

    my $ann_bed = $annotation_id;
    my $input = $chr.":".$start."-".$end;

    ##print "ANN BED = ".$ann_bed."\n";
    #print "ANN INPUT = ".$input."\t";
    #print "ANN REF/ALT = ".$ref."\/".$alt."\t";
	
    @var = split /\n/, `tabix $ann_bed $input | awk ' { print \$0 } '`;
    
    #print "tabix $ann_bed $input | awk ' { print \$0 } '";

    #print "tabix $ann_bed $input | awk ' { print \$0 } ";
    my $polyphen = "";
    my $arrSize = @var;
    #if($arrSize > 0) {
    #    print $arrSize."\n";
    #}
    
    if (0+@var <= 0) {
        return ".";
    } else {
        # split tabix returned string - return string with highest cell count for the transcription factor (ENCODE - wgEncodeRegTfbsClusteredV2)
        foreach my $snp (@var) {
            @var_array = split('\t', $snp);
            #print $ref." = ".$var_array[3]."\t".$alt." = ".$var_array[4]."\t";
            #print join(", ", @var_array); print "\t"; print $ref."|".$alt."\t"; print "\n";
            #print "$ref=$var_array[3]\t$alt=$var_array[4]\t";
            if ($ref eq $var_array[3] && $alt eq $var_array[4]) {
                #print "$ref=$var_array[3]\t$alt=$var_array[4]\t";                
                #$polyphen = $var_array[5]."(".$var_array[6].")";
                $polyphen = "del";
                last;
            }
            else {
                $polyphen = ".";
            }
        }
    }
    #print "\n";
    return $polyphen;
}

# get polyphen scores
sub bed_annotate_cadd {
    # my $annotation_tabix = shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $annotation_id = shift;
    my $ref = shift;
    my $alt = shift;
    #my @tabix_current = shift;
    my @var = ();
    my @var_array = ();
    my @annotation_array = ();

    my $ann_bed = $annotation_id;
    my $input = $chr.":".$start."-".$end;

    ##print "ANN BED = ".$ann_bed."\n";
    #print "ANN INPUT = ".$input."\t";
    #print "ANN REF/ALT = ".$ref."\/".$alt."\t";
	
    @var = split /\n/, `tabix $ann_bed $input | awk ' { print \$0 } '`;
    
    #print "tabix $ann_bed $input | awk ' { print \$0 } '";

    #print "tabix $ann_bed $input | awk ' { print \$0 } ";
    my $polyphen = "";
    my $arrSize = @var;
    #if($arrSize > 0) {
    #    print $arrSize."\n";
    #}
    
    if (0+@var <= 0) {
        return ".";
    } else {
        # split tabix returned string - return string with highest cell count for the transcription factor (ENCODE - wgEncodeRegTfbsClusteredV2)
        foreach my $snp (@var) {
            @var_array = split('\t', $snp);
            #print $ref." = ".$var_array[3]."\t".$alt." = ".$var_array[4]."\t";
            #print join(", ", @var_array); print "\t"; print $ref."|".$alt."\t"; print "\n";
            #print "$ref=$var_array[3]\t$alt=$var_array[4]\t";
            if ($ref eq $var_array[3] && $alt eq $var_array[4]) {
                #print "$ref=$var_array[3]\t$alt=$var_array[4]\t";                
                #$polyphen = $var_array[5]."(".$var_array[6].")";
                $polyphen = $var_array[5];
                last;
            }
            else {
                $polyphen = ".";
            }
        }
    }
    #print "\n";
    return $polyphen;
}

# get gwascatalog data
sub bed_annotate_gwascatalog {
    # my $annotation_tabix = shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $annotation_id = shift;
    my @var = ();
    my @var_array = ();
    my @gwascatalog_temp = ();
    my $gwascatalog = ".";
    
    @var = split('\n', $main::tabix[$annotation_id]->read($main::tabix[$annotation_id]->query( $chr, $start, $end)));

    foreach my $snp (@var) {
        @var_array = split('\t', $snp);
        if (0+@var_array > 0) {
            @gwascatalog_temp = split(';', $var_array[3]);
            $gwascatalog = $gwascatalog_temp[0];        
        }    
    }

    return $gwascatalog;

    #if (0+@var > 0) {
    #    return "T";
    #} else {
    #    return ".";
    #}
}

# get gwascatalog data
sub bed_annotate_gwascatalog {
    # my $annotation_tabix = shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $annotation_id = shift;
    my @var = ();
    my @var_array = ();
    my @gwascatalog_temp = ();
    my $gwascatalog = ".";
    
    @var = split('\n', $main::tabix[$annotation_id]->read($main::tabix[$annotation_id]->query( $chr, $start, $end)));

    foreach my $snp (@var) {
        @var_array = split('\t', $snp);
        if (0+@var_array > 0) {
            @gwascatalog_temp = split(';', $var_array[3]);
            $gwascatalog = $gwascatalog_temp[0];        
        }    
    }

    return $gwascatalog;

    #if (0+@var > 0) {
    #    return "T";
    #} else {
    #    return ".";
    #}
}


# get hapmap data
sub bed_annotate_hapmap {
    # my $annotation_tabix = shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $annotation_id = shift;
    my $ref = shift;
    my $alt = shift;
    my $maf_filter = shift;
    my @var = ();
    my @var_array = ();
    my @hapmap_temp = ();
    my $hapmap = 0;
    
    @var = split('\n', $main::tabix[$annotation_id]->read($main::tabix[$annotation_id]->query( $chr, $start, $end)));

    foreach my $snp (@var) {
        @var_array = split('\t', $snp);
        if (0+@var_array > 0 && $ref eq $var_array[3] && $alt eq $var_array[5] && $var_array[6] <= $maf_filter) {
            # @hapmap_temp = split(';', $var_array[3]);
            # $hapmap = $hapmap_temp[0];        
            $hapmap = 1;
        }
    }
    return $hapmap;
}

my %region_hash = (
    NONE => [ 'NONE', 'CHROMOSOME', 'CUSTOM', 'CDS', ],
    INTERGENIC => [ 'INTERGENIC', 'INTERGENIC_CONSERVED', ],
    UPSTREAM => [ 'UPSTREAM', ],
    UTR_5_PRIME => [ 'UTR_5_PRIME', 'UTR_5_DELETED', 'UTR_5_GAINED' ],
    SPLICE_SITE_ACCEPTOR => [ 'SPLICE_SITE_ACCEPTOR', ],    
    SPLICE_SITE_DONOR => [ 'SPLICE_SITE_DONOR', ],    
    SPLICE_SITE_REGION => [ 'SPLICE_SITE_REGION', ],
    EXON => [ 'INTRAGENIC', 'START_LOST', 'SYNONYMOUS_START', 'NON_SYNONYMOUS_START', 'GENE', 'TRANSCRIPT', 'EXON', 'EXON_DELETED', 'NON_SYNONYMOUS_CODING', 'SYNONYMOUS_CODING', 'FRAME_SHIFT', 'CODON_CHANGE', 'CODON_INSERTION', 'CODON_CHANGE_PLUS_CODON_INSERTION', 'CODON_DELETION', 'CODON_CHANGE_PLUS_CODON_DELETION', 'STOP_GAINED', 'SYNONYMOUS_STOP', 'STOP_LOST', 'RARE_AMINO_ACID', ],
    INTRON => [ 'INTRON', 'INTRON_CONSERVED', ],
    UTR_3_PRIME => [ 'UTR_3_PRIME', 'UTR_3_DELETED', ],
    DOWNSTREAM => [ 'DOWNSTREAM', ],
    REGULATION => [ 'REGULATION', ],
);


# get snpEff Region from effect type
sub get_effect_type_region {
    my $type = shift;
    my $region = ".";

    foreach my $key (keys %region_hash) {
        #print $key."\t";
        foreach (@{$region_hash{$key}}) {
            if($type eq $_){
                $region = $key;
            }
        }
    }
    return $region;
}

1;

