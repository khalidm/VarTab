package Var::Tab;
use strict;
#use warnings;
use Carp;
use Data::Dumper;
use Vcf;

use Exporter qw(import);
 
our @EXPORT_OK = qw(get_gene_counts get_ann bed_annotate_cpg bed_annotate_tfbs bed_annotate_polyphen bed_annotate_cadd bed_annotate_gwascatalog);

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


1;

