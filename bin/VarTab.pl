#!/usr/bin/env perl

use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::FeatureIO;
use Bio::DB::EMBL;
use Bio::DB::GenBank;
use Bio::Perl;
#use Bio::EnsEMBL::Registry;
use Bio::DB::Fasta;
use JSON;
use Data::Dumper;
use Carp;
use Tabix;
use File::Spec;

use Exporter qw(import);
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use lib dirname(dirname abs_path $0) . '/lib';
use Var::Tab qw( get_gene_counts get_ann bed_annotate_cpg bed_annotate_tfbs 
    bed_annotate_polyphen bed_annotate_cadd bed_annotate_gwascatalog 
    bed_annotate_hapmap get_effect_type_region get_variant_contingency_table);
use Misc::Calc qw( add multiply percent_to_count_threshold maf_filter );
use Misc::PrintFile qw(print_to_file print_to_var_type);

use strict;
#use warnings;
use Carp;
use Vcf;

my $opts = parse_params();
our %bedhash = {};
our %vargenecounts = {};
our @tabix;
our @annotation_beds;
my $dbconfig = "lib/db.config";
open (DBCONFIG, $dbconfig);
my $index = 0;

# bed_to_hash($opts);

my @annotation_beds = ();
if ( exists($$opts{annotate}) )
{
    while (my $line = <DBCONFIG>) {
        chomp($line);
        if($line =~ m/^\d/) {
            #print $line."\n";
            my @dbconfig_temp = split(' ', $line);
            $annotation_beds[$index] = @dbconfig_temp[2];
            $tabix[$index] = Tabix->new('-data' => $annotation_beds[$index]);
            $index++;
            #print "$_\n";
        }
    }
}
else 
{
    # set to default
    # $annotation_beds = "False";
}

#convert_to_tab($opts, $slice_adaptor);

#convert_to_tab($opts);
#get_gene_counts($opts);
get_variant_contingency_table($opts);

#add(2,2);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg )
    {
        croak @msg;
    }
    die
        "Usage: vcf-to-tab [OPTIONS] -v in.vcf > out.tab\n",
        "Options:\n",
        "   -h, -?, --help                   This help message.\n",
        "   -i, --iupac                      Use one-letter IUPAC codes\n",
        "   -f, --frequency                  Sample frequency threshold\n",
        "   -s, --sequence                   Print flanking sequence\n",
        "   -b, --bedfile                    Annotation file (bed file)\n",
        "   -a, --annotate                   Annotation file (bgzip bed + tabix indexed file)\n",
        "   -n, --nondbsnp                   Filter dbsnp variants (assumes the ID tag is populated in the vcf)\n",
        #"   -hapmap,                         Filter HapMap variants (HapMap CEU variants)\n",
        "   -v, --input                      Input vcf file.\n",
        "   -m, --maf                        Minor allele frequency threhold.\n",
        "   -k, --maf1kg                     1000 genome minor allele frequency threshold.\n",
        "   -r, --remove                     TODO: Remove variants seen in a population e.g. 1000 genome minor allele frequency threshold.\n",
        "   -o, --output                     Output prefix.\n",
        "\n";
}

sub testhash
{
    print "size of the hash:  " . keys( %bedhash ) . ".\n";
    my $a = $bedhash{"10:8097700"};
    print "$a\n";
}

sub parse_params
{
    if(@ARGV==0){
        error();
    }
    my $opts = { iupac=>0 };
    while (my $arg=shift(@ARGV))
    {
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        if ( $arg eq '-i' || $arg eq '--iupac' ) { $$opts{iupac}=1; next; }
        if ( $arg eq '-f' || $arg eq '--frequency' ) { $$opts{frequency}=shift(@ARGV); next; }
        if ( $arg eq '-s' || $arg eq '--sequence' ) { $$opts{sequence}=shift(@ARGV); next; }
        if ( $arg eq '-b' || $arg eq '--bedfile' ) { $$opts{bedfile}=shift(@ARGV); next; }
        #if ( $arg eq '-a' || $arg eq '--annotate' ) { $$opts{annotate}=shift(@ARGV); next; }
        if ( $arg eq '-a' || $arg eq '--annotate' ) { $$opts{annotate}=1; next; }
        if ( $arg eq '-n' || $arg eq '--nondbsnp' ) { $$opts{nondbsnp}=1; next; }
        if ( $arg eq '--hapmap' ) { $$opts{hapmap}=1; next; }
        if ( $arg eq '-v' || $arg eq '--input' ) { $$opts{input}=shift(@ARGV); next; }
        if ( $arg eq '-m' || $arg eq '--maf' ) { $$opts{maf}=shift(@ARGV); next; }
        if ( $arg eq '-k' || $arg eq '--maf1kg' ) { $$opts{maf1kg}=shift(@ARGV); next; }
        if ( $arg eq '-o' || $arg eq '--output' ) { $$opts{output}=shift(@ARGV); next; }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    
    #if ( $opts eq '-f' || $opts eq '--frequency' )
    #{
    #    my $frequency = shift(@ARGV);
    #}

    if ( $$opts{iupac} )
    {
        $$opts{iupac} = 
        {
            'GG' => 'G',
            'CC' => 'C',
            'TT' => 'T',
            'AA' => 'A',
            'GT' => 'K',
            'TG' => 'K',
            'AC' => 'M',
            'CA' => 'M',
            'CG' => 'S',
            'GC' => 'S',
            'AG' => 'R',
            'GA' => 'R',
            'AT' => 'W',
            'TA' => 'W',
            'CT' => 'Y',
            'TC' => 'Y',
            '..' => '.',
        };
    }
    return $opts;
}

sub bed_to_hash
{
    #my ($opts) = @_;
    my $opts = shift;
    
    #my $bed = shift;
    
    my $annotation_bed = "FALSE";
    
    if ( exists($$opts{bedfile}) )
    {
        $annotation_bed = $$opts{bedfile};
        #my %bed = ();        
        open(IN,'<',$annotation_bed) || die "Could not open $annotation_bed: $!\n";
        my $header = <IN>;
        while(<IN>)
        {
            chomp;
            my @line = split(/\t/, $_);
            my $myKey = "$line[0]".":"."$line[2]";
            my $myValue = $line[5];
            
            #$bed{$.} = $_;            
            #$bed->{$myKey} = $myValue;
            
            $bedhash{$myKey} = $myValue;
        }
        close(IN);
    }
    else
    {
        print "Error: no BED file found.\n";
    }
    #return %bed;
}


sub convert_to_tab
{
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
    my $outputtab = $output.".tab";
    my $snpoutputtab = $output.".snp.tab";
    my $insoutputtab = $output.".ins.tab";
    my $deloutputtab = $output.".del.tab";
    my $otheroutputtab = $output.".other.tab";
    #my $outputburden = $output.".burden";
    
    # OPEN OUTPUT FILE PREFIX.tab
    # open(OUTTAB, ">$outputtab") || die "Can't open the file: $outputtab: $!\n";
    open my $OUTTAB, '>', $outputtab or die "...$!\n";
    
    # SNP INS DEL OTHER
    open my $SNPOUTTAB, '>', $snpoutputtab or die "Could not open file $!\n";
    open my $INSOUTTAB, '>', $insoutputtab or die "Could not open file $!\n";
    open my $DELOUTTAB, '>', $deloutputtab or die "Could not open file $!\n";
    open my $OTHEROUTTAB, '>', $otheroutputtab or die "Could not open file $!\n";


    #my $vcf = Vcf->new(fh=>\*STDIN);
    my $vcf = Vcf->new(file=>$input);
    $vcf->parse_header();

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
            $header = "dbSNP_ID\tChr\tPOS\tREF\tALT\tGENE\tAA\tTYPE\tREGION\tCADD\tFUNSEQ\tPOLY_SIFT\tCONS\tTFBS\tDNASE\tCpG\tGWAS\tRMSK\tncRNA_p-gene\tCLINICAL\tMAF\tSampleFreq.HOM_REF\tSampleFreq.HET\tSampleFreq.HOM_ALT";
            # ADD FUNCTION TO PRINT TO FILE
            print_to_file($OUTTAB, $header); print_to_file($SNPOUTTAB, $header); print_to_file($INSOUTTAB, $header);
            print_to_file($DELOUTTAB, $header); print_to_file($OTHEROUTTAB, $header);
            
            if($getseq eq "T") { print OUTTAB "\tFASTA"; }
            for my $col (sort keys %{$$x{gtypes}})
            {
                print_to_file($OUTTAB, "\t".$col); print_to_file($SNPOUTTAB, "\t".$col); print_to_file($INSOUTTAB, "\t".$col);
                print_to_file($DELOUTTAB, "\t".$col); print_to_file($OTHEROUTTAB, "\t".$col);
            }
            print_to_file($OUTTAB, "\n"); print_to_file($SNPOUTTAB, "\n"); print_to_file($INSOUTTAB, "\n");
            print_to_file($DELOUTTAB, "\n"); print_to_file($OTHEROUTTAB, "\n");

            $header_printed = 1;
        }

        if( $$opts{nondbsnp} )
        {
            if($$x{ID} =~ m/^\./)
            {         
                # print_info($x, $freq_threshold, $getseq, $vcf, $output);
                if( $$opts{maf1kg} ) {
                    # get the ALT
                    my $alt;
                    for my $alt (@{$$x{ALT}}) {
                        if ( $alt eq '.' ) { $alt=$$x{REF}; }
                    }
                    @maf_output = maf_filter($$x{CHROM},$$x{POS},$$x{POS}+1,7,$$opts{maf1kg});
                    my $hapmap = bed_annotate_hapmap($$x{CHROM},$$x{POS},$$x{POS}+1, 13, $$x{REF}, $alt, $$opts{maf1kg});
                    # ONLY PRINT VARIANT TO OUTPUT IF 1KG MAF IS BELOW THRESHOLD OR NOT REPORTED IN THE BED FILE
                    # ALSO ONLY PRINT THOSE VARIANTS WITH HAPMAP MAF > 0.05
                    if ( $maf_output[0] == 1 && $hapmap == 0 ) {
                        $print_string = print_info($x, $freq_threshold, $getseq, $vcf, $output);
                        print_to_var_type($r, $a, $OUTTAB, $SNPOUTTAB, $INSOUTTAB, $DELOUTTAB, $OTHEROUTTAB, $print_string);
                        
                    } elsif ( $maf_output[0] == 2 && $hapmap == 0 ) {
                        $print_string = print_info($x, $freq_threshold, $getseq, $vcf, $output);
                        print_to_var_type($r, $a, $OUTTAB, $SNPOUTTAB, $INSOUTTAB, $DELOUTTAB, $OTHEROUTTAB, $print_string);
                    }
                } else {
                    $print_string = print_info($x, $freq_threshold, $getseq, $vcf, $output);
                    print_to_var_type($r, $a, $OUTTAB, $SNPOUTTAB, $INSOUTTAB, $DELOUTTAB, $OTHEROUTTAB, $print_string);
                }
            }
        }
        else
        {            
            if( $$opts{maf1kg} ) {
                @maf_output = maf_filter($$x{CHROM},$$x{POS},$$x{POS}+1,7,$$opts{maf1kg});
                if ( $maf_output[0] == 1 ) {
                    $print_string = print_info($x, $freq_threshold, $getseq, $vcf, $output);
                    #print OUTTAB $print_string;
                    print_to_var_type($r, $a, $OUTTAB, $SNPOUTTAB, $INSOUTTAB, $DELOUTTAB, $OTHEROUTTAB, $print_string);
                } elsif ( $maf_output[0] == 2 ) {
                    $print_string = print_info($x, $freq_threshold, $getseq, $vcf, $output);
                    #print OUTTAB $print_string;
                    print_to_var_type($r, $a, $OUTTAB, $SNPOUTTAB, $INSOUTTAB, $DELOUTTAB, $OTHEROUTTAB, $print_string);
                }
            } else {
                $print_string = print_info($x, $freq_threshold, $getseq, $vcf, $output);
                #print OUTTAB $print_string;
                print_to_var_type($r, $a, $OUTTAB, $SNPOUTTAB, $INSOUTTAB, $DELOUTTAB, $OTHEROUTTAB, $print_string);
            }
        }
    }
    #print OUTTAB "Total = $total\n";
    close $OUTTAB;
    close $SNPOUTTAB;
    close $INSOUTTAB;
    close $DELOUTTAB;
    close $OTHEROUTTAB;
    $vcf->close();
}

sub print_info
{
    my $x = shift;
    my $freq_threshold = shift;
    my $getseq = shift;
    my $vcf = shift;
    my $bed = shift;
    my $tabix = shift;
    my $output = shift;
    my $slice_adaptor;
    my $fasta_str;
    my $return_info;

    # output tab filename
    #my $outputtab = $output.".tab";
    #print $outputtab."--\n";
    # OPEN OUTPUT FILE PREFIX.tab
    #open(OUTTAB, ">$outputtab") || die "Can't open the file: $outputtab: $!\n";    

    #if (exists($$opts{sequence}))
    if ($getseq eq "T")
    {
        #my $reg = 'Bio::EnsEMBL::Registry';
        #$reg->load_registry_from_db(
        #    -host => 'ensembldb.ensembl.org',
        #    -user => 'anonymous'
        #);
        #$slice_adaptor = $reg->get_adaptor( 'human', 'core', 'slice');
    }
    
    my @pairs = $$x{ALT};
    my $alt_str = join(', ', @pairs);
    my $alt_len = scalar @{$$x{ALT}};
    my $gt_size = keys(%{$$x{gtypes}});
    
    my %types;
    for my $alt (@{$$x{ALT}})
    {
        if ( $alt eq '.' ) { $alt=$$x{REF}; }
        
        my $gt_counts = 0;
        my $gt_het = 0;
        my $gt_hom = 0;
        my $gt_string = "";
        my @gt_array = (0, 0, 0, 0); # hom_ref - het - hom_alt - sample count 
        
        for my $col (sort keys %{$$x{gtypes}})
        {
            $gt_counts++;
            my ($al1,$sep,$al2) = exists($$x{gtypes}{$col}{GT}) ? $vcf->parse_alleles($x,$col) : ('.','/','.');
            my $gt = $al1.'/'.$al2;
            my ($current_gt, $gt_index) = get_gt_type($$x{gtypes}{$col}{GT});
            
            $gt_string = "$gt_string\t$current_gt";
            
            $gt_array[$gt_index]++;
                             
            if($alt eq $al2)
            {
                if($$x{gtypes}{$col}{GT} eq "0/1")
                {
                     $gt_het++;
                }
                elsif($$x{gtypes}{$col}{GT} eq "1/1")
                {
                    $gt_hom++;
                }
            } else 
            {
                #$gt_string = "$gt_string\t*$alt/$al2";
                #$gt_string = "$gt_string\t.";
            }
        }

        my $prop_het = ( $gt_array[1]/$gt_counts ) * 1.0;
        my $prop_hom_ref = ( $gt_array[0]/$gt_counts ) * 1.0;
        my $prop_hom_alt = ( $gt_array[2]/$gt_counts ) * 1.0;        
        my $maf = ((($gt_array[0] * 0.0) + ($gt_array[1] * 1.0) + ($gt_array[2] * 2.0))/$gt_counts ) * 1.0;
        
        my $prop_homr_str = sprintf("%.3f", $prop_hom_ref);
        my $prop_het_str = sprintf("%.3f", $prop_het);
        my $prop_homa_str = sprintf("%.3f", $prop_hom_alt);
        my $maf_str = sprintf("%.3f", $maf);
        
        my $temp_count = ( $freq_threshold/100.0 ) * $gt_size * 1.0;
        #$freq_threshold = ( $temp_count/$gt_size ) * 1.0;
        my $gt_freq_temp = ( ( $gt_counts/$gt_size ) * 1.0 );
        my $gt_freq = sprintf("%.3f", $gt_freq_temp);
        #$gt_freq = percent_to_count_threshold($freq_threshold, $gt_size, $gt_counts);

        #print " $gt_array[0] \t $gt_array[1] \t $gt_array[2] \t $gt_array[3] \t $gt_counts \t $temp_count \n";
        
        #if( $gt_counts <= $temp_count )        
        if( ($gt_counts-$gt_array[3]) <= $temp_count )
        {
            # TODO: make the feature of getting flanking fasta sequences optional
            if ($getseq eq "T")
            {   
                my $db = Bio::DB::Fasta->new('genome.fa');
                #my $slice = $db->seq($$x{CHROM}, $$x{POS}-60 => $$x{POS}+60);
                my $slice = $db->seq($$x{CHROM}, $$x{POS}-60 => $$x{POS}+60);
                #my $seq_up = $slice->subseq(1 => 60);
                #my $seq_down = $slice->subseq(62 => 121);
                my $seq_up = substr $slice, 0, 60;
                my $seq_down = substr $slice, 61, 121;
                $fasta_str = $seq_up."[".$$x{REF}."/".$alt."]".$seq_down;

                #$fasta_str = $seq_up."[".$$x{REF}."/".$alt."]".$seq_down;
            }

            my $info_string = "";
            my ($fun_gene, $fun_type, $fun_network, $fun_tfp, $fun_dnase, $highlight) = parse_fun_bed($x);
            
            #other annotations
            my @gene_info = split(':',get_annotations($x, "GENEINFO"));
            my $snpeff_gene = get_annotations($x, "SNPEFF_GENE_NAME");
            my $gene_clndbn = get_annotations($x, "CLNDBN");
            my $aa_change = get_annotations($x, "SNPEFF_AMINO_ACID_CHANGE");
            my $phastcons = get_annotations($x, "PhastCons");
            # TO-DO
            my ($snpeff_type, $snpeff_gene, $snpeff_aa_change) = get_snpEffannotations($x);

            # 
            my $snpeff_region = get_effect_type_region($snpeff_type);

            # TO-DO ADD 
            my $rmsk = bed_annotate($$x{CHROM},$$x{POS}-1,$$x{POS}, 0); # RMSK
            $rmsk =~ s/\;$//g;
            #my $rmsk = ".";
            #my $gwas = bed_annotate($$x{CHROM},$$x{POS}-1,$$x{POS}, 1);
            my $cpg = bed_annotate_cpg($$x{CHROM},$$x{POS}-1,$$x{POS}, 1); # CPG
            my $clinvar = bed_annotate($$x{CHROM},$$x{POS},$$x{POS}+1, 2); # CLINVAR
            my $gwascatalog = bed_annotate_gwascatalog($$x{CHROM},$$x{POS},$$x{POS}+1, 3); # GWAS
            
            # ENCODE 0-based
            #my $encode_tfbs = bed_annotate_tfbs($$x{CHROM},$$x{POS}-1,$$x{POS}, 5);
            my $encode_tfbs = bed_annotate_tfbs($$x{CHROM},$$x{POS}-1,$$x{POS}, File::Spec->rel2abs($annotation_beds[4])); # TFBS
            #my $encode_tfbs = bed_annotate_tfbs($$x{CHROM},$$x{POS},$$x{POS}+1, File::Spec->rel2abs($annotation_beds[4])); # TFBS
            $encode_tfbs =~ s/^\;//g;
            #my $encode_dnase = bed_annotate($$x{CHROM},$$x{POS},$$x{POS}+1, 5);
            my $encode_dnase = bed_annotate($$x{CHROM},$$x{POS}-1,$$x{POS}, 5);
           
            #my $polyphen_chr = $$x{CHROM};
            #$polyphen_chr =~ s/^chr//g;
            my $polyphen = bed_annotate_polyphen($$x{CHROM},$$x{POS},$$x{POS}+1, File::Spec->rel2abs($annotation_beds[6]), $$x{REF}, $alt);
            
            #my $gerpelement = bed_annotate($$x{CHROM},$$x{POS},$$x{POS}+1, 8); # GERP ELEMENTS
            my $gerpelement = bed_annotate($$x{CHROM},$$x{POS}-1,$$x{POS}, 8); # GERP ELEMENTS
            
            my $ncrna = bed_annotate($$x{CHROM},$$x{POS},$$x{POS}+1, 9); # ncRNA + p_genes
            my $ncsensitive = bed_annotate($$x{CHROM},$$x{POS},$$x{POS}+1, 10); # ncRNA + p_genes
            # CADD
            #my $cadd = bed_annotate($$x{CHROM},$$x{POS},$$x{POS}+1, 11); # cadd score
            my $cadd = bed_annotate_cadd($$x{CHROM},$$x{POS},$$x{POS}+1, File::Spec->rel2abs($annotation_beds[11]), $$x{REF}, $alt);
            
            # HapMap
            #my $hapmap = bed_annotate_hapmap($$x{CHROM},$$x{POS}-1,$$x{POS}, 13, $$x{REF}, $alt);
            #print "TEST = ".$hapmap."\n";
            
            # data/hg19.funseq.nc_sensitive.bed.gz
            
            # check gene name from several sources
            if($fun_gene eq ".")
            {
                $info_string = $snpeff_gene;
            }
            else
            {
                $info_string = $fun_gene;
            }

            if ($info_string eq "."){
                $info_string = bed_annotate($$x{CHROM},$$x{POS},$$x{POS}+1, 12);
                if($info_string ne "."){
                    $info_string = $info_string."(p)";
                }
            }
            
            $return_info = "$$x{ID}\t$$x{CHROM}\t$$x{POS}\t$$x{REF}\t$alt";
            $return_info = $return_info."\t$info_string";
            $return_info = $return_info."\t$snpeff_aa_change";

		    $return_info = $return_info."\t$snpeff_type";
		    $return_info = $return_info."\t$snpeff_region";
		    
            $return_info = $return_info."\t$cadd";
		    $return_info = $return_info."\t$ncsensitive";

            $return_info = $return_info."\t$polyphen";
            $return_info = $return_info."\t$gerpelement";
		    
            # $return_info = $return_info."\t$fun_network";
            
		    $return_info = $return_info."\t$encode_tfbs";
            $return_info = $return_info."\t$encode_dnase";
		    
            if($getseq eq "T") { $return_info = $return_info."\t$fasta_str"; }
		    
            $return_info = $return_info."\t$cpg";
            $return_info = $return_info."\t$gwascatalog";
		    $return_info = $return_info."\t$rmsk";
            #$return_info = $return_info."\t$gwas:$gwascatalog";
            #$return_info = $return_info."\t$polyphen";
            $return_info = $return_info."\t$ncrna";
            #$return_info = $return_info."\t$gene_clndbn:$clinvar";
            $return_info = $return_info."\t$clinvar";
            
            $return_info = $return_info."\t$maf_str";
            $return_info = $return_info."\t$prop_homr_str";
            $return_info = $return_info."\t$prop_het_str";
            $return_info = $return_info."\t$prop_homa_str";
            
            #if($highlight eq "T") { $return_info = $return_info."\t*"; } else { $return_info = $return_info."\t."; }
            $return_info = $return_info."$gt_string";
            $return_info = $return_info."\n";
            
        }
    }
    #close OUTTAB;
    return $return_info;
}

# PRINT HTML
sub print_info_html
{
    exit;
} 

# PARSE FUNSEQ OUTPUT IN BED FORMAT
# TO-DO CHANGE TO TABIX 
sub parse_fun_bed
{
    #print "size of the hash:  " . keys( %bedhash ) . ".\n";
    #my $a = $bedhash{"10:8097700"};
    #print "$a\n";

    my $x = shift;
    my $output = "";

    my $chr_temp = $$x{CHROM};
    my $pos = $$x{POS};
    #$chr = ~ s/^chr/""/g;
    my $chr = substr $chr_temp, 3;
    my $k = $chr.":".$pos;

    my $info = $bedhash{$k};
    my @values = split(';', $info);

    #print $$x{CHROM}.":".$$x{POS}."->".$k."->".$info;
    #print $k."->".$info."->".$values[8].",".$values[0];
    #exit;
    
    my $type = ".";
    my $tfp = ".";
    my $tfm = ".";
    my $network = ".";
    my $gene = ".";
    my $dnase = ".";

    # coding
    if($values[0] eq "Yes")
    {
        $type = "coding";
    }
    elsif($values[0] eq "No")
    {
        $type = "non-coding";
    }
    else
    {
        $type = ".";
    }

    # network
    if ($values[2] ne "")
    {
        $network = $values[2];
    }
    # gene
    if ($values[8] ne "")
    {
        $gene = $values[8];
    }

    # encode
    #if(index($values[4], "TFP(") != -1)
    my @encode_values = split(',', $values[4]);
    my @tfp_values;
    my @dhs_values;
    if($encode_values[0] =~ /^TFP/)
    {
        @tfp_values = split('\|', $encode_values[0]);
        #$tfp = "$tfp_values[0]";
        $tfp = $tfp_values[0];
        $tfp =~ s/^TFP\(//g;
        #$tfp = "transcription factor binding peak";
    }
    elsif($encode_values[0] =~ /^DHS/)
    {
        #$dnase = "DNase1 hypersensitive sites";
        @dhs_values = split('\|', $encode_values[0]);
        $dnase = $dhs_values[0];
        $dnase =~ s/^DHS\(//g;
        #$tfp = "transcription factor binding peak";

    }
    else
    {
        $tfp = ".";
        $dnase = ".";
    }

    my $highlight = "F";
    if($values[9] > 1 && $values[9] ne ".")
    {
        $highlight = "T";
    }
    elsif($values[10] > 1 && $values[10] ne ".")
    {
        $highlight = "T";
    }

    #$output = $gene."\t".$type."\t".$network."\t".$tfp."\t".$dnase;
    # BUG: TYPE ISSUE
    return ($gene, $type, $network, $tfp, $dnase, $highlight);    
    #return $output;
}

sub get_annotations
{
    my $x = shift;
    my $annotation_name = shift;

    my $str = "";
    #$str = $$x{"$annotation_name"};
    $str = $$x{INFO}{$annotation_name};
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

sub get_snpEffannotations
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
    return ($eff_type, $eff_gene, $eff_aa_change);
}

sub get_gt_type
{
    my $gt = shift;
    if($gt eq "0/0")
    {
        return ("HOM_REF",0);
    }
    elsif($gt eq "0/1")
    {
        return ("HET",1);
    }
    elsif($gt =~ /^1\//)
    {
         return ("HOM_ALT",2);
    }
    else
    {
        return (".",3);
    }
}

sub bed_annotate
{
    # my $annotation_tabix = shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $annotation_id = shift;
    my @var = ();
    my @var_array = ();
    my $out_str = "";
        
    #@var =split('\t', $tabix[$annotation_id]->read($tabix[$annotation_id]->query( $chr, $start, $end)));
    #@snp_tabix = split /\n/, `tabix $tgp_snp -B $input | awk '{FS="\t";OFS="\t"} \$4 >= $maf_cut_off'`;	
    @var = split('\n', $tabix[$annotation_id]->read($tabix[$annotation_id]->query( $chr, $start, $end)));
    
    # return "." if empty tabix return empty string
    if (0+@var == 0) {
        return ".";
    } else {
        # split tabix returned string
        foreach my $snp (@var){
            @var_array = split('\t', $snp);
            if (0+@var_array > 0) {
                if($var_array[3] ne ""){
                    return ($var_array[3]);
                }
                else{
                    return (".");
                }
            }
            else
            {
                return (".");
            }
        }
    } 
}

# # CALL PRINT FUNCTION BASED ON REF AND ALT FORMS
# sub print_to_var_type
# {
#     #print_to_var_type($r, $a, $OUTTAB, $SNPOUTTAB, $INSOUTTAB, $DELOUTTAB, $OTHEROUTTAB, $print_string);
#     my $ref = shift;
#     my $alt = shift;
#     my $OUTTAB = shift;
#     my $SNPOUTTAB = shift;
#     my $INSOUTTAB = shift;
#     my $DELOUTTAB = shift;
#     my $OTHEROUTTAB = shift;
#     my $print_string = shift;
#     # print_to_file($OUTTAB, $print_string); if() { print_to_file($SNPOUTTAB, $print_string)}; print_to_file($INSOUTTAB, $print_string);
#     # print_to_file($DELOUTTAB, $print_string); print_to_file($OTHEROUTTAB, $print_string);

#     if( length($ref) == 1 && length($alt) == 1 ) {
#         print_to_file($OUTTAB, $print_string); print_to_file($SNPOUTTAB, $print_string);
#     } elsif( length($ref) == 1 && length($alt) > 1 ) {
#         print_to_file($OUTTAB, $print_string); print_to_file($INSOUTTAB, $print_string);
#     } elsif (length($ref) > 1 && length($alt) == 1 ) {
#         print_to_file($OUTTAB, $print_string); print_to_file($DELOUTTAB, $print_string);
#     } else {
#         print_to_file($OUTTAB, $print_string); print_to_file($OTHEROUTTAB, $print_string);
#     }
# }

# # PRINT STRING TO OUTPUT FILE
# sub print_to_file
# {
#     my $fh = shift;
#     my $str = shift;

#     print $fh $str;
# }



