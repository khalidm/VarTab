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
use Var::Tab qw( get_gene_counts get_ann bed_annotate_cpg bed_annotate_tfbs bed_annotate_polyphen bed_annotate_gwascatalog );
use Misc::Calc qw( percent_to_count_threshold );

use strict;
#use warnings;
use Carp;
use Vcf;

my $opts = parse_params();
our %bedhash = {};
our %vargenecounts = {};
our @tabix;
our @annotation_beds;

bed_to_hash($opts);

my @annotation_beds = ();
if ( exists($$opts{annotate}) )
{
    @annotation_beds = split(',', $$opts{annotate});
    if ( -e $annotation_beds[0] && -e $annotation_beds[1] && -e $annotation_beds[2])
    {
        $tabix[0] = Tabix->new('-data' => $annotation_beds[0]); # rmsk
        #$tabix[1] = Tabix->new('-data' => $annotation_beds[1]); # gwas
        $tabix[1] = Tabix->new('-data' => $annotation_beds[1]); # cpg
        $tabix[2] = Tabix->new('-data' => $annotation_beds[2]); # clinvar
        $tabix[3] = Tabix->new('-data' => $annotation_beds[3]); # gwas catalog
        $tabix[4] = Tabix->new('-data' => $annotation_beds[4]); # tfbs
        $tabix[5] = Tabix->new('-data' => $annotation_beds[5]); # dnase
        $tabix[6] = Tabix->new('-data' => $annotation_beds[6]); # polyphen
    }
    else
    {
        print "Error: annotaion file $annotation_beds[0] not found. Please check if the file is bgzipped and tabix index.\n";
    }
    # bed_annotate($annotation_bed);
}
else 
{
    # set to default
    # $annotation_beds = "False";
}

#convert_to_tab($opts, $slice_adaptor);
convert_to_tab($opts);
get_gene_counts($opts);
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
        "   -n, --nondbsnp                   Keep only non-dbsnp variants (assumes the ID tag is populated in the vcf)\n",
        "   -v, --input                      Input vcf file.\n",
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
    my $opts = { iupac=>0 };
    while (my $arg=shift(@ARGV))
    {
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        if ( $arg eq '-i' || $arg eq '--iupac' ) { $$opts{iupac}=1; next; }
        if ( $arg eq '-f' || $arg eq '--frequency' ) { $$opts{frequency}=shift(@ARGV); next; }
        if ( $arg eq '-s' || $arg eq '--sequence' ) { $$opts{sequence}=shift(@ARGV); next; }
        if ( $arg eq '-b' || $arg eq '--bedfile' ) { $$opts{bedfile}=shift(@ARGV); next; }
        if ( $arg eq '-a' || $arg eq '--annotate' ) { $$opts{annotate}=shift(@ARGV); next; }
        if ( $arg eq '-n' || $arg eq '--nondbsnp' ) { $$opts{nondbsnp}=1; next; }
        if ( $arg eq '-v' || $arg eq '--input' ) { $$opts{input}=shift(@ARGV); next; }
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
    #my $outputburden = $output.".burden";
    # OPEN OUTPUT FILE PREFIX.tab
    open(OUTTAB, ">$outputtab") || die "Can't open the file: $outputtab: $!\n";

    #my $vcf = Vcf->new(fh=>\*STDIN);
    my $vcf = Vcf->new(file=>$input);
    $vcf->parse_header();

    my $header_printed=0;
    my $total = 0;
    my $print_string;

    while (my $x=$vcf->next_data_hash())
    {
        $total++;
        # print Dumper($x);
        if ( !$header_printed ) 
        {
            #print "#CHROM,POS,REF,ALT,ID,FREQ";
            #print "dbSNP_ID\tChr\tPOS\tREF\tALT\tCount\tFreq\tGENE\tTYPE\tNETWORK\tTFP\tDNASE\tp_HET\tp_HOM";
            #NO HET HOM percentage
            #print "dbSNP_ID\tChr\tPOS\tREF\tALT\tCount\tFreq\tGENE\tTYPE\tNETWORK\tTFP\tDNASE";
            #NO EARLY COUNTS
            print OUTTAB "dbSNP_ID\tChr\tPOS\tREF\tALT\tGENE\tTYPE\tNETWORK\tTF_binding_peak\tDNASE";
            print OUTTAB "\tCLINICAL\tAA_CHANGE\tConservation\tRMSK\tCpG\tGWAS\tPOLYPHEN\tSampleFreq.HOM_REF\tSampleFreq.HET\tSampleFreq.HOM_ALT\tFS";
            if($getseq eq "T") { print OUTTAB "\tFASTA"; }
            for my $col (sort keys %{$$x{gtypes}})
            {
                print OUTTAB "\t$col";
            }
            print OUTTAB "\n";

            $header_printed = 1;
        }

        if( $$opts{nondbsnp} )
        {
            if($$x{ID} =~ m/^\./)
            {         
                print_info($x, $freq_threshold, $getseq, $vcf, $output);
            }
        }
        else
        {            
            #print_info($x, $freq_threshold, $getseq, $vcf, $output);
            $print_string = print_info($x, $freq_threshold, $getseq, $vcf, $output);
            print OUTTAB $print_string;
        }
    }
    print OUTTAB "Total = $total\n";
    close OUTTAB;
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
        my $prop_homr_str = sprintf("%.3f", $prop_hom_ref);
        my $prop_het_str = sprintf("%.3f", $prop_het);
        my $prop_homa_str = sprintf("%.3f", $prop_hom_alt);
        
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
            # TO-DO ADD 
            my $rmsk = bed_annotate($$x{CHROM},$$x{POS}-1,$$x{POS}, 0);
            $rmsk =~ s/\;$//g;
            #my $rmsk = ".";
            #my $gwas = bed_annotate($$x{CHROM},$$x{POS}-1,$$x{POS}, 1);
            my $cpg = bed_annotate_cpg($$x{CHROM},$$x{POS}-1,$$x{POS}, 1);
            my $clinvar = bed_annotate($$x{CHROM},$$x{POS},$$x{POS}+1, 2);
            my $gwascatalog = bed_annotate_gwascatalog($$x{CHROM},$$x{POS},$$x{POS}+1, 3);
            #my $encode_tfbs = bed_annotate_tfbs($$x{CHROM},$$x{POS}-1,$$x{POS}, 5);
            my $encode_tfbs = bed_annotate_tfbs($$x{CHROM},$$x{POS},$$x{POS}+1, File::Spec->rel2abs($annotation_beds[4]));
            $encode_tfbs =~ s/^\;//g;
            my $encode_dnase = bed_annotate($$x{CHROM},$$x{POS},$$x{POS}+1, 5);
            my $polyphen_chr = $$x{CHROM};
            $polyphen_chr =~ s/^chr//g;
            my $polyphen = bed_annotate_polyphen($polyphen_chr,$$x{POS},$$x{POS}+1, File::Spec->rel2abs($annotation_beds[6]), $$x{REF}, $alt);
            
            # check gene name from several sources
            if($fun_gene eq ".")
            {
                $info_string = $snpeff_gene;
            }
            else
            {
                $info_string = $fun_gene;
            }
            
            #my $info_string = "$gene_info[0]|$snpeff_gene|$fun_gene";
            #
            #
		    $return_info = "$$x{ID}\t$$x{CHROM}\t$$x{POS}\t$$x{REF}\t$alt";
            $return_info = $return_info."\t$info_string";
		    $return_info = $return_info."\t$snpeff_type";
		    $return_info = $return_info."\t$fun_network";
            
            #$return_info = $return_info."\t$fun_tfp";
		    $return_info = $return_info."\t$encode_tfbs";
		    
            #$return_info = $return_info."\t$fun_dnase";
            $return_info = $return_info."\t$encode_dnase";
		    
            if($getseq eq "T") { $return_info = $return_info."\t$fasta_str"; }
		    $return_info = $return_info."\t$gene_clndbn:$clinvar";
		    $return_info = $return_info."\t$snpeff_aa_change";
		    $return_info = $return_info."\t$phastcons";
		    $return_info = $return_info."\t$rmsk";
		    $return_info = $return_info."\t$cpg";
            $return_info = $return_info."\t$gwascatalog";
            #$return_info = $return_info."\t$gwas:$gwascatalog";
            $return_info = $return_info."\t$polyphen";
		    $return_info = $return_info."\t$prop_homr_str";
		    $return_info = $return_info."\t$prop_het_str";
            $return_info = $return_info."\t$prop_homa_str";
            if($highlight eq "T") { $return_info = $return_info."\t*"; } else { $return_info = $return_info."\t."; }
            $return_info = $return_info."$gt_string";
            $return_info = $return_info."\n";
            #
            ##print OUTTAB "$$x{ID}\t$$x{CHROM}\t$$x{POS}\t$$x{REF}\t$alt";
            ###print OUTTAB "\t$gt_counts";
            ###print OUTTAB "\t$gt_freq";
            ##print OUTTAB "\t$info_string";
            
            ###BUG ISSUE #3 print OUTTAB "\t$fun_type:";
            ##print "\t$snpeff_type";
            
            ##print OUTTAB "\t$fun_network";
            ##print OUTTAB "\t$fun_tfp";
            ##print OUTTAB "\t$fun_dnase";
            ##if($getseq eq "T") { print OUTTAB "\t$fasta_str"; }
            ##print OUTTAB "\t$gene_clndbn:$clinvar";
            ###print OUTTAB ",$snpeff_gene";
            ##print OUTTAB "\t$snpeff_aa_change";
            ##print OUTTAB "\t$phastcons";
            ##print OUTTAB "\t$rmsk";
            ##print OUTTAB "\t$cpg";
            
            ###print OUTTAB "\t$gwas";
            ##print OUTTAB "\t$gwas:$gwascatalog";
            
            ##print OUTTAB "\t$prop_homr_str";
            ##print OUTTAB "\t$prop_het_str";
            ##print OUTTAB "\t$prop_homa_str";
            ##if($highlight eq "T") { print OUTTAB "\t*"; } else { print OUTTAB "\t."; }
            ##print OUTTAB "$gt_string";
            ###print OUTTAB "\t$snpeff_gene";
            ###print OUTTAB "\t$gt_array[3]";
            ##print OUTTAB "\n";
        }
    }
    #close OUTTAB;
    return $return_info;
}

# PRINT HTML
sub print_info_html
{
    my $x = shift;
    my $freq_threshold = shift;
    my $getseq = shift;
    my $vcf = shift;
    my $bed = shift;
    my $slice_adaptor;
    my $fasta_str;

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
            }
        }
        
        my $prop_het = ( $gt_array[1]/$gt_counts ) * 1.0;
        my $prop_hom_ref = ( $gt_array[0]/$gt_counts ) * 1.0;
        my $prop_hom_alt = ( $gt_array[2]/$gt_counts ) * 1.0;        
        my $prop_homr_str = sprintf("%.3f", $prop_hom_ref);
        my $prop_het_str = sprintf("%.3f", $prop_het);
        my $prop_homa_str = sprintf("%.3f", $prop_hom_alt);
        
        my $temp_count = ( $freq_threshold/100.0 ) * $gt_size * 1.0;
        #$freq_threshold = ( $temp_count/$gt_size ) * 1.0;
        my $gt_freq_temp = ( ( $gt_counts/$gt_size ) * 1.0 );
        my $gt_freq = sprintf("%.3f", $gt_freq_temp);

        #print " $gt_array[0] \t $gt_array[1] \t $gt_array[2] \t $gt_array[3] \t $gt_counts \t $temp_count \n";
        
        if( ($gt_counts-$gt_array[3]) <= $temp_count )
        {
            # TODO: make the feature of getting flanking fasta sequences optional
            if ($getseq eq "T")
            {
                my $db = Bio::DB::Fasta->new('genome.fa');
                my $slice = $db->seq($$x{CHROM}, $$x{POS}-60 => $$x{POS}+60);
                my $seq_up = substr $slice, 0, 60;
                my $seq_down = substr $slice, 61, 121;
                $fasta_str = $seq_up."[".$$x{REF}."/".$alt."]".$seq_down;
            }

            my $info_string = "";
            my ($fun_gene, $fun_type, $fun_network, $fun_tfp, $fun_dnase, $highlight) = parse_fun_bed($x);
            
            #other annotations
            my @gene_info = split(':',get_annotations($x, "GENEINFO"));
            my $snpeff_gene = get_annotations($x, "SNPEFF_GENE_NAME");
            my $gene_clndbn = get_annotations($x, "CLNDBN");
            my $aa_change = get_annotations($x, "SNPEFF_AMINO_ACID_CHANGE");
            my $phastcons = get_annotations($x, "PhastCons");
            my ($snpeff_type, $snpeff_gene, $snpeff_aa_change) = get_snpEffannotations($x);

            # check gene name from several sources
            if($fun_gene eq ".")
            {
                $info_string = $snpeff_gene;
            }
            else
            {
                $info_string = $fun_gene;
            }
            

            print "$$x{ID}\t$$x{CHROM}\t$$x{POS}\t$$x{REF}\t$alt";
            #print "\t$gt_counts";
            #print "\t$gt_freq";
            print "\t$info_string";
            
            print "\t$fun_type";
            print ":$snpeff_type";
            
            print "\t$fun_network";
            print "\t$fun_tfp";
            print "\t$fun_dnase";
            if($getseq eq "T") { print "\t$fasta_str"; }
            print "\t$gene_clndbn";
            #print ",$snpeff_gene";
            print "\t$snpeff_aa_change";
            print "\t$phastcons";
            print "\t$prop_homr_str";
            print "\t$prop_het_str";
            print "\t$prop_homa_str";
            if($highlight eq "T") { print "\t*"; } else { print "\t."; }
            print "$gt_string";
            #print "\t$snpeff_gene";
            #print "\t$gt_array[3]";
            print "\n";
        }
    }
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

