package Var::Tab;
use strict;
#use warnings;
use Carp;
use Data::Dumper;
use Vcf;
 
use Exporter qw(import);
 
our @EXPORT_OK = qw(add multiply percent_to_count_threshold maf_filter);
 
sub add {
    my ($x, $y) = @_;
    return $x + $y;
}
 
sub multiply {
    my ($x, $y) = @_;
    return $x * $y;
}

sub percent_to_count_threshold
{
    my ( $freq_threshold, $gt_size, $gt_counts ) = @_;
    #my $temp_count = shift;
    #my $gt_size = shift;
    #my $gt_counts = shift;
    my $temp_count = ( $freq_threshold/100.0 ) * $gt_size * 1.0;
    #$freq_threshold = ( $temp_count/$gt_size ) * 1.0;
    my $gt_freq_temp = ( ( $gt_counts/$gt_size ) * 1.0 );
    my $gt_freq = sprintf("%.3f", $gt_freq_temp);

    return $gt_freq; 

}

# filter 1kgp variants with maf
sub maf_filter {
    
    my $input = shift;   #Input File
	my $input_format = shift;  #Input File Format
	my $tgp_snp = shift;    #1KG snp file
	my $maf_cut_off  = shift;  #MAF cut-off
	my $output = shift;   #output

    #########
    #########
    #########
## Filter variants against 1000 genomes phase 1 with minor allele frequency...
sub snv_filter{
# Arguments taken : input file ,input format ,1000 genome file, minor allele frequency input ,output
	
	my $self = shift;
	my @ARGV = @_;
	unless (scalar @ARGV == 5){
	print <<EOF;
	This script is to filter SNVs in 1000 Phase1 project.
	snv_filter : Incorrect number of arguments ... 
EOF
	exit;	
	}
	
	my $input = $ARGV[0];   #Input File
	my $input_format = $ARGV[1];  #Input File Format
	my $tgp_snp = $ARGV[2];    #1KG snp file
	my $maf_cut_off  = $ARGV[3];  #MAF cut-off
	my $output = $ARGV[4];   #output
	
	
	my $ErrorHead = "ERROR: snv_filter"; 
	my $input_line; 
	my @snp_tabix;
	my $snp;
	my $id;
	my %saw;
	
	open(OUT,">$output")||die;
	`sed -i 's/^chr//gI' $input`;
	`sed -i -e '/^#/! s/^/chr/' $input`;

	@snp_tabix = split /\n/, `tabix $tgp_snp -B $input | awk '{FS="\t";OFS="\t"} \$4 >= $maf_cut_off'`;
	foreach $snp (@snp_tabix){
		$id = join("\t",(split /\s+/,$snp)[0..1]);
		$id =~ s/chr//;
		$saw{$id} =1;
	}
	if ($input_format =~ /vcf/i){
		open(IN,$input)|| die "Input file not found ...\n";
		while(<IN>){
		   	if (!/^#/){
		   		my @tmp = split /\s+/,$_;
		   		# only SNVs
		   		if (length($tmp[3])==1 && length($tmp[4])==1 && $tmp[3]=~ /[ATCG]/i && $tmp[4]=~ /[ATCG]/i){
		   			$tmp[0]=~ s/chr//;
		   			my $id = join("\t",$tmp[0],$tmp[1]-1);
		   			if(defined $saw{$id}){
		   			}else{
		   				$self->{DES} ->{$id} = join("\t",@tmp[0 ..7]);
		   				print OUT "chr",$tmp[0],"\t",join("\t",@tmp[1 ..6]),"\t.\n";
		   			}
		   		}
		   	}
		}
		close IN;
	}elsif($input_format =~ /bed/i){
		open(IN,$input)|| die "Input file not found ...\n";
		while(<IN>){
		   	my @tmp = split /\s+/,$_;
		   	# only SNVs
		   	if (length($tmp[3])==1 && length($tmp[4])==1 && $tmp[3]=~ /[ATCG]/i && $tmp[4]=~ /[ATCG]/i){
		   		my $id = join("\t",$tmp[0],$tmp[1]);
		   		$id =~ s/chr//g;
		   		if(defined $saw{$id}){
		   		}else{
		   			$self->{DES} -> {$id} = join("\t",@tmp[0 ..4]);
		   			print OUT join("\t",@tmp[0 ..4]),"\n";
		   		}
		   	}
		}
		close IN;
	}else{
		die "$ErrorHead : unknown format, please use vcf or bed format";	
	}
	close OUT;
}
    
    #########
    #########
    

}

