package Misc::Calc;
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

    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $annotation_id = shift;
    my $maf_threshold = shift;
    my @var = ();
    my @var_array = ();
    my @maf = (0,0);
    
    @var = split('\n', $main::tabix[$annotation_id]->read($main::tabix[$annotation_id]->query( $chr, $start, $end)));
    
    if ( 0+@var > 0 ) {
        foreach my $snp (@var) {
            @var_array = split('\t', $snp);
            if($var_array[3] <= $maf_threshold) {
                $maf[0] = 1;
                $maf[1] = $var_array[3];
                #return $maf; 
            } else {
                $maf[0] = 0;
            }
        }
    } else {
        $maf[0] = 2;
    }   

    return @maf;
}

1;
