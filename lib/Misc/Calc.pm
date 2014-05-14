package Var::Tab;
use strict;
#use warnings;
use Carp;
use Data::Dumper;
use Vcf;
 
use Exporter qw(import);
 
our @EXPORT_OK = qw(add multiply percent_to_count_threshold);
 
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
