package Misc::PrintFile;
use strict;
#use warnings;
use Carp;
use Data::Dumper;
use Vcf;
 
use Exporter qw(import);
 
our @EXPORT_OK = qw(print_to_file print_to_var_type_single print_to_var_type);

# PRINT STRING TO OUTPUT FILE
sub print_to_file
{
    my $fh = shift;
    my $str = shift;

    print $fh $str;
}

# CALL PRINT FUNCTION BASED ON REF AND ALT FORMS
sub print_to_var_type_single
{
    #print_to_var_type($r, $a, $OUTTAB, $SNPOUTTAB, $INSOUTTAB, $DELOUTTAB, $OTHEROUTTAB, $print_string);
    my $ref = shift;
    my $alt = shift;
    my $OUTTAB = shift;
    # my $SNPOUTTAB = shift;
    # my $INSOUTTAB = shift;
    # my $DELOUTTAB = shift;
    # my $OTHEROUTTAB = shift;
    my $print_string = shift;
    # print_to_file($OUTTAB, $print_string); if() { print_to_file($SNPOUTTAB, $print_string)}; print_to_file($INSOUTTAB, $print_string);
    # print_to_file($DELOUTTAB, $print_string); print_to_file($OTHEROUTTAB, $print_string);

    if( length($ref) == 1 && length($alt) == 1 ) {
        print_to_file($OUTTAB, $print_string);
    }
    # } elsif( length($ref) == 1 && length($alt) > 1 ) {
    #     print_to_file($OUTTAB, $print_string); print_to_file($INSOUTTAB, $print_string);
    # } elsif (length($ref) > 1 && length($alt) == 1 ) {
    #     print_to_file($OUTTAB, $print_string); print_to_file($DELOUTTAB, $print_string);
    # } else {
    #     print_to_file($OUTTAB, $print_string); print_to_file($OTHEROUTTAB, $print_string);
    # }
}

# CALL PRINT FUNCTION BASED ON REF AND ALT FORMS
sub print_to_var_type
{
    #print_to_var_type($r, $a, $OUTTAB, $SNPOUTTAB, $INSOUTTAB, $DELOUTTAB, $OTHEROUTTAB, $print_string);
    my $ref = shift;
    my $alt = shift;
    my $OUTTAB = shift;
    my $SNPOUTTAB = shift;
    my $INSOUTTAB = shift;
    my $DELOUTTAB = shift;
    my $OTHEROUTTAB = shift;
    my $print_string = shift;
    # print_to_file($OUTTAB, $print_string); if() { print_to_file($SNPOUTTAB, $print_string)}; print_to_file($INSOUTTAB, $print_string);
    # print_to_file($DELOUTTAB, $print_string); print_to_file($OTHEROUTTAB, $print_string);

    if( length($ref) == 1 && length($alt) == 1 ) {
        print_to_file($OUTTAB, $print_string); print_to_file($SNPOUTTAB, $print_string);
    } elsif( length($ref) == 1 && length($alt) > 1 ) {
        print_to_file($OUTTAB, $print_string); print_to_file($INSOUTTAB, $print_string);
    } elsif (length($ref) > 1 && length($alt) == 1 ) {
        print_to_file($OUTTAB, $print_string); print_to_file($DELOUTTAB, $print_string);
    } else {
        print_to_file($OUTTAB, $print_string); print_to_file($OTHEROUTTAB, $print_string);
    }
}

1;
