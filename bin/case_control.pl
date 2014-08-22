#!/usr/bin/env perl

use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::FeatureIO;
use Bio::DB::EMBL;
use Bio::DB::GenBank;
use Bio::Perl;
use Bio::DB::Fasta;
use JSON;
use Data::Dumper;
use Carp;
use Tabix;

use strict;

our %caseshash = {};
our %controlhash = {};

my $cases = $ARGV[0];
my $control = $ARGV[1];

open(IN,'<',$cases) || die "Could not open $cases: $!\n";
#my $header = <IN>;
while(<IN>) {
    chomp;
    my @line = split(/\t/, $_);
    my $myKey = $line[0];
    my $myValue = $line[1];
    $caseshash{$myKey} = $myValue;
}
close(IN);

open(IN2,'<',$control) || die "Could not open $control: $!\n";
#my $header = <IN2>;
while(<IN2>) {
    chomp;
    my @line = split(/\t/, $_);
    my $myKey = $line[0];
    my $myValue = $line[1];
    $controlhash{$myKey} = $myValue;
}
close(IN2);

foreach my $k ( keys %caseshash )
{
    if (exists $controlhash{$k}) {
        print "$k\t$caseshash{$k}\t$controlhash{$k}\n";
    } else {
        #print "$k\t$caseshash{$k}\t0.0\n";
    }    
}


