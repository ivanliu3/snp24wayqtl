#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
$Data::Dumper::Sortkeys=1;

my $infile = shift;
open my $fh, "<$infile" or die "Cannot open the file!\n";

while (<$fh>) {
    chomp;
    my @row = split /\t/, $_;
    print join ("\t", @row[0,2,3,4,5,6,7]), "\n";
    my @f2 = @row[7..$#row];
    my %counts;
    $counts{$_}++ for @f2;
    print Dumper(\%counts);
}
