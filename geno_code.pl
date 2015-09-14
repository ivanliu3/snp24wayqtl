#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $infile = shift;

open my $fh, "<$infile" or die $!;

my $header = <$fh>;
chomp $header;
my @header = split /\t/, $header;

print join ("\t", @header[0..3]), "\t", join ("\t", @header[8..$#header]) , "\n";


while ( <$fh> ) {
    chomp;
    my %hash;
    my @row =split (/\t/, $_);
    my ($f1_female, $f1_male) = @row[2..3];
    my @f1_female_all = map lc ,(split //, $f1_female);
    my @f1_male_all = map lc, (split //, $f1_male);
    $hash{'A'} = $f1_female_all[0];
    $hash{'B'} = $f1_female_all[1];
    $hash{'C'} = $f1_male_all[0];
    $hash{'D'} = $f1_male_all[1];
    print join ("\t", @row[0..3]) ,"\t";

#    print Dumper \%hash;

    foreach my $genotype (@row[8..$#row]) {
	if ($genotype eq "NA") {
	    print "$genotype";
	} elsif ($genotype eq $hash{'A'} . $hash{'D'} && $genotype eq $hash{'B'} . $hash{'C'}) {
	    print "10";
	} elsif ($genotype eq $hash{'A'} . $hash{'C'} && $genotype eq $hash{'B'} . $hash{'D'}) {
	    print "9";
	} elsif ($genotype eq $hash{'A'} . $hash{'D'} && $genotype eq $hash{'B'} . $hash{'D'}) {
	    print "8";
	} elsif ($genotype eq $hash{'A'} . $hash{'C'} && $genotype eq $hash{'B'} . $hash{'C'}) {
	    print "7";
	} elsif ($genotype eq $hash{'B'} . $hash{'C'} && $genotype eq $hash{'B'} . $hash{'D'}) {
	    print "6";
	} elsif ($genotype eq $hash{'A'} . $hash{'C'} && $genotype eq $hash{'A'} . $hash{'D'}) {
	    print "5";
	} elsif ($genotype eq $hash{'B'} . $hash{'D'}) {
	    print "4";
	} elsif ($genotype eq $hash{'A'} . $hash{'D'}) {
	    print "3";
	} elsif ($genotype eq $hash{'B'} . $hash{'C'}) {
	    print "2";
	} elsif ($genotype eq $hash{'A'} . $hash{'C'}) {
	    print "1";
	}  else {
	    print "unknown\t"
	}
	print "\t";
    } # end of foreach
    print "\n";

}
