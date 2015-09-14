#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $infile = shift;

open my $fh, "<$infile" or die $!;

my $header = <$fh>;
print $header;
my @cross_type;
while ( <$fh> ) {
    chomp;
    my @row = split /\t/, $_;
    my $no_miss = scalar ( grep /^\-+$/, @row);
    ### filter out data missingness
    if ( $no_miss <= 10 ) {
	my ($f1_female, $f1_male) = @row[2..3];
	push @cross_type, "$f1_female.$f1_male";
	my ($f1f_female, $f1f_male, $f1m_female, $f1m_male) = @row[4..7];
	## deal with alleles
	my @f1_female_all = map lc ,(split //, $f1_female);
	my @f1f_female_all = map lc, (split //, $f1f_female);
	my @f1f_male_all = map lc , (split //, $f1f_male);

	my @f1_male_all = map lc, (split //, $f1_male);
	my @f1m_female_all = map lc, (split //, $f1m_female);
	my @f1m_male_all = map lc, (split //, $f1m_male);
	

	if ( keys %{ { map {$_, 1} @f1_male_all } } == 1 or join (",", @f1m_female_all) eq join(",", @f1m_male_all) ) { ## Cannot phase: f1 homozygetes or f0 same genotype
	    print $row[0], "\t",$row[1], "\t", join ("\t", @row[2..$#row]), "\n";
#	    print join (",", @f1m_female_all) ,"\n";
#	    print join (",", @f1m_male_all) , "\n"; 
	} else  { ## Can be pahsed: f1 heterzygetes
	    print "mp$row[0]", "\t",$row[1];
	    my %pos_all;
	    $pos_all{'0'} = $f1_male_all[0];
	    $pos_all{'1'} = $f1_male_all[1];

	    ## go through each allele below
            foreach my $key ( '0','1' ) {
		if ( (grep /$pos_all{$key}/, @f1m_female_all) && !(grep /$pos_all{$key}/, @f1m_male_all) ) {
		    print "bingo3-$pos_all{$key}\t";
		    $pos_all{'0'} = $pos_all{$key};
		    $pos_all{'1'} = join "", grep {$_ ne $pos_all{$key}} @f1_male_all;
		    last;
		} elsif ( (grep /$pos_all{$key}/, @f1m_male_all) && !(grep /$pos_all{$key}/, @f1m_female_all) ) {
		    print "bingo4-$pos_all{$key}\t";
		    $pos_all{'1'} = $pos_all{$key};
		    $pos_all{'0'} = join "", grep {$_ ne $pos_all{$key}} @f1_male_all;
		    last;
		} else {
		    print "violate\t" if ($key eq '1');
		    next;
		}
	    } # end of first foreach

            ## print genotype
	    print $row[2], "\t";
	    print $pos_all{'0'};
	    print $pos_all{'1'};
	    print "\t";
	    print  join ("\t",@row[4..$#row]), "\n";
	   # print Dumper %pos_all ;

	}# end of else (f1 heterozygotes)
    }# end of if (data missingness)
}
 print join ("\t", uniq(@cross_type) ), "\n";

sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}	


