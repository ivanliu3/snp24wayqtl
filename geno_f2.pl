#!/usr/bin/perl

use strict;
use warnings;

my $infile = shift;

open my $fh, "<$infile" or die $!;

my $header = <$fh>;
print $header;
my @cross_type;

while ( <$fh> ) {
    my @possible_geno;
    chomp;
    my @row = split "\t", $_;
    
    my ($f1_female, $f1_male) = @row[2..3];
    push @cross_type, "$f1_female.$f1_male";
    my ($f1f_female, $f1f_male, $f1m_female, $f1m_male) = @row[4..7];
                                                                                                                                 
    
    my @f1_female_all = map lc, (split //, $f1_female);
    my @f1_male_all = map lc, (split //, $f1_male);
    my @uniq_f1_female_all = uniq(@f1_female_all);
    my @uniq_f1_male_all = uniq(@f1_male_all);

    foreach my $i ( @uniq_f1_female_all ) {
	foreach my $j ( @uniq_f1_male_all ) {
	    my $combination = $i.$j;
	    push @possible_geno, $combination;
	}
    }

    my @sorted_geno = map { join "", (sort split ("", $_)) } @possible_geno;
    my @uniq_f2 =  uniq (@row[8..$#row]);
   # print $row[0],"\t", join ("\t", @possible_geno) , "\n";
   # print $row[0],"\t", join ("\t", @sorted_geno) , "\n";
 #   print join ("\t", @row[0..3]), "\t", join ("\t",@uniq_f2), "\n";


    ### parse genotypes now
    if ( join (",", sort @f1_female_all) eq join (",",sort @f1_male_all) ) { # *1 f1 same genotype
	my @new_row = map {lc} @row;
	@new_row = map { $_ =~ /^( )*(\-)+/ ? "NA" : $_ } @new_row;
	print join ("\t", @new_row), "\n";

    } else { # *1 f1 different genotype
	## deal with f2 alleles
	#print "f2-$row[0]", "\t", join ("\t", @row[1..3]), "\t";
	my @new_row = map {lc} @row;
	@new_row = map { $_ =~ /^( )*(\-)+/ ? "NA" : $_ } @new_row;
	print join ("\t", @new_row[0..7]) , "\t";
	foreach my $f2 (@new_row[8..$#new_row]) {
	    my @f2_all = split //, $f2;
	    if ( keys %{ { map {$_, 1} @f2_all } } == 1 && (grep /$f2/, @sorted_geno) ) { ## *2 f2 homozygotes and f2 allele is mapped
		print $f2 , "\t";
	    } else { ## *2 f2 heterzgyotes or f2 homozygotes not mapped
		my $sorted_f2 = join "", sort   (split "", $f2);
		### *3 if genotype is expected
		if ( (grep /$sorted_f2/ , @sorted_geno) ) {
		    my (@index) = grep {$sorted_geno[$_] eq $sorted_f2}  0..$#sorted_geno;
		    print @possible_geno[@index], "\t" ;
		} else { ### *3 if genotype is not expected
		    print "NA\t";
		}
	    }
	
	} # end of foreach f2 samples
	print "\n";
    } # end of else ( f1 different genotypes )
}


sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}
