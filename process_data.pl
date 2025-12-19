#!/usr/bin/perl

# ./process_data.pl amino_acid_genotypes_to_brightness.tsv 2 > gfp.txt
my $shift = $ARGV[1];
open(IN, $ARGV[0]) or die "Can't open file $ARGV[0] for reading: $!\n";
$_ = <IN>;
print "Genotype\tPhenotype\n";
while(<IN>) {
    s/[\r\n]//g;
    ($gen, $foo, $fitness) = split(/\t/);
#    if($gen !~ /\*/) {
	if($gen eq "") {
	    $new_gen = "wt";
	} else {
	    @singles = split(/\:/, $gen);
	    for my $i(0 .. $#singles) {
		$length = length($singles[$i]);
		$wt_letter = substr($singles[$i], 1, 1);
		# Correction for the numeration shift of two residues
		$num = substr($singles[$i], 2, $length - 3) + $shift;
		$letter = substr($singles[$i], -1);
		$singles[$i] = $wt_letter.$num.$letter;
	    }
	    $new_gen = join(':', @singles);
	}
	print "$new_gen\t$fitness\n";
#    }
}

