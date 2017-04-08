#!usr/bin/perl -w
use strict;
#use Statistics::ChiSquare;
use Statistics::Distributions;
use Statistics::Multtest;

use Statistics::Multtest qw(bonferroni holm hommel hochberg BH BY qvalue);
use Statistics::Multtest qw(:all);

die "Usage: perl $0 [1.aln.fa (DNA or amino acid)] [2.aln.fa] [outgroup] [output]\n" unless (@ARGV == 4);
open (FA1, $ARGV[0]) or die "$ARGV[0] $!\n";
open (FA2, $ARGV[1]) or die "$ARGV[1] $!\n";
open (FA3, $ARGV[2]) or die "$ARGV[2] $!\n";

open OUT, ">$ARGV[3]" or die "$ARGV[3] $!\n";
$/=">";
<FA1>;
<FA2>;
<FA3>;

my $outgroup = <FA3>;
chomp $outgroup;
my @line = split/\n+/, $outgroup; 
my $outname = shift @line;
my $outseq = join "", @line;
my @outseq = split /\s*/, $outseq;

my (@group1, @group2);
my (@group1_spp, @group2_spp);
while(<FA1>){
	my $sppa = $_;
	chomp $sppa;
	@line = split/\n+/, $sppa;
	my $sppaname = shift @line;
	my $sppaseq = join "", @line;
	push @group1, $sppaseq;
	push @group1_spp, $sppaname;
}
close FA1;

while(<FA2>){
	my $sppb = $_;
	chomp $sppb;
	@line = split /\n+/, $sppb;
	my $sppbname = shift @line;
	my $sppbseq = join "", @line;
	push @group2, $sppbseq;
	push @group2_spp, $sppbname;
}
close FA2;

my $sig = 0;
my $total = 0;
my $pval; 
foreach my $seq1 (0..$#group1){
	foreach my $seq2 (0..$#group2){

		$total ++;
		my @seq1 = split/\s*/, $group1[$seq1];
		my @seq2 = split/\s*/, $group2[$seq2];
		my $diffA = 0;
		my $diffB = 0; 

#print "$seq1\n$seq2\n$outseq\n";
		foreach my $k (0..$#outseq){
#print "***$outseq[$k]\t$seq1[$k]\t$seq2[$k]\n";
			if(($outseq[$k] ne $seq1[$k])&&($outseq[$k] eq $seq2[$k])){
				$diffA ++;
#				print "\n###$k\n";
			}

			if(($outseq[$k] ne $seq2[$k])&&($outseq[$k] eq $seq1[$k])){
				$diffB ++;
#print "";
			}
		}
		my $chi = (($diffA - $diffB)**2) / ($diffA + $diffB);
#		print "$chi\n";

		my $chisprob = Statistics::Distributions::chisqrprob (1,$chi);
		my $hash_key = $group1_spp[$seq1]." ".$group2_spp[$seq2];
		$$pval{$hash_key} = $chisprob;
#print OUT "$group1_spp[$seq1]\t$group2_spp[$seq2]\t$chisprob\n";
#		if($chi >= 3.841 ){
#		if($chisprob <= 0.05 ){
#			$sig ++;
#		}
		
	}
}

#print "%$pval\n";
my $holm_pval = holm($pval);
#print "@$holm_pval\n";

foreach my $p (sort keys %$holm_pval){
#	print OUT "### $p\t$$pval{$p}\n";
	print OUT "$p\t$$holm_pval{$p}\n";
	if($$holm_pval{$p} <= 0.05 ){
		$sig ++;
	}
}
print "Group1: $ARGV[0]\nGroup2: $ARGV[1]\nOutgroup: $ARGV[2]\n $sig - Significant comparisions\n $total - Total number of comparisions\n";
#print OUT "Group1: $ARGV[0]\nGroup2: $ARGV[1]\nOutgroup: $ARGV[2]\n $sig - Significant comparisions\n $total - Total number of comparisions\n"; 
close OUT;
print "DONE!";
