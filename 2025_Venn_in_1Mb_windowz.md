# Overlapping in 1 Mb windowz

I'm interested in checking how many positions that are heterozygous in each parent are also heterozygous in sons and daughters. I'd like to tabulate this for clivii in 1 Mb windowz.

Here is a script that reads in a list of heterozygous positions in one parent for one chromosome and another list for the same chromosome that has positions that are heterozygous in some number of offspring. 

For the X. clivii family, we have a weird sex ratio skew where 17 of 18 offspring are male. One explanation is that the mother was WZ and the father was ZY. This would generate 25% WZ daughters (4.5 out of 18 expected; 1 observed) and 75% WY, ZZ, and ZY sons (13.5 expected; 17 observed). We'd expect lots (~all discounting genotype errors) of the maternal heterozygous positions to be heterozygous in the daughter and many to also be heterozygous in sons. We'd expect some paternal het positions to also be heteroz in sons.

Here's a script to do this:
```
#!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
use List::MoreUtils qw/ uniq /;

#  This program reads in two files that have lists of positions
#  and finds the intersection between them in 1Mb windows

# it is useful to figure out how many son-specific heteroz positions
# in a window match paternal and maternal het sites

# the output file will be three columns:
# the number of sites in each 1Mb window in the first file
# the number of sites in each 1Mb window in the second file
# the overlap of these two numbers

# ./Intersect_two_lists_1Mb_windows_VENN.pl ChrJBJN_1_pat.out fivehetsons_pat_chrJBJN_1.out pat_sonz_chrJBJN_1.out

my $inputfile = $ARGV[0]; # assume no header
my $inputfile2 = $ARGV[1]; # assume no header
my $outputfile1 = $ARGV[2];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open DATAINPUT2, $inputfile2) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE1, ">$outputfile1"))  {
	print "I can\'t write to $outputfile1\n";
	exit;
}
print "Creating output file: $outputfile1\n";

print OUTFILE1 "Pos\tMatpat_het_sites\tSonDau_het_sites\tOverlap\n";

my @temp;
my @overlap;
my %file1;
my %file2;
my %overlap_hash;


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split /[\t\/]/,$line;
	$file1{$temp[1]}=1; # this position is in the first file
}

while ( my $line = <DATAINPUT2>) {
	chomp($line);
	@temp=split /[\t\/]/,$line;
	$file2{$temp[1]}=1; # this position is in the second file
	if(exists($file1{$temp[1]})){
		push(@overlap, $temp[1]);
	}
}	

print "hello @overlap\n";

my $countr = 0;
my $bin=1000000;

# load up a hash
foreach my $key (sort { $a <=> $b } keys %file1){
	if($key ne ""){
		if($key < $bin){
			$countr+=1;
		}
		else{
			$overlap_hash{$bin}[0]=$countr; # this is the number of sites in first file
			$overlap_hash{$bin}[1]=0; # just initialize the other counts
			$overlap_hash{$bin}[2]=0;
			$bin+=1000000;
			$countr=1;
		}
	}
}
# record last one
$overlap_hash{$bin}[0]=$countr; # this is the number of sites in first file
$overlap_hash{$bin}[1]=0; # just initialize the other counts
$overlap_hash{$bin}[2]=0;


$countr = 0;
$bin=1000000;

foreach my $key (sort { $a <=> $b } keys %file2){
	if($key ne ""){
		if($key < $bin){
			$countr+=1;
		}
		else{
			if(exists($overlap_hash{$bin}[0])){ # check if this has been initialized by the counts in file1
				$overlap_hash{$bin}[1]=$countr; # this is the number of sites in second file
			}
			else{
				$overlap_hash{$bin}[1]=$countr; # this is the number of sites in second file
				$overlap_hash{$bin}[0]=0; # just initialize the other counts
				$overlap_hash{$bin}[2]=0; # just initialize the other counts
			}	
			$bin+=1000000;
			$countr=1;
		}
	}
}

$countr = 0;
$bin=1000000;


foreach(@overlap){
	    print $bin," ",$countr," ",$_,"\n";
		if($_ < $bin){
			$countr+=1;
		}
		elsif($_ > $bin){
			$overlap_hash{$bin}[2]=$countr; # this is the number of sites that overlap
			while($bin < $_){
				print "hello ",$bin,"\n";
				$bin+=1000000;
			}
			$countr=1;
		}
}

foreach my $key (sort { $a <=> $b } keys %overlap_hash){
	if($key ne ""){
		print OUTFILE1 $key,"\t",$overlap_hash{$key}[0],"\t",$overlap_hash{$key}[1],"\t",$overlap_hash{$key}[2],"\n";
	}
}

close OUTFILE1;


```
