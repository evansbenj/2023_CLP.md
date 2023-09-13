# 2023_clivii_largeni_pygmaeus

This is a repo with Martin Knytl for the X. clivii, X. largeni, X. pygmaeus genomics project.

The goal of this project is to identify sex-linked regions in these species based on genomic data. 

# Searching for sex-specific heterozygous positions

For sex-linked regions, we expect divergence in the heterogametic sex and not the homogametic sex. We can identify these sites and map them.  

Here is a perl script that will look for these sites in a tab delimited file:
```perl #!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
use List::MoreUtils qw/ uniq /;

# Prepare input file
# module load tabix
# bgzip -c file.vcf > file.vcf.gz
# tabix -p vcf file.vcf.gz
#  Now use vcftools to make a tab delimited file:
# module load StdEnv/2020 vcftools/0.1.16
# zcat file.vcf.gz | vcf-to-tab > out.tab

# on computecanada:
# module load perl/5.30.2
# module load gcc/9.3.0
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited genotype file generated
#  by vcftools (vcf2tab) and searches for sites that are homozygous
#  in one sex for one SNP and at least partially heterozygous in the other sex

# to execute type Parse_tab.pl inputfile.tab 1111100110000111100011100110010100002200 interesting_sites.out proportion
# where 1111100110000111100011100110010100002200 refers to whether or not each individual in the ingroup 
# in the vcf file is (0) male, (1) female, and or (2) skipped

# proportion is the proportion of genotyped alleles in the heterogametic sex that are required to be
# different from the homogametic sex in order for the position to be reported.  This is a way to reduce reporting
# of low frequency polymorphisms (which are unlikely to be sex-linked but likely to have one sex all homozygous).
# the proportion parameter should be less than or equal to 0.5 

# if it is 0.5, this means all females are heterozygous and all males are homozygous (for positions with only 2 variants)

# we will also use this proportion to be a requirement for male-specific or female-specific SNPs, meaning at least
# this proportion of the individuals within each sex is required to have a genotype.

# het_sites.out sex_specific_sites.out diverged_sites.out are the output files that have the positions and chr of interesting sites

# example for 2023_clivii
# clivii
# perl Parse_tab.pl tabfile 1111100000222222222222222222 interesting_sites.out 0.5
# largeni
# perl Parse_tab.pl tabfile 2222222222111110000022222222 interesting_sites.out 0.5
# pyg
# perl Parse_tab.pl tabfile 2222222222222222222211110000 interesting_sites.out 0.5
    
# perl Parse_tab.pl clivii_unfiltered_removed_allchrs.vcf.tab 111111111111111111111111110000000000000000000 interesting_sites.out 0.35

# Example for XB_WGS
# perl Parse_tab.pl XB_WGS_not_filtered_allchrs.vcf.gz.tab 100110011101010000102222 interesting_sites.out 0.5

my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $outputfile1 = $ARGV[2];
my $proportion = $ARGV[3];

print "hello ",$proportion,"\n";

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE1, ">$outputfile1"))  {
	print "I can\'t write to $outputfile1\n";
	exit;
}
print "Creating output file: $outputfile1\n";
print OUTFILE1 "CHR\tPOS\tTYPE\tCATEGORY\tn_FEMs\tn_MALS\n";




my @sexes = split("",$ARGV[1]);

my @males=();
my @females=();
my @temp;
my @unique_male_nucleotides;
my @unique_female_nucleotides;
my $y;
my $x;
my $counter=0;
my $diverged=0;
my $number_of_male_individuals_genotyped=0;
my $number_of_female_individuals_genotyped=0;

for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if($sexes[$y] == 0){
		$number_of_male_individuals_genotyped +=1;
	}	
}	
for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if($sexes[$y] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}
print "This includes ",$number_of_female_individuals_genotyped," female(s) and  ", $number_of_male_individuals_genotyped," males\n";

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split /[\t\/]/,$line;
	if($temp[0] ne '#CHROM'){
		if($#temp ne (($#sexes+1)*2)+2){
			print "The number of individuals in the input line does not match the number of individuals genotyped ",
			$temp[0],"\t",$temp[1],"\t",$#temp," ",(($#sexes+1)*2)+2,"\n";
		}

		# parse the bases in all genotypes in each sex
		@males=();
		@females=();
		$counter=0;
		for ($y = 3 ; $y <= $#temp; $y=$y+2 ) {
			if(($temp[$y] ne ".")&&($temp[$y+1] ne ".")){
				if($sexes[$counter] == 0){
						push(@males, $temp[$y]);
						push(@males, $temp[$y+1]);
				}
				elsif($sexes[$counter] == 1){
					push(@females, $temp[$y]);
					push(@females, $temp[$y+1]);
				}	
			}
			$counter+=1;
		}	
		# OK I should have all the bases loaded for non-missing genotypes for each male and each female
		
		@unique_male_nucleotides = uniq @males;
		@unique_female_nucleotides = uniq @females;
		#print @females," ",@males,"\n";
		#print $#unique_male_nucleotides," ",$#unique_female_nucleotides,"\n";
		# looks fine
		if(($#unique_male_nucleotides != -1)&&($#unique_female_nucleotides != -1)){
			# we can compare homoz and het genotypes because both sexes have data
			if(($#unique_male_nucleotides == 0)&&($#unique_female_nucleotides > 0)){
				# all males are homoz but at least some females are hets or homoz for another SNP
				# check if the proportion of divergent positions in females is high enough
				$diverged=0;
				for ($x = 0 ; $x <= $#females ; $x++ ) {
					if($females[$x] ne $males[0]){
						$diverged+=1;
					}
				}
				if($diverged > $proportion*($#females+1)){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_SNP\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_SNP
					# Category of 1 means ZW or a female-specific SNP
				}
				# now check for the extreme case where all females are heterozygous and all males are homoz
				# this is rare because we expect some genotypes in females to be undercalled, even in sex-linked regions
				$diverged=0;
				for ($x = 0 ; $x <= $#females ; $x=$x+2 ) {
					if($females[$x] ne $females[$x+1]){
						$diverged+=1;
					}
				}
				if($diverged == ($#females+1)/2){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_heterozygosity\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_heterozygosity
					# 1 means all females are hets and all males are homoz
				}					
			}
			elsif(($#unique_male_nucleotides > 0)&&($#unique_female_nucleotides == 0)){
				# all females are homoz but at least some males are hets or homoz for another SNP
				# check if the proportion of divergent positions in males is high enough
				$diverged=0;
				for ($x = 0 ; $x <= $#males ; $x++ ) {
					if($males[$x] ne $females[0]){
						$diverged+=1;
					}
				}
				if($diverged > $proportion*($#males+1)){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_SNP\t-1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_SNP
					# -1 means XY or a male-specific SNP
				}	
				# now check for the extreme case where all females are heterozygous and all males are homoz
				# this is rare because we expect some genotypes in females to be undercalled, even in sex-linked regions
				$diverged=0;
				for ($x = 0 ; $x <= $#males ; $x=$x+2 ) {
					if($males[$x] ne $males[$x+1]){
						$diverged+=1;
					}
				}
				if($diverged == ($#males+1)/2){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_heterozygosity\t-1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_heterozygosity
					# -1 means all males are hets and all females are homoz
				}
			}
			elsif(($#unique_male_nucleotides == 0)&&
			($#unique_female_nucleotides == 0)&&
			($unique_male_nucleotides[0] ne $unique_female_nucleotides[0])){
				# males are homoz, females are homoz, but fixed for a different diverged nucleotide
				print OUTFILE1 $temp[0],"\t",$temp[1],"\tFixed_divergence\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
				# Fixed_divergence
				# 1 means diverged
			}
		}
		elsif(($#unique_male_nucleotides != -1)&&($#unique_female_nucleotides == -1)){
			# females have no data
			# could be male-specific
			if((($#males +1)/2) > $proportion*$number_of_male_individuals_genotyped){
				print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_nucleotides\t-1\t",($#females+1)/2,"\t",($#males+1)/2,"\n";
				# Sex_specific_nucleotides
				# -1 means male specific or male specific
			}	
		}
		elsif(($#unique_male_nucleotides == -1)&&($#unique_female_nucleotides != -1)){
			# males have no data
			# could be female-specific
			if((($#females +1)/2) > $proportion*$number_of_female_individuals_genotyped){
				print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_nucleotides\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
				# Sex_specific_nucleotides
				# 1 means fem specific or female specific
			}	
		}
	}
} # end while	
close OUTFILE1;
```


