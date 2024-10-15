# 2023_clivii_largeni_pygmaeus

This is a repo with Martin Knytl for the X. clivii, X. largeni, X. pygmaeus genomics project.

The goal of this project is to identify sex-linked regions in these species based on genomic data. 

# Searching for sex-specific heterozygous positions

First make tab delimited files out of the genotype files (gvcf files)

```
module load StdEnv/2020 vcftools/0.1.16
zcat file.vcf.gz | vcf-to-tab > out.tab
```

For sex-linked regions, we expect divergence in the heterogametic sex and not the homogametic sex. We can identify these sites and map them.  

Here is a perl script that will look for these sites in a tab delimited file:
```perl #!/usr/bin/env perl
#!/usr/bin/env perl
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

# example for clivii
# perl Parse_tab.pl clivii_unfiltered_removed_allchrs.vcf.tab 111111111111111111111111110000000000000000000 interesting_sites.out 0.35
# include only Eritrea:
# 222222221111111111111112222222222200000000222

# exclude Eritrea:
# 111111112222222222222221110000000022222222000

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

# Plot Output with R

```R
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2024_cliv_allo_WGS/allo_Parsetab")
library(ggplot2)
library(plyr)
library(viridis)
options(scipen=999)
# Parsetab collects stats on these site patterns:
# [1] "Fixed_divergence"  # each sex is homoz for a different nucleotide       
# [2] "Sex_specific_heterozygosity" # this is the most stringent and should be used - all individuals
        # from one sex are heteroz and all individuals of the other sex are homoz
        # 1 means female-specific heterozygosity; -1 means male-specific heterozygosity
# [3] "Sex_specific_SNP"  # a heterozygous genotype is found only in one sex          
# [4] "Sex_specific_nucleotides" # I think this menas that a nucleotide is found in only one sex, 
                                 # and that it can be homozygous or heterozygous


# # I want to make an overlay histogram that shows "Sex_specific_heterozygosity" and "Sex_specific_SNP"
dat <-read.table("Chr7L_ZW_parsetab.txt",header=T)
# dat <-read.table("all_larg.txt",header=F)

colnames(dat) <- c("CHR","POS","TYPE","CATEGORY","n_FEMs","n_MALS")
head(dat)

# Make a column that summarizes the number of fems and male
dat$FandM <- paste(dat$n_FEMs,dat$n_MALS,sep="_")

dat$POS <- as.numeric(dat$POS)
dat$CATEGORY <- as.numeric(dat$CATEGORY)
dat$n_FEMs <- as.numeric(dat$n_FEMs)
dat$n_MALS <- as.numeric(dat$n_MALS)

unique(dat$FandM)

summary <- count(dat, "FandM");summary

# remove rows where FandM is "1_1" - this means only 2 individuals are genotyped
dat <- dat[dat$FandM != "1_1", ]
summary <- count(dat, "FandM");summary

# make separate dataframes for female-specific and male-specific heterozygosity
fem_specific <- dat[dat$CATEGORY == 1, ]
mal_specific <- dat[dat$CATEGORY == -1, ]
summary <- count(fem_specific, "FandM");summary
summary <- count(mal_specific, "FandM");summary

#all_SL_ZW <- fem_specific[!is.na(fem_specific$CHR) &
#                     ((fem_specific$TYPE == "Sex_specific_heterozygosity")|(fem_specific$TYPE == "Sex_specific_SNP")) &
#                     (fem_specific$n_FEMs >= 8) &  
#                     (fem_specific$n_MALS >= 8) &
#                     (fem_specific$CATEGORY == "1")   ,]
all_SL_ZW <- fem_specific[!is.na(fem_specific$CHR) &
                              (fem_specific$CHR == "Chr7L") &
                              ((fem_specific$TYPE == "Sex_specific_heterozygosity")|(fem_specific$TYPE == "Sex_specific_SNP")) &
                              (fem_specific$n_FEMs >= 5) &  
                              (fem_specific$n_MALS >= 5) &
                              (fem_specific$CATEGORY == "1")   ,]

all_SL_XY <- mal_specific[!is.na(mal_specific$CHR) &
                              (mal_specific$CHR == "Chr7L") &
                              ((mal_specific$TYPE == "Sex_specific_heterozygosity")|(mal_specific$TYPE == "Sex_specific_SNP")) &
                              (mal_specific$n_FEMs >= 5) &  
                              (mal_specific$n_MALS >= 5) &
                              (mal_specific$CATEGORY == "-1")   ,]

summary <- count(all_SL_ZW, "FandM");summary
summary <- count(all_SL_XY, "FandM");summary

# females
png(filename = "allo_ZW_55.png",w=1200, h=500,units = "px", bg="transparent")
p <- ggplot(all_SL_ZW, aes(x=POS/1000000, fill=TYPE)) +
    scale_fill_manual(values=c("red","blue"))+
    geom_histogram(binwidth = 0.3) + #, position = 'identity')+
    xlab("Position(Mb)") + ylab("Count") +
    xlim(0,230)+
    ylim(0,35)+
    facet_wrap(~CHR, ncol = 2, scales = "free_x") + 
    theme_classic() +
    theme(text = element_text(size = 20))
p
dev.off()

# males
png(filename = "allo_XY_55.png",w=1200, h=500,units = "px", bg="transparent")
p <- ggplot(all_SL_XY, aes(x=POS/1000000, fill=TYPE)) +
    scale_fill_manual(values=c("red","blue"))+
    geom_histogram(binwidth = 0.3) + #, position = 'identity')+
    xlab("Position(Mb)") + ylab("Count") +
    xlim(0,230)+
    ylim(0,35)+
    facet_wrap(~CHR, ncol = 2, scales = "free_x") + 
    theme_classic() +
    theme(text = element_text(size = 20))
p
dev.off()



# get the maximum bins on Chr7Lonly
oneChronly <- all_SL_XY[(all_SL_XY$CHR == "Chr7L"),]
Chr1onlyhist <- ggplot(oneChronly,aes(x=POS/1000000)) + geom_histogram( binwidth = 0.3)+
    stat_bin(aes(y=after_stat(count), label=after_stat(count)), binwidth = 0.3,geom="text", vjust=-.5) 
max_count <- max(ggplot_build(Chr1onlyhist)$data[[1]]$count);max_count
second_best <- sort(ggplot_build(Chr1onlyhist)$data[[1]]$count, decreasing = T)[2];second_best
third_best <- sort(ggplot_build(Chr1onlyhist)$data[[1]]$count, decreasing = T)[3];third_best
fourth_best <- sort(ggplot_build(Chr1onlyhist)$data[[1]]$count, decreasing = T)[4];fourth_best

over <- layer_data(Chr1onlyhist)
head(over)
# this is the best bin; we need to multiply it by 1,000,000 to get the coordinates
over[(over$count == max_count),]$x * 1000000
# It ends here
over[(over$count == max_count),]$x * 1000000 + (1000000*0.1)
# allo
# 17700000, 17800000  
# mki67
# upstream of foxi2

# this is the second best bin; we need to multiply it by 1,000,000 to get the coordinates
over[(over$count == second_best),]$x * 1000000
# It ends here
over[(over$count == second_best),]$x * 1000000 + (1000000*0.1)
# 42600000 42700000
# LOC108699156 uncharacterized LOC108699156, 
# eif3f.L eukaryotic translation initiation factor 3 subunit F, 
# LOC108699158 eukaryotic translation initiation factor 3 subunit F
# cfap126 cilia and flagella associated protein 126
# sdhc succinate dehydrogenase complex subunit C

# this is the third best bin; we need to multiply it by 1,000,000 to get the coordinates
over[(over$count == third_best),]$x * 1000000
# It ends here
over[(over$count == third_best),]$x * 1000000 + (1000000*0.1)
# 130800000, 130900000
# LOC108699447 erythroid membrane-associated protein
# LOC121397250 butyrophilin subfamily 3 member A2-like 


# this is the fourth best bin; we need to multiply it by 1,000,000 to get the coordinates
over[(over$count == fourth_best),]$x * 1000000
# It ends here
over[(over$count == fourth_best),]$x * 1000000 + (1000000*0.1)
# 129900000, 130000000
# LOC108699144 low affinity immunoglobulin gamma Fc region receptor II-like
# LOC108699145 Fc receptor-like protein 5 
# LOC121397354 low affinity immunoglobulin gamma Fc region receptor III-B-like
# LOC121397355 high affinity immunoglobulin gamma Fc receptor I-like
# LOC108699146 Fc receptor-like protein 4
# LOC121397239 high affinity immunoglobulin gamma Fc receptor I-like



View(oneChronly)









# less stringent
all_SL_ZW <- fem_specific[(fem_specific$TYPE == "Sex_specific_SNP") &
                              (fem_specific$n_FEMs >= 8) &  
                              (fem_specific$n_MALS >= 8) &
                              (fem_specific$CATEGORY == "1")   ,]
summary <- count(all_SL_ZW, "FandM");summary


all_SL_XY <- dat[(dat$TYPE == "Sex_specific_heterozygosity") &
                     (dat$n_FEMs >= 3) &  
                     (dat$n_MALS >= 3) &
                     (dat$CATEGORY == "-1")   ,]
summary <- count(all_SL_XY, "FandM");summary




all_SL_ZW_55 <- fem_specific[(fem_specific$TYPE == "Sex_specific_heterozygosity") &
                              (fem_specific$n_FEMs == 5) &  
                              (fem_specific$n_MALS == 5) &
                              (fem_specific$CATEGORY == "1")   ,]
summary <- count(all_SL_ZW_55, "FandM");summary
all_SL_XY_55 <- dat[(dat$TYPE == "Sex_specific_heterozygosity") &
                     (dat$n_FEMs == 5) &  
                     (dat$n_MALS == 5) &
                     (dat$CATEGORY == "-1")   ,]
summary <- count(all_SL_XY_55, "FandM");summary


# below not used

only_55_dat <- dat[dat$FandM == "5_5", ]

fem_specific_noscaf <- fem_specific[fem_specific$CHR != "Scaffolds", ]


# make the combined plot more stringent
fem_specific_noscaf_specific <- fem_specific_noscaf[(fem_specific_noscaf$FandM == "5_5"), ]
mal_specific_specific <- mal_specific[(mal_specific$FandM == "5_5"), ]
dim(fem_specific_noscaf_specific)
dim(mal_specific_specific)
# females
overlay_hist <- ggplot(fem_specific_noscaf,aes(x=POS/1000000, fill="FandM")) + 
    geom_histogram( binwidth = 0.1, fill = "black") +
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    #geom_histogram(data=subset(subset_locations,TYPE == 'W_linked'),fill = "red", alpha = 0.5, binwidth = 20000) +
    #geom_histogram(data=subset(subset_locations,TYPE == 'WY_linked'),fill = "blue", alpha = 0.5, binwidth = 20000) +
    #geom_histogram(data=subset(subset_locations,TYPE == 'Z_linked'),fill = "black", alpha = 0.5, binwidth = 20000) +
    #  scale_x_continuous(breaks = round(seq(9577608, 9593589, by = 50000))) +
    #  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    #geom_rect(aes(xmin=139294668,xmax=139350317,ymin=-2,ymax=Inf),color="blue",alpha=0) + # dmrt1L
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr1L"), aes(xintercept = 139294668/1000000),color="blue",alpha=0.3) + # dmrt1L
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr1L"), aes(xintercept = 111921700/1000000),color="blue",alpha=0.3) + # amh.L 111921700..111936408
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr1S"), aes(xintercept = 118799066/1000000),color="blue",alpha=0.3) + # dmrt1S 118799066..118842776
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr1S"), aes(xintercept = 93431716/1000000),color="blue",alpha=0.3) + # amh.S 93431716..93458671 
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr4S"), aes(xintercept = 344155/1000000),color="blue",alpha=0.3) + # sf1.S 344155..361085
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr5L"), aes(xintercept = 137927930/1000000),color="blue",alpha=0.3) + # foxL2.L 137927930..137931728
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr5L"), aes(xintercept = 49625129/1000000),color="blue",alpha=0.3) + # estrogen receoptor 1L 49625129..49743175
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr5L"), aes(xintercept = 137927930/1000000),color="blue",alpha=0.3) + # foxL2.L 137927930..137931728
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr5S"), aes(xintercept = 40136271/1000000),color="blue",alpha=0.3) + # estrogen receoptor 1S 40136271..40247358
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr5S"), aes(xintercept = 115779807/1000000),color="blue",alpha=0.3) + # foxl2.S 115779807..115781431
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr7L"), aes(xintercept = 98176992/1000000),color="blue",alpha=0.3) + # wnt4.l 98176992..98220522
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr7S"), aes(xintercept = 82417208/1000000),color="blue",alpha=0.3) + # wnt4.S 82417208..82467306
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr8L"), aes(xintercept = 30541025/1000000),color="blue",alpha=0.3) + #androgen receptor 30541025..30689464
    #geom_vline(data = subset(fem_specific_noscaf_specific, CHR == "Chr8L"), aes(xintercept = 92341707/1000000),color="blue",alpha=0.3) + #estrogen receptor 2L 92341707..92374149
    xlim(0,240)+
    xlab("Position(MB)") + ylab("Count") +
    #geom_text(stat= "count", aes(label=after_stat(count)), vjust=-1, size=3) +
    #stat_bin(aes(y=after_stat(count), label=after_stat(count)), binwidth = 0.05,geom="text", vjust=-.5) +
    theme_classic();overlay_hist

# males
moverlay_hist <- ggplot(mal_specific_specific,aes(x=POS/1000000, fill="FandM")) + 
    geom_histogram( binwidth = 0.1, fill = "black") +
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    #geom_histogram(data=subset(subset_locations,TYPE == 'W_linked'),fill = "red", alpha = 0.5, binwidth = 20000) +
    #geom_histogram(data=subset(subset_locations,TYPE == 'WY_linked'),fill = "blue", alpha = 0.5, binwidth = 20000) +
    #geom_histogram(data=subset(subset_locations,TYPE == 'Z_linked'),fill = "black", alpha = 0.5, binwidth = 20000) +
    #  scale_x_continuous(breaks = round(seq(9577608, 9593589, by = 50000))) +
    #  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    geom_rect(aes(xmin=139294668,xmax=139350317,ymin=-2,ymax=Inf),color="blue",alpha=0) + # dmrt1L
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr1L"), aes(xintercept = 139294668/1000000),color="blue",alpha=0.3) + # dmrt1L
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr1L"), aes(xintercept = 111921700/1000000),color="blue",alpha=0.3) + # amh.L 111921700..111936408
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr1S"), aes(xintercept = 118799066/1000000),color="blue",alpha=0.3) + # dmrt1S 118799066..118842776
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr1S"), aes(xintercept = 93431716/1000000),color="blue",alpha=0.3) + # amh.S 93431716..93458671 
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr4S"), aes(xintercept = 344155/1000000),color="blue",alpha=0.3) + # sf1.S 344155..361085
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr5L"), aes(xintercept = 137927930/1000000),color="blue",alpha=0.3) + # foxL2.L 137927930..137931728
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr5L"), aes(xintercept = 49625129/1000000),color="blue",alpha=0.3) + # estrogen receoptor 1L 49625129..49743175
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr5L"), aes(xintercept = 137927930/1000000),color="blue",alpha=0.3) + # foxL2.L 137927930..137931728
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr5S"), aes(xintercept = 40136271/1000000),color="blue",alpha=0.3) + # estrogen receoptor 1S 40136271..40247358
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr5S"), aes(xintercept = 115779807/1000000),color="blue",alpha=0.3) + # foxl2.S 115779807..115781431
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr7L"), aes(xintercept = 98176992/1000000),color="blue",alpha=0.3) + # wnt4.l 98176992..98220522
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr7S"), aes(xintercept = 82417208/1000000),color="blue",alpha=0.3) + # wnt4.S 82417208..82467306
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr8L"), aes(xintercept = 30541025/1000000),color="blue",alpha=0.3) + #androgen receptor 30541025..30689464
    geom_vline(data = subset(mal_specific_specific, CHR == "Chr8L"), aes(xintercept = 92341707/1000000),color="blue",alpha=0.3) + #estrogen receptor 2L 92341707..92374149
    xlim(0,240)+
    xlab("Position(MB)") + ylab("Count") +
    #geom_text(stat= "count", aes(label=after_stat(count)), vjust=-1, size=3) +
    #stat_bin(aes(y=after_stat(count), label=after_stat(count)), binwidth = 0.05,geom="text", vjust=-.5) +
    theme_classic();overlay_hist

pdf("./larg_fem_noscaf_fixed_divergence_hist_5_5.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    overlay_hist
dev.off()

pdf("./larg_male_noscaf_fixed_divergence_hist_5_5.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    moverlay_hist
dev.off()

# cliv peaks: Chr3L
# larg peaks: Chr3S, Chr6L
# get the maximum bins on Chr1Lonly
oneChronly <- fem_specific_noscaf_specific[(fem_specific_noscaf_specific$CHR == "Chr8L"),]
Chr1onlyhist <- ggplot(oneChronly,aes(x=POS/1000000)) + geom_histogram( binwidth = 0.05)+
    stat_bin(aes(y=after_stat(count), label=after_stat(count)), binwidth = 0.05,geom="text", vjust=-.5) 
max_count <- max(ggplot_build(Chr1onlyhist)$data[[1]]$count);max_count
second_best <- sort(ggplot_build(Chr1onlyhist)$data[[1]]$count, decreasing = T)[2];second_best
over <- layer_data(Chr1onlyhist)
head(over)
# this is the bin; we need to multiply it by 1,000,000 to get the coordinates
over[(over$count == max_count),]$x * 1000000
# It ends here
over[(over$count == max_count),]$x * 1000000 + (1000000*0.05)
# this is the second best bin; we need to multiply it by 1,000,000 to get the coordinates
over[(over$count == second_best),]$x * 1000000
# It ends here
over[(over$count == second_best),]$x * 1000000 + (1000000*0.05)
View(oneChronly)
write.table(oneChronly, file = "pygm_Chr8L_SL_4_4.csv", sep = ",")


fplot_noscaf<-ggplot(fem_specific_noscaf_specific, aes(x = POS/1000000)) +
    geom_density(linewidth=0.15,adjust = 0.01)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

mplot_noscaf<-ggplot(mal_specific, aes(x = POS/1000000)) +
    geom_density(linewidth=0.15,adjust = 0.01)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

pdf("./larg_female_noscaf_specific_heterozygosity_5.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    fplot_noscaf
dev.off()
pdf("./larg_male_noscaf_specific_heterozygosity_5.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    mplot_noscaf
dev.off()


pdf("./cliv_female_specific_heterozygosity.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    fplot
dev.off()

pdf("./cliv_male_specific_heterozygosity.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    mplot
dev.off()


pdf("./cliv_female_specific_allcombined.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    fplot_allcombined
dev.off()



# make separate 55 dataframes for female-specific and male-specific heterozygosity
fem_55_specific <- dat[only_55_dat$CATEGORY == 1, ]
mal_55_specific <- dat[only_55_dat$CATEGORY == -1, ]
summary <- count(fem_55_specific, "FandM");summary
summary <- count(mal_55_specific, "FandM");summary


f55plot<-ggplot(fem_55_specific, aes(x = POS/1000000, fill=FandM)) +
    geom_density(linewidth=0.15,adjust = 0.05)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

pdf("./larg_55_female_specific.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    f55plot
dev.off()

# need to add endpoints to prevent densities from tanking
dat[nrow(dat) + 1,] <- c("Chr1L","232529967","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr1S","196169796","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr2L","184566229","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr2S","167897111","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr3L","145564449","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr3S","127416162","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr4L","156120765","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr4S","131359388","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr5L","174499024","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr5S","139053354","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr6L","157843502","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr6S","137668413","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr7L","136892544","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr7S","105895006","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr8L","123836259","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr8S","105436522","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr9_10L","135078614","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Chr9_10S","110702964","Sex_specific_heterozygosity","1","5","5","5_5")
dat[nrow(dat) + 1,] <- c("Scaffolds","42147407","Sex_specific_heterozygosity","1","5","5","5_5")
#View(all_SL)






View(all_SL)
# write.csv(all_SL,"all_SL_mappedtoXL.csv", row.names = FALSE)
other_types_of_SL <- dat[((dat$TYPE == "Sex_specific_SNP")|(dat$TYPE == "Fixed_divergence")|(dat$TYPE == "Sex_specific_nucleotides")) &
                  (as.numeric(dat$n_FEMs) >= 10) &  
                  (as.numeric(dat$n_MALS) >= 10),]
View(other_types_of_SL)
high_het_SL <- all_SL[(as.numeric(all_SL$POS) >=140372057)&(as.numeric(all_SL$POS) <=140433397),]
View(high_het_SL) # lots!
the_rest_of_dmrt1L <- all_SL[(as.numeric(all_SL$POS) >140449273)&(as.numeric(all_SL$POS) <=140462368),]
View(the_rest_of_dmrt1L) # nothing!

high_het_other_SL <- other_types_of_SL[(as.numeric(other_types_of_SL$POS) >=140372057)&(as.numeric(other_types_of_SL$POS) <=140433397),]
View(high_het_other_SL) # nothing!
the_other_rest_of_dmrt1L <- other_types_of_SL[(as.numeric(other_types_of_SL$POS) >140449273)&(as.numeric(other_types_of_SL$POS) <=140462368),]
View(the_other_rest_of_dmrt1L) # nothing!

#    geom_rect(aes(xmin=140372057,xmax=140433397,ymin=-2,ymax=Inf),color="blue",alpha=0) + # high heteroz
#    geom_rect(aes(xmin=140462145,xmax=140462368,ymin=-2,ymax=Inf),color="red",alpha=0) + # XB dmrt1L ex2
#    geom_rect(aes(xmin=140449273,xmax=140449553,ymin=-2,ymax=Inf),color="red",alpha=0) + # XB dmrt1L ex3
#    geom_rect(aes(xmin=140395623,xmax=140395753,ymin=-2,ymax=Inf),color="red",alpha=0) + # XB dmrt1L ex4
#    ex5 missing
#    geom_rect(aes(xmin=140449273,xmax=140449319,ymin=-2,ymax=Inf),color="red",alpha=0) + # XB dmrt1L ex6    
    
all_SL_ZW <- dat[(dat$TYPE == "Sex_specific_heterozygosity") &
              (dat$n_FEMs == "5") &  
              (dat$n_MALS == "5") &
              (dat$CATEGORY == "1")   ,]

all_SL_XY <- dat[(dat$TYPE == "Sex_specific_heterozygosity") &
                     (dat$n_FEMs == "5") &  
                     (dat$n_MALS == "5") &
                     (dat$CATEGORY == "-1")   ,]

dim(all_SL_ZW)
dim(all_SL_XY)
View(all_SL_ZW)
View(all_SL_XY)


all_SL_ZW_Chr1L <- dat[(dat$TYPE == "Sex_specific_heterozygosity") &
                     (dat$n_FEMs == "5") &  
                     (dat$n_MALS == "5") &
                     (dat$CATEGORY == "1") &
                     (dat$CHR == "Chr1L"),]

ggplot(all_SL_ZW_Chr1L, aes(x = POS)) +
        geom_density(adjust = 0.05)

library(plyr)
all_SL_ordered_ZW <- all_SL_ZW[order(all_SL_ZW$CHR, all_SL_ZW$POS),]
summary <- count(all_SL_ordered_ZW, "CHR")
chr_sizes <- c(232529967,196169796,184566229,167897111,145564449,127416162,156120765,
               131359388,174499024,139053354,157843502,137668413,136892544,105895006,123836259,
               105436522,135078614,110702964,42147407)
summary <- cbind(summary,chr_sizes)
summary$proportions <- summary$freq/summary$chr_sizes
summary[summary$proportions == max(summary$proportions),]




#View(temp)
temp <-c("red","blue")

p<-ggplot(all_SL_ordered_ZW, aes(x = as.numeric(POS)/1000000)) +
    geom_density(adjust = 0.05)+
    #geom_vline(xintercept=49000000)+ 
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlab("Position") + ylab("Density") +
    ylim(0,1) +
    theme_bw()


pdf("./larg_female_specific_heterozygosity.pdf",w=10, h=30.0, version="1.4", bg="transparent")
p
dev.off()


# filter data to retain only sites with a reasonable amount of
# genotypes in both sexes

subset_data <- dat[((as.numeric(dat$n_FEMs) + as.numeric(dat$n_MALS)) > 6)&(dat$CATEGORY==1),]
dim(subset_data)

fplot<-ggplot(dat, aes(x = as.numeric(POS)/1000000, fill=FandM)) +
    geom_density(linewidth=0.15, adjust = 0.001)+
    facet_wrap(~CHR, ncol = 2, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

pdf("./larg_gt8_femspecific.pdf",w=18, h=18.0, version="1.4", bg="transparent")
    fplot
dev.off()

p<-ggplot(dat, aes(x=POS, y=CATEGORY*(n_FEMs+n_MALS))) + 
    # add points
    geom_point(size=2, alpha = 0.7 ) +
    # add loess line
    # geom_smooth() +
    # color the stuff the way I want
    facet_wrap(~CHR, ncol = 1) +
    # get rid of gray background
    theme_bw()
pdf("./out_0.5.pdf",w=18, h=18.0, version="1.4", bg="transparent")
    p
dev.off()





fplot<-ggplot(all_SL_ZW, aes(x = POS/1000000, fill=FandM)) +
    geom_density(adjust=0.01)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()
pdf("./larg_femspecific.pdf",w=6, h=28.0, version="1.4", bg="transparent")
fplot
dev.off()

mplot<-ggplot(all_SL_XY, aes(x = POS/1000000, fill=FandM)) +
    geom_density(adjust=0.01)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

pdf("./larg_malspecific.pdf",w=6, h=28.0, version="1.4", bg="transparent")
mplot
dev.off()


fplot<-ggplot(all_SL_ZW_55, aes(x = POS/1000000, fill=FandM)) +
    geom_density(adjust=0.01)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()
pdf("./larg_femspecific_55.pdf",w=6, h=28.0, version="1.4", bg="transparent")
fplot
dev.off()

mplot<-ggplot(all_SL_XY_55, aes(x = POS/1000000, fill=FandM)) +
    geom_density(adjust=0.01)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

pdf("./larg_malspecific_55.pdf",w=6, h=28.0, version="1.4", bg="transparent")
mplot
dev.off()

all_SL_ZW_54 <- fem_specific[(fem_specific$TYPE == "Sex_specific_heterozygosity") &
                                 (fem_specific$n_FEMs == 5) &  
                                 (fem_specific$n_MALS == 4) &
                                 (fem_specific$CATEGORY == "1")   ,]
summary <- count(all_SL_ZW_54, "FandM");summary
all_SL_XY_54 <- dat[(dat$TYPE == "Sex_specific_heterozygosity") &
                        (dat$n_FEMs == 5) &  
                        (dat$n_MALS == 4) &
                        (dat$CATEGORY == "-1")   ,]
summary <- count(all_SL_XY_54, "FandM");summary

fplot<-ggplot(all_SL_ZW_54, aes(x = POS/1000000, fill=FandM)) +
    geom_density(adjust=0.01)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()
pdf("./larg_femspecific_54.pdf",w=6, h=28.0, version="1.4", bg="transparent")
fplot
dev.off()

mplot<-ggplot(all_SL_XY_54, aes(x = POS/1000000, fill=FandM)) +
    geom_density(adjust=0.01)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

pdf("./larg_malspecific_54.pdf",w=6, h=28.0, version="1.4", bg="transparent")
mplot
dev.off()


all_SL_ZW_45 <- fem_specific[(fem_specific$TYPE == "Sex_specific_heterozygosity") &
                                 (fem_specific$n_FEMs == 4) &  
                                 (fem_specific$n_MALS == 5) &
                                 (fem_specific$CATEGORY == "1")   ,]
summary <- count(all_SL_ZW_45, "FandM");summary
all_SL_XY_45 <- dat[(dat$TYPE == "Sex_specific_heterozygosity") &
                        (dat$n_FEMs == 4) &  
                        (dat$n_MALS == 5) &
                        (dat$CATEGORY == "-1")   ,]
summary <- count(all_SL_XY_45, "FandM");summary

fplot<-ggplot(all_SL_ZW_45, aes(x = POS/1000000, fill=FandM)) +
    geom_density(adjust=0.01)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()
pdf("./larg_femspecific_45.pdf",w=6, h=28.0, version="1.4", bg="transparent")
fplot
dev.off()

mplot<-ggplot(all_SL_XY_45, aes(x = POS/1000000, fill=FandM)) +
    geom_density(adjust=0.01)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

pdf("./larg_malspecific_45.pdf",w=6, h=28.0, version="1.4", bg="transparent")
mplot
dev.off()


```


