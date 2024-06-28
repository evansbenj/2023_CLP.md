# Kmers!

path for 2024_clp:
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/2024_raw_data_larg_pyg/2024_pygm
```

I'm adapting a pipeline I developed earlier for XT_WZY:
```
https://github.com/evansbenj/XT_WW_WZ_WY/blob/main/kmers_unique_to_Y_and_W.md
```
# Kmer size 29 for new pygmaeus, 21 for original clivii

# Using trimmed fq files!

# Install meryl for counting and intersecting kmer dbs

I installed meryl (https://github.com/marbl/meryl) here on graham:
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/bin/meryl/build/bin
```
in case the command doesn't work:
```
export PATH=/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/bin/meryl/build/bin:$PATH
```

Make meryl db for forward and reverse reads (separately) like this:
```
#!/bin/sh
#SBATCH --job-name=meryl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=meryl.%J.out
#SBATCH --error=meryl.%J.err
#SBATCH --account=def-ben


/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/bin/meryl/build/bin/meryl count ${1} threads=4 memory=128 k=29 output ${1}_meryldb.out
```

# Make a union-sum for each  sample
This makes a new kmer db of kmers that are in the for or rev read, or in both. The count is the sum over the for and rev db.

```
#!/bin/sh
#SBATCH --job-name=meryl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=132gb
#SBATCH --output=meryl.%J.out
#SBATCH --error=meryl.%J.err
#SBATCH --account=def-ben

# sbatch 2020_meryl_union_kmer_dbs.sh db1 db2 out

/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/bin/meryl/build/bin/meryl union-sum ${1} ${2} threads=4 memory=128 k=29 output ${3}
```

# Make intersection sum for all samples within each sex
This requires a kmer to be present in all samples from a given sex. This should cut down on background stemming from sample-specific SNPs

```
#!/bin/sh
#SBATCH --job-name=meryl_intersect
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=128gb
#SBATCH --output=meryl_intersect.%J.out
#SBATCH --error=meryl_intersect.%J.err
#SBATCH --account=def-ben

# intersect-sum is makes the sum of counts that are in both
# this is not the union

/home/ben/scratch/2023_clp_for_real/bin/meryl/build/bin/meryl intersect-sum ${1} ${2} output ${1}_${2}_intersect_sum.db
```

# Make a union-sum of all samples within each sex (this will be substracted from the intersect-sum for each sex)

```
#!/bin/sh
#SBATCH --job-name=meryl_unionsum
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=128gb
#SBATCH --output=meryl_unionsum.%J.out
#SBATCH --error=meryl_unionsum.%J.err
#SBATCH --account=rrg-ben

# make symbolic links so that all dbs are in one directory
# and launch this in that directory
# sbatch 2020_meryl_union_kmer_dbs.sh out


/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/bin/meryl/build/bin/meryl union-sum *_R1R2 threads=4 memory=128 k=29 output ${1}
```

# Subtract union-sum from the intersect-sum of the other sex. For example, do this:
## intersect-sum_females - union-sum_males

This will give kmers that are present in all females and no males

This should be done in each way (F-M and M-F)
```
#!/bin/sh
#SBATCH --job-name=meryl_difference
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=128gb
#SBATCH --output=meryl_difference.%J.out
#SBATCH --error=meryl_difference.%J.err
#SBATCH --account=def-ben


/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/bin/meryl/build/bin/meryl difference ${1} ${2} output ${3}_differnece.db
```
This is the meryl db of the female-specifc kmers for pygm:
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/2024_raw_data_larg_pyg/2024_pygm/fem_pygm/in_all_fem_notinanyof_8malez_differnece.db
```
This directory has meryl dbs of female-specific and male-specific 21mer databases for clivii:
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/raw_data/cliv_pe_trim
```

# print this output
```
/home/ben/scratch/2023_clp_for_real/bin/meryl/build/bin/meryl print in_all_fems_Z23338_Z23340_Z23341_Z23342intersectsum.db_not_all_males_Z23337_Z23349_Z23339_Z23350_intersect_sum.db_differnece.db > fems_only_kmers.txt
```

# Extract paired reads that have sex-specific kmers
I'm going to use cookie cutter to extract reads with these kmers (https://github.com/NikoLichi/Cookiecutter)

According to this example, the format of the kmer file is almost the same as the output of meryl, except there are spaces instead of tabs between the kmer and the counts:
```
https://github.com/NikoLichi/Cookiecutter/blob/master/data/alpha.dat
```

```
#!/bin/sh
#SBATCH --job-name=cookie_extract
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=24gb
#SBATCH --output=cookie_extract.%J.out
#SBATCH --error=cookie_extract.%J.err
#SBATCH --account=def-ben

/home/ben/scratch/2023_clp_for_real/bin/Cookiecutter/bin/extract -1 ${1} -2 ${2} -o ${4} --fragments ${3}
```

# Combine se reads with R1 reads



Need to change the header of the SE reverse reads after concatenating:
```
sed -i 's/2:N:0:/1:N:0:/g'  male_specific_concat_all_left.fastq
```

# Add '/1' after each read
```
awk '{ if (NR%4==1) { print $1"_"$2"/1" } else { print } }' 2024_pygm_femspecific_goodkmers.fq > 2024_pygm_femspecific_goodkmers_rename.fq 
```

# Assemble sex-specific reads

```
#!/bin/sh
#SBATCH --job-name=trinity
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=166:00:00
#SBATCH --mem=128gb
#SBATCH --output=trinity.%J.out
#SBATCH --error=trinity.%J.err
#SBATCH --account=def-ben

module purge
module load StdEnv/2020 gcc/9.3.0 salmon/1.7.0 samtools/1.17 jellyfish/2.3.0  trinity/2.14.0 python scipy-stack

Trinity --seqType fq --max_memory 50G --left ${1}  --right ${2} --CPU 6
```

This either did not work or did not finish on computecanada. Instead it did work on info2020 in this directory:
```
/home/ben/2024_cliv_larg_pygm/raw_data/
```
using commands like this one:
```
/usr/local/trinity/Trinity --seqType fq --max_memory 120G --single 2024_pygm_femspecific_goodkmers.fq --CPU 12 --normalize_reads
```
as suggested here, I combined the paired and single end reads:
```
https://github.com/trinityrnaseq/trinityrnaseq/wiki/How-do-I-combine-reads%3F
```
Here is the assembly for X. pygmaeus with 29mers and 8 females, 8 males:
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/2024_raw_data_larg_pyg/2024_pygm/fem_pygm
```

# De novo assembled sex-specific contigs

Are here on graham:
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/raw_data/larg_pe_trim/combined_male_reads
```
# Map to XL genome using minimap2

```
module load StdEnv/2023 minimap2/2.26
minimap2 -x asm10 -a --secondary=no -t8 ../../../../2021_XL_v10_refgenome/XENLA_10.1_genome.fa larg_mal_only_trinity_denovo.fasta >alignments.sam 
```
Using this script:
```
#!/bin/sh
#SBATCH --job-name=minimap2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=16gb
#SBATCH --output=minimap2.%J.out
#SBATCH --error=minimap2.%J.err
#SBATCH --account=rrg-ben

# sbatch 2024_minimap2.sh ref.fasta query.fasta

module load StdEnv/2023 minimap2/2.26
minimap2 -x asm10 --secondary=no -t8 ${1} ${2} > ${2}_alignments.paf
```
which is here:
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/ben_scripts
```

Now pull out the reads with that map to a chromosome like this:
```
egrep 'Chr1L|Chr2L|Chr3L|Chr4L|Chr5L|Chr6L|Chr7L|Chr8L|Chr9_10L|Chr9_10S|Chr8S|Chr7S|Chr6S|Chr5S|Chr4S|Chr3S|Chr2S|Chr1S' larg_fem_only_trinity_denovo.fasta_alignments.sam >larg_fem_only_mappings.txt
```

# Plotting mappings using histograms and differences between histograms:
```
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2023_clivii_largeni_pygmaeus/kmers")
library(ggplot2)
library(plyr)
library(viridis)
library(dplyr)
options(scipen=999)


dat <-read.table("larg_fem_only_trinity_denovo.fasta_alignments.paf",header=F)
dat <-read.table("larg_mal_only_trinity_denovo.fasta_alignments.paf",header=F)
dat <-read.table("cliv_fem_only_trinity_denovo.fasta_alignments.paf",header=F)
dat <-read.table("cliv_mal_only_trinity_denovo.fasta_alignments.paf",header=F)
dat <-read.table("pygm_fem_only_trinity_denovo.fasta_alignments.paf",header=F)
dat <-read.table("pygm_mal_only_trinity_denovo.fasta_alignments.paf",header=F)

colnames(dat) <- c("query","query_len","query_start","query_end","strand","target","target_len",
                   "target_start","target_end","n_matches","n_bp","map_qual")
head(dat)


# Get rid of the scaffold data
my_df_chrsonly <- dat[(dat$target == "Chr1L")|(dat$target == "Chr2L")|(dat$target == "Chr3L")|
                          (dat$target == "Chr4L")|(dat$target == "Chr5L")|(dat$target == "Chr6L")|
                          (dat$target == "Chr7L")|(dat$target == "Chr8L")|(dat$target == "Chr9_10L")|
                          (dat$target == "Chr1S")|(dat$target == "Chr2S")|(dat$target == "Chr3S")|
                          (dat$target == "Chr4S")|(dat$target == "Chr5S")|(dat$target == "Chr6S")|
                          (dat$target == "Chr7S")|(dat$target == "Chr8S")|(dat$target == "Chr9_10S"),]

# save only mappings with unique matches
my_df_chrsonly_unique <- my_df_chrsonly %>% 
    distinct(query, .keep_all = T)
# now subset to include only mappings with map_qual>=60
my_df_chrsonly_unique_mq60 <- my_df_chrsonly_unique[(my_df_chrsonly_unique$map_qual >=60),]


my_df_chrsonly$target_start <- as.numeric(my_df_chrsonly$target_start)

png(filename = "larg_femonly_kmer_mapping_histo_.png",w=1200, h=1800,units = "px", bg="transparent")
    ggplot(my_df_chrsonly_unique_mq60, aes(x=target_start/1000000)) +
        #scale_fill_manual(values=c("red","blue"))+
        geom_histogram(binwidth = 0.1)+
        xlab("Position(Mb)") + ylab("Count") +
        facet_wrap(~ target, ncol=1) + 
        theme_classic() +
        theme(text = element_text(size = 20))
dev.off()


# make a histogram of the difference between the fem and mal histograms...
library(data.table)

fem_dat <-read.table("larg_fem_only_trinity_denovo.fasta_alignments.paf",header=F)
mal_dat <-read.table("larg_mal_only_trinity_denovo.fasta_alignments.paf",header=F)
fem_dat <-read.table("cliv_fem_only_trinity_denovo.fasta_alignments.paf",header=F)
mal_dat <-read.table("cliv_mal_only_trinity_denovo.fasta_alignments.paf",header=F)

colnames(fem_dat) <- c("query","query_len","query_start","query_end","strand","target","target_len",
                   "target_start","target_end","n_matches","n_bp","map_qual")
colnames(mal_dat) <- c("query","query_len","query_start","query_end","strand","target","target_len",
                       "target_start","target_end","n_matches","n_bp","map_qual")

fem_dat$sex <- "fem"
mal_dat$sex <- "mal"

# save only mappings with unique matches
fem_dat_unique <- fem_dat %>% 
    distinct(query, .keep_all = T)
# now subset to include only mappings with map_qual>=60
fem_dat_unique_mq60 <- fem_dat_unique[(fem_dat_unique$map_qual >=60),]

# save only mappings with unique matches
mal_dat_unique <- mal_dat %>% 
    distinct(query, .keep_all = T)
# now subset to include only mappings with map_qual>=60
mal_dat_unique_mq60 <- mal_dat_unique[(mal_dat_unique$map_qual >=60),]



all_dat <- rbind(fem_dat_unique_mq60,mal_dat_unique_mq60)

# Get rid of the scaffold data
all_dat_chrsonly <- all_dat[(all_dat$target == "Chr1L")|(all_dat$target == "Chr2L")|(all_dat$target == "Chr3L")|
                          (all_dat$target == "Chr4L")|(all_dat$target == "Chr5L")|(all_dat$target == "Chr6L")|
                          (all_dat$target == "Chr7L")|(all_dat$target == "Chr8L")|(all_dat$target == "Chr9_10L")|
                          (all_dat$target == "Chr1S")|(all_dat$target == "Chr2S")|(all_dat$target == "Chr3S")|
                          (all_dat$target == "Chr4S")|(all_dat$target == "Chr5S")|(all_dat$target == "Chr6S")|
                          (all_dat$target == "Chr7S")|(all_dat$target == "Chr8S")|(all_dat$target == "Chr9_10S"),]



all_dat_chrsonly$target_start <- as.numeric(all_dat_chrsonly$target_start)

# define an empty df
all_diffs <- data.frame(xmin = c("NA"),
                 xmax = c("NA"),
                 variable = c("NA"),
                 value = c("NA"),
                 Chr = c("NA")
                )


# plotting difference; modified from here:
# https://stackoverflow.com/questions/36049729/r-ggplot2-get-histogram-of-difference-between-two-groups
for (chromosome in unique(all_dat_chrsonly$target)) {
   # chromosome <- "Chr1L"
    all_dat_chrsonly <- all_dat[(all_dat$target == chromosome),]
    p <- ggplot(all_dat_chrsonly, aes(x=target_start/1000000, group=sex, color=sex, fill=sex)) + 
        geom_histogram(binwidth = 0.1, position="identity") +
        facet_wrap(~ target, ncol=1)
    
    p_data <- as.data.table(ggplot_build(p)$data[1])[,.(count,xmin,xmax,group)]
    p1_data <- p_data[group==1]
    p2_data <- p_data[group==2]

    newplot_data <- merge(p1_data, p2_data, by=c('xmin','xmax'), suffixes = c('.p1','.p2'),allow.cartesian=TRUE)
    newplot_data <- newplot_data[,diff:=count.p1 - count.p2]
    setnames(newplot_data, old=c('count.p1','count.p2'), new=c('k1','k2'),skip_absent=TRUE)

    df2 <- melt(newplot_data,id.vars =c('xmin','xmax'),measure.vars=c('k1','diff','k2'))
    only_diff <- df2[(df2$variable == "diff"),]
    only_diff$Chr <- chromosome
    all_diffs <- rbind(all_diffs,only_diff)
}

# get rid of first row
all_diffs = all_diffs[-1,]

all_diffs$xmin <- as.numeric(all_diffs$xmin)
all_diffs$xmax <- as.numeric(all_diffs$xmax)
all_diffs$value <- as.numeric(all_diffs$value)


difference_plot <- ggplot(all_diffs, aes(xmin=xmin,xmax=xmax,ymax=value,ymin=0)) + 
    geom_rect() +
    facet_wrap(~ Chr, ncol=1)

png(filename = "larg_kmer_mapping_histo_difference.png",w=1200, h=1800,units = "px", bg="transparent")
    difference_plot
dev.off()
```
