# Kmers!

path for 2024_clp:
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/2024_raw_data_larg_pyg/2024_pygm
```

I'm redoing this with new meryl and no cookiecutter:

# Kmer size 29 for pygmaeus, 21 for original clivii 

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

Make meryl db for forward and reverse reads (together) like this:
```
#!/bin/sh
#SBATCH --job-name=makemeryldb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=166:00:00
#SBATCH --mem=128gb
#SBATCH --output=makemeryldb.%J.out
#SBATCH --error=makemeryldb.%J.err
#SBATCH --account=rrg-ben

/home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl count ${1}*R[1,2].fq.gz threads=4 memory=128 k=29 output ${1}_meryldb.out threads=16
```


# Make intersection sum for female samples 
This requires a kmer to be present in all samples from a given sex.

```
#!/bin/sh
#SBATCH --job-name=meryl_intersect
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=2:00:00
#SBATCH --mem=128gb
#SBATCH --output=meryl_intersect.%J.out
#SBATCH --error=meryl_intersect.%J.err
#SBATCH --account=rrg-ben

# sbatch 2026_meryl_intersect.sh fq_meryl kmermeryl
/home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl intersect-sum \
output allfemz_intersect_sum.meryl \
fem_pyg_15_meryldb.out \
fem_pygm_ELI1682_meryldb.out \
fem_pygm_ELI2081_meryldb.out \
fem_pygm_ELI2372_meryldb.out \
fem_pygm_ELI3012_meryldb.out \
Z23338_female_meryldb.out \
Z23340_female_meryldb.out \
Z23341_female_meryldb.out \
Z23342_female_meryldb.out 

```

# Make a union-sum of all male samples  (this will be substracted from the intersect-sum for females)

This is required for males because we want to remove kmers that are present in all females but only some males. If we only subtract female-fixed kmers from male-fixed kmers, then some of the resulting kmers will be present in one sex and some individuals of the other sex.

```
#!/bin/sh
#SBATCH --job-name=meryl_unionsum
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=2:00:00
#SBATCH --mem=128gb
#SBATCH --output=meryl_intersect.%J.out
#SBATCH --error=meryl_intersect.%J.err
#SBATCH --account=rrg-ben

# sbatch 2026_meryl_intersect.sh fq_meryl kmermeryl
/home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl union-sum \
output allmalez_unionsum.meryl \
mal_pygm_ELI1681_meryldb.out \
mal_pygm_ELI2347_meryldb.out \
mal_pygm_ELI2370_meryldb.out \
mal_pygm_ELI2545_meryldb.out \
Z23337_male_meryldb.out \
Z23339_male_meryldb.out \
Z23349_male_meryldb.out \
Z23350_male_meryldb.out

```

# Subtract union-sum from males the intersect-sum of females

This will give kmers that are present in all females and no males

```
#!/bin/sh
#SBATCH --job-name=meryl_difference
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --mem=128gb
#SBATCH --output=meryl_intersect.%J.out
#SBATCH --error=meryl_intersect.%J.err
#SBATCH --account=rrg-ben

# sbatch 2026_meryl_intersect.sh fq_meryl kmermeryl
/home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl difference ${1} ${2} output in_${1}_but_not_${2}.meryl threads=16

```
This is the meryl db of the female-specifc kmers for pygm:
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/2024_raw_data_larg_pyg/2024_pygm/in_allfemz_intersect_sum.meryl_but_not_allmalez_unionsum.meryl.meryl
```

Can do the same for males, but RADseq data indicates that females are the heterogametic sex...


# print this output (not needed)
```
/home/ben/scratch/2023_clp_for_real/bin/meryl/build/bin/meryl print in_all_fems_Z23338_Z23340_Z23341_Z23342intersectsum.db_not_all_males_Z23337_Z23349_Z23339_Z23350_intersect_sum.db_differnece.db > fems_only_kmers.txt
```

# Extract paired reads that have sex-specific kmers


```
#!/bin/sh
#SBATCH --job-name=extractreadz
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --mem=128gb
#SBATCH --output=extractreadz.%J.out
#SBATCH --error=extractreadz.%J.err
#SBATCH --account=rrg-ben


/home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl-lookup -include \
  -sequence ${1}*R1.fq.gz ${1}*R2.fq.gz \
  -mers /home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/2024_raw_data_larg_pyg/2024_pygm/in_allfemz_intersect_sum.meryl_but_not_allmalez_unionsum.meryl.meryl \
  -output ${1}_fixedfemonly_R1.fq.gz ${1}_fixedfemonly_R2.fq.gz

```

# Assemble sex-specific reads using spades

```
#!/bin/sh
#SBATCH --job-name=spades
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=125gb
#SBATCH --output=spades.%J.out
#SBATCH --error=spades.%J.err
#SBATCH --account=rrg-ben


module load StdEnv/2023 spades/4.2.0
spades.py -1 ${1} -2 ${2} --isolate -o ${3}

```

# Blast XL CDS to contigs
```
 blastn -query ../../../../2021_XL_v10_refgenome/XL_CDS_only_nospaces.fasta -db /project/6019307/ben/2023_cliv_larg_pyg/2024_raw_data_larg_pyg/2024_pygm/fem_pygm/2024_pygm_femspecific_goodkmers_trinity_out_dir.Trinity.fasta_blastable -outfmt "6 qseqid sseqid pident length evalue bitscore qcovhsp qcovs" -qcov_hsp_perc 80 -out XL_CDS_to_pygmfemspecific.txt
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

# Blast XL transcripts to trinity assembly
in this directory:
```
/home/ben/projects/rrg-ben/ben/2024_cliv_allo_WGS/fq/2024_allo
```

```
module load StdEnv/2023 blast/2.2.26
module load blast+/2.14.1  StdEnv/2023 gcc/12.3
```
```
makeblastdb -dbtype nucl -in allo_only5mals_trinity_out_dir.Trinity.fasta -out allo_only5mals_trinity_out_dir.Trinity.fasta_blastable
```
```
blastn -query ../../../2021_XL_v10_refgenome/XENLA_10.1_Xenbase.transcripts.fa -db allo_only5mals_trinity_out_dir.Trinity.fasta_blastable -outfmt 6 -out XLtranscrips_to_allo_only5mals_trinity -evalue 1e-20 -task megablast
```
# focus on only the matches that are > 500 bp
```
cat XLtranscrips_to_allo_only5mals_trinity | awk '$4 < 500 { next } { print }'> XLtranscrips_to_allo_only5mals_trinity_match_atleast_500bp.txt
```
# Get a list of transcripts with a match
```
cat XLtranscrips_to_allo_only5mals_trinity_match_atleast_500bp.txt | cut -f1 | uniq > unique_XL_500bp_matches_to_allo5males.txt
```
# Now search for the full entry for these transcripts; or for a particular transcript
```
grep -f unique_XL_500bp_matches_to_allo5males.txt ../../../2021_XL_v10_refgenome/XENLA_10.1_Xenbase.transcripts.fa | grep 'notch'
```
