# Install meryl for counting and intersecting kmer dbs

I installed meryl (https://github.com/marbl/meryl) here on graham:
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/bin/meryl/build/bin
```
in case the command doesn't work:
```
export PATH=/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/bin/meryl/build/bin:$PATH
```

Make meryl db like this:
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


# sbatch 2020_meryl_make_kmerdb.sh fastqfile
# ../raw_data/NS.1462.004.IDT_i7_102---IDT_i5_102.XT7_WY_R1.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_102---IDT_i5_102.XT7_WY_R2.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_trim.R1_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_trim.R2_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_WZ_R1.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_WZ_R2.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_trim.R1_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_trim.R2_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_WW_R1.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_WW_R2.fastq.gz


/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/bin/meryl/build/bin/meryl count ${1} threads=4 memory=128 k=29 ou
tput ${1}_meryldb.out
```
