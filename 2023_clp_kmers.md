# Kmers!
I'm adapting a pipeline I developed earlier for XT_WZY:
```
https://github.com/evansbenj/XT_WW_WZ_WY/blob/main/kmers_unique_to_Y_and_W.md
```

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

# Make intersection sum for all samples within each sex
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

# Subtract these databases from each other
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


/home/ben/scratch/2023_clp_for_real/bin/meryl/build/bin/meryl difference ${1} ${2} output in_${1}_not_${2}_differnece.db
```

# print this output
```
/home/ben/scratch/2023_clp_for_real/bin/meryl/build/bin/meryl print in_all_fems_Z23338_Z23340_Z23341_Z23342intersectsum.db_not_all_males_Z23337_Z23349_Z23339_Z23350_intersect_sum.db_differnece.db > fems_only_kmers.txt
```

I'm going to use cookie cutter to extract reads with these kmers (https://github.com/NikoLichi/Cookiecutter)
