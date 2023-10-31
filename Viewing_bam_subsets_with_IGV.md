# Subset the bam files
```#!/bin/sh                                                                                                           
#SBATCH --nodes=1                                                                                                   
#SBATCH --ntasks-per-node=1                                                                                         
#SBATCH --time=4:00:00                                                                                              
#SBATCH --mem=8gb                                                                                                   
#SBATCH --output=samtools_subset.%J.out                                                                             
#SBATCH --error=samtools_subset.%J.err                                                                              
#SBATCH --account=rrg-ben                                                                                           

# run by passing the path to the sorted bam files like this                                                         
# sbatch ./2021_samtools_subset_bamfiles.sh bamfile_prefix region                                                   

module load StdEnv/2020 samtools/1.12
samtools view ${1}.bam ${2} -b -o ${1}_${2}.bam
samtools index ${1}_${2}.bam
```
