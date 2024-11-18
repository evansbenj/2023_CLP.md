# Merge bam files

For the 2024 allo and cliv WGS data, we received two pairs of fq files for each sample. We aligned each pair to the XL ref and then we beed to merge these bam files for each sample. Here is an sbatch script that does exactly that:
```
#!/bin/sh
#SBATCH --job-name=merge2
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=32:00:00
#SBATCH --mem=16gb
#SBATCH --output=merge2.%J.out
#SBATCH --error=merge2.%J.err
#SBATCH --account=rrg-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2024_samtools_merge2.sh merged_path_and_file in1_path_and_file in2_path_and_file 
module load StdEnv/2023  gcc/12.3 samtools/1.20

samtools merge ${1} ${2} ${3}
samtools index ${1}
```
