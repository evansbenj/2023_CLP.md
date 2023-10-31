# Subset the bam files
```#!/bin/sh                                                                                                           
#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=8gb
#SBATCH --output=samtools_subset.%J.out
#SBATCH --error=samtools_subset.%J.err
#SBATCH --account=rrg-ben

# run by passing the path to the sorted bam files like this
# sbatch ./2021_samtools_subset_bamfiles.sh directory region

module load StdEnv/2020 samtools/1.12
for file in ${1}*_rg.bam
do
    samtools view ${file} ${2} -b -o ${file}_${2}.bam
    samtools index ${file}_${2}.bam
done
```
