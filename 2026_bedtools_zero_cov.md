# Identify regions with zero coverage in a set of bam files
First get regions with zero coverage that are at least 50bp long in each bam file:
```
#!/bin/sh
#SBATCH --job-name=bedtools_zero_regionz
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --output=bedtools_zero_regionz.%J.out
#SBATCH --error=bedtools_zero_regionz.%J.err
#SBATCH --account=rrg-ben

# sbatch 2026_bedtools_regions_with_zerocoverage_bamfile.sh bamfile

module load bedtools
bedtools genomecov -ibam ${1} -bga | awk '$4==0 && ($3-$2)>50' > ${1}_zero_regionz.txt
```
Now find overlapping regions in multiple bam files- do this separately for males and females (need to check this):
```
bedtools multiinter -i *_zero_regionz.txt | awk -v N=$(ls *.zero.bed | wc -l) '$4==N'
```
