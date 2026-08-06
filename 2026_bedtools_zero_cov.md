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
need to sort these regions too:
```
sort -k1,1 -k2,2n pygm_Z23342_female_to_iridian_sorted_rg.bam_zero_regionz_all.txt > pygm_Z23342_female_to_iridian_sorted_rg.bam_zero_regionz_all_sorted.txt
```

Now find overlapping regions in multiple bam files- do this separately for males and females (need to check this):
```
bedtools multiinter -i *_zero_regionz_all.txt -g ../../../../../cliv_ref/cliv_male_GCA_046118455.1_ASM4611845v1_genomic.fa.bed -empty | grep '1,2,3,4,5,6,7,8,9' > 9femalez_regions_with_no_coverage.txt
```
```
bedtools multiinter -i *_zero_regionz_all.txt -g ../../../Z23349_male_spades_assembly/Z23349_male_spades_assembly/scaffolds.bed -empty | grep '1,2,3,4,5,6,7,8,9' > 9femalez_regions_with_no_coverage.txt
```

Extract regions that are bigger than 1000bp
```
awk '{
diff = $1 - $2
if (diff < 0) diff = -diff
if (diff > 1000) print}' 8malez_regions_with_no_coverage.txt > 8malez_regions_with_no_coverage_gt_1000.txt

```
