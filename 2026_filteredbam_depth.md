# Directory
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/2026_pygm_mapped_to_boum
```

# Filter bams to include CDS only then calculate depth in windows

There seems to be some background/noise due to repetitive regions. I'm going to try to calculate depth in windows using bams that only have reads mapped to CDS

# First blast XL CDS against the Xboum genome

Save hits that have at least 150 bp and at least 80% of the query length:
```
blastn \
  -query XL_CDS_only_nospaces.fasta \
  -db boumb.wbubble.fa_blastable \
  -outfmt "6 qseqid sseqid pident length qlen sstart send evalue bitscore" \
| awk '$4 >= 150 && $4/$5 >= 0.8' \
> XL_CDS_to_boum_gt150bp_gt80percentquery.txt
```
This will retain multiple hits for each XL query, but this is sensible since the Xboum genome is octoploid.

Make a bedfile out of this.
```
awk 'BEGIN{OFS="\t"}{
    print $2, ($6<$7?$6:$7), ($6<$7?$7:$6), $1
}' XL_CDS_to_boum_gt150bp_gt80percentquery.txt > XL_CDS_to_boum_gt150bp_gt80percentquery.bed
```

# Filter bamz

```
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
for file in ${1}*rg.bam
do
     samtools view -h -b -L /home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/2026_pygm_mapped_to_boum/XL_CDS_to_boum_gt150bp_gt80percentquery.bed ${file} > ${file}_CDS_only.bam
    samtools index ${file}_CDS_only.bam
done
```
