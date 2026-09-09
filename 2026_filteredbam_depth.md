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

# Filter bamz

```
samtools view -h -b -L Regions.bed alignments.bam > alignments_in_regions.bam
```
