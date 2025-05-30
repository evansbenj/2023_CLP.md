# Genotyping cliv to cliv ref
There is a surprisingly high quality genome assembly for Xenopus clivii now. I am using this to genotype the clivii samples.

We have bam files now and need to make vcf files. This takes a while and did not complete. So we need to patch together separate runs.

# Compress and index:
```
#!/bin/sh
#SBATCH --job-name=bgzip
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=6:00:00
#SBATCH --mem=2gb
#SBATCH --output=bgzip.%J.out
#SBATCH --error=bgzip.%J.err
#SBATCH --account=def-ben

module load tabix

for file in ${1}*vcf
do
    bgzip -c ${file} > ${file}.gz
    tabix -p vcf ${file}.gz 
done
```

# Now extract chromosomes that finished:
Check with `tail -n1 vcf`. Then extract the finished ones like this:

```
module load StdEnv/2023  gcc/12.3 bcftools/1.19
bcftools view fem_CAS260390_all_sorted.bam_rg.bam.g.vcf.gz --regions JBJNIE010000001.1,JBJNIE010000002.1,JBJNIE010000003.1,JBJNIE010000004.1,JBJNIE010000005.1,JBJNIE010000006.1,JBJNIE010000007.1,JBJNIE010000008.1,JBJNIE010000009.1,JBJNIE010000010.1,JBJNIE010000011.1,JBJNIE010000012.1,JBJNIE010000013.1 > fem_CAS260390_all_sorted.bam_rg.bam_upto_13.g.vcf.gz
```
from this list:
```
JBJNIE010000001.1,JBJNIE010000002.1,JBJNIE010000003.1,JBJNIE010000004.1,JBJNIE010000005.1,JBJNIE010000006.1,JBJNIE010000007.1,JBJNIE010000008.1,JBJNIE010000009.1,JBJNIE010000010.1,JBJNIE010000011.1,JBJNIE010000012.1,JBJNIE010000013.1,JBJNIE010000014.1,JBJNIE010000015.1,JBJNIE010000016.1,JBJNIE010000017.1,JBJNIE010000018.1,JBJNIE010000019.1
```
Now genotype the other large chromsomes
```
sbatch 2021_HaplotypeCaller_onebam_onechr.sh ../cliv_ref/GCA_046118455.1_ASM4611845v1_genomic.fa ../2024_cliv/rg_bamz/fem_CAS260390_all_sorted.bam_rg.bam JBJNIE010000014.1,JBJNIE010000015.1,JBJNIE010000016.1,JBJNIE010000017.1,JBJNIE010000018.1,JBJNIE010000019.1
```
