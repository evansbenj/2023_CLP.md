# Association with 16 genomes

```
#!/bin/sh
#SBATCH --job-name=angsd_allo
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=256gb
#SBATCH --output=angsd_allo.%J.out
#SBATCH --error=angsd_allo.%J.err
#SBATCH --account=rrg-ben


module load StdEnv/2023 angsd/0.940

angsd -yBin bin_sex.ybin -doAsso 1 -GL 1 -out out_additive_F1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minInd 4 -ba
m bam.filelist -P 5 -doCounts 1 -setMinDepthInd 2 -setMaxDepthInd 100 -Pvalue 1
```
