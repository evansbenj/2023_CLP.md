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
or for independent chromosomes:
```
/home/ben/projects/rrg-ben/ben/2023_cliv_larg_pyg/ben_scripts/2023_angsd_doassoc.sh
```
```
#!/bin/sh
#SBATCH --job-name=angsd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=312gb
#SBATCH --output=angsd.%J.out
#SBATCH --error=angsd.%J.err
#SBATCH --account=rrg-ben

# http://popgen.dk/angsd/index.php/Input
# this describes how to specify regions
# here we need to add the chromosome after the sbatch command

module load StdEnv/2023 angsd/0.940
angsd -yBin bin_sex.ybin -doAsso 1 -GL 1 -out alloWGS_out_${1}_additive_F1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minInd 4 -bam bam.filelist -r ${1}: -P 5 -doCounts 1 -setMinDepthInd 2 -setMaxDepthInd 100 -Pvalue 1
```
