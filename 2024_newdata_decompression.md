# Decompressing ora formatting
We collected new data from 6 clivii and 10 allofraseri (for Jade). These came in a compression format called 'ora'. I had to download and install the latest version of this software from Illumina (DRAGEN ORA Helper Suite v2.0 Installer) here: https://support.illumina.com/downloads.html?filters=web%253Asoftware%252Fdragen-decompression-software

Then I ran this script to decompress all the ora files in a directory:

```
#!/bin/sh
#SBATCH --job-name=ora
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=2gb
#SBATCH --output=ora.%J.out
#SBATCH --error=ora.%J.err
#SBATCH --account=rrg-ben



for file in ${1}/*ora ; do         # Use ./* ... NEVER bare *    
    if [ -e "$file" ] ; then   # Check whether file exists.
	/home/ben/projects/rrg-ben/ben/2024_cliv_allo_WGS/0FFSI2G/EVA28870.20241001/20240927_LH00600_0021_B22CYYGLT4/oraHelperSuite/orad --gz ${file}
  fi
done
```
