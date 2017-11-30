#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH -p main
#SBATCH --export=ALL
#SBATCH -D /home/ag1349/SNPTB/code 
#SBATCH -o SNPTB%j.out
#SBATCH -e SNPTB%j.err

module load python/2.7.11
module load bowtie2/2.2.6
module load samtools/1.2
module load bcftools/1.2
module load java/1.8.0_66

srun python organize_data.py <full path to data directory>
srun python qcmap.py <full path to data directory>
srun python get_avgdepth_genomecov.py <full path to data directory>
srun python snp_annotation.py <full path to data directory>
