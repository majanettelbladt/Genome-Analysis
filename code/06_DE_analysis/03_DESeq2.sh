#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -J DEseq2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maja.nettelbladt.9823@student.uu.se
#SBATCH --output=%x.%j.out

module load bioinfo-tools
module load R/4.3.1
module load R_packages/4.3.1

Rscript /home/mane9823/Genome-Analysis/code/06_DE_analysis/03_DESeq2.R
