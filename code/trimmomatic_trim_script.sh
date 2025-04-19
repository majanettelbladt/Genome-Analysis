#!/bin/bash
#SBATCH -A uppmax2025-3-3         # Projekt-ID
#SBATCH -M snowy                  # Kluster Snowy
#SBATCH -p core                   # Typ av noder
#SBATCH -n 2                      # 2 cores
#SBATCH -t 02:00:00               # Max tid (2 timmar)
#SBATCH -J trimmomatic            # Jobbnamn
#SBATCH --mail-user=maja.nettelbladt.9823@student.uu.se  # E-post vid jobbets slut
#SBATCH --mail-type=ALL     # Skicka mail vid slut eller fail
#SBATCH --output=%x.%j.out        # Outputfilens namn

# Ladda moduler
module load bioinfo-tools
module load trimmomatic

# 7. Ladda trimmomatic för data preprocessing
module load bioinfo-tools trimmomatic
mkdir -p /home/mane9823/Genome-Analysis/data/trimmed_data

# 8. Trimma Illumina reads
trimmomatic PE -threads 2 \  /home/mane9823/Genome-Analysis/data/raw_data/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz \  /home/mane9823/Genome-Analysis/data/raw_data/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz \  /home/mane9823/Genome-Analysis/data/trimmed_data/E745-1.L500_SZAXPI015146-56_1_paired.fq.gz \  /home/mane9823/Genome-Analysis/data/trimmed_data/E745-1.L500_SZAXPI015146-56_1_unpaired.fq.gz \  /home/mane9823/Genome-Analysis/data/trimmed_data/E745-1.L500_SZAXPI015146-56_2_paired.fq.gz \  /home/mane9823/Genome-Analysis/data/trimmed_data/E745-1.L500_SZAXPI015146-56_2_unpaired.fq.gz \
  ILLUMINACLIP:/sw/bioinfo/trimmomatic/0.39/snowy/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Trimming är klart!"
