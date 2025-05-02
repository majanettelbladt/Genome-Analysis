#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2  
#SBATCH -t 10:00:00  
#SBATCH -J trimmomatic_Serum
#SBATCH --mail-type=ALL
#SBATCH --mail-user maja.nettelbladt.9823@student.uu.se 
#SBATCH --output=%x.%j.out   

module load bioinfo-tools
module load trimmomatic/0.39

# Sökvägar
RAW_DIR="/home/mane9823/Genome-Analysis/data/raw_data/transcriptomics_data/RNA-Seq_Serum"
OUTPUT_DIR="/home/mane9823/Genome-Analysis/data/trimmed_data/RNA-Seq_Serum"
ADAPTERS="/sw/apps/bioinfo/trimmomatic/0.39/rackham/adapters/TruSeq3-PE.fa"

# Skapa output-mapp om den inte finns
mkdir -p $OUTPUT_DIR

# Kör Trimmomatic på varje par
for R1 in $RAW_DIR/trim_paired_*_pass_1.fastq.gz; do
    # Matchande R2-fil
    R2="${R1/_pass_1.fastq.gz/_pass_2.fastq.gz}"
    
    # Extrahera ID
    ID=$(basename "$R1" | cut -d'_' -f3)

    # Definiera output-filer
    OUT_P1="$OUTPUT_DIR/${ID}_trimmed_R1.fastq.gz"
    OUT_UP1="$OUTPUT_DIR/${ID}_trimmed_R1_unpaired.fastq.gz"
    OUT_P2="$OUTPUT_DIR/${ID}_trimmed_R2.fastq.gz"
    OUT_UP2="$OUTPUT_DIR/${ID}_trimmed_R2_unpaired.fastq.gz"

    # Kör Trimmomatic
    trimmomatic PE -threads 2 -phred33 \
      "$R1" "$R2" \
      "$OUT_P1" "$OUT_UP1" \
      "$OUT_P2" "$OUT_UP2" \
      ILLUMINACLIP:$ADAPTERS:2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

echo "Trimmomatic är klar för alla RNA-Seq_Serum-filer."


