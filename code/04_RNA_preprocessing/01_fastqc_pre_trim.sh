#!/bin/bash
#SBATCH -A uppmax2025-3-3         # Projekt-ID
#SBATCH -M snowy                  # Kluster Snowy
#SBATCH -p core                   # Typ av noder
#SBATCH -n 2                      # 2 cores
#SBATCH -t 02:00:00               # Max tid (2 timmar)
#SBATCH -J illumina_qc            # Jobbnamn
#SBATCH --mail-user=maja.nettelbladt.9823@student.uu.se  # E-post vid jobbets slut
#SBATCH --mail-type=END,FAIL      # Skicka mail vid slut eller fail
#SBATCH --output=%x.%j.out        # Outputfilens namn

# Ladda moduler
module load bioinfo-tools
module load FastQC/0.11.9

# Definiera input och output
INPUT_DIR_RNA_Serum="/home/mane9823/Genome-Analysis/data/raw_data/transcriptomics_data/RNA-Seq_Serum"
INPUT_DIR_RNA_BH="/home/mane9823/Genome-Analysis/data/raw_data/transcriptomics_data/RNA-Seq_BH"

OUTPUT_DIR_Serum="/home/mane9823/Genome-Analysis/analyses/RNA_analyses/01_preprocessing/fastqc_raw/fastqc_raw_serum"
OUTPUT_DIR_BH="/home/mane9823/Genome-Analysis/analyses/RNA_analyses/01_preprocessing/fastqc_raw/fastqc_raw_BH"


# Skapa output-mapp om den inte finns
mkdir -p $OUTPUT_DIR_Serum
mkdir -p $OUTPUT_DIR_BH

# Kör FastQC på Illumina-sekvenser
fastqc -o $OUTPUT_DIR_Serum $INPUT_DIR_RNA_Serum/trim_paired_*.fastq.gz
fastqc -o $OUTPUT_DIR_BH $INPUT_DIR_RNA_BH/trim_paired_*.fastq.gz

echo "FastQC på Illumina-sekvenser klar!"

