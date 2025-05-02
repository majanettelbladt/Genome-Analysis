#!/bin/bash
#SBATCH -A uppmax2025-3-3         
#SBATCH -M snowy                  
#SBATCH -p core                   
#SBATCH -n 2                      
#SBATCH -t 02:00:00               
#SBATCH -J 03_fastqc_post_trim
#SBATCH --mail-user=maja.nettelbladt.9823@student.uu.se
#SBATCH --mail-type=ALL    
#SBATCH --output=%x.%j.out        

# Ladda moduler
module load bioinfo-tools
module load FastQC/0.11.9

# Definiera input och output
INPUT_DIR_RNA_Serum="/home/mane9823/Genome-Analysis/data/trimmed_data/RNA-Seq_Serum"
INPUT_DIR_RNA_BH="/home/mane9823/Genome-Analysis/data/trimmed_data/RNA-Seq_BH"

OUTPUT_DIR_Serum="/home/mane9823/Genome-Analysis/analyses/RNA_analyses/01_preprocessing/fastqc_trimmed/fastqc_trimmed_serum"
OUTPUT_DIR_BH="/home/mane9823/Genome-Analysis/analyses/RNA_analyses/01_preprocessing/fastqc_trimmed/fastqc_trimmed_BH"

mkdir -p $OUTPUT_DIR_Serum
mkdir -p $OUTPUT_DIR_BH

# Kör FastQC på Illumina-sekvenser
fastqc -o $OUTPUT_DIR_Serum $INPUT_DIR_RNA_Serum/*.fastq.gz

fastqc -o $OUTPUT_DIR_BH $INPUT_DIR_RNA_BH/*.fastq.gz

echo "FastQC på Illumina-sekvenser klar!"

