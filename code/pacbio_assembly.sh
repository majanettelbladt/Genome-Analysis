#!/bin/bash
#SBATCH -A uppmax2025-3-3         # Projekt-ID
#SBATCH -M snowy                  # Kluster Snowy
#SBATCH -p core                   # Typ av noder
#SBATCH -n 4                      # 4 cores
#SBATCH -t 06:00:00               # Max tid (6 timmar)
#SBATCH -J pacbio_assembly          # Jobbnamn
#SBATCH --mail-user=maja.nettelbladt.9823@student.uu.se  # E-post vid jobbets slut
#SBATCH --mail-type=END,FAIL      # Skicka mail vid slut eller fail
#SBATCH --output=%x.%j.out

# Ladda moduler
module load bioinfo-tools
module load canu

# Definiera input och output
INPUT_DIR="/home/mane9823/Genome-Analysis/data/raw_data/genomics_data/PacBio"
OUTPUT_DIR="/home/mane9823/Genome-Analysis/analyses/02_genome_assembly"

# Skapa output-mapp om den inte finns
mkdir -p $OUTPUT_DIR

# Kör Canu med PacBio data
canu -p pacbio_assembly -d $OUTPUT_DIR \
     genomeSize=3m \
     -pacbio-raw $INPUT_DIR/*.fastq.gz \

echo "Canu assembly done!"

