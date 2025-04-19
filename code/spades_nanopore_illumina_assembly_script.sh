#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 03:00:00
#SBATCH -J hybrid_spades
#SBATCH --mail-user=maja.nettelbladt.9823@student.uu.se
#SBATCH --mail-type=ALL
#SBATCH --output=hybrid_spades.%j.out

# Ladda moduler
module load bioinfo-tools
module load spades/4.0.0

# Skapa output-katalog
mkdir -p /home/mane9823/Genome-Analysis/analyses/02_genome_assembly/02_nanopore_illumina_raw

# Kör SPAdes hybrid assembly
spades.py \
-1 /home/mane9823/Genome-Analysis/data/raw_data/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz \
-2 /home/mane9823/Genome-Analysis/data/raw_data/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz \
--nanopore /home/mane9823/Genome-Analysis/data/raw_data/genomics_data/Nanopore/E745_all.fasta.gz \
--isolate -o /home/mane9823/Genome-Analysis/analyses/02_genome_assembly/02_nanopore_illumina_raw

