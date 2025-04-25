#!/bin/bash
#SBATCH -A uppmax2025-3-3        # Projekt-ID
#SBATCH -M snowy                  # Kluster-miljö
#SBATCH -p core                   # Typ av resurser (core)
#SBATCH -n 2                       # Antal kärnor
#SBATCH -t 02:00:00               # Maximal kö-tid
#SBATCH -J prokka_annotation      # Jobbnamn
#SBATCH --mail-user=maja.nettelbladt.9823@student.uu.se  # Din e-postadress
#SBATCH --mail-type=ALL           # När du vill få e-post om jobbstatus
#SBATCH --output=prokka.%j.out    # Utkörningslogg

# Ladda moduler (om det behövs, beroende på din miljö)
module load bioinfo-tools
module load prokka/1.45-5b58020

# Skapa output-mapp (om den inte finns)
mkdir -p /home/mane9823/Genome-Analysis/analyses/DNA_analyses/03_structural_annotation/prokka_output

# Kör Prokka med Enterococcus faecium som genus och art
prokka --outdir /home/mane9823/Genome-Analysis/analyses/03_structural_annotation/prokka_output \
       --prefix pacbio_annotation \
       --genus Enterococcus \
       --species faecium \
       --force \
       --cpus 2 \
       /home/mane9823/Genome-Analysis/analyses/DNA_analyses/02_genome_assembly/01_pacbio/assembly_results/pacbio_assembly.contigs.fasta

