#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J mummerplot_pacbio_vs_ref
#SBATCH --output=mummerplot.%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maja.nettelbladt.9823@student.uu.se  # E-post vid jobbets slut

# Ladda moduler
module load bioinfo-tools
module load MUMmer/3.23

# Definiera sökvägar
REF="/home/mane9823/Genome-Analysis/data/raw_data/reference_data/reference_E745.fasta"
QUERY="/home/mane9823/Genome-Analysis/analyses/DNA_analyses/02_genome_assembly/01_pacbio/assembly_results/pacbio_assembly.contigs.fasta"
OUTDIR="/home/mane9823/Genome-Analysis/analyses/DNA_analyses/02_genome_assembly/01_pacbio/mummer_output"

# Skapa output-mapp
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Kör alignment
nucmer --prefix=pacbio_vs_E745 "$REF" "$QUERY"

# Skapa plot

mummerplot --png --large --fat --layout -p pacbio_vs_E745 \
  -R "$REF" -Q "$QUERY" pacbio_vs_E745.delta

