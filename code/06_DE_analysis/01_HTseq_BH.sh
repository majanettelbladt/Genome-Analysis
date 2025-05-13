#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -J read_counting_BH
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maja.nettelbladt.9823@student.uu.se
#SBATCH --output=%x.%j.out

module load bioinfo-tools
module load htseq/2.0.2

# Stigar
GFF_FILE="/home/mane9823/Genome-Analysis/analyses/DNA_analyses/03_structural_annotation/prokka_output/pacbio_annotation_clean.gff"

OUTPUT_DIR="/home/mane9823/Genome-Analysis/analyses/RNA_analyses/03_DE_analysis/01_HTseq/RNA-seq_BH"

# Loopa igenom alla .bam-filer i BH-mappen
for BAM in /home/mane9823/Genome-Analysis/analyses/RNA_analyses/02_RNA_mapping/mapped_RNA-seq_BH/*.bam; do
           
    BASENAME=$(basename "$BAM" _mapped.bam)
    htseq-count \
        --format=bam \
        --order=pos \
        --type=CDS \
        --idattr=locus_tag \
        --mode=union \
        --stranded=reverse \
        "$BAM" "$GFF_FILE" > "$OUTPUT_DIR/${BASENAME}_counts.txt"
done

