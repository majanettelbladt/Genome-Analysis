#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -J mapping_BH
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maja.nettelbladt.9823@student.uu.se
#SBATCH --output=%x.%j.out

module load bioinfo-tools
module load bwa
module load samtools

# Genom
GENOME="/home/mane9823/Genome-Analysis/analyses/RNA_analyses/02_RNA_mapping/BWA_indexed_genome_files/pacbio_annotation.fna"

# Output-mapp
OUTDIR="/home/mane9823/Genome-Analysis/analyses/RNA_analyses/02_RNA_mapping/mapped_RNA-seq_BH"
mkdir -p "$OUTDIR"

# Lista över BH-prover
for SAMPLE in ERR1797972 ERR1797973 ERR1797974
do
    R1="/home/mane9823/Genome-Analysis/data/raw_data/transcriptomics_data/RNA-Seq_BH/trim_paired_${SAMPLE}_pass_1.fastq.gz"
    R2="/home/mane9823/Genome-Analysis/data/raw_data/transcriptomics_data/RNA-Seq_BH/trim_paired_${SAMPLE}_pass_2.fastq.gz"
    BAM="${OUTDIR}/${SAMPLE}_mapped.bam"

    # Kör BWA + samtools direkt
    bwa mem "$GENOME" "$R1" "$R2" | \
    samtools view -@ 2 -Sb - | \
    samtools sort -@ 2 -o "$BAM" -

    # Indexera bam-filen
    samtools index "$BAM"
done

