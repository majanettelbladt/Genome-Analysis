# Genome-Analysis
# Genome and Transcriptome Analysis of *Enterococcus faecium* E745

This project contains genome and RNA-seq analyses for the bacterial strain *Enterococcus faecium* E745, based on long-read PacBio data and short-read Illumina data. The work follows the pipeline from Zhang et al. (2017) and includes preprocessing, genome assembly, annotation, RNA mapping, and DE analysis.

## Project Overview

- **Genome assembly** was performed with:
  - Canu (PacBio)
  - SPAdes (hybrid assembly using Illumina + Nanopore)
- **Annotation** was done with Prokka.
- **RNA-Seq** data were mapped using BWA.
- **Gene expression quantification** was done using HTSeq.
- **Differential expression** was analyzed with DESeq2.
- **Quality assessments** included FastQC, QUAST, and MUMmerplot.

## Folder Structure

The folder tree is documented in [`folders_structure.txt`](./folders_structure.txt). Main folders include:

- `analyses/`: All analysis results (DNA and RNA)
- `code/`: Scripts for running jobs and analyses
- `data/`: Raw and intermediate data

## Reproducibility

The major parts of the analysis are made using SLURM job scripts located. These are available in the `code/` folder. 

To reproduce the analysis:
1. Follow the preprocessing steps for DNA and RNA respectively in `code/01_DNA_preprocessing/` and `code/04_RNA_preprocessing/`.
2. Genome assemblies using `code/02_Genome_Assembly/`.
3. Annotate assemblies with Prokka.
4. Run BWA mapping scripts using `code/05_RNA_mapping/`.
5. Count reads with HTSeq and DE analysis using `code/06_DE_analysis/`.

## Reference

Original article: Zhang, X. et al. (2017). "Comparative genomics and functional analysis of *Enterococcus faecium* E745."
