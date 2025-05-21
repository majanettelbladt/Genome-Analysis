# Genome-Analysis
# Genome and Transcriptome Analysis of *Enterococcus faecium* E745

This project contains genome and RNA-seq analyses for the bacterial strain *Enterococcus faecium* E745, based on long-read PacBio data and short-read Illumina data. The work follows the pipeline from Zhang et al. (2017) and includes preprocessing, genome assembly, annotation, RNA mapping, and differential expression analysis.

## Project Overview

- **Genome assembly** was performed with:
  - Canu (PacBio-only)
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
- `data/`: Input and intermediate data

## Reproducibility

All major steps are scripted using SLURM job scripts located in the `code/` folder. Filenames and folders follow a consistent naming convention.

To reproduce the analysis:
1. Follow the preprocessing steps in `code/01_DNA_preprocessing/` and `code/04_RNA_preprocessing/`.
2. Run genome assemblies from `code/02_Genome_Assembly/`.
3. Annotate assemblies with Prokka.
4. Run BWA mapping scripts in `code/05_RNA_mapping/`.
5. Count reads with HTSeq and analyze DEGs in `code/06_DE_analysis/`.

## Reference

Original data and experimental design are based on:
> Zhang, X. et al. (2017). "Comparative genomics and functional analysis of *Enterococcus faecium* E745." *Journal XYZ*. [NCBI BioProject: PRJNA224116]

---


