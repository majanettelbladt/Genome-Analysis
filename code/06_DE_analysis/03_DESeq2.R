# Differential Expression script 

# Load all necessary libraries
library("DESeq2")
library("ggplot2")
library("apeglm")
library("pheatmap")
library("dplyr")

# ------------------------ 1. Read data --------------------------


# List reads counts:
serum_files <- list.files("/home/mane9823/Genome-Analysis/analyses/RNA_analyses/03_DE_analysis/01_HTseq/RNA-seq_Serum", pattern="*.txt$", full.names=TRUE)
bh_files <- list.files("/home/mane9823/Genome-Analysis/analyses/RNA_analyses/03_DE_analysis/01_HTseq/RNA-seq_BH", pattern="*.txt$", full.names=TRUE)

count_files <- c(bh_files, serum_files)
print(length(count_files))  

sample_names <- c("BH1", "BH2", "BH3", "Serum1", "Serum2", "Serum3")

# Function to read the HTSeq-files
read.file <- function(file) {
  read.delim(file, col.names=c("gene", "count"), sep="\t", colClasses=c("character", "numeric"))
}

# Read and combine all files
all_data <- read.file(count_files[1])
for (i in 2:length(count_files)) {
  temp <- read.file(count_files[i])
  all_data <- cbind(all_data, temp$count)
}
print(dim(all_data))  


# Remove techincal HTSeq-rows such as "__no_feature" etc.
dsq_data <- all_data[!grepl("^__", all_data$gene), ]
rownames(dsq_data) <- dsq_data$gene
dsq_data <- dsq_data[, -1]  # removing the "gene"-column

# Naming columns, need to match with sample names
colnames(dsq_data) <- sample_names

# ------------------------ 2. Meta data --------------------------

# Change the order to: BH = reference levels
metadata <- data.frame(
  row.names = sample_names,
  condition = factor(c("BH", "BH", "BH", "Serum", "Serum", "Serum"), levels = c("BH", "Serum")),
  libType = rep("paired-end", 6)
)

# ------------------------ 3. DESeq2 analysis --------------------------

dds <- DESeqDataSetFromMatrix(countData = dsq_data, colData = metadata, design = ~ condition)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds)

# ------------------------ 4. Visualisation --------------------------

plotMA(res, main = "MA-plot (default)")
summary(res)

# Shrinked LFC
resLFC <- lfcShrink(dds, coef="condition_Serum_vs_BH", type="apeglm", lfcThreshold=1)
plotMA(resLFC, main = "Shrinked LFC MA-plot")
abline(h=c(-1,1), col="dodgerblue", lwd=2)

# ------------------------ 5. Significant genes --------------------------

sign_genes <- subset(res, padj < 0.001 & abs(log2FoldChange) > 1)
summary(sign_genes)

# ------------------------ 6. Heatmap --------------------------

ntd <- normTransform(dds)

top_upregulated <- rownames(sign_genes[order(-sign_genes$log2FoldChange), ])[1:20]
top_downregulated <- rownames(sign_genes[order(sign_genes$log2FoldChange), ])[1:20]
select_top20 <- c(top_upregulated, top_downregulated)

pheatmap(assay(ntd)[select_top20, ],
         cluster_rows=TRUE,
         show_rownames=TRUE,
         cluster_cols=FALSE,
         annotation_col=metadata,
         main="Top 20 DE genes")

# ------------------------ 7. Prokka-annotation--------------------------

# Fetching Prokka annotation:
prokka_file <- "/home/mane9823/Genome-Analysis/analyses/DNA_analyses/03_structural_annotation/prokka_output/pacbio_annotation.tsv"
prokka_data <- read.delim(prokka_file, header=TRUE, sep="\t")

# Matching significant genes with annotations 
filtered_upregulated <- prokka_data %>% filter(locus_tag %in% top_upregulated)
filtered_downregulated <- prokka_data %>% filter(locus_tag %in% top_downregulated)

write.table(filtered_upregulated, "./filtered_upregulated_annotations.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(filtered_downregulated, "./filtered_downregulated_annotations.txt", sep="\t", row.names=FALSE, quote=FALSE)

