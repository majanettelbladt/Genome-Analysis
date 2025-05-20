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

# Filter significant genes based on the article's criteria: q < 0.001 and fold change > 2 or < 0.5
sign_genes <- res[which(res$padj < 0.001 & abs(res$log2FoldChange) > 1), ]

# Print total number of significant genes according to the criteria
cat("Number of significantly differentially expressed genes (padj < 0.001 and |log2FC| > 1):", nrow(sign_genes), "\n")

# Separate into up- and downregulated genes
upregulated <- sign_genes[sign_genes$log2FoldChange > 1, ]
downregulated <- sign_genes[sign_genes$log2FoldChange < -1, ]

cat("Number of upregulated genes:", nrow(upregulated), "\n")
cat("Number of downregulated genes:", nrow(downregulated), "\n")

# Save results to files
write.table(upregulated, "./significant_upregulated_genes.txt", sep="\t", row.names=TRUE, quote=FALSE)
write.table(downregulated, "./significant_downregulated_genes.txt", sep="\t", row.names=TRUE, quote=FALSE)

# Optionally save the full list of significant genes as well
write.table(sign_genes, "./significant_genes_all.txt", sep="\t", row.names=TRUE, quote=FALSE)


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

