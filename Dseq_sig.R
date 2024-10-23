library(DESeq2)
library(tximport)

# Define the directory containing Salmon outputs
salmon_dir <- ".../Work/RNA"
salmon_folders <- list.dirs(salmon_dir, full.names = TRUE, recursive = FALSE)
salmon_folders <- salmon_folders[grepl("_sal_out$", salmon_folders)]

# Construct file paths for quant.sf files
files <- file.path(salmon_folders, "quant.sf")
names(files) <- basename(salmon_folders)

# Import transcript quantification data
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Read metadata correctly
metadata <- read.csv("metadata.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(metadata) <- c("sample", "condition")
rownames(metadata) <- metadata$sample

# Check the dimensions and row names of metadata and txi
print(dim(metadata))
print(rownames(metadata))
print(dim(txi$counts))
print(colnames(txi$counts))

# Ensure the sample names match between txi and metadata
if (!all(colnames(txi$counts) %in% rownames(metadata))) {
  stop("Sample names in txi and metadata do not match")
}

# Subset txi and metadata to match the sample names
common_samples <- intersect(colnames(txi$counts), rownames(metadata))
txi_subset <- txi
txi_subset$counts <- txi$counts[, common_samples]
txi_subset$abundance <- txi$abundance[, common_samples]
txi_subset$length <- txi$length[, common_samples]
metadata_subset <- metadata[common_samples, ]

# Create DESeqDataSet object
dds <- DESeqDataSetFromTximport(txi_subset, colData = metadata_subset, design = ~ condition)

# Ensure the first level of the condition factor is the Control
dds$condition <- relevel(dds$condition, ref = "Control")

# Run the DESeq function
dds <- DESeq(dds)

# Get results for Control vs. Mutant 1
res_mutant1 <- results(dds, contrast = c("condition", "Mutant 1", "Control"))
res_mutant1 <- res_mutant1[order(res_mutant1$padj), ]
sig_genes_mutant1 <- subset(res_mutant1, padj < 0.05)
write.csv(as.data.frame(res_mutant1), file = "deseq2_results_mutant1.csv")
write.csv(as.data.frame(sig_genes_mutant1), file = "significant_genes_mutant1.csv")

# Get results for Control vs. Mutant 2-5
res_mutant2_5 <- results(dds, contrast = c("condition", "Mutant 2-5", "Control"))
res_mutant2_5 <- res_mutant2_5[order(res_mutant2_5$padj), ]
sig_genes_mutant2_5 <- subset(res_mutant2_5, padj < 0.05)
write.csv(as.data.frame(res_mutant2_5), file = "deseq2_results_mutant2_5.csv")
write.csv(as.data.frame(sig_genes_mutant2_5), file = "significant_genes_mutant2_5.csv")

# Get results for Control vs. Mutant 3-4
res_mutant3_4 <- results(dds, contrast = c("condition", "Mutant 3-4", "Control"))
res_mutant3_4 <- res_mutant3_4[order(res_mutant3_4$padj), ]
sig_genes_mutant3_4 <- subset(res_mutant3_4, padj < 0.05)
write.csv(as.data.frame(res_mutant3_4), file = "deseq2_results_mutant3_4.csv")
write.csv(as.data.frame(sig_genes_mutant3_4), file = "significant_genes_mutant3_4.csv")

res_mutant1_vs_mutant2_5 <- results(dds, contrast = c("condition", "Mutant 2-5", "Mutant 1"))
res_mutant1_vs_mutant2_5 <- res_mutant1_vs_mutant2_5[order(res_mutant1_vs_mutant2_5$padj), ]
sig_genes_mutant1_vs_mutant2_5 <- subset(res_mutant1_vs_mutant2_5, padj < 0.05)
write.csv(as.data.frame(res_mutant1_vs_mutant2_5), file = "deseq2_results_mutant1_vs_mutant2_5.csv")
write.csv(as.data.frame(sig_genes_mutant1_vs_mutant2_5), file = "significant_genes_mutant1_vs_mutant2_5.csv")

# Mutant 1 vs. Mutant 3-4
res_mutant1_vs_mutant3_4 <- results(dds, contrast = c("condition", "Mutant 3-4", "Mutant 1"))
res_mutant1_vs_mutant3_4 <- res_mutant1_vs_mutant3_4[order(res_mutant1_vs_mutant3_4$padj), ]
sig_genes_mutant1_vs_mutant3_4 <- subset(res_mutant1_vs_mutant3_4, padj < 0.05)
write.csv(as.data.frame(res_mutant1_vs_mutant3_4), file = "deseq2_results_mutant1_vs_mutant3_4.csv")
write.csv(as.data.frame(sig_genes_mutant1_vs_mutant3_4), file = "significant_genes_mutant1_vs_mutant3_4.csv")

# Mutant 2-5 vs. Mutant 3-4
res_mutant2_5_vs_mutant3_4 <- results(dds, contrast = c("condition", "Mutant 3-4", "Mutant 2-5"))
res_mutant2_5_vs_mutant3_4 <- res_mutant2_5_vs_mutant3_4[order(res_mutant2_5_vs_mutant3_4$padj), ]
sig_genes_mutant2_5_vs_mutant3_4 <- subset(res_mutant2_5_vs_mutant3_4, padj < 0.05)
write.csv(as.data.frame(res_mutant2_5_vs_mutant3_4), file = "deseq2_results_mutant2_5_vs_mutant3_4.csv")
write.csv(as.data.frame(sig_genes_mutant2_5_vs_mutant3_4), file = "significant_genes_mutant2_5_vs_mutant3_4.csv")













