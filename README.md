# RNA-sequence-analysis

## To find the differetial gene abundnace between wild type and mutant strains using RNA-seq data post SALMON quantification

The script performs differential gene abundance analysis using DESeq2, a package designed for analyzing count data from RNA-seq and other high-throughput sequencing assays. Here's a breakdown of how the differential gene abundance analysis is carried out in the script:

1. Data Import and Preparation:

Transcript Quantification Data: The script uses tximport to read the transcript-level quantifications from Salmon (quant.sf files) and aggregates them for gene-level analysis.
Metadata Import: The metadata file (metadata.txt) provides information about the samples and their conditions (e.g., Control, Mutant 1, etc.). It is used to map the samples to their respective conditions.
Data Consistency Check: The script ensures that the sample names in the metadata match those in the imported transcript quantifications to maintain consistency.
Data Subsetting:

The script subsets both the transcript quantification (txi) data and metadata to include only the samples that are common between them. This ensures that only relevant data is used for downstream analysis.
DESeqDataSet Creation:

A DESeqDataSet object is created using DESeqDataSetFromTximport(), combining the count data (txi) and metadata (colData). The design formula (~ condition) specifies that the analysis will be based on the experimental conditions provided in the metadata.
Model Fitting with DESeq:

The DESeq() function is run on the DESeqDataSet object to normalize the counts and fit a negative binomial generalized linear model (GLM) for each gene. This model accounts for variability in the data and computes statistical tests to identify differentially abundant genes between conditions.
Comparisons Between Conditions:

The script defines multiple comparisons between different experimental conditions using the results() function:
Control vs. Mutant 1
Control vs. Mutant 2-5
Control vs. Mutant 3-4
Mutant 1 vs. Mutant 2-5
Mutant 1 vs. Mutant 3-4
Mutant 2-5 vs. Mutant 3-4
For each comparison, the function calculates the log fold change (LFC) for each gene, estimating the change in gene abundance between the two conditions.

2. Statistical Testing:

For each gene, DESeq2 performs a Wald test or a likelihood ratio test (LRT) to determine if the difference in abundance between conditions is statistically significant. The script uses the default Wald test.
The script orders the results by the adjusted p-value (padj), applying the Benjamini-Hochberg correction to control for false discovery rate (FDR).
Filtering for Significant Genes:

Genes with an adjusted p-value (padj) less than 0.05 are considered statistically significant. The script subsets these genes for each comparison and outputs them as CSV files.
