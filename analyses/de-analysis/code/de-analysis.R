# Differential Expression Analysis of “HGG, H3 wildtype” vs “HGG, H3 wildtype, TP53” using DESeq2
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia
# Script adapted from ALSF's DE Analysis Tutorial (https://alexslemonade.github.io/refinebio-examples/03-rnaseq/differential-expression_rnaseq_01.html)

## LOAD LIBRARIES
# Library for DE Analysis
library(DESeq2)
# Library data manipulation
#library(dplyr)
# Library for Plotting
library(ggplot2)
# We will need this so we can use the pipe: %>%
library(magrittr)

# Set seed because jitter plot function involves randomness
set.seed(1234)


## SET DIRECTORIES
analysis_dir <- "/Users/shehbeel/Documents/PBTA-Mito/analyses/de-analysis"
data_dir <- "/Users/shehbeel/Documents/PBTA-Mito/data"

# Set output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

# Make output directories if they don't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Declare input file paths
metadata_file <- file.path(data_dir, "openpbta_v12_metadata.csv")
data_file <- file.path(data_dir, "openpbta_v12_counts.csv")

#######
# Import metadata and data 
metadata <- readr::read_csv(metadata_file)
expression_df <- readr::read_csv(data_file)

#######
## PREPROCESS THE DATA
# Select only the HGAT Sample_IDs
hgg_metadata <- metadata %>%
  dplyr::filter(short_histology == "HGAT") %>%
  dplyr::filter(molecular_subtype == "HGG, H3 wildtype" | molecular_subtype == "HGG, H3 wildtype, TP53") %>%
  dplyr::filter(tumor_descriptor == "Initial CNS Tumor") %>%
  dplyr::filter(CNS_region == "Hemispheric")

expression_df <- expression_df %>% 
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(hgg_metadata$Kids_First_Biospecimen_ID)

# Check if this is in the same order
all.equal(colnames(expression_df), hgg_metadata$Kids_First_Biospecimen_ID)


# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

## Create DESeq2Dataset
# round all expression counts
gene_matrix <- round(filtered_expression_df)

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = hgg_metadata,
  # Supply our experimental variable to `design`
  design = ~molecular_subtype
)

## Run Differential Expression Analysis
deseq_object <- DESeq(ddset)
#resultsNames(deseq_object)

# Export normalized read counts
normCounts <- counts(deseq_object, normalized=TRUE)
# Save results as CSV
write.csv(normCounts, file.path(results_dir, "HGG_H3wildtype_TP53_vs_HGG_H3wildtype_gsea_normalized_counts.csv"))

# Extract results table showing "LOW" as reference
deseq_results <- results(deseq_object, contrast=c("molecular_subtype", "HGG, H3 wildtype, TP53", "HGG, H3 wildtype"), cooksCutoff=FALSE)

# Use lfcShrink() function to obtain shrunken log fold change estimates based on 
# negative binomial distribution. This will add the estimates to your results table. 
# Using lfcShrink() can help decrease noise and preserve large differences between 
# groups (it requires that apeglm package be installed) (Zhu et al., Bioinformatics 2018).
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)

# Sort and filter DESeq2 results table and convert to dataframe
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(deseq_df, file.path(results_dir, "HGG_H3wildtype_TP53_vs_HGG_H3wildtype_DEG.csv"))


# View results sorted by adjusted p-value
deseq_df %>%
  dplyr::arrange(padj)

# Check results by plotting one gene
plotCounts(ddset, gene = "TP53", intgroup = "molecular_subtype")


## Create volcano plot
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
volcano_plot
# Save volcano plot
ggsave(
  plot = volcano_plot,
  file.path(plots_dir, "HGG_H3wildtype_vs_TP53_volcano_plot.png"),
  width = 10,
  height = 8
)


