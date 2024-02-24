# Gene Set Enrichment Analysis of “HGG, H3 wildtype” vs “HGG, H3 wildtype, TP53” using DESeq2 using ClusterProfiler
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia
# Script adapted from ALSF's DE Analysis Tutorial (https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html)

## LOAD LIBRARIES
# Attach the library
library(clusterProfiler)
# Package that contains MSigDB gene sets in tidy format
library(msigdbr)
# Human annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)
# Visualization library
library(ggplot2)
# We will need this so we can use the pipe: %>%
library(magrittr)


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
dge_results_file <- file.path(results_dir, "HGG_H3wildtype_vs_TP53_DEG.csv")

# Read in the contents of the differential expression results file
dge_df <- readr::read_csv(dge_results_file)

#######
# Specifying MSigDB gene sets of interest
hs_hallmark_sets <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "H"
)

# Choose other genesets here:https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H

#######
## PERFORM GSEA

# 1. Determine our pre-ranked genes list
# Check if there are any duplicate genes present
any(duplicated(dge_df$Gene)) # None

# Create a named vector ranked based on the log2 fold change values
lfc_vector <- dge_df$log2FoldChange
names(lfc_vector) <- dge_df$Gene

# Sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# 2. Run GSEA using the GSEA() function
# Set the seed so our results are reproducible:
set.seed(2020)
# Run GSEA
gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_hallmark_sets,
    gs_name,
    gene_symbol
  )
)

# We can access the results from our `gsea_results` object using `@result`
head(gsea_results@result)

# Convert GSEA results object to dataframe
gsea_result_df <- data.frame(gsea_results@result)

# Save GSEA results
readr::write_csv(
  gsea_result_df,
  file.path(
    results_dir,
    "HGG_H3wildtype_TP53_vs_HGG_H3wildtype_gsea_results.csv"
  )
)

# 3. Visualize GSEA
# Look at the 3 gene sets with the most positive NES
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(NES, n = 3)

# Make GSEA plot
nes_plot <- enrichplot::gseaplot2(
  gsea_results,
  geneSetID = "HALLMARK_P53_PATHWAY",
  title = "HALLMARK_P53_PATHWAY"
)
nes_plot

# Save GSEA enrichment plot as tiff
# Plot path
nes_plot_path <- file.path(plots_dir, "HGG_H3wildtype_vs_TP53_P53_pathway_down_plot.png")
# Save plot
ggplot2::ggsave(nes_plot_path,
                width=8,
                height=5,
                device="png"
)

###############################

# Make dotplot of enriched gene sets
require(DOSE)
dotplot(gsea_results, showCategory=10, split=".sign") + facet_grid(.~.sign)

# Save GSEA dotplot as tiff
# Plot path
gsea_dotplot_path <- file.path(plots_dir, "HGG_H3wildtype_TP53_vs_HGG_H3wildtype_gsea_dotplot.png")
# Save plot
ggplot2::ggsave(gsea_dotplot_path,
                width=10,
                height=6,
                device="png"
)

# Make ridgeplot of enriched gene sets
ridgeplot(gsea_results) + labs(x = "enrichment distribution")

# Save GSEA dotplot as tiff
# Plot path
gsea_ridgeplot_path <- file.path(plots_dir, "HGG_H3wildtype_TP53_vs_HGG_H3wildtype_gsea_ridgeplot.png")
# Save plot
ggplot2::ggsave(gsea_ridgeplot_path,
                width=13,
                height=10,
                device="png"
)

###############################################
# Customs Genesets for GSEA
# mito_genes_file <- file.path(data_dir, "Human_Genes_MitoCarta3.csv")
mito_pathways_file <- file.path(data_dir, "MitoPathways3.gmt")
mito_pathways_gs <- read.gmt(mito_pathways_file)
# mito_genes <- readr::read_csv(mito_genes_file)

# Select only mitochondrial genes
# mito_dge_df <- dge_df %>%
#   dplyr::filter(Gene %in% mito_genes$Symbol)

## PERFORM GSEA

# 1. Determine our pre-ranked genes list
# Check if there are any duplicate genes present
# any(duplicated(mito_dge_df$Gene)) # None

# Create a named vector ranked based on the log2 fold change values
# lfc_vector <- mito_dge_df$log2FoldChange
# names(lfc_vector) <- mito_dge_df$Gene

# Sort the log2 fold change values in descending order here
# lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# 2. Run GSEA using the GSEA() function
# Set the seed so our results are reproducible:
set.seed(2020)
# Run GSEA
mito_pathways_gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 10, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mito_pathways_gs,
    term,
    gene
  )
)

# We can access the results from our `gsea_results` object using `@result`
head(mito_pathways_gsea_results@result)

# Convert GSEA results object to dataframe
mito_gsea_result_df <- data.frame(mito_pathways_gsea_results@result)

# Save GSEA results
readr::write_csv(
  gsea_result_df,
  file.path(
    results_dir,
    "HGG_H3wildtype_TP53_vs_HGG_H3wildtype_mito_gsea_results.csv"
  )
)

# Make GSEA plot
nes_plot <- enrichplot::gseaplot2(
  mito_pathways_gsea_results,
  geneSetID = "OXPHOS",
  title = "OXIDATIVE PHOSPHORYLATION"
)
nes_plot

# Save GSEA enrichment plot as tiff
# Plot path
nes_plot_path <- file.path(plots_dir, "HGG_H3wildtype_TP53_vs_HGG_H3wildtype_OXPHOS_down_plot.png")
# Save plot
ggplot2::ggsave(nes_plot_path,
                width=8,
                height=5,
                device="png"
)
