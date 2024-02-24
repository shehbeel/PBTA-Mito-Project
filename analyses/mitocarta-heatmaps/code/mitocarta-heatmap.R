# Mitochondrial Gene Expression Heatmaps for Open Pediatric Brain Tumor Atlas (OpenPBTA Dataset)
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia

## LOAD LIBRARIES
library(dplyr)
library(ComplexHeatmap)

################################################################################
## SET DIRECTORIES
analysis_dir <- "/Users/shehbeel/Documents/PBTA-Mito/analyses/mitocarta-heatmaps"
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
#mito_genes_data_file <- file.path(data_dir, "Human_Genes_MitoCarta3.csv")
mito_pathways_file <- file.path(data_dir, "MitoPathways3.gmt")
#data_file <- file.path(data_dir, "gene-counts-rsem-expected_count-collapsed.rds")
#expression_df <- readr::read_rds(data_file)
#######
# Import metadata and data 
metadata <- readr::read_csv(metadata_file)
expression_df <- readr::read_csv(data_file)
#mito_genes_df <- readr::read_csv(mito_genes_data_file)
mito_pathways_gs <- read.gmt(mito_pathways_file)

# metadata <- metadata %>%
#   dplyr::filter(cohort == "PBTA" & experimental_strategy == "RNA-Seq")
# expression_df <- expression_df %>%
#   dplyr::select(metadata$Kids_First_Biospecimen_ID)

# Select only mitochondrial-related genes
# expression_df <- expression_df %>%
#   tibble::rownames_to_column(var = "Gene") 
# expression_df <- expression_df %>%
#   dplyr::filter(Gene %in% mito_genes_df$Symbol)

# write.csv(metadata, "/Users/shehbeel/Documents/PBTA-Mito/data/openpbta_v12_metadata.csv")
# write.csv(expression_df, "/Users/shehbeel/Documents/PBTA-Mito/data/openpbta_v12_mito_counts.csv")


#######
## PREPROCESS THE DATA
# Select all HGAT histology samples in Metadata
hgat_metadata <- metadata %>%
  dplyr::filter(short_histology == "HGAT") %>%
  dplyr::filter(molecular_subtype == "HGG, H3 wildtype" | molecular_subtype == "HGG, H3 wildtype, TP53") %>%
  dplyr::filter(tumor_descriptor == "Initial CNS Tumor") %>%
  dplyr::filter(CNS_region == "Hemispheric")
  
# Select all counts in Expression counts data
hgat_expression_df <- expression_df %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(hgat_metadata$Kids_First_Biospecimen_ID)

# Select only mitochondrial-related genes
hgat_expression_df <- hgat_expression_df %>%
  tibble::rownames_to_column(var = "Gene") 
# Select OXPHOS genes
oxphos_genes <- mito_pathways_gs %>%
  dplyr::filter(term == "OXPHOS")
hgat_expression_df <- hgat_expression_df %>%
  dplyr::filter(Gene %in% oxphos_genes$gene) %>%
  tibble::column_to_rownames("Gene")



hgat_expression_mat <- hgat_expression_df %>%
  scale() %>%
  as.matrix()
  

#################
## HEATMAP ##

# Set up column annotation from metadata
col_annot_df <- hgat_metadata %>%
  # Only select the treatment and sample ID columns
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype) %>%
  # Add on the eigengene expression by joining with sample IDs
  # dplyr::inner_join(module_eigengene, by = "refinebio_accession_code") %>%
  # Arrange by cluster and short_histology
  dplyr::arrange(molecular_subtype) %>%
  # Store sample
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")

# Create the ComplexHeatmap column annotation object
col_annot <- ComplexHeatmap::HeatmapAnnotation(
  # Supply Cluster labels
  molecular_subtype = col_annot_df$molecular_subtype,
  # Pick colors for each miRNA cluster in cluster
  col = list(molecular_subtype = c("HGG, H3 wildtype, TP53" = "#0073C2FF", "HGG, H3 wildtype" = "#CD534CFF")
  )
)

# Create a color function based on standardized scale
color_func <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#67a9cf", "#f7f7f7", "#ef8a62")
  # c("blue", "white", "red")
)

ComplexHeatmap::Heatmap(hgat_expression_mat, 
                        name = "z-score",
                        show_row_names = FALSE,
                        row_names_side = "left",
                        show_column_names = FALSE,
                        cluster_rows = TRUE, 
                        cluster_columns = FALSE,
                        top_annotation = col_annot,
                        col = color_func#,
                        #width = unit(20, "cm"),
                        #height = unit(20, "cm")
)

