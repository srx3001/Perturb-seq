library(Seurat)
library(Matrix)
library(ggplot2)
library(argparse)
library(dplyr)

# Set up argparse to accept command line arguments
parser <- ArgumentParser(description = "Seurat analysis for Perturb-seq data")
parser$add_argument('--data_dir', type = "character", required = TRUE, help = "Directory containing filtered feature-barcode matrix (barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz)")
parser$add_argument('--output_dir', type = "character", default = ".", help = "Directory to save output plots and results")
args <- parser$parse_args()

# args <- list()  # Create an empty list to simulate the argparse object
# args$data_dir <- "/home/srdjan/Loka/perturbseq_analysis/outs/per_sample_outs/perturbseq_analysis/count/sample_filtered_feature_bc_matrix"
# args$output_dir <- "results"


#create out_dir so ggplot doesnt complain
if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE)
}

# Load the filtered feature-barcode matrix
data <- Read10X(data.dir = args$data_dir)  # Adjust the path to your dataset

# Step 2: Create Seurat Object for Gene Expression Data
gene_expression_data <- data[["Gene Expression"]]
seurat_obj <- CreateSeuratObject(counts = gene_expression_data, project = "PerturbSeq_GeneExpression")

# Step 3: Add CRISPR Guide Data as Metadata
# Get the CRISPR guide capture data and transpose it to match cells
crispr_guides <- as.data.frame(t(data[["CRISPR Guide Capture"]]))

# Dynamically detect the CRISPR guide columns (non-RNA metadata columns)
guide_columns <- colnames(crispr_guides)

# Step 4: Create a single CRISPR guide identity column
# For each cell, determine which guide was assigned by checking for non-zero values
seurat_obj$CRISPR_Guide <- apply(crispr_guides, 1, function(x) {
  guide_name <- guide_columns[which(x != 0)]
  if (length(guide_name) == 0) return("None")  # For cells with no guide
  return(guide_name[1])  # If more than one guide, pick the first (you can customize this)
})

# Step 5: Normalize the Gene Expression Data
seurat_obj <- NormalizeData(seurat_obj)

# Step 6: Identify Variable Features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Step 7: Scale the Data
seurat_obj <- ScaleData(seurat_obj)

# Step 8: Run PCA
seurat_obj <- RunPCA(seurat_obj)

# Step 9: Cluster the Cells
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Step 10: Perform Generalized Differential Expression Analysis Across CRISPR Guides
# Now that we have the CRISPR_Guide identity column, we can use it for DE analysis
seurat_obj_filtered <- subset(seurat_obj, subset = CRISPR_Guide != "None")

# Since we arent certain on experimental design we can just run on all objects
de_markers_all <- FindAllMarkers(seurat_obj_filtered)

# View the top differentially expressed markers
print(head(de_markers_all))

# Step 11: Visualize the Top Differentially Expressed Genes

top_genes <- de_markers_all %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.table(top_genes, paste0(args$output_dir, "/top_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
de_markers_all$significant <- ifelse(de_markers_all$p_val < 0.05 & abs(de_markers_all$avg_log2FC) > 1, "Significant", "Not Significant")


heatmap_plot  <- DoHeatmap(seurat_obj_filtered, features = top_genes$gene)
ggsave(paste0(args$output_dir, "/seurat_heatmap.png"), plot = heatmap_plot, width = 8, height = 6, dpi = 300)

# Generate the volcano plot
ggplot(de_markers_all, aes(x = avg_log2FC, y = -log10(p_val), color = significant)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.position = "top")
ggsave(paste0(args$output_dir, "/volcano_plot.png"), plot = last_plot(), width = 8, height = 6, dpi = 300)



# CRISPR_Guide based DE

Idents(seurat_obj_filtered) <- "CRISPR_Guide"
de_markers_all_guide <- FindAllMarkers(seurat_obj_filtered)


top_genes_guide <- de_markers_all_guide %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.table(top_genes_guide, paste0(args$output_dir, "/top_genes_guide.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

de_markers_all_guide$significant <- ifelse(de_markers_all_guide$p_val < 0.05 & abs(de_markers_all_guide$avg_log2FC) > 1, "Significant", "Not Significant")


# Generate the volcano plot
ggplot(de_markers_all_guide, aes(x = avg_log2FC, y = -log10(p_val), color = significant)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(paste0(args$output_dir, "/volcano_plot_guide.png"), plot = last_plot(), width = 8, height = 6, dpi = 300)
