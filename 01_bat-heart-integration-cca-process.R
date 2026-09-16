# ==============================================================================
# Script: Bat Heart scRNA-seq Full Pipeline (Integration to Visualization)
# Description: SCTransform integration, UMAP projection, Annotation, and Plotting
# ==============================================================================

# --- 1. Environment Setup & Dependencies ---
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(khroma)
library(ggrepel)
library(qs)
library(tibble)

options(future.globals.maxSize = 300000 * 1024^2)
data.path <- "/your/data/path/" # Replace with your actual path if needed

# ==============================================================================
# PART 0: DATA LOADING & DOUBLET REMOVAL
# ==============================================================================
cat("Starting Part 0: Data Loading and QC...\n")

# Define the files to load (sample files stored in qs format after DoubletFinder)
# Automatically find all files ending in "unrmd.qs" in your directory
file_paths <- list.files(path = data.path, pattern = "unrmd\\.qs$")

# Automatically generate clean names for the barcodes (e.g., "CSP-1")
# This extracts the species and number (e.g., "CSP.1") and replaces the dot with a dash
clean_names <- stringr::str_extract(file_paths, "(CSP|MFU|MUS|RSI)\\.\\d+")
clean_names <- stringr::str_replace(clean_names, "\\.", "-")

# Create the final named list for the lapply function
sample_files <- as.list(file_paths)
names(sample_files) <- clean_names

# Load each file, isolate Singlets, and store in a list
seurat_list <- lapply(sample_files, function(file) {
  obj <- qs::qread(paste0(data.path, file))
  Idents(obj) <- obj$DoubletFinder
  obj <- subset(obj, idents = "Singlet")
  return(obj)
})

# Merge all Singlet objects into one combined Seurat object
x.merge <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list)
)

# Extract Sample and Species IDs from cell barcodes
x.merge$sample <- stringr::str_split(colnames(x.merge), '_', simplify = TRUE)[, 1]
x.merge$species <- stringr::str_split(colnames(x.merge), '-', simplify = TRUE)[, 1]

# Calculate Mitochondrial Percentage
x.merge[["percent.mt"]] <- PercentageFeatureSet(x.merge, pattern = "^MT-")

# Apply uniform filtration criteria
cat("Applying uniform filtration criteria...\n")
x.merge <- subset(x.merge, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 5)

# Save the full filtered & merged dataset
qs::qsave(x.merge, "bat-heart-merged.qs")

# ==============================================================================
# PART 1: DATA FILTERING & INTEGRATION
# ==============================================================================
cat("Starting Part 1: Data Integration...\n")

# Load merged unintegrated data
x.merge <- qs::qread("bat-heart-merged.qs")
Idents(x.merge) <- x.merge$sample

# Subset to target samples
target_samples <- c("CSP-2", "CSP-4", "MFU-3", "MFU-4", "MUS-1", "MUS-2", "RSI-1", "RSI-4")
x.merge <- subset(x.merge, idents = target_samples)

# Split and perform SCTransform independently per sample
obj_list <- SplitObject(x.merge, split.by = 'sample')
for (i in 1:length(obj_list)) {
  obj_list[[i]] <- SCTransform(obj_list[[i]], verbose = FALSE)
}

# Integrate using CCA + SCT
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features, verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT", anchor.features = features, verbose = FALSE)
x <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

# Dimensionality Reduction
x <- RunPCA(object = x, verbose = FALSE)
x <- RunUMAP(object = x, dims = 1:30, verbose = FALSE)
x <- FindNeighbors(object = x, dims = 1:30, verbose = FALSE)

# ==============================================================================
# PART 2: CLUSTERING & ANNOTATION
# ==============================================================================
cat("Starting Part 2: Clustering and Annotation...\n")

# Cluster at resolution 0.1 to obtain the 12 base clusters
x <- FindClusters(object = x, verbose = FALSE, resolution = 0.1)

# Prepare RNA assay for downstream marker testing
DefaultAssay(x) <- "RNA"
x <- NormalizeData(x)
x <- ScaleData(x, features = rownames(x))

# Apply manual annotations based on resolution 0.1 clusters (0 through 11)
cell.type2 <- c("Cardiomyocyte", "Fibroblast", "Lym-Endothelia", "Schwann", 
                "Vas-Endothelia", "Macrophage", "Pericyte", "Cardiomyocyte", 
                "T", "Adipocyte", "Meso", "Leukocyte")
names(cell.type2) <- levels(x)

x <- RenameIdents(x, cell.type2)
x$cell.type2 <- Idents(x)

# Save the fully integrated and annotated object
qs::qsave(x, "bat-heart-cca-annotated.qs")

# ==============================================================================
# PART 3: COLOR PALETTE SETUP
# ==============================================================================
# Define target clusters and lock factor levels
original_clusters <- c("Cardiomyocyte", "Vas-Endothelia", "Lym-Endothelia", 
                       "Fibroblast", "Pericyte", "Adipocyte", "Mesothelia", "Schwann", 
                       "Macrophage", "T", "Leukocyte")
					   
# Create numeric IDs for cleaner plotting
numeric_ids <- as.character(seq_along(original_clusters))
names(numeric_ids) <- original_clusters
x$cluster_id_numeric <- factor(x$cell.type2, levels = original_clusters, labels = numeric_ids)

# Define robust color palette via khroma
nClust <- length(original_clusters)
smooth_rainbow <- khroma::color("smooth rainbow", reverse = TRUE)
my_colors_named <- smooth_rainbow(nClust, range = c(0.1, 0.7))
names(my_colors_named) <- original_clusters # Bind directly to cell type names for universal mapping

# Legend strings (e.g., "1 - Cardiomyocyte-1")
combined_legend_labels <- paste0(numeric_ids, " - ", original_clusters)
names(combined_legend_labels) <- original_clusters

# ==============================================================================
# PART 4: UMAP VISUALIZATIONS
# ==============================================================================
cat("Starting Part 4: Plotting UMAPs...\n")

x$species <- factor(x$species, levels = c("MUS", "CSP", "MFU", "RSI"))

# A. Global UMAP with numbered labels
centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cell.type2, 
                     data = FetchData(x, vars = c("UMAP_1", "UMAP_2", "cell.type2")), 
                     FUN = median)
centers$label_num <- numeric_ids[as.character(centers$cell.type2)]

p_global <- DimPlot(x, reduction = "umap", group.by = "cell.type2", cols = my_colors_named, raster = TRUE) + 
  geom_label_repel(data = centers, 
                   aes(x = UMAP_1, y = UMAP_2, label = label_num, fill = cell.type2),
                   color = "white", size = 5, box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"), segment.color = 'grey50',
                   show.legend = FALSE, label.r = unit(0.5, "lines")) +
  scale_fill_manual(values = my_colors_named) +
  scale_color_manual(values = my_colors_named, labels = combined_legend_labels, name = "Cell Type") +
  theme_classic() +
  ggtitle("Bat Heart scRNA-seq Atlas")

pdf("bat-heart-cca-global-umap.pdf", width = 10, height = 8)
print(p_global)
dev.off()

# B. Split UMAPs (Species & Sample)
pdf("bat-heart-cca-split-species.pdf", width = 20, height = 5)
DimPlot(x, split.by = 'species', cols = my_colors_named, label = TRUE, raster = TRUE) + 
  theme_classic() + NoLegend()
dev.off()

pdf("bat-heart-cca-split-sample.pdf", width = 20, height = 10)
DimPlot(x, split.by = 'sample', cols = my_colors_named, label = TRUE, raster = TRUE, ncol = 4) + 
  theme_classic() + NoLegend()
dev.off()

# --- C. EXPORT METADATA, UMAP COORDS, COLORS, AND MATRIX ---
cat("Exporting metadata, coordinates, and plain text matrix...\n")

# 1. Extract Metadata and UMAP coordinates
meta_export <- x@meta.data
umap_coords <- Embeddings(x, "umap")
export_df <- cbind(meta_export, umap_coords)

# 2. Map exact plotting colors to each cell based on cell.type2
export_df$plot_color <- my_colors_named[as.character(export_df$cell.type2)]

# 3. Save combined metadata + coords + colors to CSV
write.csv(export_df, "bat-heart-metadata-umap-colors.csv", row.names = TRUE)

# 4. Export the normalized expression matrix as a plain table
# Note: as.matrix() expands the sparse matrix to dense format. Ensure your server has sufficient RAM.
exp_matrix <- GetAssayData(x, assay = "RNA", slot = "data") 
write.csv(as.matrix(exp_matrix), file = "bat-heart-rna-normalized-matrix.csv", quote = FALSE)

# ==============================================================================
# PART 5: MARKER ANALYSIS & DOTPLOT
# ==============================================================================
cat("Starting Part 5: Marker Analysis...\n")

# Calculate All Markers
x.markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(x.markers, "bat-heart-cca-markers.csv", row.names = FALSE)

# Plot Canonical Markers
canonical_markers <- c(
  "RYR2", "ACTN2",      # Cardiomyocyte
  "FLT1", "ADGRL4",     # Vas-Endothelia
  "LYVE1", "CCL21",     # Lym-Endothelia
  "LUM", "DCN",         # Fibroblast
  "MYO1B", "KCNJ8",     # Pericyte
  "PPARG", "PNPLA2",    # Adipocyte
  "PARVA", "IGFBP5",    # Meso
  "ZNF536", "CADM2",    # Schwann
  "CD163", "C1QA",      # Macrophage
  "IL7R", "PTPRC",      # T
  "CD74", "HLA-DRA"     # Leukocyte
)            

pdf("bat-heart-canonical-markers-dotplot.pdf", height = 6, width = 9)
DotPlot(x, features = canonical_markers, dot.scale = 8, cols = "RdYlBu") + 
  RotatedAxis() + 
  scale_y_discrete(limits = rev(levels(x))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Marker Genes", y = "Cell Types")
dev.off()

# ==============================================================================
# PART 6: CELL COMPOSITION STACKED BARPLOT
# ==============================================================================
cat("Starting Part 6: Composition Barplot...\n")

# Extract metadata and calculate cell counts per sample
meta_df <- x@meta.data
comp_data <- meta_df %>%
  group_by(sample, species, cell.type2) %>%
  summarise(CellCount = n(), .groups = 'drop')

# Plot stacked barplot of cell counts
p_comp <- ggplot(comp_data, aes(x = sample, y = CellCount, fill = cell.type2)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = my_colors_named, name = "Cell Type") +
  theme_bw() +
  labs(title = "Cell Type Composition per Sample",
       x = "Sample",
       y = "Number of Cells") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(. ~ species, scales = "free_x", space = "free_x") # Group bars visually by species

pdf("bat-heart-cell-composition-barplot.pdf", width = 8, height = 6)
print(p_comp)
dev.off()

cat("Pipeline completed successfully!\n")