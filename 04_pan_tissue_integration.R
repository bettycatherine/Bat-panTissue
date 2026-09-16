# ==============================================================================
# Script Name: Global_Integration_and_Annotation.R
# Description: Cross-species integration (Harmony), cell type annotation, 
#              and validation (MetaNeighbor) for the Bat Multi-Organ Atlas.
# Species: MUS (Mouse), CSP, MFU, RSI
# ==============================================================================

# ==============================================================================
# 1. ENVIRONMENT SETUP
# ==============================================================================
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(qs)
library(RColorBrewer)
library(MetaNeighbor)
library(SingleCellExperiment)
library(ROGUE)

# Set working environment and options
options(future.globals.maxSize = 300000 * 1024^2)

# ==============================================================================
# 2. DATA LOADING & PREPROCESSING
# ==============================================================================
# Load the pre-merged dataset (Brain, Heart, Kidney, Liver, Lung, Spleen)
x <- qs::qread("bat-merged-v5.qs")

# Format RNA Assay and calculate mitochondrial percentage
x[["RNA"]] <- as(object = x[["RNA"]], Class = "Assay5")
x[["RNA"]] <- JoinLayers(x[["RNA"]])
x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")

# Standard Seurat normalization and scaling
x <- NormalizeData(x)
x <- FindVariableFeatures(x, verbose = FALSE)
x <- ScaleData(x, verbose = FALSE)
x <- RunPCA(x, verbose = FALSE)

# ==============================================================================
# 3. CROSS-SPECIES & TISSUE INTEGRATION (HARMONY)
# ==============================================================================
# Run Harmony using 'spts' (species + tissue) to correct complex batch effects
x <- RunHarmony(
  object = x,
  group.by.vars = "spts",
  reduction.use = "pca",
  dims = 1:30,
  verbose = TRUE
)

# Clustering and UMAP based on Harmony reduction
x <- FindNeighbors(x, reduction = "harmony", dims = 1:30)
x <- FindClusters(x, resolution = 1, algorithm = 4)
x <- RunUMAP(x, reduction = "harmony", dims = 1:30, reduction.name = "umap")

# Save the post-integration object
# qs::qsave(x, "bat-4sp-harmony.qs")

# ==============================================================================
# 4. UNIFIED CELL TYPE ANNOTATION
# ==============================================================================
# Dictionary to map fine-grained clusters into unified global lineages
cell_type_dict <- c(
  "Hepatocyte", "Stellate", "Lymphocyte", "Endothelia", "Cholangiocyte", 
  "Kupffer", "Hepatoblasts", "Micro", "ExN", "InN", "Oligodendrocyte", "ExN", "Astrocyte", 
  "MSN", "ExN", "ExN", "Endothelia", "ExN", "ExN", "ExN", "OPC", "ExN", "ExN", 
  "ExN", "Granule", "Purkinje", "Monocyte", "T", "B", "Fibroblast", "NKT", 
  "Macrophage", "Mesothelia", "DC", "SM", "Hematopoietic_Stem", "T", "Erythroblast", 
  "Macrophage", "cDC", "Plasma", "pDC", "Leukocyte", "Proximal_Tubule", 
  "Adipocyte", "Fibroblast", "Collecting_Duct_Principal", 
  "Collecting_Duct_Intercalated", "Collecting_Duct_Principal", "Podocyte", 
  "Stromal", "Distal_Convoluted_Tubular", "Neuron", "Non_PT_epithelia", "T", 
  "AT1", "Secretory", "Myeloid", "Macrophage", "Ciliated", "AT2", "Cycling", 
  "Cardiomyocyte", "Endothelia", "Pericyte", "Endothelia", "Cardiomyocyte", "Schwann"
)

# Map dictionary names to current active idents
names(cell_type_dict) <- levels(x$cell.type)
x <- RenameIdents(x, cell_type_dict)
x$cell.type2 <- Idents(x)

# Reorder factor levels for logical downstream plotting (Immune -> Epithelial -> Neural -> Stroma)
unified_levels <- c(
  "Leukocyte", "Lymphocyte", "B", "Plasma", "NKT", "T", "Myeloid", "Monocyte", 
  "Kupffer", "Microglia", "Macrophage", "DC", "cDC", "pDC", "Cycling", 
  "Hematopoietic_Stem", "Erythroblast", "Collecting_Duct_Intercalated", 
  "Proximal_Tubule", "Collecting_Duct_Principal", "Distal_Convoluted_Tubular", 
  "Non_PT_epithelia", "AT1", "AT2", "Secretory", "Ciliated", "Hepatocyte", 
  "Hepatoblasts", "Cholangiocyte", "Podocyte", "ExN", "InN", "Purkinje", "MSN", 
  "Granule", "Neuron", "Astrocyte", "Oligodentrocyte", "OPC", "Schwann", "Stellate", "Mesothelia", 
  "Fibroblast", "Adipocyte", "Stromal", "Endothelia", "Pericyte", 
  "Cardiomyocyte", "SM"
)
x$cell.type2 <- factor(x$cell.type2, levels = unified_levels)
Idents(x) <- x$cell.type2

# Define standardized species order
x$species <- factor(x$species, levels = c("MUS", "CSP", "MFU", "RSI"))

# ==============================================================================
# 5. GLOBAL VISUALIZATIONS
# ==============================================================================
# Setup global color palette matching the unified levels
nClust <- length(levels(x$cell.type2))
my_colors <- colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "Paired"))(nClust)
names(my_colors) <- levels(x$cell.type2)
species_colors <- c("MUS" = "#727171", "CSP" = "#0000ff", "MFU" = "#ff3828", "RSI" = "#2ca02c")

# Export Global UMAP (Colored by Cell Type)
png("bat-4sp-harmony-ct2-umap.png", width = 2000, height = 2000, bg = "transparent", res = 600)
DimPlot(x, reduction = "umap", group.by = "cell.type2", cols = my_colors, raster = TRUE) +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        panel.grid = element_blank(), legend.position = "none")
dev.off()

# Export Global UMAP (Group by Species and tissue)
png("bat-4sp-harmony-ct2-umap-SP.png", width = 8000, height = 2000, bg = "transparent", res = 600)
DimPlot(x, reduction = "umap", group.by = "species", cols = species_colors, raster = TRUE) +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        panel.grid = element_blank(), legend.position = "none")
dev.off()

tissue_colors <- c("brain" = "#8491b4", "heart" = "#00a087", "lung" = "#3c5488", 
                   "liver" = "#4dbbd5", "kidney" = "#ef7b9b", "spleen" = "#f39b7f")
png("bat-4sp-harmony-ct2-umap-Tissue.png", width = 2000, height = 2000, bg = "transparent", res = 600)
DimPlot(x, reduction = "umap", group.by = "tissue", cols = tissue_colors, raster = TRUE) +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        panel.grid = element_blank(), legend.position = "none")
dev.off()

# Export distribution tables
write.csv(table(x$cell.type2, x$spts), "cell-dis-4sp-tissue.csv")
write.csv(table(x$cell.type2, x$species), "cell-dis-4sp-species.csv")

# ==============================================================================
# 6. CROSS-SPECIES SIMILARITY (METANEIGHBOR)
# ==============================================================================
# Convert to SingleCellExperiment for MetaNeighbor
DefaultAssay(x) <- "RNA"
x.sc <- as.SingleCellExperiment(x)

# Identify global highly variable genes
global_hvgs <- variableGenes(dat = x.sc, exp_labels = x.sc$cell.type2)

# Compute AUROC scores across species-tissue identities
aurocs <- MetaNeighborUS(var_genes = global_hvgs,
                         dat = x.sc,
                         study_id = x.sc$spts,
                         cell_type = x.sc$cell.type2,
                         fast_version = TRUE)

write.csv(aurocs, "bat-4sp-harmony-mn-CT2-aurocs.csv")

# Save final annotated global object
qs::qsave(x, "bat-4sp-harmony-annotated.qs")