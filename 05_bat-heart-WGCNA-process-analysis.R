library(dplyr)
library(Seurat)
library(hdWGCNA)
library(UCell)
library(patchwork)
library(qs)

options(future.globals.maxSize = 300000 * 1024^2)

# 1. Load Data
x <- qs::qread('bat-heart-cca.qs')
DefaultAssay(x) <- "RNA"
Idents(x) <- x$cell.type2
x$species <- factor(x$species, levels = c("MUS", "CSP", "MFU", "RSI"))

# ===================================================================
# LINEAGE 1: MYO (Cardiomyocytes)
# ===================================================================
x <- SetupForWGCNA(x, gene_select = "fraction", fraction = 0.05, wgcna_name = "myo")
x <- MetacellsByGroups(x, group.by = c("cell.type2", "species"), reduction = 'pca', k = 25, max_shared = 10, ident.group = 'cell.type2')
x <- NormalizeMetacells(x)
x <- SetDatExpr(x, group_name = c("Cardiomyocyte"), group.by = 'cell.type2', assay = 'RNA', layer = 'data')

# 1A. Test and Plot Soft Powers
x <- TestSoftPowers(x, networkType = 'signed')
plot_list_myo <- PlotSoftPowers(x)

pdf("bat-heart-wgcna-myo-soft-power.pdf")
print(wrap_plots(plot_list_myo, ncol = 2))
dev.off()

# 1B. Construct Network (hdWGCNA auto-selects soft power when argument is omitted)
x <- ConstructNetwork(x, tom_name = 'heart_myo', overwrite_tom = TRUE)

pdf("bat-heart-wgcna-myo-dendro.pdf")
PlotDendrogram(x, main='bat heart myo hdWGCNA Dendrogram')
dev.off()

x <- ScaleData(x, features = VariableFeatures(x))
x <- ModuleEigengenes(x, group.by.vars = "species")
x <- ModuleConnectivity(x, group.by = 'cell.type2', group_name = c("Cardiomyocyte"))
x <- ResetModuleNames(x, new_name = "myo-M")
x <- ModuleExprScore(x, n_genes = 25, method = 'UCell')

# 1C. Calculate DMEs for Myo
x <- SetActiveWGCNA(x, 'myo')
desired_cell_types_myo <- c("Cardiomyocyte")
group1_myo <- rownames(subset(x@meta.data, cell.type2 %in% desired_cell_types_myo & species == "MUS"))
group2_myo <- rownames(subset(x@meta.data, cell.type2 %in% desired_cell_types_myo & species == "RSI"))
group3_myo <- rownames(subset(x@meta.data, cell.type2 %in% desired_cell_types_myo & species == "MFU"))
group4_myo <- rownames(subset(x@meta.data, cell.type2 %in% desired_cell_types_myo & species == "CSP"))

dme_rsi_myo <- FindDMEs(x, barcodes1 = group1_myo, barcodes2 = group2_myo, test.use = 'wilcox', wgcna_name = 'myo') %>% mutate(species_comp = "Mouse vs RSI")
dme_mfu_myo <- FindDMEs(x, barcodes1 = group1_myo, barcodes2 = group3_myo, test.use = 'wilcox', wgcna_name = 'myo') %>% mutate(species_comp = "Mouse vs MFU")
dme_csp_myo <- FindDMEs(x, barcodes1 = group1_myo, barcodes2 = group4_myo, test.use = 'wilcox', wgcna_name = 'myo') %>% mutate(species_comp = "Mouse vs CSP")
combined_DMEs_myo <- rbind(dme_rsi_myo, dme_mfu_myo, dme_csp_myo)
write.csv(combined_DMEs_myo, "bat-heart-myo-DMEs.csv", row.names = FALSE)

# ===================================================================
# 3. Save Fully Analyzed Object
# ===================================================================
qs::qsave(x, 'bat-heart-hdWGCNA-analyzed.qs')