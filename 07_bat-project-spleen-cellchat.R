# ==============================================================================
# 1. SETUP AND DATA PREPARATION
# ==============================================================================
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(forcats)
library(CellChat)
library(ComplexHeatmap)

options(future.globals.maxSize = 300000 * 1024^2)

# Load and subset Seurat object
x <- qs::qread("bat-spleen-cca.qs")
Idents(x) <- x$cell.type2
x <- subset(x, idents = c("B", "T-1", "NKT", "T-2", "MAC-1", "MAC-2", "DC", "Plasma", "Monocyte", "pDC", "cDC"))
levels(x) <- c("B", "Plasma", "T-1", "T-2", "NKT", "MAC-1", "MAC-2", "DC", "pDC", "cDC", "Monocyte")

species_list <- c("MFU", "RSI", "MUS", "CSP")

# ==============================================================================
# 2. CREATE CELLCHAT OBJECTS FOR EACH SPECIES (Skip if already generated)
# ==============================================================================
for(species_use in species_list) {
  y <- subset(x, subset = species == species_use)
  data.input <- GetAssayData(y, assay = "RNA", slot = "data")
  labels <- Idents(y)
  meta <- data.frame(group = labels, row.names = names(labels)) 
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  cellchat@DB <- CellChatDB.human # using human DB for bat/mouse homology
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  saveRDS(cellchat, paste0("cellchat-spleen-", species_use, ".rds"))
}

# Load processed objects
cellchat.m <- readRDS("cellchat-spleen-MUS.rds")
cellchat.c <- readRDS("cellchat-spleen-CSP.rds")
cellchat.f <- readRDS("cellchat-spleen-MFU.rds")
cellchat.r <- readRDS("cellchat-spleen-RSI.rds")

# Merge objects for cross-species comparison
object.list <- list(MUS = cellchat.m, CSP = cellchat.c, MFU = cellchat.f, RSI = cellchat.r)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
cellchat_merged <- liftCellChat(cellchat_merged, group.new = levels(cellchat_merged@idents))

# ==============================================================================
# 3. FIGURE 7G: GLOBAL SIGNALING ROLE SCATTER PLOTS (Macrophage Hubs)
# ==============================================================================
gg <- list()
species_labels <- c("MUS", "CSP", "MFU", "RSI")
objects <- list(cellchat.m, cellchat.c, cellchat.f, cellchat.r)

for (i in seq_along(objects)) {
  obj <- objects[[i]]
  num.link <- rowSums(obj@net$count) + colSums(obj@net$count) - diag(obj@net$count)
  weight.MinMax <- c(min(num.link), max(num.link))
  gg[[i]] <- netAnalysis_signalingRole_scatter(obj, 
                                               title = species_labels[i], 
                                               weight.MinMax = weight.MinMax)
}

# Display Fig 7g
patchwork::wrap_plots(plots = gg)

# ---> EXPORT SOURCE DATA FOR FIG 7G <---
fig7g_data <- lapply(seq_along(gg), function(i) {
  plot_data <- gg[[i]]$data
  plot_data$Species <- species_labels[i]
  return(plot_data)
}) %>% bind_rows() %>%
  rename(Outgoing_Signaling = x, Incoming_Signaling = y)

write.csv(fig7g_data, "SourceData_Fig7g_Global_Signaling_Roles.csv", row.names = FALSE)

# ==============================================================================
# 4. FIGURE 7H: CONSISTENTLY REGULATED PATHWAYS & VCAM ELEVATION
# ==============================================================================
# Generate pairwise comparisons vs MUS
gg1 <- rankNet(cellchat_merged, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,2)) # MUS vs CSP
gg2 <- rankNet(cellchat_merged, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,3)) # MUS vs MFU
gg3 <- rankNet(cellchat_merged, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,4)) # MUS vs RSI

# Helper function to extract robust pathways
get_robust_directions <- function(df, p_val_thresh = 0.05, min_flow_thresh = 0.005) {
  df %>%
    group_by(name) %>%
    filter(pvalues < p_val_thresh & max(contribution) >= min_flow_thresh) %>%
    summarize(
      enriched_in = as.character(group[which.max(contribution)]),
      direction = ifelse(enriched_in == "MUS", "Down", "Up"),
      .groups = "drop"
    )
}

# Apply thresholds and find consistent overlaps across all three bat species
sig_csp <- get_robust_directions(gg1$data)
sig_mfu <- get_robust_directions(gg2$data)
sig_rsi <- get_robust_directions(gg3$data)

consistent_overlap <- sig_csp %>%
  inner_join(sig_mfu, by = "name", suffix = c("_csp", "_mfu")) %>%
  inner_join(sig_rsi, by = "name") %>%
  rename(direction_rsi = direction) %>%
  filter(direction_csp == direction_mfu & direction_mfu == direction_rsi)

overlap_pathways <- consistent_overlap$name

# Plot 4-group Information Flow (Fig 7h)
if(length(overlap_pathways) > 0) {
  gg_overlap_4groups <- rankNet(cellchat_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE, 
                                signaling = overlap_pathways, comparison = c(1, 2, 3, 4))
  data_4_groups <- gg_overlap_4groups$data
  data_4_groups$group <- factor(data_4_groups$group, levels = c("RSI", "MFU", "CSP", "MUS"))
  
  fig_7h <- ggplot(data_4_groups, aes(x = contribution, y = fct_rev(reorder(name, contribution)), fill = group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("MUS" = "#727171", "CSP" = "#0000ff", "MFU" = "#ff3828", "RSI" = "#2ca02c")) +
    theme_minimal() +
    labs(title = "Information Flow of Consistently Regulated Pathways", x = "Information Flow (Contribution)", y = "Signaling Pathway", fill = "Species") +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(axis.text.y = element_text(size = 11, face = "bold"), panel.grid.major.y = element_blank())
  
  print(fig_7h)
  
  # ---> EXPORT SOURCE DATA FOR FIG 7H <---
  write.csv(data_4_groups, "SourceData_Fig7h_InformationFlow_ConservedPathways.csv", row.names = FALSE)
}

# ==============================================================================
# 5. FIGURE 7I: VCAM SENDER-RECEIVER INTERACTION HEATMAPS
# ==============================================================================
pathways.show <- "VCAM"

# Ensure names are correctly applied to the object list for titles
names(object.list) <- c("MUS", "CSP", "MFU", "RSI")
ht <- list()

# Generate heatmaps for each species to show macrophage self-regulation
for (i in seq_along(object.list)) {
  obj <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
  
  ht[[i]] <- netVisual_heatmap(obj, 
                               signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "-", names(object.list)[i])) 
}

# Display Fig 7i side-by-side
ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]] + ht[[4]], ht_gap = unit(0.5, "cm"))

# ---> EXPORT SOURCE DATA FOR FIG 7I <---
# Extract the underlying communication probability matrices for the heatmaps
fig7i_data <- lapply(seq_along(object.list), function(i) {
  obj <- object.list[[i]]
  # Check if VCAM exists in this species to prevent errors
  if (pathways.show %in% obj@netP$pathways) {
    # Extract the 2D matrix (Sender x Receiver)
    prob_mat <- obj@netP$prob[, , pathways.show]
    # Melt into a flat dataframe
    df <- as.data.frame(as.table(prob_mat))
    colnames(df) <- c("Sender", "Receiver", "Communication_Probability")
    df$Species <- names(object.list)[i]
    return(df)
  }
}) %>% bind_rows() %>% 
  filter(Communication_Probability > 0) # Optionally remove zero-links to keep file size clean

write.csv(fig7i_data, paste0("SourceData_Fig7i_", pathways.show, "_Sender_Receiver_Matrix.csv"), row.names = FALSE)