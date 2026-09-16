library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(hdWGCNA)
library(stringr)
library(scales)
library(tidyr)
library(qs)

# 1. Load Data
x <- qs::qread("bat-heart-hdWGCNA-analyzed.qs")
species_order <- c("MUS", "CSP", "MFU", "RSI")
x$species <- factor(x$species, levels = species_order)
spp_colors <- c("MUS" = "#727171", "CSP" = "#0000ff", "MFU" = "#ff3828", "RSI" = "#2ca02c")

# Define target lineages
combinations <- list(
  "myo" = list(
    wgcna_name = "myo",
    modules = c("myo-M4"), # Add more modules here if needed
    cell_types = c("Cardiomyocyte"),
    dme_file = "bat-heart-myo-DMEs.csv"
  )
)

# 2. Execute Plotting Loop
for (combo_name in names(combinations)) {
  
  cat("=========================================\n")
  cat("Plotting combination:", combo_name, "\n")
  
  current_wgcna <- combinations[[combo_name]]$wgcna_name
  target_CTs <- combinations[[combo_name]]$cell_types
  selected_modules <- combinations[[combo_name]]$modules
  
  # Set WGCNA active and get MEs
  x <- SetActiveWGCNA(x, current_wgcna)
  MEs <- GetMEs(x, harmonized = TRUE)
  mods <- levels(GetModules(x)$module); mods <- mods[mods != 'grey']
  
  # --- A. Global Module Plots (KME, DotPlot, VlnPlot) ---
  pdf(sprintf("bat-heart-wgcna-%s-kme.pdf", combo_name))
  print(PlotKMEs(x, ncol = 3))
  dev.off()
  
  # Add MEs to metadata temporarily for Seurat plotting
  x@meta.data <- cbind(x@meta.data, MEs)
  
  pdf(sprintf("bat-heart-wgcna-%s-module-dot-split-sp.pdf", combo_name), width = 12, height = 12)
  print(DotPlot(x, features = mods, group.by = 'cell.type2', split.by = "species", cols = "RdYlBu") + 
          RotatedAxis() + scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue'))
  dev.off()
  
  pdf(sprintf("bat-heart-wgcna-%s-module-vln-split-sp.pdf", combo_name), width = 12, height = 12)
  print(VlnPlot(x, features = mods, group.by = 'cell.type2', split.by = "species", stack = TRUE))
  dev.off()
  
  # Clean up metadata so columns don't clash on next loop
  x@meta.data <- x@meta.data[, !colnames(x@meta.data) %in% colnames(MEs)]
  
  # --- B. RAW ME BAR Plot ---
  # Define all modules (excluding grey) right after you get MEs
  all_network_modules <- colnames(MEs)
  all_network_modules <- all_network_modules[!grepl("grey", all_network_modules, ignore.case = TRUE)]
  
  # --- C. RAW ME BAR Plot (Now plots ALL modules like the Liver/Brain) ---
  cat("  -> Generating RAW Barplots for ALL modules...\n")
  cells_found <- sum(x$cell_type_merged %in% target_CTs)
  
  if(length(all_network_modules) > 0 & cells_found > 0) {
    plot_data_long_raw <- MEs %>%
      select(all_of(all_network_modules)) %>%
      bind_cols(x@meta.data %>% select(species, cell_type = cell_type_merged)) %>% 
      filter(cell_type %in% target_CTs) %>%
      mutate(cell_type = as.character(cell_type)) %>% 
      group_by(cell_type, species) %>%
      summarise(across(all_of(all_network_modules), \(val) mean(val, na.rm = TRUE)), .groups = 'drop') %>%
      pivot_longer(cols = all_of(all_network_modules), names_to = "Module", values_to = "Raw_Expression") %>%
      mutate(Module = str_replace(Module, sprintf("^%s-", current_wgcna), ""), 
             cell_type = str_replace_all(cell_type, "_", " "))
    
    # Enforce strict ordering
    plot_data_long_raw$species <- factor(plot_data_long_raw$species, levels = species_order)
    plot_data_long_raw$cell_type <- factor(plot_data_long_raw$cell_type, levels = str_replace_all(target_CTs, "_", " "))
    
    # Plotting using Raw Expression and advanced theme
    p_bar_raw <- ggplot(plot_data_long_raw, aes(x = cell_type, y = Raw_Expression, fill = species)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9, color = "black", linewidth = 0.2) +
      scale_fill_manual(values = spp_colors) + 
      facet_wrap(~ Module, ncol = 3, scales = "free_y") +
      theme_minimal() +
      theme(
        text = element_text(family = "sans"),
        axis.title.x = element_blank(),          
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
        axis.text.y = element_text(size = 10, color = "black"),          
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(size = 12, face = "bold", color = "black", 
                                   angle = 45, hjust = 1, vjust = 1), 
        strip.background = element_rect(fill = "gray20", color = "black"),
        strip.text = element_text(color = "white", face = "bold", size = 14, margin = margin(t=6, b=6)),  
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 12),
        panel.grid.major.y = element_line(color = "gray80", linetype = "dashed"), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), 
        axis.line.x = element_line(color = "grey30", linewidth = 1) 
      ) + 
      labs(y = "Raw Eigengene Expression", fill = "Species")
    
    ggsave(sprintf("bat-heart-Barplot_%s_RAW_Global.pdf", combo_name), plot = p_bar_raw, width = 12, height = 10)
    write.csv(plot_data_long_raw, sprintf("bat-heart-plotDF-%s_RAW_Global.csv", combo_name), row.names = FALSE)
  }
  
  # --- C. Custom SCALED Barplot ---
  me_cols <- colnames(MEs)[colnames(MEs) %in% selected_modules]
  if(length(me_cols) > 0) {
    plot_data_long <- MEs %>%
      select(all_of(me_cols)) %>%
      # 🌟 Pull from the merged column
      bind_cols(x@meta.data %>% select(species, cell_type = cell_type_merged)) %>% 
      filter(cell_type %in% target_CTs) %>%
      mutate(cell_type = as.character(cell_type)) %>% 
      group_by(cell_type, species) %>%
      summarise(across(all_of(me_cols), \(val) mean(val, na.rm = TRUE)), .groups = 'drop') %>%
      pivot_longer(cols = all_of(me_cols), names_to = "Module", values_to = "Raw_Expression") %>%
      
      group_by(Module) %>% 
      mutate(Scaled_Expression = rescale(Raw_Expression, to = c(0.1, 1))) %>% 
      ungroup() %>%
      
      mutate(Module = str_replace(Module, sprintf("^%s-", current_wgcna), ""), 
             cell_type = str_replace_all(cell_type, "_", " "))
    
    plot_data_long$species <- factor(plot_data_long$species, levels = species_order)
    
    p_bar <- ggplot(plot_data_long, aes(x = cell_type, y = Scaled_Expression, fill = species)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9, color = "black", linewidth = 0.2) +
      scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0), expand = c(0, 0)) +
      scale_fill_manual(values = spp_colors) + facet_wrap(~ Module, ncol = 3, scales = "free_x") +
      theme_minimal() + theme(legend.position = "bottom") + labs(y = "Scaled Eigengene Expression", x = NULL)
    
    ggsave(sprintf("bat-heart-Barplot_%s_SCALED.pdf", combo_name), plot = p_bar, width = 12, height = 6)
    write.csv(plot_data_long, sprintf("bat-heart-plotDF-%s_SCALED.csv", combo_name), row.names = FALSE) 
  }
  
  # --- D. Hub Gene Networks & REVERSED DotPlots ---
  tom_path <- sprintf("path_to_your_WGCNA_analysis/TOM/heart_%s_TOM.rda", current_wgcna)
  if (file.exists(tom_path)) {
    x@misc[[current_wgcna]]$wgcna_net$TOMFiles <- tom_path
    pdf(sprintf("./ModuleNetworks/bat-heart-wgcna-%s-hubs-network.pdf", combo_name), width = 12)
    HubGeneNetworkPlot(x, n_hubs = 10, n_other = 20, edge_prop = 0.75, mods = selected_modules)
    dev.off()
  }
  
  Idents(x) <- "cell.type2" 
  x_sub <- subset(x, idents = target_CTs)
  Idents(x_sub) <- "species"
  hub_df <- GetHubGenes(x_sub, n_hubs = 30)
  
  for (mod in selected_modules) {
    target_hubs <- hub_df %>% filter(module == mod)
    if(nrow(target_hubs) == 0) next 
    
    top_10_genes <- target_hubs %>%
      arrange(desc(kME)) %>%
      slice_head(n = 10) %>%
      pull(gene_name)
      
    # CRITICAL: Reverse order so the top hub gene appears at the top of the coord_flip() plot
    top_10_genes <- rev(top_10_genes)
    
    p_dot <- DotPlot(x_sub, features = top_10_genes, group.by = "species", dot.scale = 8, cols = "RdYlBu") + 
      coord_flip() + theme_bw() + labs(title = sprintf("Top 10 Hub Genes: %s", mod), x = "Gene", y = "Species") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(sprintf("./ModuleNetworks/bat-heart-%s-%s-hubs-dot.pdf", combo_name, mod), plot = p_dot, width = 6, height = 5)
  }
}