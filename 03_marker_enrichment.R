library(clusterProfiler)
library(dplyr)
library(stringr)

# 1. Parse function for local Enrichr databases
parse_enrichr_txt <- function(file_path) {
  lines <- readLines(file_path)
  
  t2g_list <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    parts <- parts[parts != "" & !is.na(parts)]
    
    if (length(parts) > 1) {
      term <- parts[1]
      genes <- parts[2:length(parts)]
      return(data.frame(term = term, gene = genes, stringsAsFactors = FALSE))
    } else {
      return(NULL)
    }
  })
  
  return(bind_rows(t2g_list))
}

# 2. Load the databases
db_files <- c("CellMarker_2024.txt", "Tabula_Muris.txt", "Tabula_Sapiens.txt")
message("Parsing databases...")
databases <- lapply(db_files, parse_enrichr_txt)
names(databases) <- gsub(".txt", "", db_files)

# 3. Find all tissue marker files
marker_files <- list.files(pattern = "_filtered_markers.csv$")

# 4. Initialize a master list to store ALL results
master_results <- list()

# 5. Batch Enrichment Loop
for (mf in marker_files) {
  tissue_name <- gsub("_filtered_markers.csv", "", mf)
  message(paste("Processing tissue:", tissue_name))
  
  markers <- read.csv(mf)
  markers <- markers %>% filter(p_val_adj < 0.05 & avg_log2FC > 0)
  clusters <- unique(markers$cluster)
  
  for (db_name in names(databases)) {
    t2g <- databases[[db_name]]
    
    for (clust in clusters) {
      cluster_genes <- markers %>% filter(cluster == clust) %>% pull(gene)
      cluster_genes <- toupper(cluster_genes) 
      
      # Run enrichment
      res <- enricher(gene = cluster_genes, 
                      TERM2GENE = t2g, 
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2)
      
      if (!is.null(res) && nrow(res@result) > 0) {
        res_df <- res@result %>% filter(p.adjust < 0.05)
        
        if(nrow(res_df) > 0) {
          # Add metadata columns to track where this came from
          res_df$Tissue <- tissue_name
          res_df$Database <- db_name
          res_df$Cluster <- clust 
          
          # Append to our master list
          master_results[[length(master_results) + 1]] <- res_df
        }
      }
    }
  }
}

# 6. Combine all results and write to a single Master CSV
if (length(master_results) > 0) {
  final_master_df <- bind_rows(master_results)
  
  # Reorder columns to put Tissue, Database, and Cluster at the very beginning
  final_master_df <- final_master_df %>% 
    relocate(Tissue, Database, Cluster, .before = ID)
  
  write.csv(final_master_df, "All_Tissues_Master_Enrichment.csv", row.names = FALSE)
  message("\nSuccess! All results saved to 'All_Tissues_Master_Enrichment.csv'")
} else {
  message("\nNo significant enrichment found across any tissues.")
}

####keep top5
library(dplyr)

# 1. Load the huge CSV you just generated
df <- read.csv("All_Tissues_Master_Enrichment.csv")

# 2. Group by the categories and keep only the top 5 smallest p-values
df_top5 <- df %>%
  group_by(Tissue, Database, Cluster) %>%
  arrange(p.adjust) %>%           # Sort by significance
  slice_head(n = 5) %>%           # Keep strictly the top 5
  ungroup()

# 3. Save the clean, trimmed version
write.csv(df_top5, "All_Tissues_Top5_Cleaned.csv", row.names = FALSE)

message("Done! Your CSV is now trimmed to exactly 5 terms per cluster.")
