##### create_DEG_matrix #####
create_DEG_matrix <- function(p.value=0.005,org=c("Spis","Smic")){
  cat("Creating DEG Matrix for :", org, "\n")
  # counts results_dir
  results_dir= paste("output/", org,"/DEG_tables_svalues/", sep = "")
  # get sample names from all directories for Stylophora
  sample_id <- conditions
  conditions_clean <- vector()
  for(sample in sample_id){
    cond0 <- strsplit(paste(sample, sep=""), split = "condition_")[[1]][2]
    cond1 <- strsplit(cond0, split = "_vs_")[[1]][1]
    cond2 <- sub(".txt", "", strsplit(cond0, split = "_vs_")[[1]][2])

    conditions_clean <- append(conditions_clean, cond1)
    conditions_clean <- append(conditions_clean, cond2)
  }
  conditions_clean <- unique(conditions_clean)

  DEG.matrix <- matrix(nrow = length(conditions_clean), ncol = length(conditions_clean), dimnames = list(c(conditions_clean), c(conditions_clean)))

  updown <- data.frame()
  deg_list <- list()
  for(sample in sample_id){
    cat("\t Processing ", sample, "\n")
    cond0 <- strsplit(paste(sample, sep=""), split = "condition_")[[1]][2]
    cond1 <- strsplit(cond0, split = "_vs_")[[1]][1]
    cond2 <- sub(".txt", "", strsplit(cond0, split = "_vs_")[[1]][2])

    df <- read.table(paste(results_dir, sample, ".txt",sep=""), sep="\t", header = T)
    deg <- dim(df[which(abs(df$log2FoldChange) > 2 & df$svalue < p.value),])[1]
    up <- dim(df[which(df$log2FoldChange > 2 & df$svalue < p.value),])[1]
    up_genes <- rownames(df[which(df$log2FoldChange > 2 & df$svalue < p.value),][1])
    down <- dim(df[which(df$log2FoldChange < (-2) & df$svalue < p.value),])[1]
    down_genes <- rownames(df[which(df$log2FoldChange < (-2) & df$svalue < p.value),][1])

    alter = paste("condition_", cond2, "_vs_", cond1, ".txt", sep="")
    # create names for the lost for listing diferentially expressed genes.
    list_name_up = paste(cond1, "_vs_", cond2, "_UP", sep="")
    list_name_down = paste(cond1, "_vs_", cond2, "_DOWN", sep="")

    if(length(up_genes) > 0){
      deg_list[[list_name_up]] <- up_genes
    }
    if(length(down_genes) > 0){
      deg_list[[list_name_down]] <- down_genes
    }

    updown <- rbind(updown, cbind(comparison=sample, upregulated=up, downregulated=down, cond1=cond1, cond2=cond2, alt=alter))
    DEG.matrix[cond1,cond2] <- deg
  }
  return(list(matrix=DEG.matrix, updown=updown,  deg_genes=deg_list))
}
