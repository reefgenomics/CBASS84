library('pheatmap')
##### 1. plot heatmap for top 50/100 genes #####
plot_topn_heatmap <- function(org=c("Spis","Smic"), dds=deseq_data$dds.Spis, n=50){
  cat("Creating Heatmap for :", org, "\n")
  # 1. heatmap
  #dds.Spis <- DESeq(dds.Spis, parallel=F) # run DESeq2
  vsd <- vst(dds, blind=FALSE) # data transformation for visuals
  select <- order(rowMeans(counts(dds,normalized=TRUE)), # select top 100 most diff expressed
                  decreasing=T)[1:n]

  df.annotate <- as.data.frame(colData(dds)[,c("Temp","Reef.Site.Name")])
  df.assay <- assay(vsd)[select,]
  ## do not plot initially to get row ordering and reordering
  phm <- pheatmap(df.assay, scale = "none", cluster_rows=T, silent = T)
  phmr <- phm$tree_row$order
  df.order <- df.assay[rev(phmr),rownames(df.annotate[order(factor(df.annotate$Reef.Site.Name, levels=c("Protected","Exposed","AF","ICN")),df.annotate$Temp),])]

  heatmap.plot <- pheatmap(df.order, scale = "none", cluster_rows=F, border_color = NA, show_rownames=TRUE, # heatmap
           cluster_cols=F, annotation_col=df.annotate, main = paste(org, " Top ",n, "DEGs", sep = ""), fontsize_col = 7, fontsize_row = 6, color = viridis(20) )

}
