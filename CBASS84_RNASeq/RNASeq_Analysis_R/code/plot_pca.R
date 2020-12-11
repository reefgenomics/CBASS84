library('ggplot2')
##### 2. plot_pca: CA Plot #####
plot_pca <- function(org=c("Spis","Smic"),vsd=deseq_data$vsd.Spis){
  cat("Creating PCA plot for :", org, "\n")
  vsd <- vsd
  # remove outlier ES1-36A and create PCA object
  pcaData <- plotPCA(vsd[,grep("ES1-36A", colnames(vsd), invert = T)], intgroup=c("Temp","Reef.Site.Name","condition"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  q <- ggplot(pcaData, aes(PC1, PC2, fill=Temp, shape=Reef.Site.Name) ) +
    geom_point(size=6) +
    scale_fill_manual(values=c("#ffeda0", "#feb24c", "#f03b20")) +  #for 30,33,36 use c("#ffeda0", "#feb24c", "#f03b20")
    scale_colour_manual(values=c("Black", "Black", "Black")) +
    scale_shape_manual(values=c(24,22,21,23)) + #for c(Eilat, KAUST, Exposed, Protected) use shapes or 21,24, 22, 23
    #geom_text(aes(label=condition)) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    theme_light() + coord_fixed() +
    theme(text = element_text(family = "Helvetica"),
          axis.text.y = element_text(angle = 0, size=14),
          axis.text.x = element_text(angle = 0, size=14),
          axis.title = element_text(size= 16)
    )
  print(q)
}
