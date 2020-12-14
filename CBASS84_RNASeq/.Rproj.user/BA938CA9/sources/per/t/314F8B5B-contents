library("calibrate")
library("reshape2")
library("ggplot2")
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library('dplyr')
library("BiocParallel")
library('factoextra')
library('topGO')
library('tictoc')
library('gridExtra')
library('cowplot')
library('viridis')
### Inititate parallel if needed
register(MulticoreParam(2, exportglobals = F))

##### Prepare data and meta-data #####

prepare_DESeq2_data <- function(wd="/Volumes/omics4tb2/Collaborations/Vulcan/R_Projects/Stylophora/",
                                counts_dir= "../../data/181203_K00235_0152_BHWNJ7BBXX_CBASS84_RNASeq/all_htseqcounts",
                                data_dir= "../../data/181203_K00235_0152_BHWNJ7BBXX_CBASS84_RNASeq/") {
  setwd(dir = wd)
  # counts directory
  directory <- counts_dir
  
  # get sample names from all directories for Stylophora
  sample_id <- dir(file.path(data_dir))
  
  # remove folders for samples with issues and other script folders also remove ES1-36A since it seems to have very high number of reads across samples.
  sample_id <- sample_id[!(sample_id %in% c("Lanes_metadata", "ES1-36A","S6-P-ST-36B","KS2-33A","ES3-33A","all_htseqcounts", "label_changes_README.txt"))]
  
  # create htsqcounts file names based on sample id
  sample_files <- unlist(lapply(sample_id, function(i) paste(i, "_htseqcounts.txt",sep="")))
  
  # get metadata info from metadata worksheet and select relevant columns and delete NA rows
  meta_data <- read.delim("../CBASS84_ExperimentalDesign_SampleIDs.txt", header=T, sep="\t", stringsAsFactors = F)
  meta_data <- dplyr::select(meta_data, sample = Sample.Name.CBASS84, 1:20)
  #meta_data <- meta_data[1:84,]
  
  # remove samples from metadata that were earlier filtered and order
  meta_data <- meta_data[grep("(S6_P_ST_36B|KS2_33A|ES3_33A|ES6_30A|ES1_36A)", meta_data$sample, invert=T),]
  meta_data <- meta_data[order(meta_data$sample),]
  meta_data$Temp <- as.factor(meta_data$Temp)
  
  # create a metadata column by combining temperature and Reef.site.name
  meta_data$condition <- paste(meta_data$Temp, meta_data$Reef.Site.Name, sep="_")
  
  #rebuild_sample table
  sampleTable <- cbind(sampleName=sample_id, fileName=sample_files, meta_data)
  
  # Build DSEq dataset
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~ condition, directory=directory )
  
  # Seperate Spis vs Spis
  dds.Spis <- ddsHTSeq[grep("Spis", rownames(ddsHTSeq)),]
  dds.Spis$condition <- relevel(dds.Spis$condition, ref = "30_KAUST")
  dds.Spis <- DESeq(dds.Spis)
  vsd.Spis <- vst(dds.Spis, blind = F)
  return(list(meta_data=meta_data, dds = dds.Spis, vsd= vsd.Spis, sample_table= sampleTable, directory=directory))
}

##### prepare DESEq data for KAUST only (Optional) #####
prepare_DESeq2_data_kaust <- function(wd="/Volumes/omics4tb2/Collaborations/Vulcan/R_Projects/Stylophora/",
                                counts_dir= "../../data/181203_K00235_0152_BHWNJ7BBXX_CBASS84_RNASeq/all_htseqcounts",
                                data_dir= "../../data/181203_K00235_0152_BHWNJ7BBXX_CBASS84_RNASeq/") {
  setwd(dir = wd)
  # counts directory
  directory <- counts_dir
  
  # get sample names from all directories for Stylophora
  sample_id <- dir(file.path(data_dir))
  
  # remove folders for samples with issues and other script folders also remove ES1-36A since it seems to have very high number of reads across samples.
  sample_id <- sample_id[!(sample_id %in% c("Lanes_metadata", "ES1-36A","S6-P-ST-36B","KS2-33A","ES3-33A","all_htseqcounts", "label_changes_README.txt"))]
  
  # create htsqcounts file names based on sample id
  sample_files <- unlist(lapply(sample_id, function(i) paste(i, "_htseqcounts.txt",sep="")))
  
  # get metadata info from metadata worksheet and select relevant columns and delete NA rows
  meta_data <- read.delim("../CBASS84_ExperimentalDesign_SampleIDs.txt", header=T, sep="\t", stringsAsFactors = F)
  meta_data <- dplyr::select(meta_data, sample = Sample.Name.CBASS84, 1:20)
  #meta_data <- meta_data[1:84,]
  
  # remove samples from metadata that were earlier filtered and order
  meta_data <- meta_data[grep("(S6_P_ST_36B|KS2_33A|ES3_33A|ES6_30A|ES1_36A)", meta_data$sample, invert=T),]
  meta_data <- meta_data[order(meta_data$sample),]
  meta_data$Temp <- as.factor(meta_data$Temp)
  
  # create a metadata column by combining temperature and Reef.site.name
  meta_data$condition <- paste(meta_data$Temp, meta_data$Reef.Site.Name, sep="_")
  
  #rebuild_sample table
  sampleTable <- cbind(sampleName=sample_id, fileName=sample_files, meta_data)
  
  # Build DSEq dataset
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~ condition, directory=directory )
  
  # Seperate Spis vs Spis
  dds.Spis <- ddsHTSeq[grep("Spis", rownames(ddsHTSeq)),]
  dds.Spis <- dds.Spis[,grep("Eilat", dds.Spis$Reef.Site.Name, invert=T)]
  dds.Spis$condition <- droplevels(dds.Spis$condition)
  dds.Spis$condition <- relevel(dds.Spis$condition, ref = "30_KAUST")
  dds.Spis <- DESeq(dds.Spis)
  vsd.Spis <- vst(dds.Spis, blind = F)
  return(list(meta_data=meta_data, dds = dds.Spis, vsd= vsd.Spis, sample_table= sampleTable, directory=directory))
}


##### create_conditions:  for comparison #####
create_conditions <- function(meta_data = deseq_data$meta_data){
  conds <- unique(meta_data$condition)
  combinations <- combn(conds, 2, simplify = F)
  conditions <- vector()
  count = 0
  for(comb in 1:length(combinations)){
    tmp1 <- strsplit(combinations[[comb]][1], split = "_")[[1]][2]
    tmp1.1 <- combinations[[comb]][1]
    tmp2 <- strsplit(combinations[[comb]][2], split = "_")[[1]][2]
    tmp2.1 <- combinations[[comb]][2]
    
    if(tmp1 == "KAUST" & tmp2 %in% c("Eilat", "Protected","KAUST", "Exposed")){
      cat(tmp2, tmp1, "1\n")
      name = paste("condition_", tmp2.1, "_vs_", tmp1.1, sep="")
    } else if(tmp2 == "KAUST" & tmp1 %in% c("Eilat", "Protected","KAUST", "Exposed")){
      cat(tmp1, tmp2, "2\n")
      name = paste("condition_", tmp1.1, "_vs_", tmp2.1, sep="")
    } else if(tmp1 == "Exposed" & tmp2 %in% c("Protected", "Eilat")){
      cat(tmp2, tmp1, "3.1\n")
      name = paste("condition_", tmp2.1, "_vs_", tmp1.1, sep="")
    } else if(tmp2 == "Exposed" & tmp1 %in% c("Protected", "Eilat", "Exposed")){
      cat(tmp1, tmp2, "3.2\n")
      name = paste("condition_", tmp1.1, "_vs_", tmp2.1, sep="")
    } else if(tmp2 == "Protected" & tmp1 %in% c("Eilat", "Protected")){
      cat(tmp1, tmp2, "4\n")
      name = paste("condition_", tmp1.1, "_vs_", tmp2.1, sep="")
    } else if(tmp1 == "Protected" & tmp2 %in% c("Eilat", "Protected")){
      cat(tmp2, tmp1, "4.2\n")
      name = paste("condition_", tmp2.1, "_vs_", tmp1.1, sep="")
    } else if(tmp1 == "Eilat" & tmp2 %in% c("Eilat")){
      cat(tmp2, tmp1, "5.2\n")
      name = paste("condition_", tmp2.1, "_vs_", tmp1.1, sep="")
    } else if(tmp2 == "Eilat" & tmp1 %in% c("Eilat")){
      cat(tmp2, tmp1, "5.2\n")
      name = paste("condition_", tmp2.1, "_vs_", tmp1.1, sep="")
    }
    
    conditions <- append(name, conditions)
    count <- count + 1
  }
  return(conditions)
}



##### make_deg: DESeq analysis for all conditions #####
make_deg <- function(condition=condition, dds=dds.Spis, lfc=2, write2file=T, sampleTable=sampleTable, directory=deseq_data$directory,svalue=0.005){
  dds.Spis <- dds
  # what analysis is being performed
  cat("Running DE analysis for condition:", condition, "\n")
  
  cond0 <- strsplit(paste(condition, sep=""), split = "condition_")[[1]][2]
  cond1 <- strsplit(cond0, split = "_vs_")[[1]][1]
  cond2 <- strsplit(cond0, split = "_vs_")[[1]][2]
  
  cat("Cond1:", cond1, " Cond2:", cond2, "\n")
  
  # relevel to set cond2 as our reference sample
  dds.Spis$condition <- relevel(dds.Spis$condition, ref = cond2)
  
  # Run DEseq2
  dds.Spis <- DESeq(dds.Spis, parallel=F, BPPARAM=MulticoreParam(2))
  
  # get results
  #res <- results(dds.Spis, name = condition, alpha = 0.05)
  # # get results and use lfc shrinkage for visualization
   resLFC <- lfcShrink(dds = dds.Spis,
                      coef = condition,
                      type="apeglm",
                      lfcThreshold = lfc)

  # order results and write into a file
  #res.ordered <- resLFC[order(resLFC$svalue),]
  res.ordered <- resLFC[order(resLFC$svalue),]
  

  # get summary and write into a file
  if(write2file){
    write.table(as.data.frame(res.ordered), 
                file=paste("../RNA-seq_Analysis/Spis_DEG/DEG_tables_svalues/",condition,".txt", sep=""), sep="\t")
  }
  return(as.data.frame(res.ordered))
}


##### make_deg_temp: DESeq analysis for Temp ##### 
make_deg_temp <- function(sampleTable=deseq_data$sample_table, directory=deseq_data$directory, lfc=2, svalue=0.005, write2file=T){
  
  # Build DSEq dataset
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~ Temp, directory=directory )
  
  # Seperate Spis vs Spis
  dds.Spis <- ddsHTSeq[grep("Spis", rownames(ddsHTSeq)),]
  dds.Spis$Temp <- relevel(dds.Spis$Temp, ref = "33")
  dds.Spis <- DESeq(dds.Spis)
  vsd.Spis <- vst(dds.Spis, blind = F)
  
  
  for(condition in resultsNames(dds.Spis)[2:length(resultsNames(dds.Spis))]){
    # open pdf for volcano plot
    pdf(file = paste("../RNA-seq_Analysis/Spis_DEG/volcano_plots_svalues/", condition, ".pdf", sep=""))
    
    # what analysis is being performed
    cat("Running DE analysis for condition:", condition, "\n")
    
    cond0 <- strsplit(paste(condition, sep=""), split = "condition_")[[1]][1]
    cond1 <- strsplit(cond0, split = "_")[[1]][2]
    cond2 <- strsplit(cond0, split = "_")[[1]][4]
    cat("Cond1:", cond1, " Cond2:", cond2, "\n")
    
    # get results and use lfc shrinkage for visualization
    resLFC <- lfcShrink(dds = dds.Spis,
                        coef = condition,
                        type="apeglm",
                        lfcThreshold = lfc)
    
    # order results and write into a file
    res.ordered <- resLFC[order(resLFC$svalue),]
    
    # volcano plot
    resPlot <- as.data.frame(res.ordered)
    resPlot$Gene <- rownames(resPlot)
    lfc <- lfc
    p.value <- svalue
    #sig.genes <- dim(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc))[1]
    sig.genes <- dim(resPlot[which(abs(resPlot$log2FoldChange) > 2 & resPlot$svalue < p.value),])[1]
    outliers <- sum(resPlot$baseMean > 0 & is.na(resPlot$svalue))
    low.counts <- sum(!is.na(resPlot$svalue) & is.na(resPlot$svalue))
    minim <- resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = F)][1] + (resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = F)][1]/100) * 20
    maxim <- resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = T)][1] + (resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = T)][1]/100) * 20
    
    with(resPlot, plot(log2FoldChange, -log10(svalue),
                       pch=20, main= paste(condition),
                       sub=paste( sig.genes, " significant genes (L2FC >", lfc, "& svalue <", p.value, ")"),
                       # " | outliers:", outliers,
                       # "| low counts:", low.counts ),
                       xlim=c(minim,maxim),
                       col="gray", cex.sub=0.8, col.sub="gray", cex.lab=0.8))
    # Add colored points: red if svalue< p.value, orange of log2FC > lfc, green if both)
    with(subset(resPlot, svalue < p.value ), points(log2FoldChange, -log10(svalue), pch=20, col="#2c7bb6"))
    with(subset(resPlot, abs(log2FoldChange) > lfc), points(log2FoldChange, -log10(svalue), pch=20, col="#fdae61"))
    with(subset(resPlot, svalue < p.value & (abs(log2FoldChange) > lfc)), points(log2FoldChange, -log10(svalue), pch=20, col="#d7191c"))
    # Label points with the textxy function from the calibrate plot
    if(sig.genes > 9){
      with(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc)[1:10,], textxy(log2FoldChange, -log10(svalue), labs=Gene, cex=.6, col=rgb(0,0,0, 0.5)))
    }
    if(sig.genes <= 9){
      with(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc), textxy(log2FoldChange, -log10(svalue), labs=Gene, cex=.6, col=rgb(0,0,0, 0.5)))
    }
    abline(v=c(-lfc, lfc), h=c(-log10(svalue),-log10(svalue)), col="gray", lty=2)
    dev.off()
    
    # get summary and write into a file
    if(write2file){
      write.table(as.data.frame(res.ordered), 
                  file=paste("../RNA-seq_Analysis/Spis_DEG/DEG_tables_svalues/",condition,".txt", sep=""), sep="\t")
      
  }
  
   }
  
}


##### make_deg: DESeq analysis for sites #####
make_deg_site <- function(sampleTable=deseq_data$sample_table, directory=deseq_data$directory, lfc=2, svalue=0.005, write2file=T){
  
  # Build DSEq dataset
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~ Reef.Site.Name, directory=directory )
  
  # Seperate Spis vs Spis
  dds.Spis <- ddsHTSeq[grep("Spis", rownames(ddsHTSeq)),]
  dds.Spis$Reef.Site.Name <- relevel(dds.Spis$Reef.Site.Name, ref = "KAUST")
  dds.Spis <- DESeq(dds.Spis)
  vsd.Spis <- vst(dds.Spis, blind = F)
  
  
  for(condition in resultsNames(dds.Spis)[2:length(resultsNames(dds.Spis))]){
    # open pdf for volcano plot
    pdf(file = paste("../RNA-seq_Analysis/Spis_DEG/volcano_plots_svalues/", condition, ".pdf", sep=""))
    
    # what analysis is being performed
    cat("Running DE analysis for condition:", condition, "\n")
    
    cond0 <- strsplit(paste(condition, sep=""), split = "condition_")[[1]][1]
    cond1 <- strsplit(cond0, split = "_")[[1]][2]
    cond2 <- strsplit(cond0, split = "_")[[1]][4]
    cat("Cond1:", cond1, " Cond2:", cond2, "\n")
    
    # get results and use lfc shrinkage for visualization
    resLFC <- lfcShrink(dds = dds.Spis,
                        coef = condition,
                        type="apeglm",
                        lfcThreshold = lfc)
    
    # order results and write into a file
    res.ordered <- resLFC[order(resLFC$svalue),]
    
    # volcano plot
    resPlot <- as.data.frame(res.ordered)
    resPlot$Gene <- rownames(resPlot)
    lfc <- lfc
    p.value <- svalue
    #sig.genes <- dim(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc))[1]
    sig.genes <- dim(resPlot[which(abs(resPlot$log2FoldChange) > 2 & resPlot$svalue < p.value),])[1]
    outliers <- sum(resPlot$baseMean > 0 & is.na(resPlot$svalue))
    low.counts <- sum(!is.na(resPlot$svalue) & is.na(resPlot$svalue))
    minim <- resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = F)][1] + (resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = F)][1]/100) * 20
    maxim <- resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = T)][1] + (resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = T)][1]/100) * 20
    
    with(resPlot, plot(log2FoldChange, -log10(svalue),
                       pch=20, main= paste(condition),
                       sub=paste( sig.genes, " significant genes (L2FC >", lfc, "& svalue <", p.value, ")"),
                       # " | outliers:", outliers,
                       # "| low counts:", low.counts ),
                       xlim=c(minim,maxim),
                       col="gray", cex.sub=0.8, col.sub="gray", cex.lab=0.8))
    # Add colored points: red if svalue< p.value, orange of log2FC > lfc, green if both)
    with(subset(resPlot, svalue < p.value ), points(log2FoldChange, -log10(svalue), pch=20, col="#2c7bb6"))
    with(subset(resPlot, abs(log2FoldChange) > lfc), points(log2FoldChange, -log10(svalue), pch=20, col="#fdae61"))
    with(subset(resPlot, svalue < p.value & (abs(log2FoldChange) > lfc)), points(log2FoldChange, -log10(svalue), pch=20, col="#d7191c"))
    # Label points with the textxy function from the calibrate plot
    if(sig.genes > 9){
      with(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc)[1:10,], textxy(log2FoldChange, -log10(svalue), labs=Gene, cex=.6, col=rgb(0,0,0, 0.5)))
    }
    if(sig.genes <= 9){
      with(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc), textxy(log2FoldChange, -log10(svalue), labs=Gene, cex=.6, col=rgb(0,0,0, 0.5)))
    }
    abline(v=c(-lfc, lfc), h=c(-log10(svalue),-log10(svalue)), col="gray", lty=2)
    dev.off()
    
    # get summary and write into a file
    if(write2file){
      write.table(as.data.frame(res.ordered), 
                  file=paste("../RNA-seq_Analysis/Spis_DEG/DEG_tables_svalues/",condition,".txt", sep=""), sep="\t")
      
  }
  
   }
  
}


##### volcana_plot #####
volcano_plot <- function( condition = condition, results_dir="../RNA-seq_Analysis/Spis_DEG/DEG_tables_svalues/", lfc=2, svalue=0.005){
  # what analysis is being performed
  cat("Plotting volcano plot for:", condition, "\n")
  
  cond0 <- strsplit(paste(condition, sep=""), split = "condition_")[[1]][2]
  cond1 <- strsplit(cond0, split = "_vs_")[[1]][1]
  cond2 <- strsplit(cond0, split = "_vs_")[[1]][2]
  file <- paste(results_dir, condition, ".txt", sep="")
  
  df <- read.table(file, sep="\t", header = T)
  df <- df[order(df$svalue),]
  
  # volcano plot
  resPlot <- df
  resPlot$Gene <- rownames(resPlot)
  lfc <- lfc
  p.value <- svalue
  #sig.genes <- dim(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc))[1]
  sig.genes <- sum(resPlot$svalue < p.value & abs(resPlot$log2FoldChange) > 2, na.rm=TRUE)
  outliers <- sum(resPlot$baseMean > 0 & is.na(resPlot$pvalue))
  low.counts <- sum(!is.na(resPlot$pvalue) & is.na(resPlot$svalue))
  minim <- resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = F)][1] + (resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = F)][1]/100) * 20
  maxim <- resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = T)][1] + (resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = T)][1]/100) * 20
  
  deg <- dim(df[which(abs(df$log2FoldChange) > 2 & df$svalue < p.value),])[1]
  up <- dim(df[which(df$log2FoldChange > 2 & df$svalue < p.value),])[1]
  down <- dim(df[which(df$log2FoldChange < (-2) & df$svalue < p.value),])[1]
  
  with(resPlot, plot(log2FoldChange, -log10(svalue),
                     pch=20, main= paste(condition),
                     sub=paste( sig.genes, " significant genes (L2FC >", lfc, "& svalue <", p.value, ")"),
                     # " | outliers:", outliers,
                     # "| low counts:", low.counts ),
                     xlim=c(minim,maxim),
                     col="gray", cex.sub=0.8, col.sub="gray", cex.lab=0.8))
  # Add colored points: red if svalue< p.value, orange of log2FC > lfc, green if both)
  with(subset(resPlot, svalue < p.value ), points(log2FoldChange, -log10(svalue), pch=20, col="#2c7bb6"))
  with(subset(resPlot, abs(log2FoldChange) > lfc), points(log2FoldChange, -log10(svalue), pch=20, col="#fdae61"))
  with(subset(resPlot, svalue < p.value & (abs(log2FoldChange) > lfc)), points(log2FoldChange, -log10(svalue), pch=20, col="#d7191c"))
  # Label points with the textxy function from the calibrate plot
  if(sig.genes > 9){
    with(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc)[1:10,], textxy(log2FoldChange, -log10(svalue), labs=Gene, cex=.6, col=rgb(0,0,0, 0.5)))
  }
  if(sig.genes <= 9){
    with(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc), textxy(log2FoldChange, -log10(svalue), labs=Gene, cex=.6, col=rgb(0,0,0, 0.5)))
  }
  abline(v=c(-lfc, lfc), h=c(-log10(p.value),-log10(p.value)), col="gray", lty=2)
  mtext(text = paste("UP: ", up), line = 0, adj = 1, col="gray" )
  mtext(text = paste("DOWN: ", down), line = 0, adj = 0, col="gray" )
}


##### create_DEG_matrix #####
create_DEG_matrix <- function(p.value=0.005){
  setwd("/Volumes/omics4tb2/Collaborations/Vulcan/R_Projects/Stylophora/")
  # counts directory
  directory <- "../../R_Projects/RNA-seq_Analysis/Spis_DEG/DEG_tables_svalues/"
  
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
    cat("Processing ", sample, "\n")
      cond0 <- strsplit(paste(sample, sep=""), split = "condition_")[[1]][2]
      cond1 <- strsplit(cond0, split = "_vs_")[[1]][1]
      cond2 <- sub(".txt", "", strsplit(cond0, split = "_vs_")[[1]][2])

      df <- read.table(paste(directory, sample, ".txt",sep=""), sep="\t", header = T)
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


####################### Visualizations ########################
##### 1. plot heatmap for top 50/100 genes #####
plot_heatmap <- function(hm_file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_DEG_heatmap.pdf", dds=deseq_data$dds){
  # 1. heatmap
  pdf(file = hm_file)
  dds.Spis <- dds
  #dds.Spis <- DESeq(dds.Spis, parallel=F) # run DESeq2
  vsd.Spis <- vst(dds.Spis, blind=FALSE) # data transformation for visuals
  select <- order(rowMeans(counts(dds.Spis,normalized=TRUE)), # select top 100 most diff expressed
                  decreasing=T)[1:30]
  
  df.annotate <- as.data.frame(colData(dds.Spis)[,c("Temp","Reef.Site.Name")])
  df.assay <- assay(vsd.Spis)[select,]
  ## do not plot initially to get row ordering and reordering
  phm <- pheatmap(df.assay, scale = "none", cluster_rows=T, silent = T)
  phmr <- phm$tree_row$order
  df.order <- df.assay[rev(phmr),rownames(df.annotate[order(factor(df.annotate$Reef.Site.Name, levels=c("Protected","Exposed","KAUST","Eilat")),df.annotate$Temp),])]
 
  pheatmap(df.order, scale = "none", cluster_rows=F, border_color = NA, show_rownames=TRUE, # heatmap
             cluster_cols=F, annotation_col=df.annotate, main = "Spis Top 30 DEG", fontsize_col = 7, fontsize_row = 6, color = viridis(20) )
  dev.off()
}


##### 2. plot_pca: CA Plot #####
plot_pca <- function(pca_file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_PCA.pdf", vsd=deseq_data$vsd){
  pdf(file = pca_file )
  vsd.Spis <- vsd
  # remove outlier ES1-36A and create PCA object
  pcaData <- plotPCA(vsd.Spis[,grep("ES1-36A", colnames(vsd.Spis), invert = T)], intgroup=c("Temp","Reef.Site.Name","condition"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  q <- ggplot(pcaData, aes(PC1, PC2, color=Temp, shape=Reef.Site.Name)) +
    geom_point(size=4) +
    scale_shape_manual(values=c(16,15,17,3)) + 
    #geom_text(aes(label=condition)) + 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    theme_gray() + coord_fixed()
  print(q)
  dev.off()
}


##### 2.1 plot_pca_kaust: PCA Plot with KAUST only #####
plot_pca_kaust <- function(pca_file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_PCA_KAUST.pdf", vsd=deseq_data_kaust$vsd){
  pdf(file = pca_file )
  vsd.Spis <- vsd
  # remove outlier ES1-36A and create PCA object
  pcaData <- plotPCA(vsd.Spis[,grep("Eilat", vsd.Spis$Reef.Site.Name, invert = T)], intgroup=c("Temp","Reef.Site.Name","condition"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  q <- ggplot(pcaData, aes(PC1, PC2, color=Temp, shape=Reef.Site.Name)) +
    geom_point(size=4) +
    scale_shape_manual(values=c(15,17,3)) +
    #geom_text(aes(label=condition)) + 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    theme_gray() + coord_fixed()
  print(q)
  dev.off()
}

##### 3. plot_dist_heatmap: distance heatmap #####
plot_dist_heatmap <- function(dhm_file = "../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_distance_heatmap.pdf", dds=deseq_data$dds, vsd=deseq_data$vsd){
  pdf(file=dhm_file)
  dds.Spis <- dds
  vsd.Spis <- vsd
  df <- as.data.frame(colData(dds.Spis)[,c("Temp","Reef.Site.Name","Site","condition")])
  rownames(df) <- paste(vsd.Spis$Temp,vsd.Spis$Reef.Site.Name,vsd.Spis$sample, sep = "_")
  sampleDists <- dist(t(assay(vsd.Spis)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd.Spis$Temp,vsd.Spis$Reef.Site.Name,vsd.Spis$sample, sep = "_")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors, annotation_row=df)
  dev.off()
}

##### 4. plot_barplot: DEG Barplot #####
plot_barplot <- function(barplot_file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_DEG_barplot.pdf", vsd=deseq_data$vsd, df=DEG.matrix$updown){
  pdf(file = barplot_file )
  
  df <- DEG.matrix$updown
  df$total <- sapply(rownames(df), function(i) sum(as.numeric(as.character(df[i,"upregulated"])), as.numeric(as.character(df[i,"downregulated"]))))
  df$downregulated <- -(as.numeric(as.character(df$downregulated)))
  df$upregulated <- as.numeric(as.character(df$upregulated))
  df$label <- sapply(rownames(df), function(i) paste(df[i,"cond1"], "vs", df[i,"cond2"]))
  df$site1 <- sapply(as.character(df$cond1), function(i) strsplit(i, split = "_")[[1]][2])
  df$site2 <- sapply(as.character(df$cond2), function(i) strsplit(i, split = "_")[[1]][2])
  df$temp2 <- sapply(as.character(df$cond2), function(i) strsplit(i, split = "_")[[1]][1])
  
  # order <- df[order(df$total, decreasing = T),]$label
  # df <- df[order(df$total, decreasing = T),]
  # df<-within(df, label<- factor(label, levels=order))
  df.eilat <- df[which(df$site1 == "Eilat" & df$site2 == "Eilat"),]
  df.kaust <- df[which(df$site1 == "KAUST" & df$site2 == "KAUST"),]
  df.protected <- df[which(df$site1 == "Protected" & df$site2 == "Protected"),]
  df.exposed <- df[which(df$site1 == "Exposed" & df$site2 == "Exposed"),]
  df.test <- rbind(df.eilat, df.kaust, df.protected, df.exposed)
  
  
  
  df_melted.eilat <- melt(df.eilat, measure.vars = c("upregulated", "downregulated"))
  df_melted.eilat <- df_melted.eilat[order(df_melted.eilat$temp2),]
  
  df_melted.kaust <- melt(df.kaust, measure.vars = c("upregulated", "downregulated"))
  df_melted.kaust <- df_melted.kaust[order(df_melted.kaust$temp2),]
  
  df_melted.protected <- melt(df.protected, measure.vars = c("upregulated", "downregulated"))
  df_melted.protected <- df_melted.protected[order(df_melted.protected$temp2),]
  
  df_melted.exposed <- melt(df.exposed, measure.vars = c("upregulated", "downregulated"))
  df_melted.exposed <- df_melted.exposed[order(df_melted.exposed$temp2),]
  
  p1 <- ggplot(df_melted.eilat, aes(label))
  p1 <- p1 + geom_col(aes(y=value, fill=variable), position="stack")
  #p1 <- p1 + coord_flip() + scale_y_reverse() 
  p1 <- p1 + theme(axis.text.y = element_text(angle = 0, size = 7)) + theme_gray()
  p1 <- p1 + labs(x="Comparison", y="# of DEG", title="Spis DEG for all comparisons")
  p1 <- p1 + theme(legend.position = "none")
  #p1 <- p1 + lims(y=c(-1200, 500))
  #p1 <- p1 + facet_grid(.~label)
  p1
  
  p2 <- ggplot(df_melted.kaust, aes(label))
  p2 <- p2 + geom_col(aes(y=value, fill=variable), position="stack")
  #p2 <- p2 + coord_flip() + scale_y_reverse() 
  p2 <- p2 + theme(axis.text.y = element_text(angle = 0, size = 7)) + theme_gray()
  p2 <- p2 + theme(legend.position = "none")
  p2 <- p2 + labs(x="Comparison", y="# of DEG")
  p2 <- p2 + lims(y=c(-50, 50))
  #p2 <- p2 + facet_grid(.~label)
  p2
  
  p3 <- ggplot(df_melted.protected, aes(label))
  p3 <- p3 + geom_col(aes(y=value, fill=variable), position="stack")
  #p3 <- p3 + coord_flip() + scale_y_reverse() 
  p3 <- p3 + theme(axis.text.y = element_text(angle = 0, size = 7)) + theme_gray()
  p3 <- p3 + labs(x="Comparison", y="# of DEG")
  p3 <- p3 + theme(legend.position = "none")
  p3 <- p3 + lims(y=c(-50, 50))
  #p3 <- p3 + facet_grid(.~label)
  p3
  
  p4 <- ggplot(df_melted.exposed, aes(label))
  p4 <- p4 + geom_col(aes(y=value, fill=variable), position="stack")
  #p4 <- p4 + coord_flip() + scale_y_reverse() 
  p4 <- p4 + theme(axis.text.y = element_text(angle = 0, size = 7)) + theme_gray()
  p4 <- p4 + labs(x="Comparison", y="# of DEG")
  p4 <- p4 + theme(legend.position = "none")
  p4 <- p4 + lims(y=c(-50, 50))
  #p4 <- p4 + facet_grid(.~label)
  p4
  
  plot_grid(p1,p2, p3, p4, ncol=1)
  
  
  dev.off()
}

##### Replicate test: test to see optimal number of replicates #####
replicate_test <- function(plot=T, sampleTable=deseq_data$sample_table, directory=directory){
  ### test for replicates number
  sig.df <- data.frame()
  # Build DSEq dataset
  Eilat_samples <- sampleTable[which(sampleTable$Reef.Site.Name == "KAUST" & sampleTable$Temp == "36"), "sampleName"]
  Kaust_samples <- sampleTable[which(sampleTable$Reef.Site.Name == "KAUST" & sampleTable$Temp == "30"), "sampleName"]
  
  for(replicate_no in 2:7){
    for(run in 1:500){
      cat("Running Run#:", run, " for ", replicate_no, " replicates.\n")
      Eilat_select <- sample(Eilat_samples, size = replicate_no, replace = F)
      Kaust_select <- sample(Kaust_samples, size = replicate_no, replace = F)
      merged_select <- union(Eilat_select, Kaust_select)
      sampleTable1 <- sampleTable[sapply(merged_select, function(i) sampleTable[which(sampleTable$sampleName == i),"sampleName"]),]
      
      # Build DSEq dataset
      ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable1, design = ~ condition, directory=directory )
      
      # Seperate Spis vs Spis
      dds.Spis <- ddsHTSeq[grep("Spis", rownames(ddsHTSeq)),]
      dds.Spis$condition <- relevel(dds.Spis$condition, ref = "30_KAUST")
      
      # Run DEseq2
      dds.Spis <- DESeq(dds.Spis, parallel=F)
      
      # get results and use lfc shrinkage for visualization
      resLFC <- lfcShrink(dds = dds.Spis,
                          coef = "condition_36_KAUST_vs_30_KAUST" , 
                          type="apeglm", 
                          lfcThreshold = 2, 
                          svalue = TRUE, 
                          parallel = F)
      
      # order results and write into a file
      res.ordered <- resLFC[order(resLFC$svalue),]
      sig.genes <- dim(res.ordered[which(abs(res.ordered$log2FoldChange) > 2 & res.ordered$svalue < 0.005),])[1]
      eilat_reps <- paste(Eilat_select, sep="", collapse = ":")
      kaust_reps <- paste(Kaust_select, sep="", collapse = ":")
      sig.df <- rbind(sig.df, cbind(replicates=replicate_no, run=run, sig.genes=sig.genes, eilat_reps=eilat_reps, kaust_reps=kaust_reps))
      cat("# of Significant genes:", sig.genes, "\n")
    }
  }
  if(plot==T){
    plot_data <- sig.df
    plot_data$sig.count <- as.numeric(as.character(plot_data$sig.genes))
    p <- ggplot(plot_data, aes(x=replicates, y=sig.count, group=replicates, fill=replicates))
    p <- p + geom_violin()
    p <- p + geom_point(data=plot_data, aes(x=replicates, y=sig.count, group=replicates),alpha=0.3)
    p <- p + labs(subtitle="36 KAUST vs 30 KAUST", x="Replicates", y=" Number of differentially expressed genes")
    p
  }
  return(sig.df)
}

##### Kmeans clustering: function to perform k-means clustering on DEGs #####
kmeans_clustering <- function(km_file="../RNA-seq_Analysis/Spis_DEG/Spis_STAR_DEG_kmeans_plots.pdf", vsd=deseq_data$vsd){
  
  pdf(file=km_file)
  vsd.Spis <- vsd
  # get sample names from all directories for Stylophora
  samples <- list.files("../../R_Projects/RNA-seq_Analysis/Spis_DEG/DEG_tables/", recursive = T, pattern="(*_vs_30_KAUST*)", full.names = T)
  
  all_diff <- vector()
  for(sample in samples){
    df <- read.delim(sample, sep="\t", header=T, stringsAsFactors = F)
    diff_df <- df[which(abs(df$log2FoldChange) > 2 & df$svalue < 0.01),]
    diff_genes <- row.names(diff_df)
    all_diff <- append(all_diff, diff_genes)
  }
  # collect unique differentially expressed genes
  all_diff <- unique(all_diff)
  
  # get vsd normalized count matrix from DESeq2 and filter for diff expressed genes
  vsd_matrix <- assay(vsd.Spis)[all_diff,]
  
  vsd.Spis <- as.data.frame(vsd_matrix)
  distance <- get_dist(t(vsd.Spis))
  fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
  # determine optimal number of clusters
  set.seed(123)
  ## high efficiency clustering
  res.km <- eclust(vsd.Spis, "kmeans", nboot = 500)
  fviz_silhouette(res.km, print.summary = T)
  
  fviz_nbclust(vsd.Spis, kmeans, method = "wss") # find optimal number of clusters
  # run k-means
  Spis.k <- kmeans(vsd.Spis, centers = res.km$nbclust, nstart = 25)
  fviz_cluster(Spis.k, data = vsd.Spis, geom = "point", main = "Spis kmeans Clustering")
  write.table(data.frame(Spis.k$cluster), sep="\t", file="../RNA-seq_Analysis/Spis_DEG/Spis_STAR_DEG_kmeans.txt")
  dev.off()
}


##### TopGO Analysis #####
##### create gene2go object for function enrichment #####
gene2go.object <-function(file=NULL){
  file = "/Volumes/omics4tb2/Collaborations/Vulcan/Genomics/coral_Stylophora/spis.genome.genes.annotation_clean.txt"
  #file = "~/Dropbox/Snytrophy_Portal/uniprot-mmp.txt"
  cat("Loading genome annotation file from ", file, "\n")
  Proteome <- read.delim(file, header=T, sep="\t")
  cat("creating gene2go file ", "\n")
  
  gene2go = list()
  for(gene in Proteome$target_id){
    if(length(strsplit(as.character(Proteome[which(Proteome$target_id == gene),"GO.terms"]), split = ",", fixed=T)[[1]]) > 0){
      gene2go[gene] = strsplit(as.character(Proteome[which(Proteome$target_id== gene),"GO.terms"]), split = ",", fixed=T)
    }
  }
  return(gene2go)
}


##### Get  topGO data Object  [6]    #####
get_topGO_object <- function(genes,gene2go,ontology=c("BP","MF","CC")) {
  require(topGO)
  cat("Getting topGO object","\n")
  # genes is a vector containing genes of interest
  geneList <- factor(as.integer(names(gene2go)%in%genes))
  names(geneList) <- names(gene2go)
  GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene2go, nodeSize=10)
  #GOdata can be used directly for analysis, e.g. 
  test <- runTest(GOdata,algorithm="classic",statistic="fisher")
  results <- GenTable(GOdata,test,topNodes=20,orderBy="classic")
  return(GOdata)
}


##### Perform post processing of TopGO enrichment to create a table for all groups #####
post.process <- function(resultsfile=NULL){
  cat("Running Post processing.\n")
  BP <- data.frame()
  CC <- data.frame()
  MF <- data.frame()
  goresults.table <- data.frame()
  
  for(element in names(resultsfile)){
    if(length(resultsfile[[element]]$BP$GO.ID) != 0) {
      BP <- rbind(BP, cbind(group=element,resultsfile[[element]]$BP, GO="BP"))
    }
    
    if(length(resultsfile[[element]]$CC$GO.ID) != 0) {
      CC <- rbind(CC, cbind(group=element,resultsfile[[element]]$CC, GO="CC"))
    }
    
    if(length(resultsfile[[element]]$MF$GO.ID) != 0) {
      MF <- rbind(MF, cbind(group=element, resultsfile[[element]]$MF, GO="MF"))
    }
    
    goresults.table <- rbind(goresults.table, BP, CC, MF)
  }
  goresults.table <- unique(goresults.table)
  return(goresults.table)
}


##### enrichment analysis #####
run.topGO.enrichment <-function(my.members){
  gene2go <- gene2go.object()
  # Enrichment analysis for all biclusters and all Processes [7]
  #require(multicore)
  bc1.go <- lapply(seq(1,length(my.members)), function(first){
    # do not test biclusters with less than 2 members
    if(length(my.members[[first]])<2) {    ## added DJR
      cat("\n Not analyzing group ", first, "... \n")   ## added DJR
      return(NULL) ## added DJR
    }   ## added DJR
    # skip bicluster with less than two genes mapping to go terms
    if(sum(is.element(my.members[[first]], names(gene2go))) <2){
      cat("\n Not analyzing group ", first, "... \n")
      return(NULL)
    }
    
    # Do enrichment for all GO categories for all biclusters
    o <- lapply(c("BP", "CC", "MF"), function(second){
      cat("\n Now analyzing group ", names(my.members)[[first]], "with category ", second, "... \n")
      tmp = get_topGO_object(my.members[[first]], gene2go, second)
      test <- runTest(tmp,algorithm="classic",statistic="fisher")
      results <- GenTable(tmp,test,topNodes=length(test@score))
      results <- results[results[,6]<=0.01,]
      
    })
    names(o) <- c("BP", "CC", "MF")
    
    return(o)
  })
}



################################# Pipeline steps ########################
##### 1. prepare DESeq2 data #####
deseq_data <- prepare_DESeq2_data()

         # 1.1 prepare DESeq2 data for KAUST only
         #deseq_data_kaust <- prepare_DESeq2_data_kaust()


##### 2. run create_conditions function #####
conditions <- create_conditions(meta_data= deseq_data$meta_data)

##### 3. DE analysis and make volcano plots #####
tic("Differential Expression Analysis")
for(condition in conditions){
  # DE Analysis
  make_deg(condition=condition, dds=deseq_data$dds, lfc=2, write2file=T)
  # volcano plot
  pdf(file=paste("../RNA-seq_Analysis/Spis_DEG/volcano_plots_svalues/",condition,"_volcano.pdf", sep=""))
  volcano_plot(condition = condition, results_dir="../RNA-seq_Analysis/Spis_DEG/DEG_tables_svalues/", lfc=2, svalue=0.005)
  dev.off()
}
toc()

##### 4. create DEG matrix and write results into a file #####
DEG.matrix <- create_DEG_matrix()
write.table(DEG.matrix$matrix, file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_DEG_matrix.txt", sep="\t")
write.table(DEG.matrix$updown, file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_DEG_matrix_up and_down.txt", sep="\t")

pdf(file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_DEG_matrix_heatmap.pdf")
pheatmap(mat=DEG.matrix$matrix, cluster_rows = F, cluster_cols = F, display_numbers = T, number_format = "%.0f", main ="Spis DEG Matrix")
dev.off()

plot_barplot(barplot_file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_DEG_barplot.pdf")

##### 5. heatmap #####
plot_heatmap(hm_file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_DEG_heatmap.pdf")

##### 6. PCA #####
plot_pca(pca_file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_PCA.pdf")

##### 6.1 PCA for KAUST Sites only #####
plot_pca_kaust(pca_file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_PCA_KAUST.pdf", vsd=deseq_data$vsd)

##### 7. distance heatmap #####
plot_dist_heatmap(dhm_file = "../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_distance_heatmap.pdf")

##### 8. replicate test #####
replicate_profiles <- replicate_test(sampleTable = deseq_data$sample_table, directory = deseq_data$directory)
write.table(replicate_profiles, file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_replicates_test_36KAUST_vs_30KAUST.txt", sep="\t")

##### 9. Kmeans clustering of DEGs #####
kmeans_clustering(km_file="../RNA-seq_Analysis/Spis_DEG/Overall_Plots_svalues/Spis_STAR_DEG_kmeans_plots.pdf", vsd=deseq_data$vsd)

##### 10. temperature specific DEG comparison #####
make_deg_temp(sampleTable=deseq_data$sample_table, directory=deseq_data$directory, lfc=2, svalue=0.005, write2file=T)

##### 11. Site specific DEG comparison #####
make_deg_site(sampleTable=deseq_data$sample_table, directory=deseq_data$directory, lfc=2, svalue=0.005, write2file=T)
  
##### 12. TopGO Enrichment Analysis
Spis.enrichment <- run.topGO.enrichment(DEG.matrix$deg_genes)
names(Spis.enrichment) <- names(DEG.matrix$deg_genes)
Spis.go.enrichment <- post.process(Spis.enrichment)
write.table(Spis.go.enrichment, file="../RNA-seq_Analysis/Spis_DEG/Spis_DEG_GO_Enrichment.txt", sep="\t", row.names=F)
