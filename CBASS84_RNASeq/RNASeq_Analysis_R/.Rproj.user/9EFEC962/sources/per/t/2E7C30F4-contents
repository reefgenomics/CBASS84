##### make_deg: DESeq analysis for all conditions #####
deg_analysis <- function(condition=condition, dds=dds.Spis, lfc=2, write2file=T, sampleTable=sampleTable,svalue=0.005, org=c("Spis", "Smic")){
  library("reshape2")
  library("DESeq2")

  cat("Performing DEG Analysis \n")

  # what analysis is being performed
  cat("\t Running DE analysis for condition:", condition, "\n")

  cond0 <- strsplit(paste(condition, sep=""), split = "condition_")[[1]][2]
  cond1 <- strsplit(cond0, split = "_vs_")[[1]][1]
  cond2 <- strsplit(cond0, split = "_vs_")[[1]][2]

  cat("\t\t Cond1:", cond1, " Cond2:", cond2, "\n")

  # relevel to set cond2 as our reference sample
  dds$condition <- relevel(dds$condition, ref = cond2)

  # Run DEseq2
  dds <- DESeq(dds, parallel=F, BPPARAM=MulticoreParam(8))

  # get results
  #res <- results(dds.Spis, name = condition, alpha = 0.05)
  # # get results and use lfc shrinkage for visualization
  resLFC <- lfcShrink(dds = dds,
                      coef = condition,
                      type="apeglm",
                      lfcThreshold = lfc,parallel = F,BPPARAM = MulticoreParam(8) )

  # order results and write into a file
  #res.ordered <- resLFC[order(resLFC$svalue),]
  res.ordered <- resLFC[order(resLFC$svalue),]


  # get summary and write into a file
  if(write2file){
    write.table(as.data.frame(res.ordered),
                file=paste("output/", org, "/DEG_tables_svalues/",condition,".txt", sep=""), sep="\t")
  }
  return(as.data.frame(res.ordered))
}

