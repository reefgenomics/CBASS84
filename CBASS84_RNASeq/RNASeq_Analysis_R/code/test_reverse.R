import_kallisto <- function(){
library('tximport')
  library('rhdf5')
  library("DESeq2")
  library('DT')
  library('tidyverse')
  library('stringr')
  # File directories
  my.dir <- "/mnt/omics4tb2-serdar/Collaborations/Vulcan/data/181203_K00235_0152_BHWNJ7BBXX_CBASS84_RNASeq"
  sample_id <- dir(file.path(my.dir))

  # remove folders for samples with issues and other script folders also remove ES1-36A since it seems to have very high number of reads across samples.
  sample_id <- sample_id[!(sample_id %in% c("Lanes_metadata", "ES1-36A","S6-P-ST-36B","KS2-33A","ES3-33A","all_htseqcounts", "label_changes_README.txt"))]

  # list of abundance files.
  files <- file.path(my.dir, sample_id, "results0920_kallisto_Spis_Smic", "abundance.h5")
  names(files) <- sample_id

  ## Collect sample metadata and prepare sample table
  # get metadata info from metadata worksheet and select relevant columns and delete NA rows
  meta_data <- read.delim("data/meta_data_v2.txt", header=T, sep="\t", stringsAsFactors = F)
  meta_data <- dplyr::select(meta_data, sample = Sample.Name.CBASS84, 1:3)

  # remove samples from metadata that were earlier filtered and order
  meta_data <- meta_data[grep("(S6_P_ST_36B|KS2_33A|ES3_33A|ES6_30A|ES1_36A)", meta_data$sample, invert=T),]
  meta_data <- meta_data[order(meta_data$sample),]
  meta_data$Temp <- as.factor(meta_data$Temp)
  meta_data$Reef.Site.Name <- as.factor(meta_data$Reef.Site.Name)
  meta_data <- meta_data %>%
    mutate(Reef.Site.Name = str_replace(Reef.Site.Name, "Exposed", "ExT")) %>%
    mutate(Reef.Site.Name = str_replace(Reef.Site.Name, "Protected", "PrT"))

  # create a metadata column by combining temperature and Reef.site.name
  meta_data$condition <- as.factor(paste(meta_data$Temp, meta_data$Reef.Site.Name, sep="_"))
  meta_data <- dplyr::mutate(meta_data, path = files)
  meta_data$Reef.Site.Name <- as.factor(meta_data$Reef.Site.Name)


  #rebuild_sample table
  sampleTable <- cbind(sampleName=sample_id, fileName=files, meta_data)

  ## Import kallisto abundances
  txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

  ## load data into DESEQ2
  dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)

  # Seperate Spis vs Spis
  cat("\t\t - filtering for Spis\n\n")
  dds.Spis <- dds[grep("Spis", rownames(dds)),]

  # Seperate Spis vs Spis
  cat("\t\t - filtering for Smic\n\n")
  dds.Smic <- dds[grep("Smic", rownames(dds)),]

  return(list(dds.Spis= dds.Spis, dds.Smic= dds.Smic))
}

## Run import
dds <- import_kallisto()

lfc = 2
### Condition 1
condition1 = "condition_36_ICN_vs_30_ICN"
dds1 <-dds$dds.Spis
# Run DEseq
dds1 <- DESeq(dds1)
# relevel to set cond2 as our reference sample
dds1$condition <- relevel(dds1$condition, ref = "30_ICN")
#dds1$condition <- factor(dds1$condition, levels = c("30_ICN", "33_ICN", "36_ICN","30_AF", "33_AF", "36_AF","30_ExT", "33_ExT", "36_ExT","30_PrT", "33_PrT", "36_PrT"))

# Run DEseq2
dds1 <- DESeq(dds1)


# get results and use lfc shrinkage for visualization
resLFC1 <- lfcShrink(
  dds = dds1,
  coef = condition1,
  type = "apeglm",
  lfcThreshold = lfc
)
## Order output
res1.ordered <- resLFC1[order(resLFC1$svalue),]
#########

### Condition 2
## Run import
dds <- import_kallisto()

condition2 = "condition_30_ICN_vs_36_ICN"
dds2 <-dds$dds.Spis

# Run DESEq2
dds2 <- DESeq(dds2)

# relevel to set cond2 as our reference sample
dds2$condition <- relevel(dds2$condition, ref = "36_ICN")
#dds2$condition <- factor(dds2$condition, levels = c("36_ICN", "33_ICN", "30_ICN","36_AF", "33_AF", "30_AF","36_ExT", "33_ExT", "30_ExT","36_PrT", "33_PrT", "30_PrT"))

# Run DEseq2
dds2 <- DESeq(dds2)

# get results and use lfc shrinkage for visualization
resLFC2 <- lfcShrink(
  dds = dds2,
  coef = condition2,
  type = "apeglm",
  lfcThreshold = lfc
)
## Order output
res2.ordered <- resLFC2[order(resLFC2$svalue),]


res1.ordered["Spis73", ]
res2.ordered["Spis73", ]
sum(res1.ordered$svalue < 0.0005, na.rm = TRUE)
sum(res2.ordered$svalue < 0.0005, na.rm = TRUE)

res1.ordered.0005 <- subset(res1.ordered, svalue < 0.0005)





res11 <- results(dds1, contrast = c("condition", "36_ICN", "30_ICN"), lfcThreshold = lfc, alpha = 0.05)
res11.0005 <- subset(res11, padj < 0.005)


res22 <- results(dds1, contrast = c("condition", "30_ICN", "36_ICN"), lfcThreshold = lfc, alpha = 0.05)
sum(res11$padj < 0.005, na.rm = TRUE)
sum(res22$padj < 0.005, na.rm = TRUE)

res111 <- lfcShrink(dds1, type = "normal",contrast = c("condition", "36_ICN", "30_ICN"), lfcThreshold = lfc)
res222 <- lfcShrink(dds1, type = "normal",contrast = c("condition", "30_ICN", "36_ICN"), lfcThreshold = lfc)
sum(res111$padj < 0.05, na.rm = TRUE)
sum(res222$padj < 0.05, na.rm = TRUE)

res1111 <- lfcShrink(dds1, type = "ashr",contrast = c("condition", "36_ICN", "30_ICN"), lfcThreshold = lfc)
res2222 <- lfcShrink(dds1, type = "ashr",contrast = c("condition", "30_ICN", "36_ICN"), lfcThreshold = lfc)
sum(res111$padj < 0.05, na.rm = TRUE)
sum(res222$padj < 0.05, na.rm = TRUE)
