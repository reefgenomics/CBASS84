library(DESeq2)
library(data.table)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/CBASS84/")

de.in=read.table("./outputs/CBASS84_noConta_raw.txt", header = TRUE, row.names = 1)[,1:84]
colnames(de.in)=gsub('\\.', '-', colnames(de.in))
met=read.table("./Input_files/metadata.txt", header = TRUE, row.names = 1, sep = "\t")
met$Temperature=as.factor(met$Temperature)


#################################
#### Comparing Temperatures ####
#################################


##DE analysis ICN
E.30.33.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("ES[0-7]-3[03]", names(de.in))], colData = subset(met, rownames(met) %like% "ES[0-7]-3[03]"), design= ~Temperature)
E.30.33.dss$condition = relevel(E.30.33.dss$Temperature, ref="30")  
E.30.33.dss = DESeq(E.30.33.dss) 
E.30.33.res=subset(as.data.frame(results(E.30.33.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
E.30.33.res$Comparison="ICN 30 vs 33"
E.30.33.res$OTU=rownames(E.30.33.res)

E.30.36.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("ES[0-7]-3[06]", names(de.in))], colData = subset(met, rownames(met) %like% "ES[0-7]-3[06]"), design= ~Temperature)
E.30.36.dss$condition = relevel(E.30.36.dss$Temperature, ref="30")  
E.30.36.dss = DESeq(E.30.36.dss) 
E.30.36.res=subset(as.data.frame(results(E.30.36.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
E.30.36.res$Comparison="ICN 30 vs 36"
E.30.36.res$OTU=rownames(E.30.36.res)

E.33.36.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("ES[0-7]-3[36]", names(de.in))], colData = subset(met, rownames(met) %like% "ES[0-7]-3[36]"), design= ~Temperature)
E.33.36.dss$condition = relevel(E.33.36.dss$Temperature, ref="33")  
E.33.36.dss = DESeq(E.33.36.dss) 
E.33.36.res=subset(as.data.frame(results(E.33.36.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
E.33.36.res$Comparison="ICN 33 vs 36"
E.33.36.res$OTU=rownames(E.33.36.res)

##DE analysis Al-Fahal
K.30.33.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("KS[0-7]-3[03]", names(de.in))], colData = subset(met, rownames(met) %like% "KS[0-7]-3[03]"), design= ~Temperature)
K.30.33.dss$condition = relevel(K.30.33.dss$Temperature, ref="30")  
K.30.33.dss = DESeq(K.30.33.dss) 
K.30.33.res=subset(as.data.frame(results(K.30.33.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
K.30.33.res$Comparison="Al-Fahal 30 vs 33"
K.30.33.res$OTU=rownames(K.30.33.res)

K.30.36.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("KS[0-7]-3[06]", names(de.in))], colData = subset(met, rownames(met) %like% "KS[0-7]-3[06]"), design= ~Temperature)
K.30.36.dss$condition = relevel(K.30.36.dss$Temperature, ref="30")  
K.30.36.dss = DESeq(K.30.36.dss) 
K.30.36.res=subset(as.data.frame(results(K.30.36.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
K.30.36.res$Comparison="Al-Fahal 30 vs 36"
K.30.36.res$OTU=rownames(K.30.36.res)

K.33.36.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("KS[0-7]-3[36]", names(de.in))], colData = subset(met, rownames(met) %like% "KS[0-7]-3[36]"), design= ~Temperature)
K.33.36.dss$condition = relevel(K.33.36.dss$Temperature, ref="33")  
K.33.36.dss = DESeq(K.33.36.dss) 
K.33.36.res=subset(as.data.frame(results(K.33.36.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
K.33.36.res$Comparison="Al-Fahal 33 vs 36"
K.33.36.res$OTU=rownames(K.33.36.res)

##DE analysis Exposed
X.30.33.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("S[0-7]-E-ST-3[03]", names(de.in))], colData = subset(met, rownames(met) %like% "S[0-7]-E-ST-3[03]"), design= ~Temperature)
X.30.33.dss$condition = relevel(X.30.33.dss$Temperature, ref="30")  
X.30.33.dss = DESeq(X.30.33.dss) 
X.30.33.res=subset(as.data.frame(results(X.30.33.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
X.30.33.res$Comparison="Exposed 30 vs 33"
X.30.33.res$OTU=rownames(X.30.33.res)

X.30.36.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("S[0-7]-E-ST-3[06]", names(de.in))], colData = subset(met, rownames(met) %like% "S[0-7]-E-ST-3[06]"), design= ~Temperature)
X.30.36.dss$condition = relevel(X.30.36.dss$Temperature, ref="30")  
X.30.36.dss = DESeq(X.30.36.dss) 
X.30.36.res=subset(as.data.frame(results(X.30.36.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
X.30.36.res$Comparison="Exposed 30 vs 36"
X.30.36.res$OTU=rownames(X.30.36.res)

X.33.36.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("S[0-7]-E-ST-3[36]", names(de.in))], colData = subset(met, rownames(met) %like% "S[0-7]-E-ST-3[36]"), design= ~Temperature)
X.33.36.dss$condition = relevel(X.33.36.dss$Temperature, ref="33")  
X.33.36.dss = DESeq(X.33.36.dss) 
X.33.36.res=subset(as.data.frame(results(X.33.36.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
X.33.36.res$Comparison="Exposed 33 vs 36"
X.33.36.res$OTU=rownames(X.33.36.res)

##DE analysis Protected
P.30.33.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("S[0-7]-P-ST-3[03]", names(de.in))], colData = subset(met, rownames(met) %like% "S[0-7]-P-ST-3[03]"), design= ~Temperature)
P.30.33.dss$condition = relevel(P.30.33.dss$Temperature, ref="30")  
P.30.33.dss = DESeq(P.30.33.dss) 
P.30.33.res=subset(as.data.frame(results(P.30.33.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
P.30.33.res$Comparison="Exposed 30 vs 33"
P.30.33.res$OTU=rownames(P.30.33.res)

P.30.36.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("S[0-7]-P-ST-3[06]", names(de.in))], colData = subset(met, rownames(met) %like% "S[0-7]-P-ST-3[06]"), design= ~Temperature)
P.30.36.dss$condition = relevel(P.30.36.dss$Temperature, ref="30")  
P.30.36.dss = DESeq(P.30.36.dss) 
P.30.36.res=subset(as.data.frame(results(P.30.36.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
P.30.36.res$Comparison="Exposed 30 vs 36"
P.30.36.res$OTU=rownames(P.30.36.res)

P.33.36.dss = DESeqDataSetFromMatrix(countData = de.in[,grepl("S[0-7]-P-ST-3[36]", names(de.in))], colData = subset(met, rownames(met) %like% "S[0-7]-P-ST-3[36]"), design= ~Temperature)
P.33.36.dss$condition = relevel(P.33.36.dss$Temperature, ref="33")  
P.33.36.dss = DESeq(P.33.36.dss) 
P.33.36.res=subset(as.data.frame(results(P.33.36.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
P.33.36.res$Comparison="Exposed 33 vs 36"
P.33.36.res$OTU=rownames(P.33.36.res)

DESeq=do.call(rbind, mget(ls(pattern="*.res")))
tax = read.csv("./Input_files/CBASS84.final.taxonomy", header = TRUE, sep = "\t")
tax$Taxonomy = gsub("[(0-9)]", "", tax$Taxonomy) 
tax$Taxonomy = gsub("[a-z]__", "", tax$Taxonomy) 
tax[c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')] <- colsplit(tax$Taxonomy,';',c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
tax=tax[, -c(3)]

DESeq$Family=tax$Family[match(DESeq$OTU,tax$OTU)]
DESeq$Genus=tax$Genus[match(DESeq$OTU,tax$OTU)]

#ggplot(data = DESeq,  aes(x = rownames(DESeq), y = log2FoldChange)) + geom_bar(stat = "identity", aes(fill = Family), position=position_dodge()) +  ylab("Log2 Fold Change") +  coord_flip() + facet_grid(~Comparison)

write.table(DESeq, "./outputs/DESeq_summary.txt", quote = FALSE, sep = "\t")


##### Number of DE OTUs

message("Number of DA OTUs in ICN 30 vs 33: ", nrow(E.30.33.res), ". More abundant: ", nrow(subset(E.30.33.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(E.30.33.res, log2FoldChange <0 )))
message("Number of DA OTUs in ICN 30 vs 36: ", nrow(E.30.36.res), ". More abundant: ", nrow(subset(E.30.36.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(E.30.36.res, log2FoldChange <0 )))
message("Number of DA OTUs in ICN 33 vs 36: ", nrow(E.33.36.res), ". More abundant: ", nrow(subset(E.33.36.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(E.33.36.res, log2FoldChange <0 )))

message("Number of DA OTUs in Al-Fahal 30 vs 33: ", nrow(K.30.33.res), ". More abundant: ", nrow(subset(K.30.33.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(K.30.33.res, log2FoldChange <0 )))
message("Number of DA OTUs in Al-Fahal 30 vs 36: ", nrow(K.30.36.res), ". More abundant: ", nrow(subset(K.30.36.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(K.30.36.res, log2FoldChange <0 )))
message("Number of DA OTUs in Al-Fahal 33 vs 36: ", nrow(K.33.36.res), ". More abundant: ", nrow(subset(K.33.36.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(K.33.36.res, log2FoldChange <0 )))

message("Number of DA OTUs in Exposed 30 vs 33: ", nrow(X.30.33.res), ". More abundant: ", nrow(subset(X.30.33.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(X.30.33.res, log2FoldChange <0 )))
message("Number of DA OTUs in Exposed 30 vs 36: ", nrow(X.30.36.res), ". More abundant: ", nrow(subset(X.30.36.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(X.30.36.res, log2FoldChange <0 )))
message("Number of DA OTUs in Exposed 33 vs 36: ", nrow(X.33.36.res), ". More abundant: ", nrow(subset(X.33.36.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(X.33.36.res, log2FoldChange <0 )))

message("Number of DA OTUs in Protected 30 vs 33: ", nrow(P.30.33.res), ". More abundant: ", nrow(subset(P.30.33.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(P.30.33.res, log2FoldChange <0 )))
message("Number of DA OTUs in Protected 30 vs 36: ", nrow(P.30.36.res), ". More abundant: ", nrow(subset(P.30.36.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(P.30.36.res, log2FoldChange <0 )))
message("Number of DA OTUs in Protected 33 vs 36: ", nrow(P.33.36.res), ". More abundant: ", nrow(subset(P.33.36.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(P.33.36.res, log2FoldChange <0 )))


#########################
#### Comparing Sites ####
#########################

##DE analysis T30
T30.E.K.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("[EK]S[0-7]-30[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "[EK]S[0-7]-30[AB]"), design= ~Region)
T30.E.K.dss$condition = relevel(T30.E.K.dss$Region, ref="ICN")  
T30.E.K.dss = DESeq(T30.E.K.dss) 
T30.E.K.res=subset(as.data.frame(results(T30.E.K.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T30.E.K.res$Comparison="T30 ICN vs Al-Fahal"
T30.E.K.res$OTU=rownames(T30.E.K.res)

T30.E.X.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("ES[0-7]-30[AB]|S[0-7]-E-ST-30[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "ES[0-7]-30[AB]|S[0-7]-E-ST-30[AB]"), design= ~Region)
T30.E.X.dss$condition = relevel(T30.E.X.dss$Region, ref="ICN")  
T30.E.X.dss = DESeq(T30.E.X.dss) 
T30.E.X.res=subset(as.data.frame(results(T30.E.X.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T30.E.X.res$Comparison="T30 ICN vs Exposed"
T30.E.X.res$OTU=rownames(T30.E.X.res)

T30.E.P.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("ES[0-7]-30[AB]|S[0-7]-P-ST-30[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "ES[0-7]-30[AB]|S[0-7]-P-ST-30[AB]"), design= ~Region)
T30.E.P.dss$condition = relevel(T30.E.P.dss$Region, ref="ICN")  
T30.E.P.dss = DESeq(T30.E.P.dss) 
T30.E.P.res=subset(as.data.frame(results(T30.E.P.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T30.E.P.res$Comparison="T30 ICN vs Protected"
T30.E.P.res$OTU=rownames(T30.E.P.res)

T30.K.X.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("KS[0-7]-30[AB]|S[0-7]-E-ST-30[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "KS[0-7]-30[AB]|S[0-7]-E-ST-30[AB]"), design= ~Region)
T30.K.X.dss$condition = relevel(T30.K.X.dss$Region, ref="Al-Fahal")  
T30.K.X.dss = DESeq(T30.K.X.dss) 
T30.K.X.res=subset(as.data.frame(results(T30.K.X.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T30.K.X.res$Comparison="T30 Al-Fahal vs Exposed"
T30.K.X.res$OTU=rownames(T30.K.X.res)

T30.K.P.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("KS[0-7]-30[AB]|S[0-7]-P-ST-30[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "KS[0-7]-30[AB]|S[0-7]-P-ST-30[AB]"), design= ~Region)
T30.K.P.dss$condition = relevel(T30.K.P.dss$Region, ref="Al-Fahal")  
T30.K.P.dss = DESeq(T30.K.P.dss) 
T30.K.P.res=subset(as.data.frame(results(T30.K.P.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T30.K.P.res$Comparison="T30 Al-Fahal vs Protected"
T30.K.P.res$OTU=rownames(T30.K.P.res)

T30.X.P.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("S[0-7]-[EP]-ST-30[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "S[0-7]-[EP]-ST-30[AB]"), design= ~Region)
T30.X.P.dss$condition = relevel(T30.X.P.dss$Region, ref="Exposed")  
T30.X.P.dss = DESeq(T30.X.P.dss) 
T30.X.P.res=subset(as.data.frame(results(T30.X.P.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T30.X.P.res$Comparison="T30 Exposed vs Protected"
T30.X.P.res$OTU=rownames(T30.X.P.res)

##DE analysis T33
T33.E.K.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("[EK]S[0-7]-33[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "[EK]S[0-7]-33[AB]"), design= ~Region)
T33.E.K.dss$condition = relevel(T33.E.K.dss$Region, ref="ICN")  
T33.E.K.dss = DESeq(T33.E.K.dss) 
T33.E.K.res=subset(as.data.frame(results(T33.E.K.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T33.E.K.res$Comparison="T33 ICN vs Al-Fahal"
T33.E.K.res$OTU=rownames(T33.E.K.res)

T33.E.X.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("ES[0-7]-33[AB]|S[0-7]-E-ST-33[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "ES[0-7]-33[AB]|S[0-7]-E-ST-33[AB]"), design= ~Region)
T33.E.X.dss$condition = relevel(T33.E.X.dss$Region, ref="ICN")  
T33.E.X.dss = DESeq(T33.E.X.dss) 
T33.E.X.res=subset(as.data.frame(results(T33.E.X.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T33.E.X.res$Comparison="T33 ICN vs Exposed"
T33.E.X.res$OTU=rownames(T33.E.X.res)

T33.E.P.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("ES[0-7]-33[AB]|S[0-7]-P-ST-33[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "ES[0-7]-33[AB]|S[0-7]-P-ST-33[AB]"), design= ~Region)
T33.E.P.dss$condition = relevel(T33.E.P.dss$Region, ref="ICN")  
T33.E.P.dss = DESeq(T33.E.P.dss) 
T33.E.P.res=subset(as.data.frame(results(T33.E.P.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T33.E.P.res$Comparison="T33 ICN vs Protected"
T33.E.P.res$OTU=rownames(T33.E.P.res)

T33.K.X.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("KS[0-7]-33[AB]|S[0-7]-E-ST-33[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "KS[0-7]-33[AB]|S[0-7]-E-ST-33[AB]"), design= ~Region)
T33.K.X.dss$condition = relevel(T33.K.X.dss$Region, ref="Al-Fahal")  
T33.K.X.dss = DESeq(T33.K.X.dss) 
T33.K.X.res=subset(as.data.frame(results(T33.K.X.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T33.K.X.res$Comparison="T33 Al-Fahal vs Exposed"
T33.K.X.res$OTU=rownames(T33.K.X.res)

T33.K.P.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("KS[0-7]-33[AB]|S[0-7]-P-ST-33[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "KS[0-7]-33[AB]|S[0-7]-P-ST-33[AB]"), design= ~Region)
T33.K.P.dss$condition = relevel(T33.K.P.dss$Region, ref="Al-Fahal")  
T33.K.P.dss = DESeq(T33.K.P.dss) 
T33.K.P.res=subset(as.data.frame(results(T33.K.P.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T33.K.P.res$Comparison="T33 Al-Fahal vs Protected"
T33.K.P.res$OTU=rownames(T33.K.P.res)

T33.X.P.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("S[0-7]-[EP]-ST-33[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "S[0-7]-[EP]-ST-33[AB]"), design= ~Region)
T33.X.P.dss$condition = relevel(T33.X.P.dss$Region, ref="Exposed")  
T33.X.P.dss = DESeq(T33.X.P.dss) 
T33.X.P.res=subset(as.data.frame(results(T33.X.P.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T33.X.P.res$Comparison="T33 Exposed vs Protected"
T33.X.P.res$OTU=rownames(T33.X.P.res)

##DE analysis T36
T36.E.K.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("[EK]S[0-7]-36[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "[EK]S[0-7]-36[AB]"), design= ~Region)
T36.E.K.dss$condition = relevel(T36.E.K.dss$Region, ref="ICN")  
T36.E.K.dss = DESeq(T36.E.K.dss) 
T36.E.K.res=subset(as.data.frame(results(T36.E.K.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T36.E.K.res$Comparison="T36 ICN vs Al-Fahal"
T36.E.K.res$OTU=rownames(T36.E.K.res)

T36.E.X.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("ES[0-7]-36[AB]|S[0-7]-E-ST-36[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "ES[0-7]-36[AB]|S[0-7]-E-ST-36[AB]"), design= ~Region)
T36.E.X.dss$condition = relevel(T36.E.X.dss$Region, ref="ICN")  
T36.E.X.dss = DESeq(T36.E.X.dss) 
T36.E.X.res=subset(as.data.frame(results(T36.E.X.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T36.E.X.res$Comparison="T36 ICN vs Exposed"
T36.E.X.res$OTU=rownames(T36.E.X.res)

T36.E.P.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("ES[0-7]-36[AB]|S[0-7]-P-ST-36[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "ES[0-7]-36[AB]|S[0-7]-P-ST-36[AB]"), design= ~Region)
T36.E.P.dss$condition = relevel(T36.E.P.dss$Region, ref="ICN")  
T36.E.P.dss = DESeq(T36.E.P.dss) 
T36.E.P.res=subset(as.data.frame(results(T36.E.P.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T36.E.P.res$Comparison="T36 ICN vs Protected"
T36.E.P.res$OTU=rownames(T36.E.P.res)

T36.K.X.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("KS[0-7]-36[AB]|S[0-7]-E-ST-36[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "KS[0-7]-36[AB]|S[0-7]-E-ST-36[AB]"), design= ~Region)
T36.K.X.dss$condition = relevel(T36.K.X.dss$Region, ref="Al-Fahal")  
T36.K.X.dss = DESeq(T36.K.X.dss) 
T36.K.X.res=subset(as.data.frame(results(T36.K.X.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T36.K.X.res$Comparison="T36 Al-Fahal vs Exposed"
T36.K.X.res$OTU=rownames(T36.K.X.res)

T36.K.P.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("KS[0-7]-36[AB]|S[0-7]-P-ST-36[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "KS[0-7]-36[AB]|S[0-7]-P-ST-36[AB]"), design= ~Region)
T36.K.P.dss$condition = relevel(T36.K.P.dss$Region, ref="Al-Fahal")  
T36.K.P.dss = DESeq(T36.K.P.dss) 
T36.K.P.res=subset(as.data.frame(results(T36.K.P.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T36.K.P.res$Comparison="T36 Al-Fahal vs Protected"
T36.K.P.res$OTU=rownames(T36.K.P.res)

T36.X.P.dss = DESeqDataSetFromMatrix(countData = de.in[,grep("S[0-7]-[EP]-ST-36[AB]", colnames(de.in))], colData = subset(met, rownames(met) %like% "S[0-7]-[EP]-ST-36[AB]"), design= ~Region)
T36.X.P.dss$condition = relevel(T36.X.P.dss$Region, ref="Exposed")  
T36.X.P.dss = DESeq(T36.X.P.dss) 
T36.X.P.res=subset(as.data.frame(results(T36.X.P.dss, lfcThreshold = 0.0, alpha=0.05)), padj < 0.05)
T36.X.P.res$Comparison="T36 Exposed vs Protected"
T36.X.P.res$OTU=rownames(T36.X.P.res)


DESeq2=do.call(rbind, mget(ls(pattern="*.res")))
tax = read.csv("./Input_files/CBASS84.final.taxonomy", header = TRUE, sep = "\t")
tax$Taxonomy = gsub("[(0-9)]", "", tax$Taxonomy) 
tax$Taxonomy = gsub("[a-z]__", "", tax$Taxonomy) 
tax[c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')] <- colsplit(tax$Taxonomy,';',c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
tax=tax[, -c(3)]

DESeq2$Family=tax$Family[match(DESeq2$OTU,tax$OTU)]
DESeq2$Genus=tax$Genus[match(DESeq2$OTU,tax$OTU)]

write.table(DESeq2, "./outputs/DESeq2_sites_summary.txt", quote = FALSE, sep = "\t")


##### Number of DE OTUs

message("Number of DA OTUs in T30 ICN vs Al-Fahal: ", nrow(T30.E.K.res), ". More abundant: ", nrow(subset(T30.E.K.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T30.E.K.res, log2FoldChange <0 )))
message("Number of DA OTUs in T30 ICN vs Exposed: ", nrow(T30.E.X.res), ". More abundant: ", nrow(subset(T30.E.X.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T30.E.X.res, log2FoldChange <0 )))
message("Number of DA OTUs in T30 ICN vs Protected: ", nrow(T30.E.P.res), ". More abundant: ", nrow(subset(T30.E.P.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T30.E.P.res, log2FoldChange <0 )))
message("Number of DA OTUs in T30 Al-Fahal vs Exposed: ", nrow(T30.K.X.res), ". More abundant: ", nrow(subset(T30.K.X.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T30.K.X.res, log2FoldChange <0 )))
message("Number of DA OTUs in T30 Al-Fahal vs Protected: ", nrow(T30.K.P.res), ". More abundant: ", nrow(subset(T30.K.P.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T30.K.P.res, log2FoldChange <0 )))
message("Number of DA OTUs in T30 Exposed vs Protected: ", nrow(T30.X.P.res), ". More abundant: ", nrow(subset(T30.X.P.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T30.X.P.res, log2FoldChange <0 )))

message("Number of DA OTUs in T33 ICN vs Al-Fahal: ", nrow(T33.E.K.res), ". More abundant: ", nrow(subset(T33.E.K.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T33.E.K.res, log2FoldChange <0 )))
message("Number of DA OTUs in T33 ICN vs Exposed: ", nrow(T33.E.X.res), ". More abundant: ", nrow(subset(T33.E.X.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T33.E.X.res, log2FoldChange <0 )))
message("Number of DA OTUs in T33 ICN vs Protected: ", nrow(T33.E.P.res), ". More abundant: ", nrow(subset(T33.E.P.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T33.E.P.res, log2FoldChange <0 )))
message("Number of DA OTUs in T33 Al-Fahal vs Exposed: ", nrow(T33.K.X.res), ". More abundant: ", nrow(subset(T33.K.X.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T33.K.X.res, log2FoldChange <0 )))
message("Number of DA OTUs in T33 Al-Fahal vs Protected: ", nrow(T33.K.P.res), ". More abundant: ", nrow(subset(T33.K.P.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T33.K.P.res, log2FoldChange <0 )))
message("Number of DA OTUs in T33 Exposed vs Protected: ", nrow(T33.X.P.res), ". More abundant: ", nrow(subset(T33.X.P.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T33.X.P.res, log2FoldChange <0 )))

message("Number of DA OTUs in T36 ICN vs Al-Fahal: ", nrow(T36.E.K.res), ". More abundant: ", nrow(subset(T36.E.K.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T36.E.K.res, log2FoldChange <0 )))
message("Number of DA OTUs in T36 ICN vs Exposed: ", nrow(T36.E.X.res), ". More abundant: ", nrow(subset(T36.E.X.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T36.E.X.res, log2FoldChange <0 )))
message("Number of DA OTUs in T36 ICN vs Protected: ", nrow(T36.E.P.res), ". More abundant: ", nrow(subset(T36.E.P.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T36.E.P.res, log2FoldChange <0 )))
message("Number of DA OTUs in T36 Al-Fahal vs Exposed: ", nrow(T36.K.X.res), ". More abundant: ", nrow(subset(T36.K.X.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T36.K.X.res, log2FoldChange <0 )))
message("Number of DA OTUs in T36 Al-Fahal vs Protected: ", nrow(T36.K.P.res), ". More abundant: ", nrow(subset(T36.K.P.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T36.K.P.res, log2FoldChange <0 )))
message("Number of DA OTUs in T36 Exposed vs Protected: ", nrow(T36.X.P.res), ". More abundant: ", nrow(subset(T36.X.P.res, log2FoldChange >0 )), ". Less abundant: ", nrow(subset(T36.X.P.res, log2FoldChange <0 )))
