
############################################################
####################### Stats ##############################
############################################################

library(vegan)
library(pairwiseAdonis)

otu.ad=read.table("./outputs/CBASS84_noConta_raw.txt", header = TRUE, row.names = 1)[,1:84]
otu.n=as.data.frame(t(sweep(otu.ad,2,colSums(otu.ad),`/`)))
rownames(otu.n)=gsub('\\.', '-', rownames(otu.n))
map=read.table("./Input_files/metadata.txt", header = TRUE, row.names = 1, sep = "\t")
otu.n$Temperature=map$Temperature[match(rownames(otu.n), rownames(map))] # change colnames by lookup table
otu.n$Site=map$Region[match(rownames(otu.n), rownames(map))]
adonis(otu.n[,1:2264] ~ otu.n$Site * otu.n$Temperature, permutations = 999, method = "bray")

#comparing temperatures
ICN=otu.n[grep("ES", rownames(otu.n)),]
al_fahal=otu.n[grep("KS", rownames(otu.n)),]
exp=otu.n[grep("-E-", rownames(otu.n)),]
prot=otu.n[grep("-P-", rownames(otu.n)),]

pairwise.adonis(ICN[,1:2264], ICN$Temperature , p.adjust.m ='fdr', sim.method = 'bray')
pairwise.adonis(al_fahal[,1:2264], al_fahal$Temperature , p.adjust.m ='fdr', sim.method = 'bray')
pairwise.adonis(exp[,1:2264], exp$Temperature , p.adjust.m ='fdr', sim.method = 'bray')
pairwise.adonis(prot[,1:2264], prot$Temperature , p.adjust.m ='fdr', sim.method = 'bray')

#comparing sites
T30=otu.n[grep(".30[AB]", rownames(otu.n)),]
T33=otu.n[grep(".33[AB]", rownames(otu.n)),]
T36=otu.n[grep(".36[AB]", rownames(otu.n)),]

pairwise.adonis(T30[,1:2264], T30$Site , p.adjust.m ='fdr', sim.method = 'bray')
pairwise.adonis(T33[,1:2264], T33$Site , p.adjust.m ='fdr', sim.method = 'bray')
pairwise.adonis(T36[,1:2264], T36$Site , p.adjust.m ='fdr', sim.method = 'bray')

