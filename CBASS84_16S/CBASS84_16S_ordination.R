############################################################
##################### phyloseq #############################
############################################################
library(phyloseq)
library(microbiome)

#data in
otu.2=read.table("./outputs/CBASS84_noConta_raw.txt", header = TRUE, row.names = 1)
colnames(otu.2)=gsub('\\.', '-', colnames(otu.2))
map=read.table("./Input_files/metadata.txt", header = TRUE, row.names = 1)
map$Temperature=as.factor(map$Temperature)
map$ID=rownames(map)
map$Region=gsub("Kaust", "KAUST", map$Region)
map$Region=factor(map$Region, levels = c("Eilat", "KAUST", "Exposed","Protected"))
tax=read.table("outputs/CBASS84_noConta_tax", header = T, row.names = 1)

#create phyloseq object
otu.t= otu_table(otu.2[,1:84], taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))
phy.all= phyloseq(otu.t, tax.t,  sam.t)

#plot ordination
P4=c('#2c7bb6','#fdae61','#df65b0', '#91003f') 

phy.t=microbiome::transform(phy.all, transform = "compositional", target = "OTU", shift = 0, scale = 1)

PCOA_br = ordinate(phy.t, method = "PCoA", distance = "bray")
pdf(file = "./outputs/CBASS84_ordination.pdf",  width = 7, height = 5, pointsize = 12) 
plot_ordination(phy.t, PCOA_br, color = "Region", shape = "Region")  + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P4) + labs(x="PCoA1: 47.3% variance",y="PCoA 2: 27.7% variance") + theme(plot.title = element_text(hjust = 0.5)) +  scale_shape_manual(values=c(19,17,15,18))
dev.off()
