library(reshape2)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/CBASS84/")

## remove potential contaminants
otu = t(read.table("./Input_files/CBASS84.final.OTU_table", header = FALSE))
colnames(otu) = otu[2, ]
otu=otu[-c(1,2,3), ] 
tax = read.csv("./Input_files/CBASS84.final.taxonomy", header = TRUE, sep = "\t")
tax$Taxonomy = gsub("[(0-9)]", "", tax$Taxonomy) 
tax$Taxonomy = gsub("[a-z]__", "", tax$Taxonomy) 
tax[c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')] <- colsplit(tax$Taxonomy,';',c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
tax=tax[, -c(3)]
all=merge(otu,tax, by.x= "Group", by.y= "OTU")
all.m=apply(all[,2:89], 2, as.numeric)
all.n=as.data.frame(sweep(all.m,2,colSums(all.m),`/`))
row.names(all.n)=all$Group
all.n$sumNeg=all.n$`neg-PCR`+all.n$`neg-PCR`+all.n$`neg-seawater`
all.n$Conta=(all.n$sumNeg/all.n$Size)*(100)
conta=subset(all.n, all.n$Conta>10) 
conta$Family=tax$Family[match(rownames(conta),tax$OTU)]
noConta=subset(as.data.frame(otu), !Group %in% rownames(conta))
noConta2=noConta[,c(2:85)]# not normalized
rownames(noConta2)=noConta$Group
write.table(noConta2, "./outputs/CBASS84_noConta_raw.txt", quote = FALSE, row.names = FALSE, sep="\t")


noConta.n=subset(as.data.frame(all.n), !rownames(all.n) %in% rownames(conta))
noConta2.n=noConta.n[,c(1:84)]
write.table(noConta2.n, "./outputs/CBASS84_noConta_norm.txt", quote = FALSE, row.names = TRUE, sep="\t")

taxNoConta=subset(tax, OTU %in% rownames(noConta2))
write.table(taxNoConta, "./outputs/CBASS84_noConta_tax", quote = FALSE, row.names = FALSE, sep="\t")
