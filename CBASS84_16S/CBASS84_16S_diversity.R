library(ggplot2)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/CBASS84/")

## data in
counts=read.table("./outputs/CBASS84_noConta_raw.txt", header = TRUE, row.names = 1)
tax=read.table("./outputs/CBASS84_noConta_tax", header = TRUE,sep = "\t")
counts$Family=tax$Family[match(rownames(counts),tax$OTU)]
fam.wid.agg=aggregate(counts[, 1:84], by = list(counts$Family), FUN =  sum)#defined sample range and group factor

## identify most abundant families and aggregate the rest to "Others"
topFamilies=fam.wid.agg[order(rowSums(fam.wid.agg[, 2:ncol(fam.wid.agg)]),decreasing = TRUE),][1:20,1]
fam.top=subset(fam.wid.agg, fam.wid.agg$Group.1 %in% topFamilies)
fam.bot=subset(fam.wid.agg, !fam.wid.agg$Group.1 %in% topFamilies)
fam.bot$Group.1=gsub(".*","zOthers", fam.bot$Group.1)
others=aggregate(fam.bot[, 2:ncol(fam.bot)], by = list(fam.bot[, 1]), FUN =  sum)
all.2 =rbind(fam.top, others)
all.l=melt(all.2, id.vars=c("Group.1"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("Family","Sample","Abundance")

## Add metadata to plot
met=read.table("./Input_files/metadata.txt", header = TRUE)
met$ID=gsub('-', '\\.', met$ID)
all.l$Genotype=met$Genotype[match(all.l$Sample,met$ID)]
all.l$Site=met$Region[match(all.l$Sample,met$ID)]
all.l$Comparison=met$Comparison[match(all.l$Sample,met$ID)]
all.l$Temperature=met$Temperature[match(all.l$Sample,met$ID)]
all.l$Site=factor(all.l$Site,levels=c("Eilat", "Kaust", "Exposed","Protected"))
all.l2=all.l %>% group_by(Family,  Site, Comparison, Temperature) %>% summarise(Abundance=sum(Abundance))
all.l2$Site=factor(all.l2$Site,levels=c("Eilat", "Kaust", "Exposed","Protected"))

## Plot bar plots
P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
pdf(file = "./ouCBASS84_barplots.pdf",  width = 10, height = 5, pointsize = 12)
ggplot() +geom_bar(aes(y = Abundance, x = Genotype, fill = Family), data = all.l, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA reads", x="Samples") + scale_fill_manual(values=P21) + facet_grid(Temperature~Site) + theme(legend.key = element_blank(), strip.background = element_blank()) + guides(fill=guide_legend(ncol=1))
dev.off()

P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
pdf(file = "./outputs/CBASS84_averagedBarplots.pdf",  width = 7, height = 5, pointsize = 12) 
ggplot() +geom_bar(aes(y = Abundance, x = Site, fill = Family), data = all.l2, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA reads", x="") + scale_fill_manual(values=P21) + facet_grid(Temperature~.) + theme(legend.key = element_blank(), strip.background = element_blank()) + guides(fill=guide_legend(ncol=1))
dev.off()

