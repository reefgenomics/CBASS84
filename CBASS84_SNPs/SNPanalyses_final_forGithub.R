#Script written by Daniel Barshis, adapted by Hannah Aichelman, and re-adapted by Daniel Barshis
#Last updated: 05/2021

#This script takes genepop file formats (.gen) and uses principal components analyses and distance clustering to consider population structure and clonality

#load required packages
library(adegenet)
library(poppr)
library(gplots)
library(RColorBrewer)

#set working directory
setwd("/Users/danbarshis/dansstuff/Manuscripts/2021_Voolstra-Valenzuela-Barshis_CBASS84/MolEcology/Reviews/SNPs")

datafile28<-read.genepop('36498_CoralSNPPooled_HEAFilters-maxmissing95_Spisonly_28samples.recode_genepop.gen', ncode=2)

sum(is.na(datafile28$tab))
datafile28 #shows info
YOURdata<-scaleGen(datafile28, NA.method='mean')
X<-YOURdata
##for 36498_28samples
datafile28$pop<-as.factor(gsub("AF7-30","AF",datafile28$pop))
datafile28$pop<-as.factor(gsub("PrT5-30","PrT",datafile28$pop))
datafile28$pop<-as.factor(gsub("ICN5-30","ICN",datafile28$pop))
datafile28$pop<-as.factor(gsub("ExT7-30","ExT",datafile28$pop))

pca1 <- dudi.pca(X,cent=T, scale=T, scannf=F, nf=3)
summary(pca1)
propvar<-format(pca1$eig/sum(pca1$eig)*100, digits=1)

#### PCAs ####
# individual labels
s.label(pca1$li)

# population elipses
s.class(pca1$li, pop(datafile28))

#color symbols, pop names
pdf("FigS3_36498_28samples_ColorPCA1v2.pdf")
Sytes<-c("ICN","AF","ExT","PrT")
datafile28$pop<-factor(datafile28$pop, Sytes)
Colorz<-c('#2c7bb6','#fdae61','#df65b0', '#91003f')
names(Colorz)<-Sytes
Syms<-c(19,17,15,18)
names(Syms)<-Sytes
plot(pca1$li$Axis1, pca1$li$Axis2, col=Colorz[datafile28$pop], pch=Syms[datafile28$pop], xlim=c(min(pca1$li$Axis1),max(pca1$li$Axis1)), ylim=c(min(pca1$li$Axis2)-50,max(pca1$li$Axis2)+50), cex=3,xlab=paste("PC1, ",propvar[1],"% variance explained", sep=""), ylab=paste("PC2, ",propvar[2],"% variance explained)", sep=""), cex.axis=1.5, cex.lab=1.5)
abline(v=0,h=0,col="grey", lty=2)
legend("bottomright", Sytes, col=Colorz, pch=Syms, cex=1.5, pt.cex=3)
dev.off()

datafile79<-read.genepop('36498_CoralSNPPooled_HEAFilters-maxmissing95_Spisonly_79samples.recode_genepop.gen', ncode=2)

sum(is.na(datafile79$tab))
datafile79 #shows info
##for 36498_79samples
datafile79$pop<-as.factor(gsub("AF7-36","AF",datafile79$pop))
datafile79$pop<-as.factor(gsub("PrT4-36","PrT",datafile79$pop))
datafile79$pop<-as.factor(gsub("ICN6-36","ICN",datafile79$pop))
datafile79$pop<-as.factor(gsub("ExT7-36","ExT",datafile79$pop))

####Dissimilarity Distance####
#Create a diss.dist dissimilarity distance matrix to determine a threshhold based on clones, https://www.rdocumentation.org/packages/poppr/versions/2.8.5/topics/diss.dist
#ignores missing data and counts the shared genotypes 
distgenDISS <- diss.dist(datafile79, percent = FALSE, mat = FALSE) 
#make percent = TRUE to get Prevosti distance
#By default, diss.dist() returns a distance reflecting the number of allelic differences between two individuals (Hamming's distance).

distgenDISS2 <- as.matrix(distgenDISS)


pdf("FigS2_36498_79samples_distance-dendrogram.pdf")
#Cluster with hclust and plot
clust_tree<-hclust(distgenDISS, "ave")
plot(clust_tree, cex=0.4)+abline(h=50)
dev.off()
