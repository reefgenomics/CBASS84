#original script by Daniel J Barshis
#modified by Christian R Voolstra and Daniel J Barshis
library(emmeans)
library(drc)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


#### ReadIn/QAQC ####
pamdata<-read.delim(file.choose(),header=T) 
pamdata$PAM<-pamdata$PAM/1000
pamdata$Geno<-as.factor(pamdata$Geno)
pamdata$SiteGeno<-paste0(pamdata$Site,pamdata$Geno)
pamdata$SiteTemp<-paste0(pamdata$Site,pamdata$Temp)
pamdata$SiteTempRep<-paste0(pamdata$Temp,pamdata$Site,pamdata$Replicate)
pamdata$FacTemp<-as.factor(pamdata$Temp)
aggregate(PAM ~ Site + Temp, data=pamdata, summary)
# Checking sample sizes
aggregate(PAM ~ Site + Temp, data=pamdata, length)
aggregate(PAM ~ Site + FacTemp, data=pamdata, length)
aggregate(PAM ~ Site + FacTemp + Geno, data=pamdata, length)

#### DRC Curve Fitting####
#getMeanFunctions()
DRCpam = drm(PAM ~ Temp, data = pamdata, curveid = Site,
              fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpam)
compParm(DRCpam, 'ed50')
compParm(DRCpam, 'ed50', "-")
plot(DRCpam)
points(pamdata$Temp, pamdata$PAM)
ED(DRCpam, c(50))[,1]

# fit to each site individually

#### DRC-ICN ####
DRCpamICN = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="ICN",],
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamICN)
DRCpamICN$coefficients[3]
ED(DRCpamICN, c(50))[,1]
# genotype-specific curve fits
DRCpamICNgeno = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="ICN",], curveid=Geno,
             fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamICNgeno)
DRCpamICNgeno$coefficients[15:21]
ED(DRCpamICNgeno, c(50))[,1]

#### DRC-AF ####
DRCpamAF = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="AF",],
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamAF)
DRCpamAF$coefficients[3]
# genotype-specific curve fits
DRCpamAFgeno = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="AF",], curveid=Geno,
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamAFgeno)
DRCpamAFgeno$coefficients[15:21]
ED(DRCpamAFgeno, c(50))[,1]

#### DRC-ExT ####
DRCpamExT = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="ExT",],
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamExT)
DRCpamExT$coefficients[3]
# genotype-specific curve fits
DRCpamExTgeno = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="ExT",], curveid=Geno,
                      fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamExTgeno)
DRCpamExTgeno$coefficients[15:21]
ED(DRCpamExTgeno, c(50))[,1]

#### DRC-PrT ####
DRCpamPrT = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="PrT",],
                    fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamPrT)
DRCpamPrT$coefficients[3]
# genotype-specific curve fits
DRCpamPrTgeno = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="PrT",], curveid=Geno,
                        fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamPrTgeno)
DRCpamPrTgeno$coefficients[15:21]
ED(DRCpamPrTgeno, c(50))[,1]

#### Merging Coeffs ####
Coeffs<-c(DRCpamICN$coefficients[3],DRCpamAF$coefficients[3],DRCpamExT$coefficients[3],DRCpamPrT$coefficients[3])
GenoCoeffs<-data.frame("ED50"=c(DRCpamICNgeno$coefficients[15:21],DRCpamAFgeno$coefficients[15:21],DRCpamExTgeno$coefficients[15:21],DRCpamPrTgeno$coefficients[15:21]),"Site"=c(rep("ICN",7), rep("AF",7), rep("ExT",7), rep("PrT",7)))

aggregate(ED50 ~ Site, data=GenoCoeffs, mean)

Coeffs
CoeffsED25<-c(ED(DRCpamICN,c(25))[,1],ED(DRCpamAF,c(25))[,1],ED(DRCpamExT,c(25))[,1],ED(DRCpamPrT,c(25))[,1])

#### ED50 ANoVA ####
aggregate(ED50 ~ Site, data=GenoCoeffs, FUN= function(x) shapiro.test(x)$p.value)
bartlett.test(ED50 ~ Site, data=GenoCoeffs)

ED50.aov<-aov(ED50 ~ Site, GenoCoeffs)
summary(ED50.aov)
TukeyHSD(ED50.aov)

aggregate(ED50 ~ Site, data=GenoCoeffs, summary)
GenoCoeffs
SumaryStats<-data.frame("Site"=names(tapply(GenoCoeffs$ED50,GenoCoeffs$Site, mean)), "MeanED50"=tapply(GenoCoeffs$ED50,GenoCoeffs$Site, mean), "ED50StdDev"=tapply(GenoCoeffs$ED50,GenoCoeffs$Site, sd), "ED50StdErr"=tapply(GenoCoeffs$ED50,GenoCoeffs$Site, sd)/sqrt(tapply(GenoCoeffs$ED50,GenoCoeffs$Site, length)))
SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]
write.table(data.frame("Genet"=gsub("ed50:", "", row.names(GenoCoeffs)), GenoCoeffs), file="./GenetED50s.txt", quote=F, sep="\t", row.names=F)
write.table(SumaryStats, file="./CBASS84_ED50_SummaryStats.txt", quote=F, sep="\t", row.names=F)

#### Plotting ####
Sytes<-c("ICN","AF","ExT","PrT")

Colorz<-c('#2c7bb6','#fdae61','#df65b0', '#91003f')
Syms<-c(19,17,15,18)
temp_x<- seq(30, 39, length = 100)

MMMs<-c(27.01, 30.76, 30.74, 30.75)

pdf("../CBASS84_DRCCurves_MMMs.pdf",10,7)
line_width=2
offsets<-c(0.1875,0.0625,-0.0625,-0.1875)
Denscity <- 45

i<-4 # PrT
matplot(temp_x, predict(DRCpamPrT, data.frame(Temp = temp_x), interval="confidence"),
        type="n",col=Colorz[i],lty=c(1,3,3),lwd=line_width,ylab="Fv/Fm",xlab="Temperature °C", xlim=c(29.5,39.5),ylim=c(0,0.65), cex.axis=1.5, cex.lab=1.5)
polygon(c(temp_x, rev(temp_x)),c(predict(DRCpamPrT, data.frame(Temp = temp_x), interval="confidence")[,2],rev(predict(DRCpamPrT, data.frame(Temp = temp_x), interval="confidence")[,3])), col=Colorz[i], density = Denscity)
matpoints(temp_x, predict(DRCpamPrT, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width)
with(pamdata[pamdata$Site==Sytes[i],],matpoints(Temp-offsets[i],PAM,pch=Syms[i], col=Colorz[i], cex=1.5))
i<-3 # ExT
polygon(c(temp_x, rev(temp_x)),c(predict(DRCpamExT, data.frame(Temp = temp_x), interval="confidence")[,2],rev(predict(DRCpamExT, data.frame(Temp = temp_x), interval="confidence")[,3])), col=Colorz[i], density = Denscity)
matpoints(temp_x, predict(DRCpamExT, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width)
with(pamdata[pamdata$Site==Sytes[i],],matpoints(Temp-offsets[i],PAM,pch=Syms[i], col=Colorz[i], cex=1.5))
i<-2 #AF
polygon(c(temp_x, rev(temp_x)),c(predict(DRCpamAF, data.frame(Temp = temp_x), interval="confidence")[,2],rev(predict(DRCpamAF, data.frame(Temp = temp_x), interval="confidence")[,3])), col=Colorz[i], density = Denscity)
matpoints(temp_x, predict(DRCpamAF, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width)
with(pamdata[pamdata$Site==Sytes[i],],matpoints(Temp-offsets[i],PAM,pch=Syms[i], col=Colorz[i], cex=1.5))
i<-1 #ICN
polygon(c(temp_x, rev(temp_x)),c(predict(DRCpamICN, data.frame(Temp = temp_x), interval="confidence")[,2],rev(predict(DRCpamICN, data.frame(Temp = temp_x), interval="confidence")[,3])), col=Colorz[i], density = Denscity)
matpoints(temp_x, predict(DRCpamICN, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width)
with(pamdata[pamdata$Site==Sytes[i],],matpoints(Temp-offsets[i],PAM,pch=Syms[i], col=Colorz[i], cex=1.5))

legend("bottomleft",c("ICN: A, 27.01","AF: B, 30.76","ExT: BC, 30.74","PrT: C, 30.75"),pch=Syms, col=Colorz,pt.cex=2, bty="n",cex=1.5)
title(main="Stylophora pistillata Fv/Fm")
abline(v=SumaryStats$MeanED50, col=Colorz)
text(SumaryStats$MeanED50+0.45, c(0.65,0.65,0.62,0.65),labels=as.character(round(SumaryStats$MeanED50, digits=2)),col=Colorz, cex=1.5)
dev.off()

