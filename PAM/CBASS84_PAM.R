setwd("/Users/...")
library(emmeans)
library(drc)

#### ReadIn/QAQC ####
pamdata<-read.delim("../rawdata/CBASS84_PAM_data.txt")
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
#### DRC-IUICoralNursery ####
DRCpamIUICoralNursery = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="IUICoralNursery",],
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamIUICoralNursery)
DRCpamIUICoralNursery$coefficients[3]
ED(DRCpamIUICoralNursery, c(50))[,1]
# genotype-specific curve fits
DRCpamIUICoralNurserygeno = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="IUICoralNursery",], curveid=Geno,
             fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamIUICoralNurserygeno)
DRCpamIUICoralNurserygeno$coefficients[15:21]
ED(DRCpamIUICoralNurserygeno, c(50))[,1]

#### DRC-AlFahal ####
DRCpamAlFahal = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="AlFahal",],
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamAlFahal)
DRCpamAlFahal$coefficients[3]

DRCpamAlFahalgeno = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="AlFahal",], curveid=Geno,
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamAlFahalgeno)
DRCpamAlFahalgeno$coefficients[15:21]
ED(DRCpamAlFahalgeno, c(50))[,1]

#### DRC-Exposed ####
DRCpamExposed = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="Exposed",],
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamExposed)
DRCpamExposed$coefficients[3]

DRCpamExposedgeno = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="Exposed",], curveid=Geno,
                      fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamExposedgeno)
DRCpamExposedgeno$coefficients[15:21]
ED(DRCpamExposedgeno, c(50))[,1]

#### DRC-Protected ####
DRCpamProtected = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="Protected",],
                    fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamProtected)
DRCpamProtected$coefficients[3]

DRCpamProtectedgeno = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="Protected",], curveid=Geno,
                        fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamProtectedgeno)
DRCpamProtectedgeno$coefficients[15:21]
ED(DRCpamProtectedgeno, c(50))[,1]

#### Central Red Sea Only ####
DRCpamCRSgeno = drm(PAM ~ Temp, data = pamdata[pamdata$Site=="AlFahal" | pamdata$Site=="Exposed" | pamdata$Site=="Protected",], curveid=Geno,
                          fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamCRSgeno)
ED(DRCpamCRSgeno, c(50))[,1]
plot(DRCpamCRSgeno)

#### Merging Coeffs ####
Coeffs<-c(DRCpamIUICoralNursery$coefficients[3],DRCpamAlFahal$coefficients[3],DRCpamExposed$coefficients[3],DRCpamProtected$coefficients[3])
GenoCoeffs<-data.frame("ED50"=c(DRCpamIUICoralNurserygeno$coefficients[15:21],DRCpamAlFahalgeno$coefficients[15:21],DRCpamExposedgeno$coefficients[15:21],DRCpamProtectedgeno$coefficients[15:21]),"Site"=c(rep("IUICoralNursery",7), rep("AlFahal",7), rep("Exposed",7), rep("Protected",7)))

aggregate(ED50 ~ Site, data=GenoCoeffs, mean)

Coeffs
CoeffsED25<-c(ED(DRCpamIUICoralNursery,c(25))[,1],ED(DRCpamAlFahal,c(25))[,1],ED(DRCpamExposed,c(25))[,1],ED(DRCpamProtected,c(25))[,1])

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
write.table(data.frame("Genet"=gsub("ed50:", "", row.names(GenoCoeffs)), GenoCoeffs), file="../rawdata/GenetED50s.txt", quote=F, sep="\t", row.names=F)
write.table(SumaryStats, file="../rawdata/ED50SummaryStats.txt", quote=F, sep="\t", row.names=F)

#### Plotting ####
Sytes<-c("IUICoralNursery","AlFahal","Exposed","Protected")

Colorz<-c('#2c7bb6','#fdae61','#df65b0', '#91003f')
Syms<-c(19,17,15,18)
temp_x<- seq(30, 39, length = 100)

MMMs<-c(27.01, 30.76, 30.74, 30.75)

pdf("../CBASS84_DRCCurves_MMMs.pdf",10,7)
line_width=2
offsets<-c(0.1875,0.0625,-0.0625,-0.1875)
i<-1 #IUICoralNursery
matplot(temp_x, predict(DRCpamIUICoralNursery, data.frame(Temp = temp_x), interval="confidence"),
        type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width,ylab="Fv/Fm",xlab="Temperature °C", xlim=c(29.5,39.5),ylim=c(0,0.65), cex.axis=1.5, cex.lab=1.5)
with(pamdata[pamdata$Site==Sytes[i],],matpoints(Temp-offsets[i],PAM,pch=Syms[i], col=Colorz[i], cex=1.5))
i<-2 #AlFahal
matpoints(temp_x, predict(DRCpamAlFahal, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width)
with(pamdata[pamdata$Site==Sytes[i],],matpoints(Temp-offsets[i],PAM,pch=Syms[i], col=Colorz[i], cex=1.5))
i<-3 # Exposed
matpoints(temp_x, predict(DRCpamExposed, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width)
with(pamdata[pamdata$Site==Sytes[i],],matpoints(Temp-offsets[i],PAM,pch=Syms[i], col=Colorz[i], cex=1.5))
i<-4 # Protected
matpoints(temp_x, predict(DRCpamProtected, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width)
with(pamdata[pamdata$Site==Sytes[i],],matpoints(Temp-offsets[i],PAM,pch=Syms[i], col=Colorz[i], cex=1.5))

legend("bottomleft",c("IUICoralNursery: A, 27.01","AlFahal: B, 30.76","Exposed: BC, 30.74","Protected: C, 30.75"),pch=Syms, col=Colorz,pt.cex=2, bty="n",cex=1.5)
title(main="Stylophora pistillata Fv/Fm")
abline(v=SumaryStats$MeanED50, col=Colorz)
text(SumaryStats$MeanED50+0.45, c(0.65,0.65,0.62,0.65),labels=as.character(round(SumaryStats$MeanED50, digits=2)),col=Colorz, cex=1.5)
dev.off()

