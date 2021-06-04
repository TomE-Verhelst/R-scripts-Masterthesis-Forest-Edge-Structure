library(lmerTest)
library(sjPlot)
library(ggpubr)
library(lme4)
library(matrixStats)
library(tree)
library(lmtest)
library(ggplot2)
library(Hmisc)
library(corrplot)
library(vegan)
library(ggvegan)
install.packages("devtools")
devtools::install_github("gavinsimpson/ggvegan")

detach(AIs)
AIs <- read.table("dawnmatrix_allindices.txt",header = TRUE)
attach(AIs)
detach(AIs)
AIs <- read.table("duskmatrix_allindices.txt",header = TRUE)
attach(AIs)

loc <- as.factor(Location)
reg <- as.factor(Region)
struc <- as.factor(Structure)
p <- as.factor(Plot)

## best fit ACI ##
fitACIm <- lm(ACIm ~ loc+p+Tmax+loc:p+loc:Tmax+p:Tmax)
anova(fitACIm)
fitACIm <- lm(ACIm ~ loc + p + Tmax + loc:p)
summary(fitACIm)
TukeyHSD(aov(fitACIm))
plot_model(fitACIm,show.values = TRUE,value.offset = .4,rm.terms = c('Tmax'))

fitACIe <- lm(ACIe ~ loc+p+Tmax+loc:p+loc:Tmax+p:Tmax)
anova(fitACIe)
fitACIe <- lm(ACIe ~ loc + p)
summary(fitACIe)

## best fit AEI ##
fitAEIm <- lm(AEIm ~ loc+p+Tmax+loc:p+loc:Tmax+p:Tmax)
anova(fitAEIm)
fitAEIm <- lm(AEIm ~ loc + p)
anova(fitAEIm)
summary(fitAEIm)


fitAEIe <- lm(AEIe ~ loc+p+Tmax+loc:p+loc:Tmax+p:Tmax)
anova(fitAEIe)
fitAEIe <- lm(AEIe ~ loc + loc:p)
summary(fitAEIe)

## best fit BI ##
fitBIm <- lm(BIm ~ loc+p+Tmax+loc:p+loc:Tmax+p:Tmax)
anova(fitBIm)
fitBIm <- lm(BIm ~ loc + Tmax)
anova(fitBIm)
summary(fitBIm)

fitBIe <- lm(BIe ~ loc+p+Tmax+loc:p+loc:Tmax+p:Tmax)
anova(fitBIe)
fitBIe <- lm(BIe ~ loc + Tmax + loc:p)
summary(fitBIe)

## structure effect on P1 (1 or 3) ##
detach(AIs)
AIs <- read.table("finalmatrix_P1_dawn.txt",header = TRUE)
attach(AIs)
loc <- as.factor(Location)
reg <- as.factor(Region)
struc <- as.factor(Structure)
p <- as.factor(Plot)

fitACI <- lmer(ACIm ~ struc + Tmax + (1+Tmax|reg:struc))
anova(fitACI)
summary(fitACI)

fitAEI <- lmer(AEIm ~ struc + Tmax + (1+Tmax|reg:struc))
anova(fitAEI)
summary(fitAEI)

fitBI <- lmer(BIm ~ struc + Tmax + (1+Tmax|reg:struc))
anova(fitBI)
summary(fitBI)

plot_model(fitACI)

plot(Tmax,BIm)
abline(fitBI[struc=='1'],lty=1)
abline(fitBI)

#structure effect on P5 (1 or 3)
detach(AIs)
AIs <- read.table("finalmatrix_P5_dawn.txt",header = TRUE)
attach(AIs)
loc <- as.factor(Location)
reg <- as.factor(Region)
struc <- as.factor(Structure)
p <- as.factor(Plot)

fitACI5 <- lmer(ACIm5 ~ struc + Tmax + (1+Tmax|reg:struc))
anova(fitACI5)
summary(fitACI5)

fitAEI5 <- lmer(AEIm5 ~ struc + Tmax + (1+Tmax|reg:struc))
anova(fitAEI5)
summary(fitAEI5)

fitBI5 <- lmer(BIm5 ~ struc + Tmax + (1+Tmax|reg:struc))
anova(fitBI5)
summary(fitBI5)


test <- lmer(ACIm ~ Cgap + Dist075 + Vheight + Tmax + (1|reg))
anova(test)
summary(test)
plot(Cheight,ACIm)
abline(test2)
test2 <- lm(ACIm ~ Cheight)
anova(test2)
summary(test2)

## final analysis ##
AIs <- read.table("finalmatrix_P1_dawn_noM3.txt",header = TRUE)
attach(AIs)
AIs5 <- read.table("finalmatrix_P5_dawn_noM3.txt",header = TRUE)
attach(AIs5)

reg <- lm(ACIm ~ Cs + Spla)
summary(reg)
## ACI ##
par(mfrow=c(4,3),cex.lab=1.2,font.lab=2)
plot(Cgap,ACIm,xlab='Canopy gap fraction')
points(Cgap,ACIm5,col='red')
abline(lm(ACIm ~ Cgap))
abline(lm(ACIm5 ~ Cgap),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Cgap))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Cgap))$r.squared,digits=3)))

plot(Cs,ACIm,xlab='Canopy slope (m/m)')
points(Cs,ACIm5,col='red')
abline(lm(ACIm ~ Cs))
abline(lm(ACIm5 ~ Cs),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Cs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Cs))$r.squared,digits=3)))

plot(Sshs,ACIm,xlab='High shrub total PC slope (1/m)')
points(Sshs,ACIm5,col='red')
abline(lm(ACIm ~ Sshs))
abline(lm(ACIm5 ~ Sshs),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Sshs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Sshs))$r.squared,digits=3)))

plot(Egap,ACIm,xlab='Edge gap fraction')
points(Egap,ACIm5,col='red')
abline(lm(ACIm ~ Egap))
abline(lm(ACIm5 ~ Egap),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Egap))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Egap))$r.squared,digits=3)))

plot(Vs,ACIm,xlab='Vegetation slope (m/m)')
points(Vs,ACIm5,col='red')
abline(lm(ACIm ~ Vs))
abline(lm(ACIm5 ~ Vs),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Vs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Vs))$r.squared,digits=3)))

plot(Ssla,ACIm,xlab='Low arboreal total PC slope (1/m)')
points(Ssla,ACIm5,col='red')
abline(lm(ACIm ~ Ssla))
abline(lm(ACIm5 ~ Ssla),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Ssla))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Ssla))$r.squared,digits=3)))

plot(Ch,ACIm,xlab='Canopy height (m)')
points(Ch,ACIm5,col='red')
abline(lm(ACIm ~ Ch))
abline(lm(ACIm5 ~ Ch),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Ch))$r.squared, digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Ch))$r.squared,digits=3)))

plot(Cs075,ACIm,xlab='Canopy slope 3/4 (m/m)')
points(Cs075,ACIm5,col='red')
abline(lm(ACIm ~ Cs075))
abline(lm(ACIm5 ~ Cs075),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Cs075))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Cs075))$r.squared,digits=3)))

plot(Sphs,ACIm,xlab='High shrub proportion slope (1/m)')
points(Sphs,ACIm5,col='red')
abline(lm(ACIm ~ Sphs))
abline(lm(ACIm5 ~ Sphs),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Sphs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Sphs))$r.squared,digits=3)))

plot(Vh,ACIm,xlab='Vegetation height (m)')
points(Vh,ACIm5,col='red')
abline(lm(ACIm ~ Vh))
abline(lm(ACIm5 ~ Vh),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Vh))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Vh))$r.squared,digits=3)))

plot(Vs075,ACIm,xlab='Vegetation slope 3/4 (m/m)')
points(Vs075,ACIm5,col='red')
abline(lm(ACIm ~ Vs075))
abline(lm(ACIm5 ~ Vs075),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Vs075))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Vs075))$r.squared,digits=3)))

plot(Spla,ACIm,xlab='Low arboreal proportion slope (1/m)')
points(Spla,ACIm5,col='red')
abline(lm(ACIm ~ Spla))
abline(lm(ACIm5 ~ Spla),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(ACIm ~ Spla))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(ACIm5 ~ Spla))$r.squared,digits=3)))

## AEI ##
par(mfrow=c(3,4),cex.lab=1.2,font.lab=2)
plot(Cgap,AEIm,xlab='Canopy gap fraction')
points(Cgap,AEIm5,col='red')
abline(lm(AEIm ~ Cgap))
abline(lm(AEIm5 ~ Cgap),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Cgap))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Cgap))$r.squared,digits=3)))

plot(Egap,AEIm,xlab='Edge gap fraction')
points(Egap,AEIm5,col='red')
abline(lm(AEIm ~ Egap))
abline(lm(AEIm5 ~ Egap),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Egap))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Egap))$r.squared,digits=3)))

plot(Ch,AEIm,xlab='Canopy height (m)')
points(Ch,AEIm5,col='red')
abline(lm(AEIm ~ Ch))
abline(lm(AEIm5 ~ Ch),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Ch))$r.squared, digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Ch))$r.squared,digits=3)))

plot(Vh,AEIm,xlab='Vegetation height (m)')
points(Vh,AEIm5,col='red')
abline(lm(AEIm ~ Vh))
abline(lm(AEIm5 ~ Vh),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Vh))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Vh))$r.squared,digits=3)))

plot(Cs,AEIm,xlab='Canopy slope (m/m)')
points(Cs,AEIm5,col='red')
abline(lm(AEIm ~ Cs))
abline(lm(AEIm5 ~ Cs),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Cs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Cs))$r.squared,digits=3)))

plot(Vs,AEIm,xlab='Vegetation slope (m/m)')
points(Vs,AEIm5,col='red')
abline(lm(AEIm ~ Vs))
abline(lm(AEIm5 ~ Vs),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Vs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Vs))$r.squared,digits=3)))

plot(Cs075,AEIm,xlab='Canopy slope 3/4 (m/m)')
points(Cs075,AEIm5,col='red')
abline(lm(AEIm ~ Cs075))
abline(lm(AEIm5 ~ Cs075),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Cs075))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Cs075))$r.squared,digits=3)))

plot(Vs075,AEIm,xlab='Vegetation slope 3/4 (m/m)')
points(Vs075,AEIm5,col='red')
abline(lm(AEIm ~ Vs075))
abline(lm(AEIm5 ~ Vs075),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Vs075))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Vs075))$r.squared,digits=3)))

plot(Sshs,AEIm,xlab='High shrub total PC slope (1/m)')
points(Sshs,AEIm5,col='red')
abline(lm(AEIm ~ Sshs))
abline(lm(AEIm5 ~ Sshs),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Sshs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Sshs))$r.squared,digits=3)))

plot(Ssla,AEIm,xlab='Low arboreal total PC slope (1/m)')
points(Ssla,AEIm5,col='red')
abline(lm(AEIm ~ Ssla))
abline(lm(AEIm5 ~ Ssla),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Ssla))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Ssla))$r.squared,digits=3)))

plot(Sphs,AEIm,xlab='High shrub proportion slope (1/m)')
points(Sphs,AEIm5,col='red')
abline(lm(AEIm ~ Sphs))
abline(lm(AEIm5 ~ Sphs),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Sphs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Sphs))$r.squared,digits=3)))

plot(Spla,AEIm,xlab='Low arboreal proportion slope (1/m)')
points(Spla,AEIm5,col='red')
abline(lm(AEIm ~ Spla))
abline(lm(AEIm5 ~ Spla),col='red')
legend("topright", bty="n", 
        legend=paste(format(summary(lm(AEIm ~ Spla))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(AEIm5 ~ Spla))$r.squared,digits=3)))

## BI ##
par(mfrow=c(3,4),cex.lab=1.2,font.lab=2)
plot(Cgap,BIm,xlab='Canopy gap fraction')
points(Cgap,BIm5,col='red')
abline(lm(BIm ~ Cgap))
abline(lm(BIm5 ~ Cgap),col='red')
legend("topright", bty="n",cex=1.2, 
        legend=paste(format(summary(lm(BIm ~ Cgap))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Cgap))$r.squared,digits=3)))


plot(Egap,BIm,xlab='Edge gap fraction')
points(Egap,BIm5,col='red')
abline(lm(BIm ~ Egap))
abline(lm(BIm5 ~ Egap),col='red')
legend("topright", bty="n",cex=1.2,
        legend=paste(format(summary(lm(BIm ~ Egap))$r.squared,digits=3)))
text(0.25,11,format(summary(lm(BIm ~ Spla))$r.squared,digits=3),col='red')
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Egap))$r.squared,digits=3)))

plot(Ch,BIm,xlab='Canopy height (m)')
points(Ch,BIm5,col='red')
abline(lm(BIm ~ Ch))
abline(lm(BIm5 ~ Ch),col='red')
legend("topright", bty="n",cex=1.2, 
        legend=paste(format(summary(lm(BIm ~ Ch))$r.squared, digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Ch))$r.squared,digits=3)))

plot(Vh,BIm,xlab='Vegetation height (m)')
points(Vh,BIm5,col='red')
abline(lm(BIm ~ Vh))
abline(lm(BIm5 ~ Vh),col='red')
legend("topright", bty="n",cex=1.2, 
        legend=paste(format(summary(lm(BIm ~ Vh))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Vh))$r.squared,digits=3)))

plot(Cs,BIm,xlab='Canopy slope (m/m)')
points(Cs,BIm5,col='red')
abline(lm(BIm ~ Cs))
abline(lm(BIm5 ~ Cs),col='red')
legend("topright", bty="n",cex=1.2, 
        legend=paste(format(summary(lm(BIm ~ Cs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Cs))$r.squared,digits=3)))

plot(Vs,BIm,xlab='Vegetation slope (m/m)')
points(Vs,BIm5,col='red')
abline(lm(BIm ~ Vs))
abline(lm(BIm5 ~ Vs),col='red')
legend("topright", bty="n",cex=1.2, 
        legend=paste(format(summary(lm(BIm ~ Vs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Vs))$r.squared,digits=3)))

plot(Cs075,BIm,xlab='Canopy slope 3/4 (m/m)')
points(Cs075,BIm5,col='red')
abline(lm(BIm ~ Cs075))
abline(lm(BIm5 ~ Cs075),col='red')
legend("topright", bty="n",cex=1.2, 
        legend=paste(format(summary(lm(BIm ~ Cs075))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Cs075))$r.squared,digits=3)))

plot(Vs075,BIm,xlab='Vegetation slope 3/4 (m/m)')
points(Vs075,BIm5,col='red')
abline(lm(BIm ~ Vs075))
abline(lm(BIm5 ~ Vs075),col='red')
legend("topright", bty="n",cex=1.2, 
        legend=paste(format(summary(lm(BIm ~ Vs075))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Vs075))$r.squared,digits=3)))

plot(Sshs,BIm,xlab='High shrub total PC slope (1/m)')
points(Sshs,BIm5,col='red')
abline(lm(BIm ~ Sshs))
abline(lm(BIm5 ~ Sshs),col='red')
legend("topright", bty="n",cex=1.2, 
        legend=paste(format(summary(lm(BIm ~ Sshs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Sshs))$r.squared,digits=3)))

plot(Ssla,BIm,xlab='Low arboreal total PC slope (1/m)')
points(Ssla,BIm5,col='red')
abline(lm(BIm ~ Ssla))
abline(lm(BIm5 ~ Ssla),col='red')
legend("topright", bty="n",cex=1.2, 
        legend=paste(format(summary(lm(BIm ~ Ssla))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Ssla))$r.squared,digits=3)))

plot(Sphs,BIm,xlab='High shrub proportion slope (1/m)')
points(Sphs,BIm5,col='red')
abline(lm(BIm ~ Sphs))
abline(lm(BIm5 ~ Sphs),col='red')
legend("topright", bty="n",cex=1.2, 
        legend=paste(format(summary(lm(BIm ~ Sphs))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Sphs))$r.squared,digits=3)))

plot(Spla,BIm,xlab='Low arboreal proportion slope (1/m)')
points(Spla,BIm5,col='red')
abline(lm(BIm ~ Spla))
abline(lm(BIm5 ~ Spla),col='red')
legend("topright", bty="n",cex=1.2, 
        legend=paste(format(summary(lm(BIm ~ Spla))$r.squared,digits=3)))
legend("bottomright", bty="n",col=c('red'),cex=1.2,
        legend=paste(format(summary(lm(BIm5 ~ Spla))$r.squared,digits=3)))

summary(lm(BIm5 ~ Cgap))

## boxplots ##
AIs <- read.table("finalmatrix_P1_dawn.txt",header = TRUE)
attach(AIs)
AIs5 <- read.table("finalmatrix_P5_dawn.txt",header = TRUE)
attach(AIs5)

par(mfrow=c(4,3),cex.lab=1.2,font.lab=2,mai=c(0.75,0.75,0.1,0.1))
boxplot(ACIm ~ Cgap,xlab='Canopy gap fraction',ylab='ACI')
boxplot(ACIm ~ Cs,xlab='Canopy slope (m/m)',ylab='ACI')
boxplot(ACIm ~ Sshs,xlab='TPPC hs slope (1/m)',ylab='ACI')
boxplot(ACIm ~ Egap,xlab='Edge gap fraction')
boxplot(ACIm ~ Vs,xlab='Vegetation slope (m/m)',ylab='ACI')
boxplot(ACIm ~ Ssla,xlab='TPPC la slope (1/m)',ylab='ACI')
boxplot(ACIm ~ Ch,xlab='Canopy height (m)')
boxplot(ACIm ~ Cs075,xlab='Canopy slope 3/4 (m/m)',ylab='ACI')
boxplot(ACIm ~ Sphs,xlab='PCP hs slope (1/m)')
boxplot(ACIm ~ Vh,xlab='Vegetation height (m)')
boxplot(ACIm ~ Vs075,xlab='Vegetation slope 3/4 (m/m)',ylab='ACI')
boxplot(ACIm ~ Spla,xlab='PCP la slope (1/m)',ylab='ACI')

par(mfrow=c(4,3),cex.lab=1.2,font.lab=2)
boxplot(AEIm ~ Cgap,xlab='Canopy gap fraction')
boxplot(AEIm ~ Cs,xlab='Canopy slope (m/m)')
boxplot(AEIm ~ Sshs,xlab='TPPC hs slope (1/m)')
boxplot(AEIm ~ Egap,xlab='Edge gap fraction')
boxplot(AEIm ~ Vs,xlab='Vegetation slope (m/m)')
boxplot(AEIm ~ Ssla,xlab='TPPC la slope (1/m)')
boxplot(AEIm ~ Ch,xlab='Canopy height (m)')
boxplot(AEIm ~ Cs075,xlab='Canopy slope 3/4 (m/m)')
boxplot(AEIm ~ Sphs,xlab='PCP hs slope (1/m)')
boxplot(AEIm ~ Vh,xlab='Vegetation height (m)')
boxplot(AEIm ~ Vs075,xlab='Vegetation slope 3/4 (m/m)')
boxplot(AEIm ~ Spla,xlab='PCP la slope (1/m)')

par(mfrow=c(4,3),cex.lab=1.2,font.lab=2)
boxplot(BIm ~ Cgap,xlab='Canopy gap fraction')
boxplot(BIm ~ Cs,xlab='Canopy slope (m/m)')
boxplot(BIm ~ Sshs,xlab='TPPC hs slope (1/m)')
boxplot(BIm ~ Egap,xlab='Edge gap fraction')
boxplot(BIm ~ Vs,xlab='Vegetation slope (m/m)')
boxplot(BIm ~ Ssla,xlab='TPPC la slope (1/m)')
boxplot(BIm ~ Ch,xlab='Canopy height (m)')
boxplot(BIm ~ Cs075,xlab='Canopy slope 3/4 (m/m)')
boxplot(BIm ~ Sphs,xlab='PCP hs slope (1/m)')
boxplot(BIm ~ Vh,xlab='Vegetation height (m)')
boxplot(BIm ~ Vs075,xlab='Vegetation slope 3/4 (m/m)')
boxplot(BIm ~ Spla,xlab='PCP la slope (1/m)')

par(mfrow=c(3,3),cex.lab=1.2,font.lab=2)
boxplot(ACIm ~ Location,xlab='Transect',ylab='ACI plot 1')
boxplot(ACIm5 ~ Location,xlab='Transect',ylab='ACI plot 5')
plot(Tmax,ACIm5,xlab='Maximum day temperature (°C)',ylab='ACI plot 5')
boxplot(AEIm ~ Location,xlab='Transect',ylab='AEI plot 1')
boxplot(AEIm5 ~ Location,xlab='Transect',ylab='AEI plot 5')
plot(Tmax,AEIm5,xlab='Maximum day temperature (°C)',ylab='AEI plot 5')
boxplot(BIm ~ Location,xlab='Transect',ylab='BI plot 1')
boxplot(BIm5 ~ Location,xlab='Transect',ylab='BI plot 5')
plot(Tmax,BIm,xlab='Maximum day temperature (°C)',ylab='BI plot 1')

par(mfrow=c(3,2),cex.lab=1.2,font.lab=2)
boxplot(ACIm ~ Location,xlab='Transect',ylab='ACI plot 1')
boxplot(ACIm5 ~ Location,xlab='Transect',ylab='ACI plot 5')
boxplot(AEIm ~ Location,xlab='Transect',ylab='AEI plot 1')
boxplot(AEIm5 ~ Location,xlab='Transect',ylab='AEI plot 5')
boxplot(BIm ~ Location,xlab='Transect',ylab='BI plot 1')
boxplot(BIm5 ~ Location,xlab='Transect',ylab='BI plot 5')

## boxplots P1 en 5 ##

fin <- read.table("finaal_P1en5.txt",header = TRUE)
attach(fin)
require(ggplot2)
a <- ggplot(data = fin, aes(x=as.factor(Cgap), y=ACIf)) + geom_boxplot(aes(fill=Plot))
a1 <- a + labs(x='Canopy gap fraction',y='ACI',col='green')
b <- ggplot(data = fin, aes(x=as.factor(Cs), y=ACIf)) + geom_boxplot(aes(fill=Plot))
b1 <- b + labs(x='Canopy slope (m/m)',y='ACI')
c <- ggplot(data = fin, aes(x=as.factor(Sshs), y=ACIf)) + geom_boxplot(aes(fill=Plot))
c1 <- c + labs(x='TPPC hs slope (1/m)',y='ACI')
d <- ggplot(data = fin, aes(x=as.factor(Egap), y=ACIf)) + geom_boxplot(aes(fill=Plot))
d1 <- d + labs(x='Edge gap fraction',y='ACI')
e <- ggplot(data = fin, aes(x=as.factor(Vs), y=ACIf)) + geom_boxplot(aes(fill=Plot))
e1 <- e + labs(x='Vegetation slope (m/m)',y='ACI')
f <- ggplot(data = fin, aes(x=as.factor(Ssla), y=ACIf)) + geom_boxplot(aes(fill=Plot))
f1 <- f + labs(x='TPPC la slope (1/m)',y='ACI')
g <- ggplot(data = fin, aes(x=as.factor(Ch), y=ACIf)) + geom_boxplot(aes(fill=Plot))
g1 <- g + labs(x='Canopy height (m)',y='ACI')
h <- ggplot(data = fin, aes(x=as.factor(Cs075), y=ACIf)) + geom_boxplot(aes(fill=Plot))
h1 <- h + labs(x='Canopy slope 3/4 (m/m)',y='ACI')
i <- ggplot(data = fin, aes(x=as.factor(Sphs), y=ACIf)) + geom_boxplot(aes(fill=Plot))
i1 <- i + labs(x='PCP hs slope (1/m)',y='ACI')
j <- ggplot(data = fin, aes(x=as.factor(Vh), y=ACIf)) + geom_boxplot(aes(fill=Plot))
j1 <- j + labs(x='Vegetation height (m)',y='ACI')
k <- ggplot(data = fin, aes(x=as.factor(Vs075), y=ACIf)) + geom_boxplot(aes(fill=Plot))
k1 <- k + labs(x='Vegetation slope 3/4 (m/m)',y='ACI')
l <- ggplot(data = fin, aes(x=as.factor(Spla), y=ACIf)) + geom_boxplot(aes(fill=Plot))
l1 <- l + labs(x='PCP la slope (1/m)',y='ACI')
ggarrange(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1, 
          ncol = 3, nrow = 4)

a <- ggplot(data = fin, aes(x=as.factor(Cgap), y=AEIf)) + geom_boxplot(aes(fill=Plot))
a1 <- a + labs(x='Canopy gap fraction',y='AEI',col='green')
b <- ggplot(data = fin, aes(x=as.factor(Cs), y=AEIf)) + geom_boxplot(aes(fill=Plot))
b1 <- b + labs(x='Canopy slope (m/m)',y='AEI')
c <- ggplot(data = fin, aes(x=as.factor(Sshs), y=AEIf)) + geom_boxplot(aes(fill=Plot))
c1 <- c + labs(x='TPPC hs slope (1/m)',y='AEI')
d <- ggplot(data = fin, aes(x=as.factor(Egap), y=AEIf)) + geom_boxplot(aes(fill=Plot))
d1 <- d + labs(x='Edge gap fraction',y='AEI')
e <- ggplot(data = fin, aes(x=as.factor(Vs), y=AEIf)) + geom_boxplot(aes(fill=Plot))
e1 <- e + labs(x='Vegetation slope (m/m)',y='AEI')
f <- ggplot(data = fin, aes(x=as.factor(Ssla), y=AEIf)) + geom_boxplot(aes(fill=Plot))
f1 <- f + labs(x='TPPC la slope (1/m)',y='AEI')
g <- ggplot(data = fin, aes(x=as.factor(Ch), y=AEIf)) + geom_boxplot(aes(fill=Plot))
g1 <- g + labs(x='Canopy height (m)',y='AEI')
h <- ggplot(data = fin, aes(x=as.factor(Cs075), y=AEIf)) + geom_boxplot(aes(fill=Plot))
h1 <- h + labs(x='Canopy slope 3/4 (m/m)',y='AEI')
i <- ggplot(data = fin, aes(x=as.factor(Sphs), y=AEIf)) + geom_boxplot(aes(fill=Plot))
i1 <- i + labs(x='PCP hs slope (1/m)',y='AEI')
j <- ggplot(data = fin, aes(x=as.factor(Vh), y=AEIf)) + geom_boxplot(aes(fill=Plot))
j1 <- j + labs(x='Vegetation height (m)',y='AEI')
k <- ggplot(data = fin, aes(x=as.factor(Vs075), y=AEIf)) + geom_boxplot(aes(fill=Plot))
k1 <- k + labs(x='Vegetation slope 3/4 (m/m)',y='AEI')
l <- ggplot(data = fin, aes(x=as.factor(Spla), y=AEIf)) + geom_boxplot(aes(fill=Plot))
l1 <- l + labs(x='PCP la slope (1/m)',y='AEI')
ggarrange(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1, 
          ncol = 3, nrow = 4)

a <- ggplot(data = fin, aes(x=as.factor(Cgap), y=BIf)) + geom_boxplot(aes(fill=Plot))
a1 <- a + labs(x='Canopy gap fraction',y='BI',col='green')
b <- ggplot(data = fin, aes(x=as.factor(Cs), y=BIf)) + geom_boxplot(aes(fill=Plot))
b1 <- b + labs(x='Canopy slope (m/m)',y='BI')
c <- ggplot(data = fin, aes(x=as.factor(Sshs), y=BIf)) + geom_boxplot(aes(fill=Plot))
c1 <- c + labs(x='TPPC hs slope (1/m)',y='BI')
d <- ggplot(data = fin, aes(x=as.factor(Egap), y=BIf)) + geom_boxplot(aes(fill=Plot))
d1 <- d + labs(x='Edge gap fraction',y='BI')
e <- ggplot(data = fin, aes(x=as.factor(Vs), y=BIf)) + geom_boxplot(aes(fill=Plot))
e1 <- e + labs(x='Vegetation slope (m/m)',y='BI')
f <- ggplot(data = fin, aes(x=as.factor(Ssla), y=BIf)) + geom_boxplot(aes(fill=Plot))
f1 <- f + labs(x='TPPC la slope (1/m)',y='BI')
g <- ggplot(data = fin, aes(x=as.factor(Ch), y=BIf)) + geom_boxplot(aes(fill=Plot))
g1 <- g + labs(x='Canopy height (m)',y='BI')
h <- ggplot(data = fin, aes(x=as.factor(Cs075), y=BIf)) + geom_boxplot(aes(fill=Plot))
h1 <- h + labs(x='Canopy slope 3/4 (m/m)',y='BI')
i <- ggplot(data = fin, aes(x=as.factor(Sphs), y=BIf)) + geom_boxplot(aes(fill=Plot))
i1 <- i + labs(x='PCP hs slope (1/m)',y='BI')
j <- ggplot(data = fin, aes(x=as.factor(Vh), y=BIf)) + geom_boxplot(aes(fill=Plot))
j1 <- j + labs(x='Vegetation height (m)',y='BI')
k <- ggplot(data = fin, aes(x=as.factor(Vs075), y=BIf)) + geom_boxplot(aes(fill=Plot))
k1 <- k + labs(x='Vegetation slope 3/4 (m/m)',y='BI')
l <- ggplot(data = fin, aes(x=as.factor(Spla), y=BIf)) + geom_boxplot(aes(fill=Plot))
l1 <- l + labs(x='PCP la slope (1/m)',y='BI')
ggarrange(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1, 
          ncol = 3, nrow = 4)

plot(Dist075,AEIm)
abline(lm(AEIm ~ Dist075))
plot(Vheight,AEIm)
abline(lm(AEIm ~ Vheight))
plot(Cheight,AEIm)
abline(lm(AEIm ~ Cheight))
plot(Cgap,AEIm)
abline(lm(AEIm ~ Cgap))

plot(Dist075,BIm)
abline(lm(BIm ~ Dist075))
plot(Vheight,BIm)
abline(lm(BIm ~ Vheight))
plot(Cheight,BIm)
abline(lm(BIm ~ Cheight))
plot(Cgap,BIm)
abline(lm(BIm ~ Cgap))

cor(Vheight,Cgap)
reg <- lm(ACIm ~ Slopesla + Intersla + Slopesla:Intersla)
summary(reg)
plot(Slopesla,ACIm)
plot(Intersla,ACIm)
reg2 <- lm(ACIm ~ Intershs)
plot(reg)

AIvals <- cbind(ACIm,AEIm,BIm)
rda_tree <- rda(ACIm ~ Tmax + Egap + Spla + Cs)
RsquareAdj(rda_tree)
plot(rda_tree, type='n', scaling=1)
orditorp(rda_tree, display='sp', cex=0.5, scaling=1, col='blue')
text(rda_tree, col='red')

autoplot(rda_tree,arrows=TRUE)

plot(Tmax,Hm)
ggplot(data = AIse) +
	geom_boxplot(mapping = aes(x = reorder(Location, ACIe, FUN = median), 
                   y = ACIe))

ggplot(data = AIs) +
  geom_point(mapping = aes(x = Tmax, y = Hm))

cormatrix <- cor(as.matrix(AIse[,6:11]))
round(cormatrix,digits=2)

cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(AIs[,6:11])
corrplot(cormatrix,method='color',tl.col="black",addCoef.col = "black",
         p.mat = p.mat, sig.level = 0.05)

cor(mydata)
cor(BIm,Hm)
drop1(fitADI,test='F')

lmACI <- lm(ACIm ~ Tmax+Cgap+Egap+Ch+Vh+Cs+Vs+Cs075+Vs075+Sshs+Ssla+Sphs+Spla)

lmaci <- lmer(ACIm ~ p + w + Tmax + w:Tmax + p:Tmax + (1|loc:p)+(1|loc:Tmax))
anova(lmACI)
summary(lmACI)

TukeyHSD(aov(fitACIm))

shapiro.test(Hm)
cor(ACI,H)

plot(Tmax,fitACI$residuals)
plot(Tmax,fitBI$residuals)
plot(Tmax,fitH$residuals)
bptest(fitACI)

fittest <- lm(Hm ~ Tmax)
bptest(fittest)

model <- tree(ACIm ~ Tmax+Cgap+Egap+Ch+Vh+Cs+Vs+Cs075+Vs075+
              Sshs+Ssla+Sphs+Spla)
plot(model)
text(model)

model2 <- prune.tree(model)
plot(model2)
model3 <- prune.tree(model,best=5)
plot(model3)
text(model3)

## MRF regression ##
library(MultivariateRandomForest)
trainX <- as.matrix(read.table("input_p5.txt",header = TRUE))
trainY <- as.matrix(read.table("output_p5.txt",header = TRUE))
testX <- as.matrix(read.table("test_p5.txt",header = TRUE))
n_tree <- 10
mtree <- 10
min_leaf <- 1

model <- build_forest_predict(trainX, trainY, n_tree, mtree, min_leaf, testX) 

