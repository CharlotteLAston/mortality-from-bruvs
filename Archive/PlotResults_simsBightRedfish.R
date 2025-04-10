# Alex Hesp 9 Feb 2025
# Simulations for Coulson et al Bight redfish paper
# comparing traditional vb growth model vs length and age-based catch curve (LACC) from L3Assess


rm(list=ls())
library("vioplot") # base R

# if using ggplot... currently not set up
library(ggplot2)
library(ggthemes)
library(ggh4x)

# get stats 

# *********************************
# traditional von Bertalanffy curve
# *********************************
dat = read.csv("vbResults_combined.csv",header = TRUE)

dat$Scenario = NA
dat$Scenario[which(dat$Sel==1 & dat$FMort==1 & dat$RecVar==1)] = 1
dat$Scenario[which(dat$Sel==2 & dat$FMort==1 & dat$RecVar==1)] = 2
dat$Scenario[which(dat$Sel==1 & dat$FMort==2 & dat$RecVar==1)] = 3
dat$Scenario[which(dat$Sel==2 & dat$FMort==2 & dat$RecVar==1)] = 4
dat$Scenario[which(dat$Sel==1 & dat$FMort==3 & dat$RecVar==1)] = 5
dat$Scenario[which(dat$Sel==2 & dat$FMort==3 & dat$RecVar==1)] = 6
dat$Scenario[which(dat$Sel==1 & dat$FMort==1 & dat$RecVar==2)] = 7
dat$Scenario[which(dat$Sel==2 & dat$FMort==1 & dat$RecVar==2)] = 8
dat$Scenario[which(dat$Sel==1 & dat$FMort==2 & dat$RecVar==2)] = 9
dat$Scenario[which(dat$Sel==2 & dat$FMort==2 & dat$RecVar==2)] = 10
dat$Scenario[which(dat$Sel==1 & dat$FMort==3 & dat$RecVar==2)] = 11
dat$Scenario[which(dat$Sel==2 & dat$FMort==3 & dat$RecVar==2)] = 12
table(dat$Scenario) # aiming for at least 200 for paper, per scenario
# 1   2   3   4   5   6   7   8   9  10  11  12 
# 175 175 175 175 175 175 175 175 175 175 175 175 

True_Linf = 600
True_vbK = 0.1
True_tzero = 0

# relative error
dat$LinfRelErr = (dat$EstLinf - True_Linf) / True_Linf
dat$vbKRelErr = (dat$EstvbK - True_vbK) / True_vbK
# can't do MRE for tzero, as true value is zero

# absolute relative error
dat$LinfAbsRelErr = abs((dat$EstLinf - True_Linf) / True_Linf)
dat$vbKAbsRelErr = abs((dat$EstvbK - True_vbK) / True_vbK)
# can't do MARE for tzero, as true value is zero

# with better plots, these are added as dots to violin plots x scenario
# median relative error x scenario
temp=aggregate(LinfRelErr~Scenario,data=dat,median)
TradvbLinfMRE = temp$LinfRelErr
temp=aggregate(vbKRelErr~Scenario,data=dat,median)
TradvbvbKMRE = temp$vbKRelErr

# median absolute relative error x scenario
temp=aggregate(LinfAbsRelErr~Scenario,data=dat,median)
TradvbLinfMARE = temp$LinfRelErr
temp=aggregate(vbKAbsRelErr~Scenario,data=dat,median)
TradvbvbKMARE = temp$vbKRelErr

# ********************************
# length and age based catch curve
# *********************************
dat2 = read.csv("LACCResults_combined.csv",header = TRUE)
dat2 <- dat2[which(dat2$Converge==0),] # remove model runs not converged
dat2 <- dat2[which(dat2$EstLinf<800),] # remove model runs with large outliers for param estimates,
# as these would be rejected as erroneous anyway, as above all data points

dat2$Scenario = NA
dat2$Scenario[which(dat2$Sel==1 & dat2$FMort==1 & dat2$RecVar==1)] = 1
dat2$Scenario[which(dat2$Sel==2 & dat2$FMort==1 & dat2$RecVar==1)] = 2
dat2$Scenario[which(dat2$Sel==1 & dat2$FMort==2 & dat2$RecVar==1)] = 3
dat2$Scenario[which(dat2$Sel==2 & dat2$FMort==2 & dat2$RecVar==1)] = 4
dat2$Scenario[which(dat2$Sel==1 & dat2$FMort==3 & dat2$RecVar==1)] = 5
dat2$Scenario[which(dat2$Sel==2 & dat2$FMort==3 & dat2$RecVar==1)] = 6
dat2$Scenario[which(dat2$Sel==1 & dat2$FMort==1 & dat2$RecVar==2)] = 7
dat2$Scenario[which(dat2$Sel==2 & dat2$FMort==1 & dat2$RecVar==2)] = 8
dat2$Scenario[which(dat2$Sel==1 & dat2$FMort==2 & dat2$RecVar==2)] = 9
dat2$Scenario[which(dat2$Sel==2 & dat2$FMort==2 & dat2$RecVar==2)] = 10
dat2$Scenario[which(dat2$Sel==1 & dat2$FMort==3 & dat2$RecVar==2)] = 11
dat2$Scenario[which(dat2$Sel==2 & dat2$FMort==3 & dat2$RecVar==2)] = 12
# table(dat2$Scenario) # aiming for at least 200 for paper, per scenario
# 1   2   3   4   5   6   7   8   9  10  11  12 
# 175 124 175 128 175 131 163 108 167 124 164 132 
# less than for traditional vb, due to rejection of models
# not converged, large outliers. 

# specify true values
True_Linf = 600
True_vbK = 0.1
NatMort <- 4.22 / 84
dat2$TrueF = NA
dat2$TrueF[which(dat2$FMort==1)] = 2/3 * NatMort 
dat2$TrueF[which(dat2$FMort==2)] = NatMort 
dat2$TrueF[which(dat2$FMort==3)] = 3/2 * NatMort 
dat2$TrueL50 = NA
dat2$TrueL50[which(dat2$Sel==1)] = 400
dat2$TrueL50[which(dat2$Sel==2)] = 550
dat2$TrueDelta = 50

# Median relative error
dat2$LinfRelErr = (dat2$EstLinf - True_Linf) / True_Linf
dat2$vbKRelErr = (dat2$EstvbK - True_vbK) / True_vbK
dat2$FMortRelErr = (dat2$EstFMort - dat2$TrueF) / dat2$TrueF
dat2$L50RelErr = (dat2$EstSelL50 - dat2$TrueL50) / dat2$TrueL50
dat2$DeltaRelErr = (dat2$EstDelta - dat2$TrueDelta) / dat2$TrueDelta

# median absolute relative error
dat2$LinfAbsRelErr = abs((dat2$EstLinf - True_Linf) / True_Linf)
dat2$vbKAbsRelErr = abs((dat2$EstvbK - True_vbK) / True_vbK)
dat2$FMortAbsRelErr = abs((dat2$EstFMort - dat2$TrueF) / dat2$TrueF)
dat2$L50AbsRelErr = abs((dat2$EstSelL50 - dat2$TrueL50) / dat2$TrueL50)
dat2$DeltaAbsRelErr = abs((dat2$EstDelta - dat2$TrueDelta) / dat2$TrueDelta)

# with better plots, these are added as dots to violin plots x scenario
# median relative error x scenario
temp=aggregate(LinfRelErr~Scenario,data=dat2,median)
LACCvbLinfMRE = temp$LinfRelErr
temp=aggregate(vbKRelErr~Scenario,data=dat2,median)
LACCvbKMRE = temp$vbKRelErr
temp=aggregate(FMortRelErr~Scenario,data=dat2,median)
LACCFMortMRE = temp$FMortRelErr
temp=aggregate(L50RelErr~Scenario,data=dat2,median)
LACCL50MRE = temp$L50RelErr
temp=aggregate(DeltaRelErr~Scenario,data=dat2,median)
LACCDeltaMRE = temp$DeltaRelErr

# median absolute relative error x scenario
temp=aggregate(LinfAbsRelErr~Scenario,data=dat2,median)
LACCvbLinfMARE = temp$LinfRelErr
temp=aggregate(vbKAbsRelErr~Scenario,data=dat2,median)
LACCvbvMARE = temp$vbKRelErr
temp=aggregate(FMortAbsRelErr~Scenario,data=dat2,median)
LACCFMortMARE = temp$FMortAbsRelErr
temp=aggregate(L50AbsRelErr~Scenario,data=dat2,median)
LACCL50MARE = temp$L50AbsRelErr
temp=aggregate(DeltaAbsRelErr~Scenario,data=dat2,median)
LACCDeltaMARE = temp$DeltaAbsRelErr

# *******************************
# plot estimated parameter values
# *******************************

par(mfrow=c(2,2),mar=c(5,4,2,2))

# Linf
vioplot(EstLinf ~ Scenario, data = dat, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Param est", xlab="Scenario", ylim=c(500,700), main=bquote("Trad_vb L" ~ infinity ~ "(mm)"), border=NA)
legend("bottomright", legend=c("low F", "mid F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=600)

vioplot(EstLinf ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Param est", xlab="Scenario", ylim=c(500,700), main=bquote("LACC L" ~ infinity ~ "(mm)"), border=NA)
legend("bottomright", legend=c("low F", "mid F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=600)


# vbK
vioplot(EstvbK ~ Scenario, data = dat, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                   "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Param est", xlab="Scenario", ylim=c(0,0.2), main=bquote("Trad_vb k" ~ y^-1), border=NA)
legend("bottomright", legend=c("low F", "mid F", "hi F"), lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=0.1)


# vbK
vioplot(EstvbK ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                               "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Param est", xlab="Scenario", ylim=c(0,0.2), main=bquote("LACC_vb k" ~ y^-1), border=NA)
legend("bottomright", legend=c("low F", "mid F", "hi F"), lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=0.1)


# tzero
vioplot(Esttzero ~ Scenario, data = dat, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Param est", xlab="Scenario", ylim=c(-10,2), main=bquote("Trad_vb_tzero"), border=NA)
legend("bottomright", legend=c("low F", "mid F", "hi F"), lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=0.1)


# FMort
vioplot(EstFMort ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                  "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Param est", xlab="Scenario", ylim=c(0,0.2), main=bquote("F" ~ y^-1), border=NA)
legend("topleft", legend=c("low F", "mid F", "hi F"), lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
lines(1:2, rep(2/3*NatMort,2))
lines(3:4, rep(NatMort,2))
lines(5:6, rep(3/2*NatMort,2))
lines(7:8, rep(2/3*NatMort,2))
lines(9:10, rep(NatMort,2))
lines(11:12, rep(3/2*NatMort,2))

# L50 selectivity
dat2$EstSelL50
vioplot(EstSelL50 ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                  "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Param est", xlab="Scenario", ylim=c(350,650), main=bquote("LACC_vb L50, mm"), border=NA)
legend("topleft", legend=c("low F", "mid F", "hi F"), lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=400)
abline(h=550)

# delta (L50 + delta = L95)
dat2$EstDelta
vioplot(EstDelta ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                   "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Param est", xlab="Scenario", ylim=c(0,100), main=bquote("LACC_vb L50, mm"), border=NA)
legend("topleft", legend=c("low F", "mid F", "hi F"), lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=50)


# ***************************************
# plot MRE values for parameter estimates
# ***************************************

par(mfrow=c(2,2),mar=c(5,4,2,2))
# Linf
vioplot(LinfRelErr ~ Scenario, data = dat, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                   "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), ylab="Median relative error",
        xlab="Scenario", ylim=c(-0.1,0.1), main=bquote("Trad_vb L" ~ infinity ~ "(mm)"), border=NA)
legend("bottomright", legend=c("low F", "mid F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=0)

# Linf
vioplot(LinfRelErr ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                    "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Median relative error", xlab="Scenario", ylim=c(-0.1,0.1), main=bquote("LACC L" ~ infinity ~ "(mm)"), border=NA)
legend("bottomright", legend=c("low F", "mid F", "hi F"), lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=0)


# vbK
vioplot(vbKRelErr ~ Scenario, data = dat, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                  "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), ylab="Median relative error",
        xlab="Scenario", ylim=c(-0.5,0.5), main=bquote("Trad_vb k" ~ y^-1), border=NA)
legend("bottomright", legend=c("low F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff", "#458cff"))
abline(h=0)

# vbK
vioplot(vbKRelErr ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                  "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), ylab="Median relative error",
        xlab="Scenario", ylim=c(-0.5,0.5), main=bquote("LACC  k" ~ y^-1), border=NA)
legend("bottomright", legend=c("low F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff", "#458cff"))
abline(h=0)

plot(0,0,cex=0, ylab="")

# FMort
vioplot(FMortRelErr ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                   "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), ylab="Median relative error",
        xlab="Scenario", ylim=c(-2,2), main=bquote("LACC F" ~ y^-1), border=NA)
legend("bottomright", legend=c("low F", "mid F", "hi F"), lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=0)


# L50
vioplot(L50RelErr ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                   "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), ylab="Median relative error",
        xlab="Scenario", ylim=c(-0.5,0.5), main=bquote("LACC L50"), border=NA)
legend("bottomright", legend=c("low F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff", "#458cff"))
abline(h=0)


# delta
vioplot(DeltaRelErr ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                   "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), ylab="Median relative error",
        xlab="Scenario", ylim=c(-1,1), main=bquote("LACC delta"), border=NA)
legend("bottomright", legend=c("low F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff", "#458cff"))
abline(h=0)


# ***************************************
# plot MARE values for parameter estimates
# ***************************************

par(mfrow=c(2,2),mar=c(5,4,2,2))
# Linf
vioplot(LinfAbsRelErr ~ Scenario, data = dat, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                   "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Abs median relative error", xlab="Scenario", ylim=c(0,0.2), main=bquote("Trad_vb L" ~ infinity ~ "(mm)"), border=NA)
legend("topright", legend=c("low F", "mid F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=0)

# Linf
vioplot(LinfAbsRelErr ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                    "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Abs median relative error", xlab="Scenario", ylim=c(0,0.2), main=bquote("LACC L" ~ infinity ~ "(mm)"), border=NA)
legend("topright", legend=c("low F", "mid F", "hi F"), lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=0)


# vbK
vioplot(vbKAbsRelErr ~ Scenario, data = dat, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                  "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Abs median relative error", xlab="Scenario", ylim=c(0,1), main=bquote("Trad_vb k" ~ y^-1), border=NA)
legend("bottomright", legend=c("low F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff", "#458cff"))
abline(h=0)

# vbK
vioplot(vbKAbsRelErr ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                   "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Abs median relative error",
        xlab="Scenario", ylim=c(0,1), main=bquote("LACC  k" ~ y^-1), border=NA)
legend("bottomright", legend=c("low F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff", "#458cff"))
abline(h=0)

plot(0,0,cex=0, ylab="")
# FMort
vioplot(FMortAbsRelErr ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                     "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"), 
        ylab="Abs median relative error", xlab="Scenario", ylim=c(0,5), main=bquote("LACC F" ~ y^-1), border=NA)
legend("bottomright", legend=c("low F", "mid F", "hi F"), lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff","#75b7ff","#458cff"))
abline(h=0)


# L50
vioplot(LinfAbsRelErr ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                   "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"),
        ylab="Abs median relative error", xlab="Scenario", ylim=c(0,0.2), main=bquote("LACC L50"), border=NA)
legend("topleft", legend=c("low F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff", "#458cff"))
abline(h=0)


# delta
vioplot(DeltaAbsRelErr ~ Scenario, data = dat2, col = c("#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff", 
                                                       "#bef7ff", "#bef7ff","#75b7ff","#75b7ff","#458cff", "#458cff"),
        ylab="Abs median relative error", xlab="Scenario", ylim=c(0,1), main=bquote("LACC delta"), border=NA)
legend("topleft", legend=c("low F", "hi F"),
       lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("#bef7ff", "#458cff"))
abline(h=0)

# ******************************************************************************
# ggplot code (Brett provided, and I've modified a bit, but then went to base R
# to progress more quickly)
# ******************************************************************************

# MREplotDynm <- ggplot(dat, aes(x=factor(Scenario), y=LinfMREVal,fill=factor(Scenario))) +
#   geom_violin(trim=TRUE,position=position_dodge(1), scale = "width") +
#   stat_summary(fun=median, geom="point", size=1, color="black",position=position_dodge(0.5)) + # calculates median
#   theme_few(base_size = 20, base_family = "ArialMT") + # colour themes
#   geom_hline(yintercept = 0) + # horizontal line
#   theme(legend.position= "none") + # removes legend
#   xlab("Scenario") + ylab("Median relative error") + 
#   theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5), strip.text = element_text(face = "italic")) +
#   scale_fill_manual(values=c("green", "green","pink","pink","green", "green","pink","pink"))
# # geom_text(data = std_data, aes(x=factor(scenario), y=y_pos, label=CV), # labels on top
# #           position = position_dodge(1), size = 1, inherit.aes = FALSE)
# MREplotDynm
# 
# 
# MREplotDynm <- ggplot(dat, aes(x=factor(Scenario), y=vbKMREVal,fill=factor(Scenario))) +
#   geom_violin(trim=TRUE,position=position_dodge(1), scale = "width") +
#   stat_summary(fun=median, geom="point", size=1, color="black",position=position_dodge(0.5)) + # calculates median
#   theme_few(base_size = 20, base_family = "ArialMT") + # colour themes
#   geom_hline(yintercept = 0) + # horizontal line
#   theme(legend.position= "none") + # removes legend
#   xlab("Scenario") + ylab("Median relative error") + 
#   theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5), strip.text = element_text(face = "italic")) +
#   scale_fill_manual(values=c("green", "green","pink","pink","green", "green","pink","pink"))
# # geom_text(data = std_data, aes(x=factor(scenario), y=y_pos, label=CV), # labels on top
# #           position = position_dodge(1), size = 1, inherit.aes = FALSE)
# MREplotDynm
# 
# # traditional vb curve (growth coefficient)
# MREplotDynm2 <- ggplot(dat, aes(x=factor(Scenario), y=vbKMREVal,fill=factor(Scenario))) +
#   geom_violin(trim=TRUE,position=position_dodge(1), scale = "width") +
#   stat_summary(fun=median, geom="point", size=1, color="black",position=position_dodge(0.5)) + # calculates median
#   theme_few(base_size = 20, base_family = "ArialMT") + # colour themes
#   geom_hline(yintercept = 0) + # horizontal line
#   theme(legend.position= "none") + # removes legend
#   xlab("Scenario") + ylab("Median relative error") + 
#   theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5), strip.text = element_text(face = "italic")) +
#   scale_fill_manual(values=c("#FF0000", "#0000FF","#FF0000","#0000FF","#FF0000", "#0000FF","#FF0000","#0000FF"))
# # geom_text(data = std_data, aes(x=factor(scenario), y=y_pos, label=CV), # labels on top
# #           position = position_dodge(1), size = 1, inherit.aes = FALSE)
# 
# MREplotDynm2
# 
# 
# 
