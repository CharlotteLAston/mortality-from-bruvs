###################################################
# Age-based SPR Calculations 
###################################################
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)

library(tidyr)
library(dplyr)
library(stringr)
library(sf)
library(ggplot2)
library(forcats)
library(ggtext)
library(tidyverse)
library(sf)
library(terra)
library(ggplot2)
library(patchwork)
library(ggnewscale)
library(scales)
library(ggspatial)
library(tidyterra)
library(L3Assess)

rm(list = ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial")

#### Pink Snapper ####
#* Set parameters ####
InitRecruit <- 1 # Initial recruitment
MaxModelAge <- 45 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)

# Growth, weight, maturity 
Linf <- c(1150, 1127) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.12, 0.12) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(-0.41, -0.46) # years - von Bertalanffy growth model parameters - Females, males

EstLenAtAge <- data.frame(EstFemLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters - from https://www.fish.wa.gov.au/Documents/research_reports/frr278.pdf
ln_lenwt_a <- -9.433484 # for log-log relationship
lenwt_b <- 2.7016 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
WLrel_Type <- 2 # 1=power, 2=log-log relationship
EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame

ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 0.5 # Ratio of females to males at age zero
FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
mat_A50 <- c(5.5, 3.8) # females, males - Logistic length (mm) at maturity relationship parameters
mat_A95 <- c(11.1, 9.3) # females, males - Logistic length (mm) at maturity relationship parameters
EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame

# Fishing selectivity
# Gear_sel_A50 <- c(4,4) # females, males - Logistic age selectivity relationship parameters - taken from logisic age based CC
# Gear_sel_A95 <- c(5.5,5.5) # females, males - Logistic age selectivity relationship parameters
# EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
Gear_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
Gear_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
FemGearSelAtAge=rep(1,MaxModelAge+1)
# FemGearSelAtAge[1:6] = 0
MalGearSelAtAge = FemGearSelAtAge
EstGearSelAtAge <- data.frame(FemGearSelAtAge=FemGearSelAtAge, MalGearSelAtAge=MalGearSelAtAge) # gear selectivity at age (from age 0), inputted as values in data frame

# Land_sel_A50 <- c(4, 4) # females, males - Logistic age selectivity relationship parameters
# Land_sel_A95 <- c(6, 6) # females, males - Logistic age selectivity relationship parameters
# EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
Land_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
Land_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
FemSelLandAtAge=rep(1,MaxModelAge+1)
FemSelLandAtAge[1:5]=0
MalSelLandAtAge = FemSelLandAtAge
EstLandSelAtAge <- data.frame(FemSelLandAtAge=FemSelLandAtAge, MalSelLandAtAge=MalSelLandAtAge) # gear selectivity at age (from age 0), inputted as values in data frame

ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
# FemRetProbAtAge=rep(1,MaxModelAge+1)
# FemRetProbAtAge[1:7]=0
# MalRetProbAtAge =rep(1,MaxModelAge+1)
# MalRetProbAtAge[1:6]=0
# EstRetenAtAge <- data.frame(FemRetProbAtAge=FemRetProbAtAge, MalRetProbAtAge=MalRetProbAtAge) # gear selectivity at age (from age 0), inputted as values in data frame
# DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)

# Recruitment relationship
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
NatMort = 0.12 # natural mortality  (year-1)
FMort <- 0.166 # estimate of fishing mortality, e.g. from catch curve analysis

#* AB-SPR calculations ####
Res=CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                             ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                             EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
                             EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
                             ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
Res$Fem_SPR
Res$Mal_SPR
Res$Eq_CombSexRelSpBiom

#* Plot outputs Pink Snapper ####
Current_F <- FMort # estimate of fishing mortality, e.g. from catch curve analysis
PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel_land at age, 4=mal mat/sel_land at age,
# 5=fem sel_land/sel_disc/prob reten/gear sel, 6=fem sel_land/sel_disc/prob reten/gear sel, 7=fem F at age,
# 8=mal F at age, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=fem biom at age,
# 13=catch at age, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
# 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
# 22=plot male depletion vs eq. recruitment (plots 1-22 plotted if !is.na(EggFertParam))
PlotPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                         lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                         ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam,
                         mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge,
                         Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95,
                         EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, Current_F)

##### Breaksea Cod ####
#* Set parameters ####
InitRecruit <- 1 # Initial recruitment
MaxModelAge <- 20 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)

# Growth, weight, maturity 
Linf <- c(561, 558) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.14, 0.17) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(-0.24, -0.04) # years - von Bertalanffy growth model parameters - Females, males

EstLenAtAge <- data.frame(EstFemLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters - from https://www.fish.wa.gov.au/Documents/research_reports/frr278.pdf
ln_lenwt_a <- log(0.000011) # for log-log relationship
lenwt_b <- 3.0697 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
WLrel_Type <- 2 # 1=power, 2=log-log relationship
EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame

ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 0.5 # Ratio of females to males at age zero
FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
mat_A50 <- c(2.17, 5.7) # females, males - Logistic length (mm) at maturity relationship parameters
mat_A95 <- c(10.21, 19.95) # females, males - Logistic length (mm) at maturity relationship parameters
EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame

# Fishing selectivity
Gear_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters - taken from logistic age based CC
Gear_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
# Gear_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
# Gear_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
# FemGearSelAtAge=rep(1,MaxModelAge+1)
# MalGearSelAtAge=rep(1,MaxModelAge+1)
# FemGearSelAtAge[1:6] = 0
# MalGearSelAtAge[1:5] = 0
# EstGearSelAtAge <- data.frame(FemGearSelAtAge=FemGearSelAtAge, MalGearSelAtAge=MalGearSelAtAge) # gear selectivity at age (from age 0), inputted as values in data frame

Land_sel_A50 <- c(4.5, 3.5) # females, males - Logistic age selectivity relationship parameters
Land_sel_A95 <- c(4.5, 6) # females, males - Logistic age selectivity relationship parameters
EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
# Land_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
# Land_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
# FemSelLandAtAge=rep(1,MaxModelAge+1)
# MalSelLandAtAge=rep(1,MaxModelAge+1)
# FemSelLandAtAge[1:6]= 0
# MalSelLandAtAge[1:5] = 0
# EstLandSelAtAge <- data.frame(FemSelLandAtAge=FemSelLandAtAge, MalSelLandAtAge=MalSelLandAtAge) # gear selectivity at age (from age 0), inputted as values in data frame

ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
# FemRetProbAtAge=rep(1,MaxModelAge+1)
# FemRetProbAtAge[1:5]=0
# MalRetProbAtAge =rep(1,MaxModelAge+1)
# MalRetProbAtAge[1:5]=0
EstRetenAtAge <- data.frame(FemRetProbAtAge=NA, MalRetProbAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)

# Recruitment relationship
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
NatMort = 0.3 # natural mortality  (year-1)
FMort <- 0.25 # estimate of fishing mortality, e.g. from catch curve analysis

#* AB-SPR calculations ####
Res=CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                             ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                             EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
                             EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
                             ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
Res$YPR
Res$CombSex_SPR
Res$Eq_CombSexRelSpBiom



#* Plot outputs Breaksea Cod ####
Current_F <- FMort # estimate of fishing mortality, e.g. from catch curve analysis
PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel_land at age, 4=mal mat/sel_land at age,
# 5=fem sel_land/sel_disc/prob reten/gear sel, 6=fem sel_land/sel_disc/prob reten/gear sel, 7=fem F at age,
# 8=mal F at age, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=fem biom at age,
# 13=catch at age, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
# 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
# 22=plot male depletion vs eq. recruitment (plots 1-22 plotted if !is.na(EggFertParam))
PlotPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                         lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                         ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam,
                         mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge,
                         Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95,
                         EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, Current_F)

##### Western King Wrasse Metro ####
#* Set parameters ####
InitRecruit <- 1 # Initial recruitment
MaxModelAge <- 20 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
Linf <- c(361, 473) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.24, 0.11) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(-0.11, -4.27) # years - von Bertalanffy growth model parameters - Females, males
EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
ln_lenwt_a <- -11.82 # for log-log relationship
lenwt_b <- 3.07 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
WLrel_Type <- 2 # 1=power, 2=log-log relationship
EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 1 # Ratio of females to males at age zero
FinalSex_A50 <- 5.4 # Logistic sex change relationship parameters (inflection point)
FinalSex_A95 <- 7.3 # Logistic sex change relationship parameters (95% of max probability)
EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
mat_A50 <- c(2.8, 2.8) # females, males - Logistic length (mm) at maturity relationship parameters
mat_A95 <- c(4.7, 4.7) # females, males - Logistic length (mm) at maturity relationship parameters
EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
Gear_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
Gear_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
FemGearSelAtAge=rep(1,MaxModelAge+1)
MalGearSelAtAge = FemGearSelAtAge
EstGearSelAtAge <- data.frame(FemGearSelAtAge=FemGearSelAtAge, MalGearSelAtAge=MalGearSelAtAge) # gear selectivity at age (from age 0), inputted as values in data frame

# Land_sel_A50 <- c(2,3) # females, males - Logistic age selectivity relationship parameters
# Land_sel_A95 <- c(5,4) # females, males - Logistic age selectivity relationship parameters
FemSelLandAtAge=rep(1,MaxModelAge+1)
FemSelLandAtAge[1:5]=0
MalSelLandAtAge = rep(1,MaxModelAge+1)
MalSelLandAtAge[1:4] =0
EstLandSelAtAge <- data.frame(FemSelLandAtAge=FemSelLandAtAge, MalSelLandAtAge=MalSelLandAtAge) # gear selectivity at age (from age 0), inputted as values in data frame


ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
# FemRetProbAtAge=rep(1,MaxModelAge+1)
# FemRetProbAtAge[1:6]=0
# MalRetProbAtAge = FemRetProbAtAge
EstRetenAtAge <- data.frame(FemRetProbAtAge=NA, MalRetProbAtAge=NA) # fish retention at age (from age 0), inputted as values in data frame
DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
NatMort = 0.3 # natural mortality  (year-1)
FMort <- 0.32 # estimate of fishing mortality, e.g. from catch curve analysis

#* AB-SPR Calculations ####
Res=CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                             ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                             EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
                             EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
                             ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)

Res$Fem_SPR
Res$Mal_SPR
Res$Eq_CombSexRelSpBiom

#* Plot outputs Western King Wrasse - Metro ####
Current_F <- FMort # estimate of fishing mortality, e.g. from catch curve analysis
PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel_land at age, 4=mal mat/sel_land at age,
# 5=fem sel_land/sel_disc/prob reten/gear sel, 6=fem sel_land/sel_disc/prob reten/gear sel, 7=fem F at age,
# 8=mal F at age, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=fem biom at age,
# 13=catch at age, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
# 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
# 22=plot male depletion vs eq. recruitment (plots 1-22 plotted if !is.na(EggFertParam))
PlotPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                         lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                         ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam,
                         mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge,
                         Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95,
                         EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, Current_F)



##### Western King Wrasse Jurien Bay ####
#* Set parameters ####
InitRecruit <- 1 # Initial recruitment
MaxModelAge <- 20 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
Linf <- c(333, 399) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.29, 0.14) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(0.05, -3.76) # years - von Bertalanffy growth model parameters - Females, males
EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
ln_lenwt_a <- -11.82 # for log-log relationship
lenwt_b <- 3.07 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
WLrel_Type <- 2 # 1=power, 2=log-log relationship
EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 1 # Ratio of females to males at age zero
FinalSex_A50 <- 4.7 # Logistic sex change relationship parameters (inflection point)
FinalSex_A95 <- 6.9 # Logistic sex change relationship parameters (95% of max probability)
EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
mat_A50 <- c(3.8, 3.8) # females, males - Logistic length (mm) at maturity relationship parameters
mat_A95 <- c(6.7, 6.7) # females, males - Logistic length (mm) at maturity relationship parameters
EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
Gear_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
Gear_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
FemGearSelAtAge=rep(1,MaxModelAge+1)
MalGearSelAtAge = FemGearSelAtAge
EstGearSelAtAge <- data.frame(FemGearSelAtAge=FemGearSelAtAge, MalGearSelAtAge=MalGearSelAtAge) # gear selectivity at age (from age 0), inputted as values in data frame

# Land_sel_A50 <- c(2,4) # females, males - Logistic age selectivity relationship parameters
# Land_sel_A95 <- c(3,5) # females, males - Logistic age selectivity relationship parameters
FemSelLandAtAge=rep(1,MaxModelAge+1)
FemSelLandAtAge[1:3]=0
MalSelLandAtAge = rep(1,MaxModelAge+1)
MalSelLandAtAge[1:5]=0
EstLandSelAtAge <- data.frame(FemSelLandAtAge=FemSelLandAtAge, MalSelLandAtAge=MalSelLandAtAge) # gear selectivity at age (from age 0), inputted as values in data frame


ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
# FemRetProbAtAge=rep(1,MaxModelAge+1)
# FemRetProbAtAge[1:6]=0
# MalRetProbAtAge = FemRetProbAtAge
EstRetenAtAge <- data.frame(FemRetProbAtAge=NA, MalRetProbAtAge=NA) # fish retention at age (from age 0), inputted as values in data frame
DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
NatMort = 0.3 # natural mortality  (year-1)
FMort <- 0.352 # estimate of fishing mortality, e.g. from catch curve analysis

#* AB-SPR Calculations ####
Res=CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                             ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                             EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
                             EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
                             ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)

Res$Fem_SPR
Res$Mal_SPR
Res$Eq_CombSexRelSpBiom

#* Plot Outputs for Western King Wrasse - Jurien Bay ####
Current_F <- FMort # estimate of fishing mortality, e.g. from catch curve analysis
PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel_land at age, 4=mal mat/sel_land at age,
# 5=fem sel_land/sel_disc/prob reten/gear sel, 6=fem sel_land/sel_disc/prob reten/gear sel, 7=fem F at age,
# 8=mal F at age, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=fem biom at age,
# 13=catch at age, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
# 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
# 22=plot male depletion vs eq. recruitment (plots 1-22 plotted if !is.na(EggFertParam))
PlotPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                         lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                         ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam,
                         mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge,
                         Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95,
                         EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, Current_F)


##### Southern Maori Wrasse Metro ####
#* Set parameters ####
InitRecruit <- 1 # Initial recruitment
MaxModelAge <- 15 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
Linf <- c(299, 303) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.33, 0.53) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(0.22, 1.2) # years - von Bertalanffy growth model parameters - Females, males
EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
ln_lenwt_a <- -12.9 # for log-log relationship
lenwt_b <- 2.265 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
WLrel_Type <- 2 # 1=power, 2=log-log relationship
EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 1 # Ratio of females to males at age zero
FinalSex_A50 <- 6.5 # Logistic sex change relationship parameters (inflection point)
FinalSex_A95 <- 9 # Logistic sex change relationship parameters (95% of max probability)
EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
mat_A50 <- c(2.6, 2.6) # females, males - Logistic length (mm) at maturity relationship parameters
mat_A95 <- c(5.1, 5.1) # females, males - Logistic length (mm) at maturity relationship parameters
EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
Gear_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
Gear_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
FemGearSelAtAge=rep(1,(MaxModelAge+1))
MalGearSelAtAge = FemGearSelAtAge
EstGearSelAtAge <- data.frame(FemGearSelAtAge=FemGearSelAtAge, MalGearSelAtAge=MalGearSelAtAge) # gear selectivity at age (from age 0), inputted as values in data frame

# Land_sel_A50 <- c(3.5,6) # females, males - Logistic age selectivity relationship parameters
# Land_sel_A95 <- c(5,8) # females, males - Logistic age selectivity relationship parameters
FemSelLandAtAge=rep(1,MaxModelAge+1)
FemSelLandAtAge[1:5]=0
MalSelLandAtAge = rep(1,MaxModelAge+1)
MalSelLandAtAge[1:8] = 0
EstLandSelAtAge <- data.frame(FemSelLandAtAge=FemSelLandAtAge, MalSelLandAtAge=MalSelLandAtAge) # gear selectivity at age (from age 0), inputted as values in data frame


ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
# FemRetProbAtAge=rep(1,MaxModelAge+1)
# FemRetProbAtAge[1:6]=0
# MalRetProbAtAge = FemRetProbAtAge
EstRetenAtAge <- data.frame(FemRetProbAtAge=NA, MalRetProbAtAge=NA) # fish retention at age (from age 0), inputted as values in data frame
DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
NatMort = 0.3 # natural mortality  (year-1)
FMort <- 0.367 # estimate of fishing mortality, e.g. from catch curve analysis

#* AB-SPR Calculations ####
Res=CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                             ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                             EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
                             EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
                             ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)

Res$Fem_SPR
Res$Mal_SPR
Res$Eq_CombSexRelSpBiom

#* Plot outputs Southern Maori Wrasse - Metro ####
Current_F <- FMort # estimate of fishing mortality, e.g. from catch curve analysis
PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel_land at age, 4=mal mat/sel_land at age,
# 5=fem sel_land/sel_disc/prob reten/gear sel, 6=fem sel_land/sel_disc/prob reten/gear sel, 7=fem F at age,
# 8=mal F at age, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=fem biom at age,
# 13=catch at age, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
# 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
# 22=plot male depletion vs eq. recruitment (plots 1-22 plotted if !is.na(EggFertParam))
PlotPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                         lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                         ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam,
                         mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge,
                         Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95,
                         EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, Current_F)

#### Brownspot Wrasse Metro ####
#* Set parameters ####
InitRecruit <- 1 # Initial recruitment
MaxModelAge <- 15 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
Linf <- c(357, 552) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.2, 0.05) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(-1.05, -8.98, 1.2) # years - von Bertalanffy growth model parameters - Females, males
EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
ln_lenwt_a <- -10.95 # for log-log relationship
lenwt_b <- 2.265 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
WLrel_Type <- 2 # 1=power, 2=log-log relationship
EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 1 # Ratio of females to males at age zero
FinalSex_A50 <- 6.3 # Logistic sex change relationship parameters (inflection point)
FinalSex_A95 <- 8.1 # Logistic sex change relationship parameters (95% of max probability)
EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
mat_A50 <- c(2.9, 2.9) # females, males - Logistic length (mm) at maturity relationship parameters
mat_A95 <- c(3.9, 3.9) # females, males - Logistic length (mm) at maturity relationship parameters
EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
Gear_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
Gear_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
FemGearSelAtAge=rep(1,(MaxModelAge+1))
MalGearSelAtAge = FemGearSelAtAge
EstGearSelAtAge <- data.frame(FemGearSelAtAge=FemGearSelAtAge, MalGearSelAtAge=MalGearSelAtAge) # gear selectivity at age (from age 0), inputted as values in data frame

# Land_sel_A50 <- c(2,6) # females, males - Logistic age selectivity relationship parameters
# Land_sel_A95 <- c(5,6.5) # females, males - Logistic age selectivity relationship parameters
FemSelLandAtAge=rep(1,MaxModelAge+1)
FemSelLandAtAge[1:5]=0
MalSelLandAtAge = rep(1,MaxModelAge+1)
MalSelLandAtAge[1:7]=0
EstLandSelAtAge <- data.frame(FemSelLandAtAge=FemSelLandAtAge, MalSelLandAtAge=MalSelLandAtAge) # gear selectivity at age (from age 0), inputted as values in data frame


ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
# FemRetProbAtAge=rep(1,MaxModelAge+1)
# FemRetProbAtAge[1:6]=0
# MalRetProbAtAge = FemRetProbAtAge
EstRetenAtAge <- data.frame(FemRetProbAtAge=NA, MalRetProbAtAge=NA) # fish retention at age (from age 0), inputted as values in data frame
DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
NatMort = 0.3 # natural mortality  (year-1)
FMort <- 0.423 # estimate of fishing mortality, e.g. from catch curve analysis

#* AB-SPR Calculations ####
Res=CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                             ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                             EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
                             EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
                             ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)

Res$Fem_SPR
Res$Mal_SPR
Res$Eq_CombSexRelSpBiom

#* Plot outputs Southern Maori Wrasse - Metro ####
Current_F <- FMort # estimate of fishing mortality, e.g. from catch curve analysis
PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel_land at age, 4=mal mat/sel_land at age,
# 5=fem sel_land/sel_disc/prob reten/gear sel, 6=fem sel_land/sel_disc/prob reten/gear sel, 7=fem F at age,
# 8=mal F at age, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=fem biom at age,
# 13=catch at age, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
# 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
# 22=plot male depletion vs eq. recruitment (plots 1-22 plotted if !is.na(EggFertParam))
PlotPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                         lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                         ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam,
                         mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge,
                         Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95,
                         EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, Current_F)
