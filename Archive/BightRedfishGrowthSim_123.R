
# Alex Hesp 9 Feb 2025
# Simulation for Bight redfish paper

# Simulate realistic length and age data for Bight redfish,
# with relatively 'early' vs 'late' selectivity schedule,
# i.e. mimic situation on ESC vs Capes, fit traditional vb
# curve vs length and age based catch curve (LACC), and
# assess bias/precision.

# Recruitment variability - explore none vs SigmaR=0.6 (with autocorr = 0.7)
# Fishing mortality - explore low vs high
# Selectivity - explore early vs late

# Specifying same growth params as same for females and males (to simplify, as
# not much different anyway), combining data and fitting traditional vb and
# LACC as single sex models

# ********************************************************************************
# if installing from .gz file (for Alex)
# install.packages("C:/~/L3Assess_0.1.0.tar.gz", source = TRUE, repos=NULL)

# if installing from github
# library(devtools)
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)
# ********************************************************************************


rm(list=ls())
library(L3Assess)
library(WAFishBiology)


# Simulate length data from dynamic age-length model, given specified biology and exploitation history
MaxModelAge <- 100 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
MaxLen <- 800
LenInc <- 20
lbnd <- seq(0,MaxLen - LenInc, LenInc)
ubnd <- lbnd + LenInc
midpt <- lbnd + (LenInc/2)
nLenCl <- length(midpt)
GrowthCurveType <- 1 # 1 = von Bertalanffy, 2 = Schnute
Linf <- c(600, 600) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
GrowthParams <- data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
RefnceAges <- NA
CVSizeAtAge = c(0.05,0.05)
lenwt_a <- 0.00001299 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
ln_lenwt_a <- NA # for log-log relationship
lenwt_b <- 3.000 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
WLrel_Type <- 1 # 1=power, 2=log-log relationship
EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
                         EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 0.5 # Ratio of females to males at recruitment age/length
FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
mat_L50 <- c(350, 350) # females, males - Logistic length (mm) at maturity relationship parameters
mat_L95 <- c(550, 550) # females, males - Logistic length (mm) at maturity relationship parameters
# Based on Coulson et al., but modified to consider area differences, and low sample sizes for 
# immature fish which would have led to quite biased results
EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
                          EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
                              EstMalGearSelAtLen=NA)
EstRetenAtLen <- data.frame(EstFemRetenAtLen=NA,
                            EstMalRetenAtLen=NA)
DiscMort <- 0 # discard mortality (e.g. 50% released fish die = 0.5)
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
NatMort <- 4.22 / 84 # natural mortality  (year-1)
# Specify exploitation history - initial in increase F, reduction in F then stable F period
nPeriods <- 1 # number of periods of mortality trends
InitYr <- 1926 # first year of mortality
PeriodEndYr <- 2025 # periods when mortality changes
InitYr_FMort <- NatMort # initial mortality, at beginning of first period
PeriodEndYr_FMort <-NatMort
randerr_sd <- 0 # annual fishing mortality variation
InitRec <- 1000 # initial recruitment level

# ************************************************************************
# Alternative simulation model options

InitYr_FMort <- 2/3 * NatMort # initial mortality, at beginning of first period
PeriodEndYr_FMort <- 2/3 * NatMort
# InitYr_FMort <- 3/2 * NatMort # initial mortality, at beginning of first period
# PeriodEndYr_FMort <- 3/2 * NatMort
sel_L50 <- c(400, 400) # females, males - Logistic length selectivity relationship parameters
sel_L95 <- c(450, 450) # females, males - Logistic length selectivity relationship parameters
# sel_L50 <- c(550, 550) # females, males - Logistic length selectivity relationship parameters
# sel_L95 <- c(600, 600) # females, males - Logistic length selectivity relationship parameters
# autocorr <- 0.7 # annual fishing mortality autocorrelation
# lnSigmaR <- 0.6 # annual recruitment variation
autocorr <- 0 # annual fishing mortality autocorrelation
lnSigmaR <- 0 # annual recruitment variation

# ignoring discarding to simplify analysis, and make easier to get specified sample
# size for retained fish only
ret_Pmax <- c(1,1) # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_L50 <- c(10,10) # females, males - Logistic fish retention at length parameters
ret_L95 <- c(30,30) # females, males - Logistic fish retention at length parameters
# ************************************************************************

set.seed(123)

res=SimExploitationHistory_DynMod(nPeriods, InitYr, PeriodEndYr, InitYr_FMort,
                                  PeriodEndYr_FMort, autocorr, randerr_sd)


# Get random age and length data
nYears <- res$nYears
FMortByYear <- res$rand_FMort
SimAnnSampSize <- 1000 # note, this is combined sample size for both sexes (for retained fish)
# actual sample sizes vary from just over 1000 to over 4000 fish, depending on region
res=SimLenAndAgeFreqData_DynMod_ALB(SimAnnSampSize, nYears, lnSigmaR, autocorr, InitRec,
                                    MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                    RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                    ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                    EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                    ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMortByYear)


CheckDiagnostics = TRUE
if (CheckDiagnostics) {
  
  # *****************************
  # check biology and selectivity
  # *****************************
  
  par(mfrow=c(3,2))
  plot(1:MaxModelAge, res$ModelDiag$MeanSizeAtAge[1,], "l", ylim=c(0,800), ylab="Length, mm", col="red") # females
  lines(1:MaxModelAge, res$ModelDiag$MeanSizeAtAge[2,], col="blue")
  plot(midpt, res$ModelDiag$FemWtAtLen, "l", col="red", xlim=c(0,800), ylab="Weight, kg")
  lines(midpt, res$ModelDiag$MalWtAtLen, col="blue")
  plot(midpt, res$ModelDiag$FemPropMatAtLen, "l", col="red", xlim=c(0,800), ylab="Prop mat")
  lines(midpt, res$ModelDiag$MalPropMatAtLen, col="blue")
  plot(midpt, res$ModelDiag$FemGearSelAtLen, "l", col="red", xlim=c(0,800), ylab="Gear Sel")
  lines(midpt, res$ModelDiag$MalGearSelAtLen, col="blue")
  plot(midpt, res$ModelDiag$FemRetProbAtLen, "l", col="red", xlim=c(0,800), ylab="Reten Prob")
  lines(midpt, res$ModelDiag$MalRetProbAtLen, col="blue")
  plot(midpt, res$ModelDiag$FemSelLandAtLen, "l", col="red", xlim=c(0,800), ylab="Sel Land & Disc")
  lines(midpt, res$ModelDiag$MalSelLandAtLen, col="blue")
  lines(midpt, res$ModelDiag$FemSelDiscAtLen, col="red")
  lines(midpt, res$ModelDiag$MalSelDiscAtLen, col="blue")
  
  # *************************
  # check population dynamics
  # *************************
  par(mfrow=c(2,2))
  plot(1:nYears, res$ModelDiag$FMortByYear, "l", ylim=c(0,1.2*max(res$ModelDiag$FMortByYear)),
       ylab="F") # fishing mortality
  abline(h=NatMort, col="blue")
  plot(1:nYears, res$ModelDiag$FemSpBiom, "l", ylim=c(0,1.2*max(res$ModelDiag$FemSpBiom)), ylab="Fem sp biom") # female spawning biomass
  plot(1:nYears, res$ModelDiag$MalSpBiom, "l", ylim=c(0,1.2*max(res$ModelDiag$MalSpBiom)), ylab="Mal sp biom") # male spawning biomass
  plot(1:nYears, res$ModelDiag$RelSpBiom, "l", ylim=c(0,1.0), ylab="Rel sp biom") # spawning biomass
  abline(h=0.4, col="green")
  abline(h=0.2, col="red")
  par(mfrow=c(2,2))
  plot(seq(1,nYears+1,1), res$ModelDiag$AnnRecruit, "l", ylim=c(0,1.2*max(res$ModelDiag$AnnRecruit)), ylab="Ann Rec") # annual recruitment
  plot(seq(1,nYears,1), res$ModelDiag$RetCatchBiom_Fem, "l", ylim=c(0,1.2*max(res$ModelDiag$RetCatchBiom_Fem)), col="red", ylab="Fem Ret catch biom") 
  lines(seq(1,nYears,1), res$ModelDiag$DiscCatchBiom_Fem, col="pink") 
  plot(seq(1,nYears,1), res$ModelDiag$RetCatchBiom_Mal, "l", ylim=c(0,1.2*max(res$ModelDiag$RetCatchBiom_Mal)), col="blue", ylab="Mal Ret catch biom") 
  lines(seq(1,nYears,1), res$ModelDiag$DiscCatchBiom_Mal, col="lightblue") 
  
  # Inspect randomly-generated data
  # age frequency data (final year) - females
  par(mfrow=c(2,2))
  plot(1:MaxModelAge,res$ResultsSummary$ObsRetCatchFreqAtDecAge_Fem[nYears,], "l", col="red")
  plot(1:MaxModelAge,res$ResultsSummary$ObsRetCatchFreqAtDecAge_Mal[nYears,], "l", col="blue")
  
  # check total sample size (both sexes, retained plus discarded should equal user input sample size)
  sum(sum(res$ResultsSummary$ObsRetCatchFreqAtDecAge_Fem[nYears,])+
        sum(res$ResultsSummary$ObsRetCatchFreqAtDecAge_Mal[nYears,]))
  
  # length frequency data (final year) - females
  # age frequency data (final year) - females
  plot(midpt,res$ResultsSummary$ObsRetCatchFreqAtLen_Fem[nYears,], "l", col="red")
  plot(midpt,res$ResultsSummary$ObsRetCatchFreqAtLen_Mal[nYears,], "l", col="blue")
  
  # check total sample size (both sexes, retained plus discarded should equal user input sample size)
  sum(sum(res$ResultsSummary$ObsRetCatchFreqAtLen_Fem[nYears,])+
        sum(res$ResultsSummary$ObsRetCatchFreqAtLen_Mal[nYears,]))
  
  # conditional length and age data (final year) - females
  par(mfrow=c(1,1))
  plot(res$ResultsSummary$ObsAgeClRetCatch_Fem[nYears,], res$ResultsSummary$ObsLenClRetCatchMidPt_Fem[nYears,],
       ylim=c(0,MaxLen), xlim=c(0,MaxModelAge), ylab="Length, mm", xlab="Ages, y", col="red")
  points(res$ResultsSummary$ObsAgeClRetCatch_Mal[nYears,], res$ResultsSummary$ObsLenClRetCatchMidPt_Mal[nYears,], col="blue")
  
}

# ************************
# Traditional growth curve
# ************************
DataType=1  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
ObsAge = as.vector(unlist(c(res$ResultsSummary$ObsAgeClRetCatch_Fem[nYears,], res$ResultsSummary$ObsAgeClRetCatch_Mal[nYears,])))
ObsLen = as.vector(unlist(c(res$ResultsSummary$ObsLenClRetCatchMidPt_Fem[nYears,], res$ResultsSummary$ObsLenClRetCatchMidPt_Mal[nYears,])))
ObsAge = ObsAge[which(!is.na(ObsAge))]
ObsLen = ObsLen[which(!is.na(ObsLen))]
plot(ObsAge,ObsLen)
params = c(log(600),log(0.2),0) # log(Linf), log(k), tzero
FittedRes=GetvonBertalanffyGrowthResults(params, nSexes=1, DataType=1, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLensd=NA)
FittedRes$ParamEst

# **********************************************************
# Constant recruitment -  1 realisation (sample size = 1000 fish)
# autocorr <- 0 # annual fishing mortality autocorrelation
# lnSigmaR <- 0 # annual recruitment variation

# early selectivity
# sel_L50 <- c(400, 400) 
# sel_L95 <- c(450, 450)
# F=2/3M 
# Estimate lw_95%CL up_95%CL
# Linf    611.60   598.00   625.50
# vbK       0.08     0.06     0.09
# tzero    -4.20    -6.58    -1.82

# F=3/2M 
# Estimate lw_95%CL up_95%CL
# Linf    600.20   583.10   617.80
# vbK       0.08     0.07     0.10
# tzero    -3.29    -5.51    -1.06

# late selectivity
# sel_L50 <- c(550, 550) 
# sel_L95 <- c(600, 600) 
# F=2/3M 
# Estimate lw_95%CL up_95%CL
# Linf    622.10   615.20   629.20
# vbK       0.11     0.08     0.15
# tzero    -4.70   -11.30     1.91

# F=3/2M 
# Estimate lw_95%CL up_95%CL


# **********************************************************


# **********************************************************
# Variable recruitment -  1 realisation (sample size = 1000 fish)
# autocorr <- 0.7 # annual fishing mortality autocorrelation
# lnSigmaR <- 0.6 # annual recruitment variation

# early selectivity
# sel_L50 <- c(400, 400) 
# sel_L95 <- c(450, 450)
# F=2/3M 
# Estimate lw_95%CL up_95%CL
# Linf    614.80   603.80   625.90
# vbK       0.08     0.07     0.10
# tzero    -2.93    -4.47    -1.38

# F=3/2M 
# Estimate lw_95%CL up_95%CL
# Linf    610.50   594.60   626.80
# vbK       0.08     0.07     0.09
# tzero    -3.75    -5.43    -2.08

# late selectivity
# sel_L50 <- c(550, 550) 
# sel_L95 <- c(600, 600) 
# F=2/3M 
# Estimate lw_95%CL up_95%CL
# Linf    625.60   618.10   633.20
# vbK       0.10     0.07     0.13
# tzero    -6.06   -12.33     0.22

# F=3/2M 
# Estimate lw_95%CL up_95%CL
# Linf    603.80   597.60   610.10
# vbK       0.17     0.10     0.26
# tzero     0.63    -5.71     6.97

# **********************************************************

plotages=seq(0, 100,0.1)
PlotFittedGrowthCurve(DataType,  nSexes=1, GrowthEqn=1, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLensd=NA,
                      params, Ref_ages, plotages, ymax=800, xmax=NA, yint=100, xint=NA, GraphTitle=NA,
                      xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)


# ****
# LACC
# ****

# length frequency
ObsRetCatchFreqAtLen = as.vector(unlist(res$ResultsSummary$ObsRetCatchFreqAtLen_Fem[nYears,] + res$ResultsSummary$ObsRetCatchFreqAtLen_Mal[nYears,]))
# age-length frequency
DecAges = 1:MaxModelAge
ObsRetCatchFreqAtLengthAndAge <- data.frame(matrix(nrow = MaxModelAge, ncol = length(lbnd)))
colnames(ObsRetCatchFreqAtLengthAndAge) <- lbnd
lbns = trunc(ObsLen/LenInc)*LenInc
ObsRetCatchFreqAtLengthAndAge = as.matrix(table(factor(ObsAge, levels=1:MaxModelAge),
                                                factor(lbns, levels=c(lbnd))))
MaxAge = MaxModelAge
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logist sel curve, 3=separate logist sel params, common growth
SelectivityAtLen = NA # selectivity vector
DiscMort = 0
MLL = NA

# fit model
InitFishMort = 0.15 # specify starting parameters
InitL50 = 450
InitDelta = 50 # L95-L50
InitLinf = 600
InitvbK = 0.15
InitCVSizeAtAge = 0.07
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))

FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                                lbnd, ubnd, midpt, SelectivityAtLen, DiscMort, MaxAge, NatMort, TimeStep)
FittedRes$ParamEst
# Estimate lw_95%CL up_95%CL
# FMort     0.030    0.025    0.037
# SelL50  397.830  390.199  405.610
# Delta    44.870   34.912   57.667
# Linf    604.222  594.151  614.463
# vbK       0.098    0.093    0.102
# CV        0.051    0.049    0.054



# Simulations with 1 realisation (i.e. 1 set of synthetic data)

# **********************************************************
# Constant recruitment (sample size = 1000 fish)
# autocorr <- 0 # annual fishing mortality autocorrelation
# lnSigmaR <- 0 # annual recruitment variation

# **********************************************************
# early selectivity
# sel_L50 <- c(400, 400) 
# sel_L95 <- c(450, 450)
# F=2/3M 
# Estimate lw_95%CL up_95%CL
# FMort     0.030    0.025    0.037
# SelL50  397.830  390.199  405.610
# Delta    44.870   34.912   57.667
# Linf    604.222  594.151  614.463
# vbK       0.098    0.093    0.102
# CV        0.051    0.049    0.054

# F=3/2M 
# Estimate lw_95%CL up_95%CL
# FMort     0.083    0.072    0.094
# SelL50  399.481  392.452  406.635
# Delta    47.186   38.508   57.820
# Linf    603.860  591.521  616.457
# vbK       0.097    0.092    0.101
# CV        0.050    0.047    0.053

# late selectivity
# F=2/3M
# sel_L50 <- c(550, 550) 
# sel_L95 <- c(600, 600) 
# Estimate lw_95%CL up_95%CL
# FMort     0.038    0.029    0.051
# SelL50  552.135  548.490  555.805
# Delta    52.143   49.107   55.367
# Linf    604.279  599.151  609.450
# vbK       0.099    0.094    0.104
# CV        0.052    0.048    0.056

# F=3/2M
# Estimate lw_95%CL up_95%CL
# FMort     0.085    0.067    0.107
# SelL50  551.736  548.381  555.112
# Delta    45.872   44.534   47.250
# Linf    597.608  593.851  601.389
# vbK       0.100    0.095    0.106
# CV        0.047    0.044    0.051

# **********************************************************
# Variable recruitment -  1 realisation (sample size = 1000 fish)
# autocorr <- 0.7 # annual fishing mortality autocorrelation
# lnSigmaR <- 0.6 # annual recruitment variation

# early selectivity
# sel_L50 <- c(400, 400) 
# sel_L95 <- c(450, 450)

# F=2/3M 
# Estimate lw_95%CL up_95%CL
# FMort     0.030    0.025    0.037
# SelL50  369.830  363.560  376.209
# Delta    37.053   29.236   46.960
# Linf    605.217  596.691  613.865
# vbK       0.104    0.101    0.107
# CV        0.047    0.045    0.050

# F=3/2M 
# Estimate lw_95%CL up_95%CL
# FMort     0.077    0.068    0.087
# SelL50  371.779  366.211  377.431
# Delta    37.457   30.275   46.343
# Linf    594.254  583.919  604.772
# vbK       0.107    0.103    0.111
# CV        0.048    0.046    0.051

# late selectivity
# F=2/3M
# sel_L50 <- c(550, 550) 
# sel_L95 <- c(600, 600) 
# Estimate lw_95%CL up_95%CL
# FMort     0.019    0.012    0.029
# SelL50  551.370  547.743  555.021
# Delta    55.164   53.446   56.937
# Linf    606.534  602.473  610.622
# vbK       0.097    0.092    0.103
# CV        0.049    0.045    0.053

# F=3/2M
# Estimate lw_95%CL up_95%CL
# FMort     0.084    0.063    0.112
# SelL50  548.996  545.439  552.575
# Delta    45.027   43.833   46.254
# Linf    594.023  590.116  597.956
# vbK       0.092    0.085    0.100
# CV        0.056    0.051    0.061
# **********************************************************
FittedRes$ParamEst


# plot catch curve results
# growth curve
PlotAgeLengthCatchCurve_Growth(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                               lbnd, ubnd, midpt, SelectivityAtLen, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
                               xaxis_lab=NA, yaxis_lab=NA, xmax=100, xint=20,
                               ymax=800, yint=200, PlotCLs=TRUE, FittedRes, nReps=200, ShowLegend=TRUE)
# marginal length
PlotAgeLengthCatchCurve_MargLength(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                   lbnd, ubnd, midpt, SelectivityAtLen, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
                                   ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)

# conditional age at length - pearson residuals
PlotAgeLengthCatchCurve_Pears_Resid(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                    lbnd, ubnd, midpt, SelectivityAtLen, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                    xaxis_lab=NA, yaxis_lab=NA, xmin=200, xmax=800, xint=NA, ymin=0, ymax=100, yint=NA, CircleScale=3, MinLenClFreq=10, ShowLegend=T, FittedRes)


# Multiple realisations (assess precision and bias)



# Save results
FileName = "vbResults_123.csv"
FileName2 = "LACCResults_123.csv"
SaveFiles <- T
set.seed(123)
randvals = sample(1:100000,100000)
nTrials = 100
nRecVarScenarios = 2 # none, moderate
nFMortScenarios = 3 # low, medium, high
nSelScenarios = 2 # early, late

vbResults <- data.frame(matrix(nrow = nTrials*nRecVarScenarios*nFMortScenarios*nSelScenarios, ncol = 8))
colnames(vbResults) <- c("Trial","RecVar","FMort","Sel","EstLinf","EstvbK","Esttzero","Converge")

LACCResults <- data.frame(matrix(nrow = nTrials*nRecVarScenarios*nFMortScenarios*nSelScenarios, ncol = 11))
colnames(LACCResults) <- c("Trial","RecVar","FMort","Sel","EstFMort","EstSelL50","EstDelta","EstLinf","EstvbK","EstCV","Converge")


k = 0
for (trials in 1:nTrials) {
  
  for (RecVarScenarios in 1:nRecVarScenarios) { # low vs high rec var
    if (RecVarScenarios == 1) {
      autocorr <- 0 
      lnSigmaR <- 0 
    } else {
      autocorr <- 0.7 
      lnSigmaR <- 0.6 
    }
    
    for (FMortScenarios in 1:nFMortScenarios) { # low vs high fishing mortality
      if (FMortScenarios == 1) {
        InitYr_FMort <- 2/3 * NatMort 
        PeriodEndYr_FMort <- 2/3 * NatMort
      } else if (FMortScenarios == 2) {
        InitYr_FMort <- NatMort 
        PeriodEndYr_FMort <- NatMort
      } else {
        InitYr_FMort <- 3/2 * NatMort 
        PeriodEndYr_FMort <- 3/2 * NatMort
      }
      
      # simulate data
      res=SimExploitationHistory_DynMod(nPeriods, InitYr, PeriodEndYr, InitYr_FMort,
                                        PeriodEndYr_FMort, autocorr, randerr_sd)
      FMortByYear = res$rand_FMort
      
      for (SelScenarios in 1:nSelScenarios) { # early vs late selectivity
        k = k + 1
        if (SelScenarios == 1) {
          sel_L50 <- c(400, 400) 
          sel_L95 <- c(450, 450) 
        } else {
          sel_L50 <- c(550, 550) 
          sel_L95 <- c(600, 600) 
        }
        
        cat("trials",trials,"k",k,"RecVarScenarios",RecVarScenarios,"FMortScenarios",FMortScenarios,"SelScenarios",SelScenarios,'\n')
        set.seed(randvals[k])
        res=SimLenAndAgeFreqData_DynMod_ALB(SimAnnSampSize, nYears, lnSigmaR, autocorr, InitRec,
                                            MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                            RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                            ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                            EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                            ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMortByYear)
        
        
        # fit traditional growth curve
        ObsAge = as.vector(unlist(c(res$ResultsSummary$ObsAgeClRetCatch_Fem[nYears,], res$ResultsSummary$ObsAgeClRetCatch_Mal[nYears,])))
        ObsLen = as.vector(unlist(c(res$ResultsSummary$ObsLenClRetCatchMidPt_Fem[nYears,], res$ResultsSummary$ObsLenClRetCatchMidPt_Mal[nYears,])))
        ObsAge = ObsAge[which(!is.na(ObsAge))]
        ObsLen = ObsLen[which(!is.na(ObsLen))]
        lbns = trunc(ObsLen/LenInc)*LenInc
        
        # set random initial parameter values
        Init_Linf = rnorm(1,600,10)
        Init_vbK = rnorm(1,0.1,0.01)
        Init_tzero = rnorm(1,0,0.01)
        params = c(log(Init_Linf),log(Init_vbK),Init_tzero) # log(Linf), log(k), tzero
        FittedRes=GetvonBertalanffyGrowthResults(params, nSexes=1, DataType=1, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLensd=NA)
        
        plot(ObsAge,ObsLen, ylim=c(0,800), xlim=c(0,100), main="TradVB")
        Ages=1:100
        EstLen = 600 * (1-exp(-0.1*(Ages-0)))
        lines(Ages,EstLen, col="blue")
        EstLen = FittedRes$ParamEst[1,1] * (1-exp(-FittedRes$ParamEst[2,1] * (Ages- FittedRes$ParamEst[3,1])))
        lines(Ages,EstLen, col="black")
        
        vbResults[k,1]=trials
        vbResults[k,2]=RecVarScenarios
        vbResults[k,3]=FMortScenarios
        vbResults[k,4]=SelScenarios
        vbResults[k,5]=round(exp(FittedRes$params[1]),3)
        vbResults[k,6]=round(exp(FittedRes$params[2]),3)
        vbResults[k,7]=round(FittedRes$params[3],3)
        vbResults[k,8]=FittedRes$convergence
        
        # fit length and age-based catch curve
        ObsRetCatchFreqAtLen = as.vector(unlist(res$ResultsSummary$ObsRetCatchFreqAtLen_Fem[nYears,] + res$ResultsSummary$ObsRetCatchFreqAtLen_Mal[nYears,]))
        ObsRetCatchFreqAtLengthAndAge <- data.frame(matrix(nrow = MaxModelAge, ncol = length(lbnd)))
        colnames(ObsRetCatchFreqAtLengthAndAge) <- lbnd
        ObsRetCatchFreqAtLengthAndAge = as.matrix(table(factor(ObsAge, levels=1:MaxModelAge),
                                                        factor(lbns, levels=c(lbnd))))
        
        
        # fit model
        InitFishMort = rnorm(1,0.1,0.01)
        InitL50 = rnorm(1,500,10)
        InitDelta = rnorm(1,50,1) # L95-L50
        InitLinf = rnorm(1,600,10)
        InitvbK = rnorm(1,0.1,0.01)
        InitCVSizeAtAge = rnorm(1,0.05,0.005)
        InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
        params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
        
        LACCResults[k,1]=trials
        LACCResults[k,2]=RecVarScenarios
        LACCResults[k,3]=FMortScenarios
        LACCResults[k,4]=SelScenarios
        LACCResults[k,5]=NA
        LACCResults[k,6]=NA
        LACCResults[k,7]=NA
        LACCResults[k,8]=NA
        LACCResults[k,9]=NA
        LACCResults[k,10]=NA
        LACCResults[k,11]=NA
        
        tryCatch( { 
          FittedRes = list(convergence=NA)
          FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                                          lbnd, ubnd, midpt, SelectivityAtLen, DiscMort, MaxAge, NatMort, TimeStep)        
          
          
        }, warning = function(w) { 
          
        print("Error in LACC fit") })
        
        if (!is.na(FittedRes$convergence)) {
          LACCResults[k,5]=FittedRes$ParamEst[1,1]
          LACCResults[k,6]=FittedRes$ParamEst[2,1]
          LACCResults[k,7]=FittedRes$ParamEst[3,1]
          LACCResults[k,8]=FittedRes$ParamEst[4,1]
          LACCResults[k,9]=FittedRes$ParamEst[5,1]
          LACCResults[k,10]=FittedRes$ParamEst[6,1]
          LACCResults[k,11]=FittedRes$convergence          
          
          plot(ObsAge,ObsLen, ylim=c(0,800), xlim=c(0,100), main="LACC")
          Ages=1:100
          EstLen = 600 * (1-exp(-0.1*(Ages-0)))
          lines(Ages,EstLen, col="blue")
          EstLen = FittedRes$ParamEst[4,1] * (1-exp(-FittedRes$ParamEst[5,1] * Ages))
          lines(Ages,EstLen, col="green")
        }

        
      } # SelScenarios
    } # FMortScenarios
  } # RecVarScenarios
  
  if (SaveFiles) {
    write.csv(vbResults, FileName, row.names = FALSE)
    write.csv(LACCResults, FileName2, row.names = FALSE)
  }
} # trials

# 
# # plot catch curve results
# # growth curve
# PlotAgeLengthCatchCurve_Growth(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#                                lbnd, ubnd, midpt, SelectivityAtLen, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#                                xaxis_lab=NA, yaxis_lab=NA, xmax=100, xint=20,
#                                ymax=800, yint=200, PlotCLs=TRUE, FittedRes, nReps=200, ShowLegend=TRUE)
# 
# 
# # marginal length
# PlotAgeLengthCatchCurve_MargLength(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#                                    lbnd, ubnd, midpt, SelectivityAtLen, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#                                    xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
#                                    ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)
# 
# # conditional age at length - pearson residuals
# PlotAgeLengthCatchCurve_Pears_Resid(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#                                     lbnd, ubnd, midpt, SelectivityAtLen, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#                                     xaxis_lab=NA, yaxis_lab=NA, xmin=200, xmax=800, xint=NA, ymin=0, ymax=100, yint=NA, CircleScale=3, MinLenClFreq=10, ShowLegend=T, FittedRes)
# 
