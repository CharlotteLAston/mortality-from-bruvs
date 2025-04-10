
# Alex Hesp 29 Jan 2025
# For Charlotte
# Exploration of an approach for deriving estimates of F, retention parameters and selectivity
# parameters from line and BRUV data
# 
# install.packages("C:/~/L3Assess_0.1.0.tar.gz", source = TRUE, repos=NULL)


library(L3Assess)
library(tidyverse)
library(dplyr)

rm(list=ls())


#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_data", sep="/")

setwd(working.dir)
source("X_My-Functions.R")
a4.width=160


# if installing from github
# library(devtools)
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)


# ******************************
# first, simulation test method
# ******************************

# produce some length data first for line data (for retained fish only),
# to estimate F and selectivity of landings/retention
# then some length data for bruv data (to estimate selectivity, using
# retention function estimated from line data to apply F to proportion
# of fish population exposed to fishing, and estimating 'gear' selectivity of BRUV gear)

# aiming to see if specified F, retention params (from line data) and selectivity of BRUV gear
# (from BRUV data) can be recovered.

# Simulate data
set.seed(123)
SampleSize=1000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
MaxAge = 20 # max model age
TimeStep = 1/2 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 0.15
FishMort = 0.2
MaxLen = 650
LenInc = 20
lbnd = seq(0,MaxLen-LenInc,LenInc)
midpt = lbnd + LenInc/2
ubnd = lbnd + LenInc
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(150, 50) # L50, L95-L50 for gear selectivity
RetenParams = c(300, 10) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality
# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 506
vbK = 0.253
CVSizeAtAge = 0.07
GrowthParams = c(Linf, vbK)
RefnceAges = NA
DiscMort = 0

# generate line data
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# plot selectivity assumptions
plot(Res$midpt, Res$ModelDiag$SelAtLength, "l")
lines(Res$midpt, Res$ModelDiag$RetAtLength, col="blue") # Retained fish L50 ~200
lines(Res$midpt, Res$ModelDiag$SelLandAtLength, col="green") # Landings is about the same because all fish are retained

# plot length-frequency data (retained fish - line)
lendat_line = Res$ObsRetCatchFreqAtLen
lendat_line2 = Res$ObsDiscCatchFreqAtLen
plot(Res$midpt, lendat_line, "o", ylim=c(0,250))

# generate BRUV data (with different realisation, i.e. starting with different random number)
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# compare length-frequency data - BRUV
lendat_bruv = Res$ObsRetCatchFreqAtLen # Retained by line fishing
lendat_bruv2 = Res$ObsDiscCatchFreqAtLen # Discarded fish - difference between retained and selectivity?
lendat_bruv3 = lendat_bruv + lendat_bruv2
lines(Res$midpt, lendat_bruv3, "o", col="red") # Gives us those extra fish before they are retained by the fishery


# *************************************************************
# specify biology and fishery parameters for  Redthroat Emperor
# *************************************************************
setwd(data_dir)

RT_bruv <- readRDS("australian-synthesis_complete_length_lethrinus_miniatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Abrolhos"))) %>%
  dplyr::filter(str_detect(campaign, "2021")) %>% 
  mutate(Method = "Stereo-BRUV")


set.seed(123)
MaxAge = 20 # https://fish.gov.au/report/364-Redthroat-Emperor-2023#!
NatMort = 0.15 #https://www.fish.gov.au/2020-Reports/redthroat_emperor
FishMort = 0.2
MaxLen = 650 ## https://www.fish.gov.au/report/364-Redthroat-Emperor-2023?jurisdictionId=2#!
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(75, 150) # L50, L95-L50 for gear selectivity
RetenParams = c(100, 75) # L50, L95-L50 for retention

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 506
vbK = 0.253
CVSizeAtAge = 0.07
GrowthParams = c(Linf, vbK)
RefnceAges = NA
DiscMort = 0

# ************************
# check growth assumptions
# ************************
SampleSize=1000
FishMort=0.1
SelParams = c(300,10) 
RetenParams=NA
SimRes=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                                  SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
# plot data
PlotOpt=0 # 0=all plots, 1=retained lengths at age and growth curves, 2=retained plus discard lengths at age
# and growth curves, 3=retained plus discard length frequency, 4=retained plus discard age frequency,
# 5=sex specific retained lengths at age and growth curves, 6=retained plus discard lengths at age and growth curves,
# 7=female length frequency, 8=male length frequency, 9=female age frequency, 10=male age frequency,
# 11=selectivity/retention, 12=F-at-age reten + disc, 13=F-at-age reten, 14=F-at-age disc
par(mfcol=c(1,1))
PlotSimLenAndAgeFreqData_EqMod(MaxAge, MaxLen, SimRes, PlotOpt=2)

# Check that this looks like the plots in the papers

# *******************************
# fit LBCC to simulated line data
# *******************************

ObsRetCatchFreqAtLen  = lendat_line
DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
RetenAtLen = NA
PropReleased=NA

# specify initial params
InitFishMort = 0.25 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 300
InitDelta = 10
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$convergence
FittedRes$ParamEst
#            Estimate lw_95%CL up_95%CL
# FMort        0.171    0.141    0.205
# L50_sel    122.871  109.320  138.101
# Delta_sel   86.142   64.496  115.054

# # True values
# FishMort = 0.1
# RetenParams = c(300, 10) # L50, L95-L50 for retention

# so, modeldoesn't capture F very well because selection is late 
# selectivity of line pretty well estimated as well
# see difference (for F estimate) if retention is set to be lower, to simulate data

# plot results
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                   ymax=0.50, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)
PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                      ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)

# *******************************
# fit LBCC to simulated BRUV data
# *******************************

ObsRetCatchFreqAtLen  = lendat_bruv3
DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
L50_Ret = FittedRes$ParamEst[2,1] # Put retention in because this shows us where to apply the fishing mortality
L95_Ret =  FittedRes$ParamEst[3,1] + L50_Ret
RetenAtLen = 1 / (1+exp(-log(19)*(midpt-L50_Ret)/(L95_Ret-L50_Ret)))
plot(midpt,RetenAtLen)

# specify initial params
InitFishMort = 0.25 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 100
InitDelta = 50
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$convergence
FittedRes$ParamEst
#           Estimate lw_95%CL up_95%CL
# FMort        0.202    0.138    0.286
# L50_sel     80.455   45.453  142.410
# Delta_sel  212.192   86.864  518.340

# # True values
# FishMort = 0.2 wasn't bad
# RetenParams = c(122, 87) # L50, L95-L50 for retention

# similarly, model not able to capture F well, due to very late selectivity
# selectivity of BRUVs not well estiamted 
# see difference (for F estimate) if retention is set to be lower, to simulate data

# plot results
# note, not really retained catch here (not relabelled) is selectivity of BRUVs
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                   ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)
PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                      ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)


# *****************************
# fit LBCC - fishery line data
# *****************************
setwd(data_dir)
# get data
dat <-read.csv("Length data for Charlotte.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "Redthroat Emperor") %>% 
  filter(Region %in% "Abrolhos")

# which(is.na(dat$TL))
range(RT_bruv$length)
table(RT_bruv$length)
FishLengths = RT_bruv$length
lbns = trunc(FishLengths/LenInc)*LenInc
ObsRetCatchFreqAtLen  = as.vector(table(factor(lbns, levels=lbnd)))
plot(midpt, ObsRetCatchFreqAtLen, "l") # selectivity affected by small hooks

DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
RetenAtLen = NA # proportion of fish released, vector including mean and sd (option probably now obsolete)

# specify initial params
InitFishMort = 0.25 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 250
InitDelta = 50
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$convergence
FittedRes$ParamEst
# > FittedRes$ParamEst
#            Estimate lw_95%CL up_95%CL
# FMort        0.897    0.001    1.000
# L50_sel    400.678  351.933  456.174
# Delta_sel   86.164   64.979  114.254

# Useless estimate and couldn't estimate the selectivity very well 

# plot results
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                   ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)

PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                      ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)


# save selectivity params from line data analysis.
# these assumed to be equivalent to retention parameters,
# and will be used in next analysis to fix retention

plot(midpt, ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen), "l", ylim=c(0,0.2)) # line data
sum(ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen))

# *****************************
# fit LBCC - BRUV data
# *****************************

L50_Ret = FittedRes$ParamEst[2,1] # Here saying that "retention" is the fish that get caught by the fishing gear so where mortality would be applied
L95_Ret =  FittedRes$ParamEst[3,1] + L50_Ret
RetenAtLen = 1 / (1+exp(-log(19)*(midpt-L50_Ret)/(L95_Ret-L50_Ret)))
# plot(midpt,RetenAtLen)

# get data
dat2 <-readRDS("australian-synthesis_complete_length_epinephelides-armatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Marmion|TwoRocks|Two.Rocks|Warnbro|Rottnest"))) %>% #Capes|Geographe|SwC|south-west|Ngari|
  dplyr::filter(str_detect(campaign,"2006|2007|2008|2009|2010|2011")) %>% 
  dplyr::filter(!str_detect(campaign, "Perth.Canyon")) 

# which(is.na(dat2$length))
head(dat2)
range(dat2$length)
table(dat2$length)
FishLengths = dat2$length
lbns = trunc(FishLengths/LenInc)*LenInc

ObsRetCatchFreqAtLen  = as.vector(table(factor(lbns, levels=lbnd)))

# compare line vs BRUV length data
lines(midpt,ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen), col="blue") # bruv data has very different shape to the line data
sum(ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen))


# specify initial params
InitFishMort = 0.25 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 100 # Selectivity of the BRUVs
InitDelta = 50
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$convergence
FittedRes$ParamEst



# > FittedRes$ParamEst
#            Estimate lw_95%CL up_95%CL
# FMort        0.363    0.275    0.460
# L50_sel    207.229  188.033  228.386
# Delta_sel   81.372   61.982  106.829

# Gives us a pretty reasonable estimate of F 

# plot results
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                   ymax=0.50, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)

PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                      ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)

plot(midpt, lendat_bruv3/sum(lendat_bruv3), type="l", ylim=c(0,0.2))
lines(midpt,ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen),col="blue") #Actual observed data



