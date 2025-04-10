###################################################
# Fitting catch curves to stereo-BRUV and DPIRD
# length data for Brownspot Wrasse
###################################################
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
sp_dir <- paste(working.dir, "Spatial_data", sep="/")

setwd(working.dir)
source("X_My-Functions.R")

#### Read in and format data ####
setwd(data_dir)

#* Stereo-BRUV data ####

BW_bruv_metro <- readRDS("australian-synthesis_complete_length_notolabrus-parilus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Two.Rocks|TwoRocks|Warnbro|Marmion|Rottnest"))) %>%
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009")) %>% 
  filter(length<410)

BW_bruv_jurien <- readRDS("australian-synthesis_complete_length_notolabrus-parilus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) %>%
  # dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Abrolhos_stereo-BRUVs"))) %>%
  # filter(length<=MaxLen) %>% 
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009")) 

#* Fisheries data ####

BW_dpird_metro <- read.csv("Length data for Charlotte with brownspot.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "Brownspot Wrasse") %>% 
  filter(Region %in% "Metro")

BW_dpird_jurien <- read.csv("Length data for Charlotte with brownspot.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "Brownspot Wrasse") %>% 
  filter(Region %in% "Jurien")


#### Catch curve for metropolitan ####
#* Parameters for Catch Curve that don't change ####
set.seed(123)
MaxAge = 11.8
NatMort = 4.22/MaxAge
FishMort = 0.2
MaxLen = 410
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 410
vbK = 0.15
GrowthParams = c(Linf, vbK)
RefnceAges = NA
DiscMort = 0

#### Estimating the coefficient of variation using simulated data ####
SampleSize = 495
CVSizeAtAge = 0.05
TimeStep=1


Sim_Dat=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                                   SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

LenInterval = 20
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
LenCats

HistData <- hist(Sim_Dat$ObsRandLenRetCatch, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)

#* Stereo BRUVS ####
range(BW_bruv_metro$length)
LenInterval = 15
LenCats <- seq(from=0, to=MaxLen+15, by=LenInterval)
LenCats

HistData <- hist(BW_bruv_metro$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData

#* Set selectivity and filter data
SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1) # selectivity vector
SelParams = c(NA, NA) #c(135, 50) # L50, L95-L50 for gear selectivity
RetenParams = c(NA, NA) #c(190, 75) # L50, L95-L50 for retention

BW_bruv_metro_trunc <- BW_bruv_metro %>% 
  filter(length > 285)

HistData <- hist(BW_bruv_metro_trunc$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData

TimeStep = 0.5 # model timestep (e.g. 1 = annual, 1/12 = monthly)
ObsRetCatchFreqAtLen = as.vector(HistData$counts) 
  
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA
midpt = as.vector(HistData$mids)
lbnd = midpt - (LenInterval/2)
ubnd = midpt + (LenInterval/2)
InitFishMort = 0.2 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 250
InitDelta = 50
DistnType=1 #2
# InitTheta = 0.3 # specify starting parameters
# InitTheta_logit = log(InitTheta/(1-InitTheta)) # logit transform
# params = c(InitFishMort_logit, log(InitL50), log(InitDelta),InitTheta_logit)
params = c(InitFishMort_logit, log(InitDelta))


FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$ParamEst


plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                                 SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                                 RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                                 xaxis_lab=NA, yaxis_lab="Proportion (observed)", xmax=700, xint=50,
                                                 ymax=0.5, yint=0.1, PlotCLs=TRUE, FittedRes, nReps=200, Error.Colour = "#84D6A4")

#* Pauly ####
LenInterval = 20
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
LenCats

HistData <- hist(BW_bruv_metro$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData

ObsRetCatchFreqAtLen = as.vector(HistData$counts)
MinFreq = 5 # set minimum frequency for larger lengths for analysis
# note, this needs to be high enough so that data for ln(n/dt) vs relative age for essentially straight
# line - if not, Z will be biased!!!
midpt = as.vector(HistData$mids)
lbnd = midpt - (LenInterval/2)
ubnd = midpt + (LenInterval/2)
ModelType = 1 # 1 = von Bertalanffy growth curve (Pauly), 2 = length-converted catch curve - Schnute growth curve)
res=GetLenConvCatchCurveResults(ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                MinFreq, lbnd, midpt, ubnd, MaxAge)
res$ParamEst


PlotLenConvCatchCurveResults(MaxAge, ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen, MinFreq,
                             lbnd, midpt, ubnd)



#* DPIRD data ####
range(BW_dpird_metro$TL)
LenInterval = 15
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
LenCats

HistData <- hist(BW_dpird_metro$TL, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData

TimeStep = 0.25 # model timestep (e.g. 1 = annual, 1/12 = monthly)
ObsRetCatchFreqAtLen = as.vector(HistData$counts)
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA
midpt = as.vector(HistData$mids)
lbnd = midpt - (LenInterval/2)
ubnd = midpt + (LenInterval/2)
InitFishMort = 0.2 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 250
InitDelta = 50
DistnType=1 #2
# InitTheta = 0.3 # specify starting parameters
# InitTheta_logit = log(InitTheta/(1-InitTheta)) # logit transform
# params = c(InitFishMort_logit, log(InitL50), log(InitDelta),InitTheta_logit)
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))

FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$ParamEst
FittedRes$RetCatch_SampleSize
FittedRes$RetCatch_EffSampleSize

plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                                 SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                                 RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                                 xaxis_lab=NA, yaxis_lab="Proportion (observed)", xmax=700, xint=50,
                                                 ymax=0.6, yint=0.1, PlotCLs=TRUE, FittedRes, nReps=200, Error.Colour = "#84D6A4")

ObsRetCatchFreqAtLen = as.vector(HistData$counts)
MinFreq = 5 # set minimum frequency for larger lengths for analysis
# note, this needs to be high enough so that data for ln(n/dt) vs relative age for essentially straight
# line - if not, Z will be biased!!!
midpt = as.vector(HistData$mids)
lbnd = midpt - (LenInterval/2)
ubnd = midpt + (LenInterval/2)
ModelType = 1 # 1 = von Bertalanffy growth curve (Pauly), 2 = length-converted catch curve - Schnute growth curve)
res=GetLenConvCatchCurveResults(ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                MinFreq, lbnd, midpt, ubnd, MaxAge)
res$ParamEst


PlotLenConvCatchCurveResults(MaxAge, ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen, MinFreq,
                             lbnd, midpt, ubnd)


#* Jurien ####
range(BW_dpird_jurien$TL)
LenInterval = 15
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
LenCats

HistData <- hist(BW_dpird_jurien$TL, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData

TimeStep = 1# model timestep (e.g. 1 = annual, 1/12 = monthly)
ObsRetCatchFreqAtLen = as.vector(HistData$counts)
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA
midpt = as.vector(HistData$mids)
lbnd = midpt - (LenInterval/2)
ubnd = midpt + (LenInterval/2)
InitFishMort = 0.2 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 250
InitDelta = 50
DistnType=1 #2
# InitTheta = 0.3 # specify starting parameters
# InitTheta_logit = log(InitTheta/(1-InitTheta)) # logit transform
# params = c(InitFishMort_logit, log(InitL50), log(InitDelta),InitTheta_logit)
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))

FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$ParamEst
FittedRes$RetCatch_SampleSize
FittedRes$RetCatch_EffSampleSize

plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                                 SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                                 RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                                 xaxis_lab=NA, yaxis_lab="Proportion (observed)", xmax=700, xint=50,
                                                 ymax=0.6, yint=0.1, PlotCLs=TRUE, FittedRes, nReps=200, Error.Colour = "#84D6A4")


ObsRetCatchFreqAtLen = as.vector(HistData$counts)
MinFreq = 20 # set minimum frequency for larger lengths for analysis
# note, this needs to be high enough so that data for ln(n/dt) vs relative age for essentially straight
# line - if not, Z will be biased!!!
midpt = as.vector(HistData$mids)
lbnd = midpt - (LenInterval/2)
ubnd = midpt + (LenInterval/2)
ModelType = 1 # 1 = von Bertalanffy growth curve (Pauly), 2 = length-converted catch curve - Schnute growth curve)
res=GetLenConvCatchCurveResults(ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                MinFreq, lbnd, midpt, ubnd, MaxAge)
res$ParamEst


PlotLenConvCatchCurveResults(MaxAge, ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen, MinFreq,
                             lbnd, midpt, ubnd)


