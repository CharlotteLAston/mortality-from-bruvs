###################################################
# Fitting catch curves to stereo-BRUV and DPIRD
# length data for Red Throat Emperor
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

RT_bruv <- readRDS("australian-synthesis_complete_length_lethrinus_miniatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Abrolhos"))) %>%
  dplyr::filter(str_detect(campaign, "2010|2011|2012|2013|2014|2015|2016")) %>% 
  mutate(Method = "Stereo-BRUV")

#* Fisheries data ####

RT_dpird <- read.csv("Redthroat length data for Charlotte.csv") %>% 
  dplyr::select(Species, Year, Location, Lat, Lon, TL,Sector) %>% 
  filter(Location %in% "Abrolhos") %>% 
  filter(Year %in% c(2010, 2011, 2012, 2013,2014,2015,2016)) %>% 
  mutate(Method = "Line fishing")

plot_data2 <- RT_dpird %>% 
  dplyr::select(Method, TL) %>% 
  rename("length"=TL)

plot_data <- RT_bruv %>% 
  dplyr::select(Method, length) %>% 
  rbind(., plot_data2)

Redthroat_plot <- ggplot()+
  geom_histogram(data=plot_data, aes(x=length, fill=Method), position="identity", alpha=0.5)+
  scale_fill_manual(values=c("#BBCC33", "#66CCEE"), labels=c("Line fishing","Stereo-BRUV"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()
Redthroat_plot

#### Parameters for Catch Curve that don't change ####
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

#* Stereo BRUVS ####
range(RT_bruv$length)
LenInterval = 20
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
LenCats

SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve

SelectivityAtLen = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1) # selectivity vector
RT_bruv_sel <- RT_bruv %>%
  filter(length >340)
HistData <- hist(RT_bruv_sel$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)

TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
ObsRetCatchFreqAtLen = as.vector(HistData$counts)
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA
midpt = as.vector(HistData$mids)
lbnd = midpt - (LenInterval/2)
ubnd = midpt + (LenInterval/2)
InitFishMort = 0.2 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 200
InitDelta = 100
DistnType=1 #2
CVSizeAtAge = 0.05
# InitTheta = 0.3 # specify starting parameters
# InitTheta_logit = log(InitTheta/(1-InitTheta)) # logit transform
# params = c(InitFishMort_logit, log(InitL50), log(InitDelta),InitTheta_logit)
# params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
params = c(InitFishMort_logit)

FittedRes_bruv=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                               lbnd, ubnd, midpt, SelectivityAtLen, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes_bruv$ParamEst
FittedRes_bruv$convergence

plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                                 SelectivityVec=SelectivityAtLen, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                                 RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                                 xaxis_lab=NA, yaxis_lab="Proportion (observed)", xmax=700, xint=50,
                                                 ymax=0.5, yint=0.1, PlotCLs=TRUE, FittedRes_bruv, nReps=200, Error.Colour = "#84D6A4")

plot.label = deparse1(bquote(.(plotting$Fest)))

#* DPIRD data ####
range(RT_dpird$TL)
LenInterval = 20
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
LenCats

SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve

SelectivityAtLen = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1) # selectivity vector
RT_dpird_sel <- RT_dpird %>%
  filter(TL >350)
HistData <- hist(RT_dpird$TL, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)

TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
ObsRetCatchFreqAtLen = as.vector(HistData$counts)
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA
midpt = as.vector(HistData$mids)
lbnd = midpt - (LenInterval/2)
ubnd = midpt + (LenInterval/2)
InitFishMort = 0.2 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 300
InitDelta = 100
DistnType=1 #2
CVSizeAtAge = 0.05
# InitTheta = 0.3 # specify starting parameters
# InitTheta_logit = log(InitTheta/(1-InitTheta)) # logit transform
# params = c(InitFishMort_logit, log(InitL50), log(InitDelta),InitTheta_logit)
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
# params = c(InitFishMort_logit)

FittedRes_dpird=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                               lbnd, ubnd, midpt, SelectivityAtLen, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes_dpird$ParamEst
FittedRes_dpird$convergence

plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                                 SelectivityVec=SelectivityAtLen, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                                 RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                                 xaxis_lab=NA, yaxis_lab="Proportion (observed)", xmax=700, xint=50,
                                                 ymax=0.5, yint=0.1, PlotCLs=TRUE, FittedRes_dpird, nReps=200, Error.Colour = "#84D6A4")

plot.label = deparse1(bquote(.(plotting$Fest)))

#### Length converted catch curve ####
RT_bruv_sel <- RT_bruv %>%
  filter(length >340)
HistData <- hist(RT_bruv_sel$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)

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

RT_dpird_sel <- RT_dpird %>%
  filter(TL >340)
HistData <- hist(RT_dpird_sel$TL, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)

ObsRetCatchFreqAtLen = as.vector(HistData$counts)
MinFreq = 1 # set minimum frequency for larger lengths for analysis
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

