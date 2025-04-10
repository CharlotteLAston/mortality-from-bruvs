###################################################
# Fitting catch curves to stereo-BRUV and DPIRD
# length data for Western King Wrasse
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
a4.width=160

#### Read in and format data ####
setwd(data_dir)

#* Stereo-BRUV data ####

WKW_bruv_metro <- readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Two.Rocks|TwoRocks|Warnbro|Marmion|Rottnest"))) %>%
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009")) 

WKW_bruv_jurien <- readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) %>%
  # dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Abrolhos_stereo-BRUVs"))) %>%
  # filter(length<=MaxLen) %>% 
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009")) 

#* Fisheries data ####

WKW_dpird_metro <- read.csv("Length data for Charlotte UPDATED.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "King Wrasse") %>% 
  filter(Region %in% "Metro") %>% 
  filter(Year != 2007)

WKW_dpird_jurien <- read.csv("Length data for Charlotte UPDATED.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "King Wrasse") %>% 
  filter(Region %in% "Jurien")

# HistData <- hist(WKW_dpird_metro$TL, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
# HistData

#### Catch curve for metropolitan ####
#* Parameters for Catch Curve that don't change ####
set.seed(123)
MaxAge = 10.2
NatMort = 0.3 #4.22/MaxAge
FishMort = 0.2
MaxLen = 400
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(75, 150) # L50, L95-L50 for gear selectivity
RetenParams = c(100, 75) # L50, L95-L50 for retention

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 484
vbK = 0.16
CVSizeAtAge = 0.07
GrowthParams = c(Linf, vbK)
RefnceAges = NA
DiscMort = 0

#### Estimating the coefficient of variation using simulated data ####
SampleSize = 275
CVSizeAtAge = 0.07
TimeStep=1


Sim_Dat=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                                   SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

LenInterval = 20
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
LenCats

HistData <- hist(Sim_Dat$ObsRandLenRetCatch, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)

#* Stereo BRUVS ####
range(WKW_bruv_metro$length)
LenInterval = 20
LenCats <- seq(from=0, to=MaxLen+40, by=LenInterval)
LenCats


SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1) # selectivity vector
WKW_bruv_sel <- WKW_bruv_metro %>%
  filter(length >280)
HistData <- hist(WKW_bruv_sel$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)

HistData <- hist(WKW_bruv_sel$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
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
InitL50 = 100
InitDelta = 50
DistnType=1 #2
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

good.plot <- ggplot()+
  geom_polygon(aes(x=plotting$x, y=plotting$y), colour=NA, fill="#66CCEE", alpha=0.5)+
  geom_point(aes(y=plotting$plotting.points$Prop, x=plotting$plotting.points$midpoint, 
                 shape=plotting$plotting.points$Obs.Est, group=plotting$plotting.points$Obs.Est), colour="grey40")+
  scale_shape_manual(values = c(1,19), name=NULL, label=c("Estimate", "Observed"))+
  xlab("Length (mm)")+
  ylab("Proportion in length class")+
  ylim(0,0.4)+
  ggplot2::annotate("text", x=100, y=0.4, parse=T, label=as.character(plot.label))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_C-auricularis_bruv.png", height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300 )


#* DPIRD data ####
range(WKW_dpird_metro$TL)
LenInterval = 20
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
LenCats

# HistData <- hist(WKW_dpird_metro$TL, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
# HistData

SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1) # selectivity vector
WKW_dpird_sel <- WKW_dpird_metro %>%
  filter(TL >280)
HistData <- hist(WKW_dpird_sel$TL, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)


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
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
params = c(InitFishMort_logit)

FittedRes_dpird=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes_dpird$ParamEst
FittedRes_dpird$convergence

plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                                 SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                                 RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                                 xaxis_lab=NA, yaxis_lab="Proportion (observed)", xmax=700, xint=50,
                                                 ymax=0.6, yint=0.1, PlotCLs=TRUE, FittedRes_dpird, nReps=200, Error.Colour = "#84D6A4")

plot.label = deparse1(bquote(.(plotting$Fest)))

good.plot <- ggplot()+
  geom_polygon(aes(x=plotting$x, y=plotting$y), colour=NA, fill="#66CCEE", alpha=0.5)+
  geom_point(aes(y=plotting$plotting.points$Prop, x=plotting$plotting.points$midpoint, 
                 shape=plotting$plotting.points$Obs.Est, group=plotting$plotting.points$Obs.Est), colour="grey40")+
  scale_shape_manual(values = c(1,19), name=NULL, label=c("Estimate", "Observed"))+
  xlab("Length (mm)")+
  ylab("Proportion in length class")+
  ylim(0,0.4)+
  ggplot2::annotate("text", x=100, y=0.4, parse=T, label=as.character(plot.label))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_C-auricularis_dpird.png", height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300 )

#### Plotting mortality as a function of length ####

X_PlotLengthBasedCatchCurve_Mortality(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, xmax=NA, xint=NA, ymax=NA, yint=NA, FittedRes_bruv, FittedRes_2=NA, colour_1="#66CCEE", colour_2="white")


#### Catch curve for Jurien ####
#* Parameters for Catch Curve that don't change ####
set.seed(123)
MaxAge = 10.2
NatMort = 0.3 #4.22/MaxAge
FishMort = 0.2
MaxLen = 400
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(250, 20) # L50, L95-L50 for gear selectivity
RetenParams = c(125, 75) # L50, L95-L50 for retention

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 378
vbK = 0.25
GrowthParams = c(Linf, vbK)
RefnceAges = NA
DiscMort = 0
CVSizeAtAge = 0.04

#### Estimating the coefficient of variation using simulated data ####
SampleSize = 800
CVSizeAtAge = 0.04
TimeStep=1


Sim_Dat=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                                   SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

LenInterval = 20
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
LenCats

HistData <- hist(Sim_Dat$ObsRandLenRetCatch, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)


#* Stereo BRUVS ####
range(WKW_bruv_jurien$length)
LenInterval = 15
LenCats <- seq(from=0, to=MaxLen+60, by=LenInterval)
LenCats

HistData <- hist(WKW_bruv_jurien$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData


SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1) # selectivity vector
WKW_bruv_sel <- WKW_bruv_jurien %>%
  filter(length >220)
HistData <- hist(WKW_bruv_sel$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)

# SelectivityType=2

CVSizeAtAge = 0.06
TimeStep = 0.25 # model timestep (e.g. 1 = annual, 1/12 = monthly)
ObsRetCatchFreqAtLen = as.vector(HistData$counts)
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA
midpt = as.vector(HistData$mids)
lbnd = midpt - (LenInterval/2)
ubnd = midpt + (LenInterval/2)
InitFishMort = 0.15 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 200
InitDelta = 50
DistnType=1 #2
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

good.plot <- ggplot()+
  geom_polygon(aes(x=plotting$x, y=plotting$y), colour=NA, fill="#66CCEE", alpha=0.5)+
  geom_point(aes(y=plotting$plotting.points$Prop, x=plotting$plotting.points$midpoint, 
                 shape=plotting$plotting.points$Obs.Est, group=plotting$plotting.points$Obs.Est), colour="grey40")+
  scale_shape_manual(values = c(1,19), name=NULL, label=c("Estimate", "Observed"))+
  xlab("Length (mm)")+
  ylab("Proportion in length class")+
  ylim(0,0.2)+
  ggplot2::annotate("text", x=100, y=0.2, parse=T, label=as.character(plot.label))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_C-auricularis_bruv_jurien.png", height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300 )



#* DPIRD data ####

range(WKW_dpird_jurien$TL)
LenInterval = 20
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
LenCats

HistData <- hist(WKW_dpird_jurien$TL, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData

SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1) # selectivity vector
WKW_dpird_sel <- WKW_dpird_jurien %>%
  filter(TL >240)

HistData <- hist(WKW_dpird_sel$TL, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData

CVSizeAtAge=0.03
TimeStep = 0.5 # model timestep (e.g. 1 = annual, 1/12 = monthly)
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
# InitTheta = 0.3 # specify starting parameters
# InitTheta_logit = log(InitTheta/(1-InitTheta)) # logit transform
# params = c(InitFishMort_logit, log(InitL50), log(InitDelta),InitTheta_logit)
# params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
params = c(InitFishMort_logit)

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

good.plot <- ggplot()+
  geom_polygon(aes(x=plotting$x, y=plotting$y), colour=NA, fill="#BBCC33", alpha=0.5)+
  geom_point(aes(y=plotting$plotting.points$Prop, x=plotting$plotting.points$midpoint, 
                 shape=plotting$plotting.points$Obs.Est, group=plotting$plotting.points$Obs.Est), colour="grey40")+
  scale_shape_manual(values = c(1,19), name=NULL, label=c("Estimate", "Observed"))+
  xlab("Length (mm)")+
  ylab("Proportion in length class")+
  ylim(0,0.4)+
  ggplot2::annotate("text", x=100, y=0.4, parse=T, label=as.character(plot.label))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_C-auricularis_dpird_jurien.png", height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300)


#### Plotting mortality as a function of length ####

X_PlotLengthBasedCatchCurve_Mortality(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, xmax=NA, xint=NA, ymax=NA, yint=NA, FittedRes_bruv, FittedRes_dpird, colour_1="#66CCEE", colour_2="#BBCC33")

#### Try with Pauly catch curve ####
Linf = 484 # Metro
vbK = 0.16

Linf = 378 # Jurien
vbK = 0.25

range(WKW_bruv_jurien$length)
LenInterval = 20
LenCats <- seq(from=0, to=MaxLen+60, by=LenInterval)
LenCats

WKW_bruv_metro_fished <- WKW_bruv_metro %>% 
  filter(length<240)

WKW_bruv_sel <- WKW_bruv_metro %>%
  filter(length >220)
HistData <- hist(WKW_bruv_sel$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)


HistData <- hist(WKW_bruv_sel$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData

ObsRetCatchFreqAtLen = as.vector(HistData$counts)
MinFreq = 10 # set minimum frequency for larger lengths for analysis
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

used_data.plot.bruv <- ggplot()+
  geom_point(aes(x=midpt, y=res$ObsCatchFreqAtLen2),colour="#66CCEE")+
  geom_point(aes(x=midpt, y=ObsRetCatchFreqAtLen), shape=1, colour="grey40")+
  xlab("Length (mm)")+
  ylab("Proportion in length class")+
  theme(axis.text=element_text(size=10))+
  #ylim(0,0.2)+
  # ggplot2::annotate("text", x=300, y=0.2, parse=T, label=as.character(plot.label))+
  theme_classic()
used_data.plot.bruv

res_x = c(res$Age_midptlencl,rev(res$Age_midptlencl)) # using shading for 95% CLs
res_y = c(res$Est_ln_n_dtlow, rev(res$Est_ln_n_dtup))

fit_data.plot.bruv <- ggplot()+
  geom_polygon(aes(x=res_x, y=res_y), colour=NA, fill="#66CCEE", alpha=0.25)+
  geom_point(aes(x=res$Age_midptlencl, y=res$Obs_ln_n_dt),colour="grey40")+
  geom_line(aes(x=res$Age_midptlencl, y=res$Est_ln_n_dt),colour="grey40")+
  xlab("Relative age")+
  ylab("ln(n/dt)")+
  theme(axis.text=element_text(size=10))+
  #ylim(0,0.2)+
  # ggplot2::annotate("text", x=300, y=0.2, parse=T, label=as.character(plot.label))+
  theme_classic()
fit_data.plot.bruv



range(WKW_dpird_jurien$TL)
LenInterval = 20
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
LenCats

HistData <- hist(WKW_dpird_jurien$TL, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData

ObsRetCatchFreqAtLen_dpird = as.vector(HistData$counts)
MinFreq = 10 # set minimum frequency for larger lengths for analysis
# note, this needs to be high enough so that data for ln(n/dt) vs relative age for essentially straight
# line - if not, Z will be biased!!!
midpt_dpird = as.vector(HistData$mids)
lbnd = midpt_dpird - (LenInterval/2)
ubnd = midpt_dpird + (LenInterval/2)
ModelType = 1 # 1 = von Bertalanffy growth curve (Pauly), 2 = length-converted catch curve - Schnute growth curve)
res_dpird=GetLenConvCatchCurveResults(ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen_dpird,
                                MinFreq, lbnd, midpt_dpird, ubnd, MaxAge)
res_dpird$ParamEst

PlotLenConvCatchCurveResults(MaxAge, ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen_dpird, MinFreq,
                             lbnd, midpt_dpird, ubnd)

used_data.plot.dpird <- ggplot()+
  geom_point(aes(x=midpt_dpird, y=res_dpird$ObsCatchFreqAtLen2),colour="#BBCC33")+
  geom_point(aes(x=midpt_dpird, y=ObsRetCatchFreqAtLen_dpird), shape=1, colour="grey40")+
  xlab("Length (mm)")+
  ylab("Proportion in length class")+
  theme(axis.text=element_text(size=10))+
  #ylim(0,0.2)+
  # ggplot2::annotate("text", x=300, y=0.2, parse=T, label=as.character(plot.label))+
  theme_classic()
used_data.plot.dpird

res_x_dpird = c(res_dpird$Age_midptlencl,rev(res_dpird$Age_midptlencl)) # using shading for 95% CLs
res_y_dpird = c(res_dpird$Est_ln_n_dtlow, rev(res_dpird$Est_ln_n_dtup))

fit_data.plot.dpird <- ggplot()+
  geom_polygon(aes(x=res_x_dpird, y=res_y_dpird), colour=NA, fill="#BBCC33", alpha=0.25)+
  geom_point(aes(x=res_dpird$Age_midptlencl, y=res_dpird$Obs_ln_n_dt),colour="grey40")+
  geom_line(aes(x=res_dpird$Age_midptlencl, y=res_dpird$Est_ln_n_dt),colour="grey40")+
  xlab("Relative age")+
  ylab("ln(n/dt)")+
  theme(axis.text=element_text(size=10))+
  #ylim(0,0.2)+
  # ggplot2::annotate("text", x=300, y=0.2, parse=T, label=as.character(plot.label))+
  theme_classic()
fit_data.plot.dpird

setwd(fig_dir)

Pauly_diagnostics <- grid.arrange(arrangeGrob(used_data.plot.bruv + theme(legend.position="none"),
                                              fit_data.plot.bruv + theme(legend.position="none")
                                              # used_data.plot.dpird + theme(legend.position="none"),
                                              # fit_data.plot.dpird + theme(legend.position="none")
                                              ))
ggsave(Pauly_diagnostics, filename="Pauly_diagnostics_KingWrasse_metro.png",height = a4.width*1, width = a4.width*1.1, units  ="mm", dpi = 300 )



