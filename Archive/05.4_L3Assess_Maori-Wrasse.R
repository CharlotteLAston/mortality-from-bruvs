rm(list=ls())
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)
library(L3Assess)
library(dplyr)
library(stringr)
library(sf)
library(ggplot2)

# Simulate data
SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
set.seed(123)
MaxAge = 11.8 # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1095-8649.2012.03446.x
TimeStep = 0.5 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/MaxAge
FishMort = 0.2 * NatMort
MaxLen = 400 # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1095-8649.2012.03446.x
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityVec = NA # selectivity vector
SelParams = c(140, 75) # L50, L95-L50 for gear selectivity
RetenParams = c(NA, NA) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
# Linf = 664
Linf = 328 # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1095-8649.2012.03446.x
# vbK = 0.241
vbK = 0.28 # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1095-8649.2012.03446.x
CVSizeAtAge = 0.1
GrowthParams = c(Linf, vbK)
RefnceAges = NA
Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                         SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# fit catch curve to simulated data
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
midpt=Res$midpt
lbnd=Res$lbnd
ubnd=Res$ubnd
ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
length(ObsRetCatchFreqAtLen)
length(midpt)
InitFishMort = 0.1 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 140
InitDelta = 70
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
DistnType = 1

FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$ResultsSummary
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=500, xint=50,
                                   ymax=0.15, yint=0.05, PlotCLs=TRUE, FittedRes, nReps=200)


# ***************************
# try fitting to 'real' data
# ***************************
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "Data", sep="/")

#* Metro ####

setwd(data_dir)

dat <- readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"2006|2007|2008|2009|2010|2011"))) %>% 
  # dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Rottnest|TwoRocks|Marmion|Warnbro|Two.Rocks"))) %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) #%>% 
  # filter(length<MaxLen)

dat <- readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"2006|2007|2008|2009|2010|2011"))) %>% 
  # dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Rottnest|TwoRocks|Marmion|Warnbro|Two.Rocks"))) %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Capes"))) %>% 
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009|2010|2011"))
# filter(length<MaxLen)
  
head(dat)
range(dat$length)
LenInterval = 20
LenCats <- seq(from=0, to=MaxLen+20, by=LenInterval)
LenCats


HistData <- hist(dat$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData

ObsRetCatchFreqAtLen = as.vector(HistData$counts)
midpt = as.vector(HistData$mids)
lbnd = midpt - (LenInterval/2)
ubnd = midpt + (LenInterval/2)
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
InitFishMort = 0.1 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 140
InitDelta = 60
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
DistnType = 1
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
CVSizeAtAge = 0.05 #0.025
TimeStep = 0.5

FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)

FittedRes$ResultsSummary


plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                     SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                     RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                     xaxis_lab=NA, yaxis_lab="Proportion (observed)", xmax=MaxLen, xint=50,
                                     ymax=0.5, yint=0.1, PlotCLs=TRUE, FittedRes, nReps=200, Error.Colour = "#FFED6F")

plot.label = deparse1(bquote(.(plotting$Fest)))

good.plot <- ggplot()+
  geom_polygon(aes(x=plotting$x, y=plotting$y), colour=NA, fill="#FFED6F")+
  geom_point(aes(y=plotting$plotting.points$Prop, x=plotting$plotting.points$midpoint, 
                 shape=plotting$plotting.points$Obs.Est, group=plotting$plotting.points$Obs.Est), colour="grey20")+
  scale_shape_manual(values = c(1,19), name=NULL, label=c("Estimate", "Observed"))+
  xlab("Length (mm)")+
  ylab("Proportion in length class")+
  ggplot2::annotate("text", x=100, y=0.3, parse=T, label=as.character(plot.label))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_O-lineolatus_Metro.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#* Capes ####

setwd(data_dir)

dat <- readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Capes|south-west|Ngari"))) %>%
  filter(str_detect(campaign, "2014|2015|2016|2016|2017|2018|2019")) %>% 
  filter(length<MaxLen) %>% 
  dplyr::select(campaign, length, sample) %>% 
  rename(campaignid = "campaign")

DBCA_data <- readRDS("DBCA_lengths.RDS") %>% 
  filter(scientific %in% "Ophthalmolepis lineolatus") %>% 
  filter(campaignid %in% "2019-02_NgariCapes_stereoBRUVs")%>% 
  dplyr::select(campaignid, length, sample) %>% 
  filter(length<MaxLen)

dat <- rbind(dat, DBCA_data)

head(dat)
range(dat$length)
LenInterval = 15
LenCats <- seq(from=0, to=MaxLen+15, by=LenInterval)
LenCats


HistData <- hist(dat$length, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="15 mm Length category", ylab="Frequency",las=1)
HistData

ObsRetCatchFreqAtLen = as.vector(HistData$counts)
midpt = as.vector(HistData$mids)
lbnd = midpt - (LenInterval/2)
ubnd = midpt + (LenInterval/2)
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
InitFishMort = 0.2 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 115
InitDelta = 55
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
DistnType = 1
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
CVSizeAtAge = 0.05 #0.025
TimeStep = 1/12

FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)

FittedRes$ResultsSummary


plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                     SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                     RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                     xaxis_lab=NA, yaxis_lab="Proportion (observed)", xmax=500, xint=50,
                                     ymax=0.5, yint=0.1, PlotCLs=TRUE, FittedRes, nReps=200, Error.Colour = "#FEB461")

plot.label = deparse1(bquote(.(plotting$Fest)))

good.plot <- ggplot()+
  geom_polygon(aes(x=plotting$x, y=plotting$y), colour=NA, fill="#FEB461")+
  geom_point(aes(y=plotting$plotting.points$Prop, x=plotting$plotting.points$midpoint, 
                 shape=plotting$plotting.points$Obs.Est, group=plotting$plotting.points$Obs.Est), colour="grey20")+
  scale_shape_manual(values = c(1,19), name=NULL, label=c("Estimate", "Observed"))+
  xlab("Length (mm)")+
  ylab("Proportion in length class")+
  ggplot2::annotate("text", x=100, y=0.3, parse=T, label=as.character(plot.label))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_O-lineolatus_Capes.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )
