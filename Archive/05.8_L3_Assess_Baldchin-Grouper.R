rm(list=ls())
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)
library(L3Assess)
library(dplyr)
library(stringr)
library(sf)
library(ggplot2)

# Simulate data
SampleSize=500 # sample size for retained catches (and same number for released fish, if an MLL is specified)

MaxAge = 22 # Wise et al. (2007) in https://ap-st01.ext.exlibrisgroup.com/61MUN_INST/upload/1719386817922/baldchin_groper.pdf?Expires=1719386938&Signature=AuRiuB690uPhl9E4m~YmGNobHDgyxJeAodOTwOiAr8ES0-zi~6LhhbmspqLPTTLPCQM4zH4iwJjM69G-Oce30moZIws2cUfHOAM0rvaYQaPSorJ7dkyTfzm0EW73SoUpJDq0uYjXhdJcyqRkw7qPX~2fz5Wij-~cimra61RKpunoXSZjVwWRT9xRHwbx2IvGRDod1b8XlbBwBPsKbIhTtYeZmWUTKtfyCPMLLaq-Ayl5KMPJE5683sd9ZR-MtlRYkmybW4aKK7E0~1YFOwyuBLNLXalXZse6jwPiN9dY1qVTg4KGMbYWizA9bFj39JYJWfDI7zqn9qnBMUyCFLZYzQ__&Key-Pair-Id=APKAJ72OZCZ36VGVASIA
TimeStep = 0.5 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 0.21 # Wise et al. (2007) in https://ap-st01.ext.exlibrisgroup.com/61MUN_INST/upload/1719386817922/baldchin_groper.pdf?Expires=1719386938&Signature=AuRiuB690uPhl9E4m~YmGNobHDgyxJeAodOTwOiAr8ES0-zi~6LhhbmspqLPTTLPCQM4zH4iwJjM69G-Oce30moZIws2cUfHOAM0rvaYQaPSorJ7dkyTfzm0EW73SoUpJDq0uYjXhdJcyqRkw7qPX~2fz5Wij-~cimra61RKpunoXSZjVwWRT9xRHwbx2IvGRDod1b8XlbBwBPsKbIhTtYeZmWUTKtfyCPMLLaq-Ayl5KMPJE5683sd9ZR-MtlRYkmybW4aKK7E0~1YFOwyuBLNLXalXZse6jwPiN9dY1qVTg4KGMbYWizA9bFj39JYJWfDI7zqn9qnBMUyCFLZYzQ__&Key-Pair-Id=APKAJ72OZCZ36VGVASIA
FishMort = 0.2 
MaxLen = 700 #Wise et al. (2007) in https://ap-st01.ext.exlibrisgroup.com/61MUN_INST/upload/1719386817922/baldchin_groper.pdf?Expires=1719386938&Signature=AuRiuB690uPhl9E4m~YmGNobHDgyxJeAodOTwOiAr8ES0-zi~6LhhbmspqLPTTLPCQM4zH4iwJjM69G-Oce30moZIws2cUfHOAM0rvaYQaPSorJ7dkyTfzm0EW73SoUpJDq0uYjXhdJcyqRkw7qPX~2fz5Wij-~cimra61RKpunoXSZjVwWRT9xRHwbx2IvGRDod1b8XlbBwBPsKbIhTtYeZmWUTKtfyCPMLLaq-Ayl5KMPJE5683sd9ZR-MtlRYkmybW4aKK7E0~1YFOwyuBLNLXalXZse6jwPiN9dY1qVTg4KGMbYWizA9bFj39JYJWfDI7zqn9qnBMUyCFLZYzQ__&Key-Pair-Id=APKAJ72OZCZ36VGVASIA
LenInc = 30
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityVec = NA # selectivity vector
SelParams = c(300, 150) # L50, L95-L50 for gear selectivity
RetenParams = c(NA, NA) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute

Linf = 534.7 #Wise et al. (2007) in https://ap-st01.ext.exlibrisgroup.com/61MUN_INST/upload/1719386817922/baldchin_groper.pdf?Expires=1719386938&Signature=AuRiuB690uPhl9E4m~YmGNobHDgyxJeAodOTwOiAr8ES0-zi~6LhhbmspqLPTTLPCQM4zH4iwJjM69G-Oce30moZIws2cUfHOAM0rvaYQaPSorJ7dkyTfzm0EW73SoUpJDq0uYjXhdJcyqRkw7qPX~2fz5Wij-~cimra61RKpunoXSZjVwWRT9xRHwbx2IvGRDod1b8XlbBwBPsKbIhTtYeZmWUTKtfyCPMLLaq-Ayl5KMPJE5683sd9ZR-MtlRYkmybW4aKK7E0~1YFOwyuBLNLXalXZse6jwPiN9dY1qVTg4KGMbYWizA9bFj39JYJWfDI7zqn9qnBMUyCFLZYzQ__&Key-Pair-Id=APKAJ72OZCZ36VGVASIA

vbK = 0.19 #Wise et al. (2007) in https://ap-st01.ext.exlibrisgroup.com/61MUN_INST/upload/1719386817922/baldchin_groper.pdf?Expires=1719386938&Signature=AuRiuB690uPhl9E4m~YmGNobHDgyxJeAodOTwOiAr8ES0-zi~6LhhbmspqLPTTLPCQM4zH4iwJjM69G-Oce30moZIws2cUfHOAM0rvaYQaPSorJ7dkyTfzm0EW73SoUpJDq0uYjXhdJcyqRkw7qPX~2fz5Wij-~cimra61RKpunoXSZjVwWRT9xRHwbx2IvGRDod1b8XlbBwBPsKbIhTtYeZmWUTKtfyCPMLLaq-Ayl5KMPJE5683sd9ZR-MtlRYkmybW4aKK7E0~1YFOwyuBLNLXalXZse6jwPiN9dY1qVTg4KGMbYWizA9bFj39JYJWfDI7zqn9qnBMUyCFLZYzQ__&Key-Pair-Id=APKAJ72OZCZ36VGVASIA
CVSizeAtAge = 0.05
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
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=700, xint=50,
                                   ymax=0.15, yint=0.05, PlotCLs=TRUE, FittedRes, nReps=200)



# ***************************
# try fitting to 'real' data
# ***************************
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "Data", sep="/")

setwd(data_dir)

dat <- readRDS("australian-synthesis_complete_length_choerodon-rubescens.RDS") %>% 
  #dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"TwoRocks|Marmion|Warnbro|Two.Rocks|NgariCapes|south-west|Capes|SwC|Geographe.Bay"))) %>% 
  filter(length<MaxLen) %>% 
  filter(campaign %in% c("2007-03_Capes.MF_stereoBRUVs", "2007-05_Abrolhos.MF_stereoBRUVs", "2007-05_Jurien.Bay.MF_stereoBRUVs", "2007-09_Rottnest.MF_stereoBRUVs",
                         "2008-10_Two.Rocks.WAMSI_stereoBRUVs","2008-04_Rottnest.MF_stereoBRUVs", "2008-04_Capes.Wrights.Bank_stereoBRUVs","2009-03_Capes.sanctuaries_stereoBRUVs",
                         "2009-07_Rottnest.sanctuaries_stereoBRUVs", "2010-04_Rottnest.sanctuaries_stereoBRUVs", "2010-04_Rottnest.WAMSI_stereoBRUVs", "2010-04_TwoRocks.WAMSI_stereoBRUVs",
                         "2010-05_Capes.WAMSI_stereoBRUVs", "2011-02_Rottnest.sanctuaries_stereoBRUVs", "2012-06_Marmion.Warnbro.Cross.Shelf_stereoBRUVs", "2013-01_Capes_stereoBRUVs",
                         "2013-06_Rottnest.MF_stereoBRUVs","2013-06_Warnbro_stereoBRUVs", "2014-12_Geographe.Bay_stereoBRUVs", "2019-02_NgariCapes_stereoBRUVs", "2020-06_south-west_stereo-BRUVs",
                         "2020-10_south-west_stereo-BRUVs", "2021-05_Jurien.Bay.MP.Monitoring_stereoBRUVs", "2023-03_SwC_stereo-BRUVs")) %>%
  filter_all(.vars_predicate = any_vars(str_detect(., "TwoRocks|Marmion|Warnbro|Two.Rocks|NgariCapes|south-west|Capes|SwC|Geographe.Bay|Rottnest"))) #%>%
  # filter(!campaign %in% c("2019-02_NgariCapes_stereoBRUVs", "2020-06_south-west_stereo-BRUVs","2020-10_south-west_stereo-BRUVs", "2023-03_SwC_stereo-BRUVs"))


DBCA_data <- readRDS("DBCA_lengths.RDS") %>% 
  filter(scientific %in% "Choerodon rubescens")

head(dat)
range(dat$length)
LenInterval = 25
LenCats <- seq(from=0, to=MaxLen, by=LenInterval)
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
InitL50 = 340
InitDelta = 140
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

