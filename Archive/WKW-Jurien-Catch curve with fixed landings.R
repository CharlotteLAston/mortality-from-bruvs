
# Alex Hesp 29 Jan 2025
# For Charlotte
# Exploration of an approach for deriving estimates of F, retention parameters and selectivity
# parameters from line and BRUV data
# 
install.packages("C:/~/L3Assess_0.1.0.tar.gz", source = TRUE, repos=NULL)


rm(list=ls())
library(L3Assess)
library(tidyverse)
library(dplyr)
library(scales)


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
TimeStep = 1/4 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/10.2
FishMort = 0.1
MaxLen = 500
LenInc = 10
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(100, 20) # L50, L95-L50 for gear selectivity
RetenParams = c(200, 100) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality
# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 378
vbK = 0.25
CVSizeAtAge = 0.025
GrowthParams = c(Linf, vbK)
RefnceAges = NA
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
plot(Res$midpt, lendat_line, "o", ylim=c(0,100))

# generate BRUV data (with different realisation, i.e. starting with different random number)
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# compare length-frequency data - BRUV
lendat_bruv = Res$ObsRetCatchFreqAtLen # Retained by line fishing
lendat_bruv2 = Res$ObsDiscCatchFreqAtLen # Discarded fish - difference between retained and selectivity?
lendat_bruv3 = lendat_bruv + lendat_bruv2
lines(Res$midpt, lendat_bruv3, "o", col="red") # Gives us those extra fish before they are retained by the fishery


# *********************************************************
# specify biology and fishery parameters for WKW in Jurien
# *********************************************************
MaxAge = 12
MaxObsAge = 10.2
TimeStep = 1/4 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = exp(1.551-1.066*log(10.2))
MaxLen = 500
LenInc = 10
lbnd = seq(0,MaxLen-LenInc,LenInc)
midpt = lbnd + LenInc/2
ubnd = lbnd + LenInc
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
DiscMort = 0 # proportion of fish that die due to natural mortality
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 378
vbK = 0.25
GrowthParams = c(Linf, vbK)
CVSizeAtAge = 0.025
RefnceAges = NA

# ************************
# check growth assumptions
# ************************
SampleSize=1000
FishMort=0.1
SelParams = c(100,20) 
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

age.length <- X_plot_length_at_age(SimRes, point.size = 1,MaxAge, MaxLen)
age.length.plot <- print(age.length)

age.length.plot <- age.length.plot + geom_text(aes(label="c)", y=400, x=0.5))

setwd(fig_dir)
ggsave(age.length.plot, filename="age_length_WKW.png", height = a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )


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
InitL50 = 120
InitDelta = 30
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$convergence
FittedRes$ParamEst
#           Estimate lw_95%CL up_95%CL
# FMort        0.164    0.079    0.309
# L50_sel    188.956  172.650  206.803
# Delta_sel   76.737   65.221   90.287

# # True values
# FishMort = 0.2
# RetenParams = c(200, 100) # L50, L95-L50 for retention

# so, model kind of gets F but with wide confidence intervals
# selectivity of line not that badly estimated
# see difference (for F estimate) if retention is set to be lower, to simulate data

# plot results
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=400, xint=50,
                                   ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)
PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=400, xint=50,
                                      ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)

# *******************************
# fit LBCC to simulated BRUV data
# *******************************

ObsRetCatchFreqAtLen  = lendat_bruv3
DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
L50_Ret = FittedRes$ParamEst[2,1]
L95_Ret =  FittedRes$ParamEst[3,1] + L50_Ret
RetenAtLen = 1 / (1+exp(-log(19)*(midpt-L50_Ret)/(L95_Ret-L50_Ret)))
plot(midpt,RetenAtLen)

# specify initial params
InitFishMort = 0.25 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 120
InitDelta = 30
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$convergence
FittedRes$ParamEst
#            Estimate lw_95%CL up_95%CL
# FMort        0.016    0.002    0.099
# L50_sel    104.698  100.849  108.694
# Delta_sel   23.101   17.945   29.737

# # True values
# FishMort = 0.2
# RetenParams = c(100, 20) # L50, L95-L50 for retention

# similarly, model not able to capture F well, due to very late selectivity
# selectivity of BRUVs quite well estimated
# see difference (for F estimate) if retention is set to be lower, to simulate data

# plot results
# note, not really retained catch here (not relabelled) is selectivity of BRUVs
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=400, xint=50,
                                   ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)
PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=400, xint=50,
                                      ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)


# *****************************
# fit LBCC - fishery line data
# *****************************
setwd(data_dir)
# get data
dat <- read.csv("Length data for Charlotte UPDATED.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "King Wrasse") %>% 
  filter(Region %in% "Jurien")

# which(is.na(dat$TL))
range(dat$TL)
table(dat$TL)
FishLengths = dat$TL
lbns = trunc(FishLengths/LenInc)*LenInc
ObsRetCatchFreqAtLen  = as.vector(table(factor(lbns, levels=lbnd)))
plot(midpt, ObsRetCatchFreqAtLen, "l") # selectivity affected by small hooks

DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
RetenAtLen = NA # proportion of fish released, vector including mean and sd (option probably now obsolete)

# specify initial params
InitFishMort = 0.25 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 200
InitDelta = 100
TimeStep = 1/12
CVSizeAtAge = 0.025
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$convergence
FittedRes$ParamEst
# > FittedRes$ParamEst
# Estimate lw_95%CL up_95%CL
# FMort        0.820    0.538    0.947
# L50_sel    269.433  265.567  273.356
# Delta_sel  108.567  104.764  112.508

# Again can't estimtae F because of late selectivity but can estimate selectivity params

# plot results
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=400, xint=50,
                                   ymax=0.40, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)

PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=400, xint=50,
                                      ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)

plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                                 SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                                 RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                                 xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                                 ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200, Data="line")

plot.label = deparse1(bquote(.(plotting$Fest)))

good.plot <- ggplot()+
  geom_polygon(aes(x=plotting$x, y=plotting$y), colour=NA, fill=line_colour, alpha=0.5)+
  geom_point(aes(y=plotting$plotting.points$Prop, x=plotting$plotting.points$midpoint,
                 shape=plotting$plotting.points$Obs.Est, group=plotting$plotting.points$Obs.Est), colour="grey40", size=1)+
  scale_shape_manual(values = c(1,19), guide="none")+
  xlab("Length (mm)")+
  ylab("Proportion in length class")+
  geom_segment(aes(x = 377 , y = 0, xend = 377, yend = 0.15), linetype="dashed", colour="grey40")+
  geom_segment(aes(x = Linf , y = 0, xend = Linf, yend = 0.15), linetype="twodash", colour="grey40")+
  ylim(0,0.2)+
  #ggplot2::annotate("text", x=75, y=0.3, parse=T, label=as.character(plot.label))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_C-auricularis_line.png", height = a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )


# *****************************
# fit LBCC - BRUV data
# *****************************

L50_Ret = FittedRes$ParamEst[2,1] # Here saying that "retention" is the fish that get caught by the fishing gear so where mortality would be applied
L95_Ret =  FittedRes$ParamEst[3,1] + L50_Ret
RetenAtLen = 1 / (1+exp(-log(19)*(midpt-L50_Ret)/(L95_Ret-L50_Ret)))
# plot(midpt,RetenAtLen)

# get data
dat2 <- readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) %>%
  # dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Abrolhos_stereo-BRUVs"))) %>%
  # filter(length<=MaxLen) %>% 
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009")) 

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
InitDelta = 20
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$convergence
FittedRes$ParamEst



# > FittedRes$ParamEst
# Estimate lw_95%CL up_95%CL
# FMort        0.095    0.042    0.203
# L50_sel     80.493   76.967   84.182 # selectivity of BRUVs
# Delta_sel   32.249   27.647   37.616

# very poor fit - model expects few fish around 100-150, and more fish around 200-250

# plot results
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=400, xint=50,
                                   ymax=0.40, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)

PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=400, xint=50,
                                      ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)

plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                                 SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                                 RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                                 xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                                 ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200, Data="line")

plot.label = deparse1(bquote(.(plotting$Fest)))

good.plot <- ggplot()+
  geom_polygon(aes(x=plotting$x, y=plotting$y), colour=NA, fill=line_colour, alpha=0.5)+
  geom_point(aes(y=plotting$plotting.points$Prop, x=plotting$plotting.points$midpoint,
                 shape=plotting$plotting.points$Obs.Est, group=plotting$plotting.points$Obs.Est), colour="grey40", size=1)+
  scale_shape_manual(values = c(1,19), guide="none")+
  xlab("Length (mm)")+
  ylab("Proportion in length class")+
  geom_segment(aes(x = 377 , y = 0, xend = 377, yend = 0.15), linetype="dashed", colour="grey40")+
  geom_segment(aes(x = Linf , y = 0, xend = Linf, yend = 0.15), linetype="twodash", colour="grey40")+
  ylim(0,0.2)+
  #ggplot2::annotate("text", x=75, y=0.3, parse=T, label=as.character(plot.label))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_C-auricularis_line.png", height = a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )


### Test the observed distributions ####
setwd(data_dir)
# get data
dat <- read.csv("Length data for Charlotte UPDATED.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "King Wrasse") %>% 
  filter(Region %in% "Jurien")

FishLengths = dat$TL
lbns = trunc(FishLengths/LenInc)*LenInc
ObsRetCatchFreqAtLen  = as.vector(table(factor(lbns, levels=lbnd)))

dat2 <- readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) %>%
  # dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Abrolhos_stereo-BRUVs"))) %>%
  # filter(length<=MaxLen) %>% 
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009")) 

# which(is.na(dat2$length))
FishLengths = dat2$length
lbns = trunc(FishLengths/LenInc)*LenInc

ObsRetCatchFreqAtLen2  = as.vector(table(factor(lbns, levels=lbnd)))

Test <- exactRankTests::wilcox.exact(ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLen2)


