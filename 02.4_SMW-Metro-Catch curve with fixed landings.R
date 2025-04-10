
# Alex Hesp 29 Jan 2025
# For Charlotte
# Exploration of an approach for deriving estimates of F, retention parameters and selectivity
# parameters from line and BRUV data
# 
# install.packages("C:/~/L3Assess_0.1.0.tar.gz", source = TRUE, repos=NULL)


rm(list=ls())
library(L3Assess)
library(tidyverse)
library(dplyr)


#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_data", sep="/")

setwd(working.dir)
source("X_My-Functions.R")
a4.width=160
line_colour <- "#BBCC33"
bruv_colour <- "#66CCEE"

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
MaxAge = 15 # max model age
TimeStep = 1/4 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort =  exp(1.551-1.066*log(11.8))
FishMort = 0.2
MaxLen = 470
LenInc = 10
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(100, 20) # L50, L95-L50 for gear selectivity
RetenParams = c(200, 100) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality
# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 328
vbK = 0.28
CVSizeAtAge = 0.03
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
# specify biology and fishery parameters for WKW in Metro
# *********************************************************
MaxAge = 20
MaxObsAge = 11.8 # max model age
TimeStep = 1/4 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort =  exp(1.551-1.066*log(11.8))
MaxLen = 500
LenInc = 20
lbnd = seq(0,MaxLen-LenInc,LenInc)
midpt = lbnd + LenInc/2
ubnd = lbnd + LenInc
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
DiscMort = 0 # proportion of fish that die due to natural mortality
# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 328
vbK = 0.28
CVSizeAtAge = 0.03
GrowthParams = c(Linf, vbK)
RefnceAges = NA
DiscMort = 0

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

age.length.plot <- age.length.plot + geom_text(aes(label="d)", y=400, x=1))

setwd(fig_dir)
ggsave(age.length.plot, filename="age_length_SMW.png", height = a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )

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
#            Estimate lw_95%CL up_95%CL
# FMort        0.141    0.085    0.226
# L50_sel    178.468  165.489  192.465
# Delta_sel   71.597   60.965   84.083

# # True values
# FishMort = 0.2
# RetenParams = c(200, 100) # L50, L95-L50 for retention

# Actually did a pretty good job with this 
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
L50_Ret = FittedRes$ParamEst[2,1] # Specify retention using the simulated line fishing data
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
#             Estimate lw_95%CL up_95%CL
# FMort        0.053    0.032    0.086
# L50_sel    109.857  105.379  114.526
# Delta_sel   28.448   22.476   36.007

# # True values
# FishMort = 0.2
# SelectParams = c(100, 20) # L50, L95-L50 for retention

# Model not able to capture F well, not sure why
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
  filter(Species %in% "Maori Wrasse") %>% 
  filter(Region %in% "Metro")


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
InitDelta = 50

CVSizeAtAge = 0.03
TimeStep = 1/4

params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes_line=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes_line$convergence
FittedRes_line$ParamEst
# > FittedRes$ParamEst
#            Estimate lw_95%CL up_95%CL
# FMort        0.834    0.515    0.959
# L50_sel    267.213  265.366  269.073
# Delta_sel   60.787   59.043   62.584

# Again can't estimate F because there's not enough data at the end of the length distribution

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
  xlab(NULL)+
  ylab(NULL)+
  geom_segment(aes(x = 320 , y = 0, xend = 320, yend = 0.3), linetype="dashed", colour="grey40")+
  geom_segment(aes(x = Linf , y = 0, xend = Linf, yend = 0.3), linetype="dotted", colour="grey40")+
  ylim(0,0.35)+
  geom_text(aes(label="e)", x=40, y=0.35))+
  #ggplot2::annotate("text", x=75, y=0.3, parse=T, label=as.character(plot.label))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_O-lineolatus_line.png", height = a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )



# save selectivity params from line data analysis.
# these assumed to be equivalent to retention parameters,
# and will be used in next analysis to fix retention

plot(midpt, ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen), "l", ylim=c(0,0.2)) # line data
sum(ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen))

# *****************************
# fit LBCC - BRUV data
# *****************************

L50_Ret = FittedRes_line$ParamEst[2,1] # Here saying that "retention" is the fish that get caught by the fishing gear so where mortality would be applied
L95_Ret =  FittedRes_line$ParamEst[3,1] + L50_Ret
RetenAtLen = 1 / (1+exp(-log(19)*(midpt-L50_Ret)/(L95_Ret-L50_Ret)))
# plot(midpt,RetenAtLen)

# get data
setwd(data_dir)
dat2 <- readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Rottnest|TwoRocks|Marmion|Warnbro|Two.Rocks"))) %>% 
  dplyr::filter(str_detect(campaign,"2006|2007|2008|2009|2010|2011"))

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
FittedRes_bruv=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes_bruv$convergence
FittedRes_bruv$ParamEst




# > FittedRes$ParamEst
#           Estimate lw_95%CL up_95%CL
# FMort        0.164    0.082    0.302
# L50_sel    190.663  166.770  217.979
# Delta_sel  104.560   88.000  124.236


# plot results
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=400, xint=50,
                                   ymax=0.15, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)

PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=400, xint=50,
                                      ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)

plot(midpt, lendat_bruv3/sum(lendat_bruv3), type="l", ylim=c(0,0.1))
lines(midpt,ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen),col="blue") #Actual observed data

plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                                 SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                                 RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                                 xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                                 ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200, Data="bruv")

plot.label = deparse1(bquote(.(plotting$Fest)))

good.plot <- ggplot()+
  geom_polygon(aes(x=plotting$x, y=plotting$y), colour=NA, fill=bruv_colour, alpha=0.5)+
  geom_point(aes(y=plotting$plotting.points$Prop, x=plotting$plotting.points$midpoint,
                 shape=plotting$plotting.points$Obs.Est, group=plotting$plotting.points$Obs.Est), colour="grey40", size=1)+
  scale_shape_manual(values = c(1,19), guide="none")+
  xlab(NULL)+
  ylab(NULL)+
  geom_segment(aes(x = 320 , y = 0, xend = 320, yend = 0.3), linetype="dashed", colour="grey40")+
  geom_segment(aes(x = Linf , y = 0, xend = Linf, yend = 0.3), linetype="dotted", colour="grey40")+
  ylim(0,0.35)+
  geom_text(aes(label="f)", x=40, y=0.35))+
  #ggplot2::annotate("text", x=75, y=0.3, parse=T, label=as.character(plot.label))+
  theme_classic()+
  theme(axis.text = element_text(size=))
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_O-lineolatus_bruv.png", height = a4.width*0.5 ,width = a4.width*0.5, units  ="mm", dpi = 300 )

## Mortality at length plot
mortality.at.length <-X_PlotLengthBasedCatchCurve_Mortality(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                                            PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                                            MaxAge, NatMort, TimeStep, xmax=1500, xint=50, ymax=1, yint=0.02, 
                                                            FittedRes_1=FittedRes_bruv, FittedRes_2=FittedRes_line, colour_1=bruv_colour, colour_2=line_colour) 

setwd(fig_dir)
ggsave(mortality.at.length, filename="Mortality-at-length_O-lineolatus.png", height = a4.width*0.5, width = a4.width, units  ="mm", dpi = 300 )
