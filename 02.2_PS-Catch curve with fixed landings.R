
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
MaxAge = 40 # max model age
TimeStep = 1/4 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 0.12
FishMort = 0.2
MaxLen = 1300
LenInc = 30
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(250, 50) # L50, L95-L50 for gear selectivity
RetenParams = c(400, 50) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality
# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 1136
vbK = 0.12
CVSizeAtAge = 0.05
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
plot(Res$midpt, lendat_line, "o", ylim=c(0,150))

# generate BRUV data (with different realisation, i.e. starting with different random number)
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# compare length-frequency data - BRUV
lendat_bruv = Res$ObsRetCatchFreqAtLen # Retained by line fishing
lendat_bruv2 = Res$ObsDiscCatchFreqAtLen # Discarded fish - difference between retained and selectivity?
lendat_bruv3 = lendat_bruv + lendat_bruv2
lines(Res$midpt, lendat_bruv3, "o", col="red") # Gives us those extra fish before they are retained by the fishery


# *********************************************************
# specify biology and fishery parameters for PS Snapper
# *********************************************************
MaxAge = 40
MaxObsAge = 35
TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 0.12
MaxLen = 1300
LenInc = 50
lbnd = seq(0,MaxLen-LenInc,LenInc)
midpt = lbnd + LenInc/2
ubnd = lbnd + LenInc
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
DiscMort = 0 # proportion of fish that die due to natural mortality
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 1136
vbK = 0.12
GrowthParams = c(Linf, vbK)
CVSizeAtAge = 0.05
RefnceAges = NA

# ************************
# check growth assumptions
# ************************
SampleSize=1000
FishMort=0.1
SelParams = c(250,50) 
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

age.length.plot <- age.length.plot + geom_text(aes(label="a)", y=1255, x=2))

setwd(fig_dir)
ggsave(age.length.plot, filename="age_length_PS.png", height = a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )

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
InitL50 = 400
InitDelta = 10
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$convergence
FittedRes$ParamEst
#             Estimate lw_95%CL      up_95%CL
# FMort        0.201    0.179  2.230000e-01
# L50_sel    392.370    0.000  4.145368e+08
# Delta_sel    2.214    0.000 9.987730e+198

# # True values
# FishMort = 0.2
# RetenParams = c(400, 100) # L50, L95-L50 for retention

# so, model captures F well because selectivity is early
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
MLL = NA

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
#             Estimate lw_95%CL up_95%CL
# FMort        0.207    0.182    0.233
# L50_sel    303.665  285.905  322.529
# Delta_sel  110.526   88.275  138.384

# # True values
# FishMort = 0.2 wasn't bad
# RetenParams = c(250, 50) # L50, L95-L50 for retention

# similarly, model not able to capture F well, due to very late selectivity
# selectivity of BRUVs quite well estimated
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
dat <- read.csv("Length data for Charlotte.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "Pink Snapper")

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
MLL=NA
# specify initial params
InitFishMort = 0.2 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 400
InitDelta = 50
CVSizeAtAge = 0.05
TimeStep = 1
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes_line=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes_line$convergence
FittedRes_line$ParamEst

# > FittedRes$ParamEst
#             Estimate lw_95%CL up_95%CL
# FMort        0.182    0.161    0.205
# L50_sel    441.782  432.256  451.519
# Delta_sel   45.078   31.497   64.515

# Pretty reasonable estimate, might be a little on the lower side

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


plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                                 SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                                 RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                                 xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                                ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200, Data="line")

plot.label = deparse1(bquote(.(plotting$Fest)))

line_colour <- "#BBCC33"
bruv_colour <- "#66CCEE"

good.plot <- ggplot()+
  geom_polygon(aes(x=plotting$x, y=plotting$y), colour=NA, fill=line_colour, alpha=0.5)+
  geom_point(aes(y=plotting$plotting.points$Prop, x=plotting$plotting.points$midpoint,
                 shape=plotting$plotting.points$Obs.Est, group=plotting$plotting.points$Obs.Est), colour="grey40", size=1)+
  scale_shape_manual(values = c(1,19), guide="none")+
  xlab(NULL)+
  ylab(NULL)+
  geom_segment(aes(x = 510 , y = 0, xend = 510, yend = 0.2), linetype="dashed", colour="grey40")+
  geom_segment(aes(x = 1136 , y = 0, xend = 1136, yend = 0.2), linetype="dotted", colour="grey40")+
  geom_text(aes(label="a)", x=100, y=0.25))+
  ylim(0,0.25)+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_C-auratus_line.png", height = a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )

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
dat2 <- readRDS("australian-synthesis_complete_length_chrysophrys_auratus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Marmion|TwoRocks|Two.Rocks|Capes|Geographe|SwC|south-west|Ngari|Warnbro|Rottnest"))) %>% 
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
InitL50 = 200 # Selectivity of the BRUVs
InitDelta = 50
TimeStep=1/2
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes_bruv=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes_bruv$convergence
FittedRes_bruv$ParamEst


# > FittedRes$ParamEst
#           Estimate lw_95%CL up_95%CL
# FMort        0.246    0.221    0.273
# L50_sel    256.869  248.708  265.298 # Selectivity of BRUVs
# Delta_sel   69.846   56.375   86.536

# Actually a pretty good fit!

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

plotting <- X_PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                                 SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                                 RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                                 xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                                 ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200, Data="bruv")

plot.label = deparse1(bquote(.(plotting$Fest)))

line_colour <- "#BBCC33"
bruv_colour <- "#66CCEE"

good.plot <- ggplot()+
  geom_polygon(aes(x=plotting$x, y=plotting$y), colour=NA, fill=bruv_colour, alpha=0.5)+
  geom_point(aes(y=plotting$plotting.points$Prop, x=plotting$plotting.points$midpoint,
                 shape=plotting$plotting.points$Obs.Est, group=plotting$plotting.points$Obs.Est), colour="grey40", size=1)+
  scale_shape_manual(values = c(1,19), guide="none")+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,0.25)+
  geom_segment(aes(x = 480 , y = 0, xend = 480, yend = 0.2), linetype="dashed", colour="grey40")+
  geom_segment(aes(x = 1136 , y = 0, xend = 1136, yend = 0.2), linetype="dotted", colour="grey40")+
  geom_text(aes(label="b)", x=100, y=0.25))+
  #ggplot2::annotate("text", x=200, y=0.25, parse=T, label=as.character(plot.label))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_C-auratus_bruv.png", height = a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )

# Plot mortality at length

mortality.at.length <-X_PlotLengthBasedCatchCurve_Mortality(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                                            PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                                            MaxAge, NatMort, TimeStep, xmax=1500, xint=100, ymax=1, yint=0.02, 
                                                            FittedRes_1=FittedRes_bruv, FittedRes_2=FittedRes_line, colour_1=bruv_colour, colour_2=line_colour) 
setwd(fig_dir)
ggsave(mortality.at.length, filename="Mortality-at-length_C-auratus.png", height = a4.width*0.5 ,width = a4.width, units  ="mm", dpi = 300 )
