
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
MaxAge = 30 # max model age
TimeStep = 1
NatMort = exp(1.551-1.066*log(20))
FishMort = 0.2
MaxLen = 560
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(150, 50) # L50, L95-L50 for gear selectivity
RetenParams = c(300, 50) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality
# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 561
vbK = 0.14
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
plot(Res$midpt, lendat_line, "o", ylim=c(0,250))

# generate BRUV data (with different realisation, i.e. starting with different random number)
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# compare length-frequency data - BRUV
lendat_bruv = Res$ObsRetCatchFreqAtLen # Retained by line fishing
lendat_bruv2 = Res$ObsDiscCatchFreqAtLen # Discarded fish - difference between retained and selectivity?
lendat_bruv3 = lendat_bruv + lendat_bruv2
lines(Res$midpt, lendat_bruv3, "o", col="red") # Gives us those extra fish before they are retained by the fishery


# *********************************************************
# specify biology and fishery parameters for  Breaksea Cod
# *********************************************************
MaxAge = 30
MaxObsAge = 20
TimeStep = 1/4 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = exp(1.551-1.066*log(20))
MaxLen = 600
LenInc = 20
lbnd = seq(0,MaxLen-LenInc,LenInc)
midpt = lbnd + LenInc/2
ubnd = lbnd + LenInc
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
DiscMort = 0 # proportion of fish that die due to natural mortality
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = c(561,558)
vbK = c(0.14,0.17)
GrowthParams = data.frame(Linf=Linf, vbK=vbK)
CVSizeAtAge = c(0.02,0.02)
RefnceAges = NA

# ************************
# check growth assumptions
# ************************
SampleSize=500
FishMort=0.1
SelParams = c(390,80) 
RetenParams=NA
SimRes=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                                  SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# plot data
PlotOpt=0 # 0=all plots, 1=retained lengths at age and growth curves, 2=retained plus discard lengths at age
# and growth curves, 3=retained plus discard length frequency, 4=retained plus discard age frequency,
# 5=sex specific retained lengths at age and growth curves, 6=retained plus discard lengths at age and growth curves,
# 7=female length frequency, 8=male length frequency, 9=female age frequency, 10=male age frequency,
# 11=selectivity/retention, 12=F-at-age reten + disc, 13=F-at-age reten, 14=F-at-age disc
# par(mfcol=c(1,1))
PlotSimLenAndAgeFreqData_EqMod(MaxAge, MaxLen, SimRes, PlotOpt=2)


age.length <- X_plot_length_at_age(SimRes, point.size = 1,MaxAge, MaxLen)
age.length.plot <- print(age.length)

age.length.plot <- age.length.plot + geom_text(aes(label="b)", y=605, x=1))

setwd(fig_dir)
ggsave(age.length.plot, filename="age_length_BC.png", height = a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )


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
#             Estimate lw_95%CL up_95%CL
# FMort        0.400    0.353    0.449
# L50_sel    297.378  283.103  312.374
# Delta_sel    8.188    1.160   57.790

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
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=650, xint=50,
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
#             Estimate lw_95%CL up_95%CL
# FMort        0.321    0.253    0.399
# L50_sel    243.428  220.148  269.171
# Delta_sel  142.565  122.856  165.436

# # True values
# FishMort = 0.2 wasn't bad, did L50 well but not Delta 
# RetenParams = c(299, 50) # L50, L95-L50 for retention

# similarly, model not able to capture F well, due to very late selectivity
# selectivity of BRUVs not well estiamted 
# see difference (for F estimate) if retention is set to be lower, to simulate data

# plot results
# note, not really retained catch here (not relabelled) is selectivity of BRUVs
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=650, xint=50,
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
  filter(Species %in% "Breaksea Cod") %>% 
  filter(Region %in% "Metro")

# which(is.na(dat$TL))
LenInc=20
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
InitFishMort = 0.2 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 300
InitDelta = 50

Timestep = 1/4
CVSizeAtAge = c(0.02,0.02)
MLL=NA

params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes_line=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes_line$convergence
FittedRes_line$ParamEst
# > FittedRes$ParamEst
#            Estimate lw_95%CL up_95%CL
# FMort        0.551    0.150    0.895
# L50_sel    387.308  347.055  432.230
# Delta_sel   78.858   58.959  105.474

# Useless estimate and couldn't estimate the selectivity very well 

# plot results
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=600, xint=50,
                                   ymax=0.20, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)

PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=650, xint=50,
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
  geom_segment(aes(x = 469 , y = 0, xend = 469, yend = 0.25), linetype="dashed", colour="grey40")+
  geom_segment(aes(x = Linf[1] , y = 0, xend = Linf[1], yend = 0.25), linetype="dotted", colour="grey40")+
  ylim(0,0.3)+
  geom_text(aes(label="c)", x=50, y=0.3))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_E-armatus_line.png", height = a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )



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
MLL=NA
# Timestep = 1/6
# CVSizeAtAge = 0.03

params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes_bruv=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes_bruv$convergence
FittedRes_bruv$ParamEst



# > FittedRes$ParamEst
#            Estimate lw_95%CL up_95%CL
# FMort        0.636    0.404    0.819
# L50_sel    167.870  148.806  189.376
# Delta_sel   48.511   27.630   85.172

# Gives us a pretty reasonable estimate of F 

# plot results
PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                   SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                   RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                   xaxis_lab=NA, yaxis_lab=NA, xmax=650, xint=50,
                                   ymax=0.30, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)

PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=50,
                                      ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)

plot(midpt, lendat_bruv3/sum(lendat_bruv3), type="l", ylim=c(0,0.2))
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
  geom_segment(aes(x = 469 , y = 0, xend = 469, yend = 0.25), linetype="dashed", colour="grey40")+
  geom_segment(aes(x = Linf[1] , y = 0, xend = Linf[1], yend = 0.25), linetype="dotted", colour="grey40")+
  ylim(0,0.3)+
  geom_text(aes(label="d)", x=50, y=0.3))+
  #ggplot2::annotate("text", x=100, y=0.25, parse=T, label=as.character(plot.label))+
  theme_classic()
good.plot

setwd(fig_dir)
ggsave(good.plot, filename="Catch-Curve_E-armatus_bruv.png", height = a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )

# Plot mortality at length
mortality.at.length <-X_PlotLengthBasedCatchCurve_Mortality(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                                            PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                                            MaxAge, NatMort, TimeStep, xmax=1500, xint=50, ymax=1, yint=0.02, 
                                                            FittedRes_1=FittedRes_bruv, FittedRes_2=FittedRes_line, colour_1=bruv_colour, colour_2=line_colour) 
setwd(fig_dir)
ggsave(mortality.at.length, filename="Mortality-at-length_E-armatus.png", height = a4.width*0.5,width = a4.width, units  ="mm", dpi = 300 )

