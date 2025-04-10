
# Alex Hesp 29 Jan 2025
# For Charlotte
# Exploration of an approach for deriving estimates of F, retention parameters and selectivity
# parameters from line and BRUV data


# rm(list=ls())
# library(L3Assess)
# 
# install.packages("/Users/22291142/Documents/L3Assess_0.1.0.tar.gz", source = TRUE, repos=NULL)
# # 
# # if installing from github
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
FishMort = 0.2
MaxLen = 400
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
CVSizeAtAge = 0.05
GrowthParams = c(Linf, vbK)
RefnceAges = NA
# generate line data
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# plot selectivity assumptions
plot(Res$midpt, Res$ModelDiag$SelAtLength, "l")
lines(Res$midpt, Res$ModelDiag$RetAtLength, col="blue")
lines(Res$midpt, Res$ModelDiag$SelLandAtLength, col="green")

# plot length-frequency data (retained fish - line)
lendat_line = Res$ObsRetCatchFreqAtLen
lendat_line2 = Res$ObsDiscCatchFreqAtLen
plot(Res$midpt, lendat_line, "o", ylim=c(0,100))

# generate BRUV data (with different realisation, i.e. starting with different random number)
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# compare length-frequency data - BRUV
lendat_bruv = Res$ObsRetCatchFreqAtLen
lendat_bruv2 = Res$ObsDiscCatchFreqAtLen
lendat_bruv3 = lendat_bruv + lendat_bruv2
lines(Res$midpt, lendat_bruv3, "o", col="red")


# **************************************
# specify biology and fishery parameters
# **************************************
MaxAge = 20
MaxObsAge = 10.2
TimeStep = 1/4 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/MaxObsAge
MaxLen = 400
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
CVSizeAtAge = 0.05
RefnceAges = NA

# ************************
# check growth assumptions
# ************************
SampleSize=1000
FishMort=0.1
SelParams = c(100,20); RetenParams=NA
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

# *******************************
# fit LBCC to simulated line data
# *******************************

ObsRetCatchFreqAtLen  = lendat_line
DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
RetenAtLen = NA

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
# Estimate lw_95%CL up_95%CL
# FMort        0.065    0.019    0.202
# L50_sel    171.292  157.529  186.257
# Delta_sel   66.432   55.112   80.077

# # True values
# FishMort = 0.2
# RetenParams = c(200, 100) # L50, L95-L50 for retention

# so, model not able to capture F well, due to very late selectivity
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
# Estimate lw_95%CL up_95%CL
# FMort        0.027    0.007    0.093
# L50_sel    104.773  101.092  108.587
# Delta_sel   23.222   18.117   29.766

# # True values
# FishMort = 0.2
# RetenParams = c(200, 100) # L50, L95-L50 for retention

# similarly, model not able to capture F well, due to very late selectivity
# selectivity of BRUVs quite well estimated
# see difference (for F estimate) if retention is set to be lower, to simulate data

# plot results
# note, not really retained catch here - not relabelled
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

# get data
dat <- read.csv("WKW_line_Jurien.csv", header=TRUE)
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


# save selectivity params from line data analysis.
# these assumed to be equivalent to retention parameters,
# and will be used in next analysis to fix retention

plot(midpt, ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen), "l", ylim=c(0,0.2)) # line data
sum(ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen))

# *****************************
# fit LBCC - BRUV data
# *****************************

L50_Ret = FittedRes$ParamEst[2,1]
L95_Ret =  FittedRes$ParamEst[3,1] + L50_Ret
RetenAtLen = 1 / (1+exp(-log(19)*(midpt-L50_Ret)/(L95_Ret-L50_Ret)))
# plot(midpt,RetenAtLen)

# get data
dat2 <- read.csv("WKW_BRUV_Jurien.csv", header=TRUE)
# which(is.na(dat2$length))
head(dat2)
range(dat2$length)
table(dat2$length)
FishLengths = dat2$length
lbns = trunc(FishLengths/LenInc)*LenInc

ObsRetCatchFreqAtLen  = as.vector(table(factor(lbns, levels=lbnd)))

# compare line vs BRUV length data
lines(midpt,ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen), col="blue") # bruv data
sum(ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen))


# specify initial params
InitFishMort = 0.25 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 100
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
                                   ymax=0.10, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)

PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                      RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                      MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=400, xint=50,
                                      ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)




