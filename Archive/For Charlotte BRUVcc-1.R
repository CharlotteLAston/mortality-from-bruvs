# For Charlotte
# Alex Hesp 27/2/2025

# Revised analysis - fitting model to BRUV data, whilst accounting for fishing mortality and retention curve

rm(list=ls())
library(L3Assess)

# Charlotte, will need to update package!!!
# install.packages("C:/~/L3Assess_0.1.0.tar.gz", source = TRUE, repos=NULL)

# if installing from github
# library(devtools)
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)

# simulate length data for stock
set.seed(123)
MaxAge = 40
MaxObsAge = 40 # max model age
TimeStep = 1/12 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort =  0.12 #exp(1.551-1.066*log(11.8))
MaxLen = 1300
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
Linf = 1136
vbK = 0.12
CVSizeAtAge = 0.02 # *** given short time step, reduced a bit
GrowthParams = c(Linf, vbK)
RefnceAges = NA
FishMort = 0.67*NatMort #0.67*NatMort #1.5*NatMort
PropReleased = NA 
SelParams = c(200, 50)
Ret_95 = Linf*0.5
Ret_Delta = 50
RetenParams = c(Ret_95-Ret_Delta, Ret_Delta) 
MLL=NA
SampleSize = 1000 # sample size is for retained fish

# **************************
# generate line fishing data
# **************************

Res_line=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType, SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
lendat_line1 = Res_line$ObsRetCatchFreqAtLen
lendat_line2 = Res_line$ObsDiscCatchFreqAtLen
plot(midpt, lendat_line1, "l")
lines(midpt, lendat_line2)
lines(midpt, lendat_line1 + lendat_line2, col="blue")

# check growth
par(mfrow=c(2,2))
plot(Res_line$ObsDecAgeRetCatch, Res_line$ObsRandLenRetCatch, xlim=c(0,40), ylim=c(0,1400))
points(Res_line$ObsDecAgeDiscCatch, Res_line$ObsRandLenDiscCatch, col="grey")

# **************************
# Fit to simulated line data
# **************************
ObsRetCatchFreqAtLen  = lendat_line1
DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
RetenAtLen = NA

# Set selectivity (retention paramaters) for line fishing
InitL50 = 500
InitDelta = 25
InitFishMort = 0.25 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes_line=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes_line$ParamEst
# Estimate lw_95%CL up_95%CL
# FMort        0.080    0.062    0.102
# L50_sel    517.229  505.751  528.967
# Delta_sel   64.468   50.384   82.489
# PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
#                                    SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
#                                    RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
#                                    xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=500,
#                                    ymax=0.15, yint=0.05, PlotCLs=TRUE, FittedRes_line, nReps=200)

# basic plot that could be used in chapter (to compare with one for BRUV data fit)
plot(midpt, ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen), "l")
lines(midpt, FittedRes_line$ModelDiag$ExpRetCatchPropInLenClass, col="blue")

# true retention params
# RetenParams
# [1] 518  50
# True Fmort
# FishMort = 0.67*NatMort
# [1] 0.0804


# ******************
# generate BRUV data
# ******************
set.seed(123)
Res_bruv=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType, SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

lendat_bruv1 = Res_bruv$ObsRetCatchFreqAtLen # Retained by line fishing
lendat_bruv2 = Res_bruv$ObsDiscCatchFreqAtLen # Fish that are selected but then discarded, what we would see on the stereo-BRUVs
lendat_bruv3 = lendat_bruv1 + lendat_bruv2 


# check growth
plot(Res_line$ObsDecAgeRetCatch, Res_line$ObsRandLenRetCatch, xlim=c(0,40), ylim=c(0,1400))
points(Res_line$ObsDecAgeDiscCatch, Res_line$ObsRandLenDiscCatch, col="grey")

L50_Ret = FittedRes_line$ParamEst[2,1] # Here saying that "retention" is the fish that get caught by the fishing gear so where mortality would be applied
L95_Ret =  FittedRes_line$ParamEst[3,1] + L50_Ret
RetenAtLen = 1 / (1+exp(-log(19)*(midpt-L50_Ret)/(L95_Ret-L50_Ret)))

ObsRetCatchFreqAtLen = lendat_bruv1 + lendat_bruv2
ObsDiscCatchFreqAtLen = NA

par(mfrow=c(2,2))
# overall length distn
plot(midpt, ObsRetCatchFreqAtLen, "l") # total

# fish just selected by BRUVs, vs fish of size that could be retained
plot(midpt, lendat_bruv1, "l")
lines(midpt, lendat_bruv2)
lines(midpt, lendat_bruv3, col="blue")

InitL50 = 200 # selectivity of the stereo-BRUVs
InitDelta = 25 # selectivity of the stereo-BRUVs
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)     

FittedRes$ParamEst
# Estimate lw_95%CL up_95%CL
# FMort        0.081    0.069    0.095
# L50_sel    198.212  190.531  206.203
# Delta_sel   47.163   35.702   62.303

# True Fmort
# FishMort = 0.67*NatMort
# [1] 0.0804

# True Selectivity
# SelParams = c(200, 50) 

# plot fit of model to data (which includes fish of all sizes measured by BRUVs)
par(mfrow=c(2,2))
plot(midpt, lendat_bruv3/sum(lendat_bruv3), "l")
lines(midpt, FittedRes$ModelDiag$ExpTotCatchPropInLenClass, col="blue")


