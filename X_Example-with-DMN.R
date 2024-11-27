# if installing from github
# library(devtools)
# devtools::install_github("SAlexHesp/WAFishBiologyRPackage", build_vignettes=TRUE, force=TRUE)
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)

rm(list=ls())
library(L3Assess)

# Fit length-based catch curve to length-frequency data generated from Dirchlet multinomial distribution.
# First, use SimLenAndAgeFreqData_EqMod function to calculate expected proportions at length
set.seed(123)
SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
MaxAge = 30
TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/MaxAge
FishMort = 0.2
MaxLen = 1200
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
RetenParams = c(500, 50) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality
# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 800
vbK = 0.2
CVSizeAtAge = 0.08
GrowthParams = c(Linf, vbK)
RefnceAges = NA
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# Now simulate length-frequency data from Dirchlet multinomial distribution.
set.seed(123)
nSampEvents = 50
nFishPerSampEvent = 20
theta_val = 0.3
midpt = Res$midpt
ExpPropAtLen = Res$ModelDiag$ExpRetCatchPropAtLen_Fem
# plot(midpt, Res$ModelDiag$ExpRetCatchPropAtLen_Fem, "l")
# lines(midpt, Res$ModelDiag$ExpRetCatchPropAtLen_Mal, col="blue") # check that same for both sexes (as single sex model)
res=SimLenFreqDat_DirMultDistn(nSampEvents, nFishPerSampEvent, theta_val, midpt, ExpPropAtLen)
plot(res$midpt, res$simLenFreq, "o")

# Fit length-based catch curve, assuming Dirichlet multinomial distribution
ObsRetCatchFreqAtLen = res$simLenFreq
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obsolete)
midpt=Res$midpt
lbnd=Res$lbnd
ubnd=Res$ubnd
InitFishMort = 0.25 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 400
InitDelta = 100
DistnType=2
InitTheta = 0.3 # specify starting parameters
InitTheta_logit = log(InitTheta/(1-InitTheta)) # logit transform
params = c(InitFishMort_logit, log(InitL50), log(InitDelta),InitTheta_logit)
FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
FittedRes$ParamEst
FittedRes$RetCatch_SampleSize
FittedRes$RetCatch_EffSampleSize


# Simulate data
SampleSize=1000
set.seed(123)
MaxAge = 30
TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/MaxAge
FishMort = 0.2
MaxLen = 1500
LenInc = 50
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
RetenParams = c(NA, NA) # L50, L95-L50 for retention
DiscMort = 0
# # von Bertalanffy
# GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
# Linf = 800
# vbK = 0.2
# CVSizeAtAge = 0.08
# GrowthParams = c(Linf, vbK)
# RefnceAges = NA
# Schnute
GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
t1 = 1 # growth - Schnute
t2 = 10 # growth - Schnute
y1 = 400 # growth - Schnute
y2 = 1000 # growth - Schnute
a = 0.1 # growth - Schnute
b = 2.0 # growth - Schnute
CVSizeAtAge = 0.08
GrowthParams = c(y1, y2, a, b)
RefnceAges = c(t1,t2)
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
ObsRetCatchFreqAtLen = as.vector(Res$ObsRetCatchFreqAtLen)
MinFreq = 20 # set minimum frequency for larger lengths for analysis
# note, this needs to be high enough so that data for ln(n/dt) vs relative age for essentially straight
# line - if not, Z will be biased!!!
midpt=Res$midpt
lbnd=Res$lbnd
ubnd=Res$ubnd
ModelType = 2 # 1 = von Bertalanffy growth curve (Pauly), 2 = length-converted catch curve - Schnute growth curve)
res=GetLenConvCatchCurveResults(ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                MinFreq, lbnd, midpt, ubnd, MaxAge)
PlotLenConvCatchCurveResults(MaxAge, ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen, MinFreq,
                             lbnd, midpt, ubnd)



# Simulate data
SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
set.seed(123)
MaxAge = 30
TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/MaxAge
FishMort = 0.2
MaxLen = 1200
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
RetenParams = c(NA, NA) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality
# single sex, von Bertalanffy
DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 800
vbK = 0.2
CVSizeAtAge = 0.08
GrowthParams = c(Linf, vbK)
RefnceAges = NA
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
# 2 sexes, von Bertalanffy
DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = c(700,850)
vbK = c(0.3,0.2)
CVSizeAtAge = c(0.08,0.08)
GrowthParams = data.frame(Linf=Linf, vbK=vbK)
RefnceAges = NA
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
# 1 sex, Schnute
DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
t1 = 0.5 # growth - Schnute
t2 = 25 # growth - Schnute
y1 = 100 # growth - Schnute
y2 = 1000 # growth - Schnute
a = 0.02 # growth - Schnute
b = 3.0 # growth - Schnute
GrowthParams = c(y1, y2, a, b)
RefnceAges = c(t1,t2)
CVSizeAtAge = 0.05
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
# 2 sexes, Schnute
DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
t1 = c(0.5,0.5) # growth - Schnute
t2 = c(25,25) # growth - Schnute
y1 = c(100,100) # growth - Schnute
y2 = c(1000,1000) # growth - Schnute
a = c(0.02,0.02) # growth - Schnute
b = c(3,3) # growth - Schnute
CVSizeAtAge = c(0.05, 0.05)
GrowthParams = data.frame(y1=y1, y2=y2, a=a, b=b)
RefnceAges = data.frame(t1=t1,t2=t2)
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obsolete)
midpt=Res$midpt
lbnd=Res$lbnd
ubnd=Res$ubnd
InitFishMort = 0.25 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 400
InitDelta = 100
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))

FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityAtLen, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)

PlotLengthBasedCatchCurve_Mortality(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                    PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                    MaxAge, NatMort, TimeStep, xmax=NA, xint=NA, ymax=NA, yint=NA, FittedRes)

