
rm(list=ls())
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)
library(L3Assess)
library(dplyr)
library(stringr)
library(ggplot2)

working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")

a4.width <- 160

#### SAMPLE SIZE FOR SPANGLED EMPEROR ####

# inputs that don't change
MaxAge = 30
TimeStep = 1 
NatMort = 4.22/MaxAge
FishMort = 0.2
MaxLen = 800
LenInc = 20
MLL=NA 
SelectivityType = 2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityVec = NA # selectivity vector
SelParams = c(388, 185)  # bruv curve

RetenParams = c(NA, NA) 
DiscMort = 0 
GrowthCurveType = 1
Linf = 573.1
vbK = 0.282
CVSizeAtAge = 0.05
GrowthParams = c(Linf, vbK)
RefnceAges = NA
ObsDiscCatchFreqAtLen = NA 
PropReleased = NA 
InitFishMort = 0.05
InitFishMort_logit = log(InitFishMort/(1-InitFishMort))
DistnType = 1

InitL50 = 400
InitDelta = 100

params = c(InitFishMort_logit, log(InitL50), log(InitDelta))

samplesize <-  c(100,250,500,750,1000,1500,2000,3000,4000,5000,6000)

results <- as.data.frame(array(0, dim=c(length(samplesize), 3)))

for(i in 1:length(samplesize)) {
  SampleSize = samplesize[i]
  Res_sim=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
  
  ObsRetCatchFreqAtLen = Res_sim$ObsRetCatchFreqAtLen
  midpt=Res_sim$midpt
  lbnd=Res_sim$lbnd
  ubnd=Res_sim$ubnd
  
  Results_sim=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                              lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  Est <- Results_sim$ParamEst[1,1]
  uprbnd <- Results_sim$ParamEst[1,3]
  lwrbnd <- Results_sim$ParamEst[1,2]
  Results <- data.frame(Est,uprbnd,lwrbnd)
  
  results[i, ] <- Results
  
}

results$samplesize <- samplesize

results.spango <- results %>% 
  rename(
    Fmort = V1,
    upperbnd = V2,
    lowerbnd = V3) %>% 
  mutate(mortality=Fmort) %>% 
  mutate(condition = 1-(mortality/(mortality+NatMort))) %>% 
  mutate(upr_condition = 1-(lowerbnd/(lowerbnd+NatMort))) %>% 
  mutate(lwr_condition = 1-(upperbnd/(upperbnd+NatMort))) 

# results1

#* Plot spangled emperor sample sizes ####

Spango_Sample_Size <- results.spango %>% 
  ggplot() +
  geom_rect(xmin=-275, xmax=6350, ymin=0.37, ymax=0.47, linetype="dotted", fill=NA, colour="grey50")+
  geom_point(aes(x=samplesize, y=condition), colour="grey20")+
  geom_errorbar(aes(x=samplesize, y=condition ,ymin=lwr_condition, ymax=upr_condition), colour="grey20")+
  theme_classic()+
  ylab("Condition Estimate")+
  xlab("Sample size")+
  geom_vline(xintercept = 404, linetype="dashed", colour="#88CBED", linewidth=1.25)
Spango_Sample_Size 

setwd(fig_dir)
ggsave(Spango_Sample_Size, filename="Sample-Size_L-nebulosus.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )



#### SAMPLE SIZE FOR Breaksea Cod ####

# Simulate data
MaxAge = 20
TimeStep = 0.5 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 0.3
FishMort = 0.2
MaxLen = 500
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityVec = NA # selectivity vector
SelParams = c(251, 93) # L50, L95-L50 for gear selectivity
RetenParams = c(NA, NA) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 561
# vbK = 0.241
vbK = 0.14
CVSizeAtAge = 0.025
GrowthParams = c(Linf, vbK)
RefnceAges = NA


# fit catch curve to simulated data
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
length(ObsRetCatchFreqAtLen)
length(midpt)
InitFishMort = 0.05 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 251
InitDelta = 93
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
DistnType = 1

samplesize <-  c(100,250,500,750,1000,1500,2000,3000,4000,5000,6000)

results <- as.data.frame(array(0, dim=c(length(samplesize), 3)))

for(i in 1:length(samplesize)) {
  SampleSize = samplesize[i]
  Res_sim=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
  
  ObsRetCatchFreqAtLen = Res_sim$ObsRetCatchFreqAtLen
  midpt=Res_sim$midpt
  lbnd=Res_sim$lbnd
  ubnd=Res_sim$ubnd
  
  Results_sim=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                              lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  Est <- Results_sim$ParamEst[1,1]
  uprbnd <- Results_sim$ParamEst[1,3]
  lwrbnd <- Results_sim$ParamEst[1,2]
  Results <- data.frame(Est,uprbnd,lwrbnd)
  
  results[i, ] <- Results
  
}

results$samplesize <- samplesize

results.breaksea <- results %>% 
  rename(
    Fmort = V1,
    upperbnd = V2,
    lowerbnd = V3) %>% 
  mutate(mortality=Fmort) %>% 
  mutate(condition = 1-(mortality/(mortality+NatMort))) %>% 
  mutate(upr_condition = 1-(lowerbnd/(lowerbnd+NatMort))) %>% 
  mutate(lwr_condition = 1-(upperbnd/(upperbnd+NatMort))) 


# results1

#* Plot breaksea sample sizes ####

Breaksea_Sample_Size <- results.breaksea %>% 
  ggplot() +
  geom_rect(xmin=-275, xmax=6350, ymin=0.54, ymax=0.64, linetype="dotted", fill=NA, colour="grey50")+
  geom_point(aes(x=samplesize, y=condition), colour="grey20")+
  geom_errorbar(aes(x=samplesize, y=condition ,ymin=lwr_condition, ymax=upr_condition), colour="grey20")+
  theme_classic()+
  ylab("Condition Estimate")+
  xlab("Sample size")+
  geom_vline(xintercept = 252, linetype="dashed", colour="#A9439A", linewidth=1)
Breaksea_Sample_Size 

setwd(fig_dir)
ggsave(Breaksea_Sample_Size, filename="Sample-Size_E-armatus.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )


#### SAMPLE SIZE FOR Western King Wrasse ####

# Simulate data
MaxAge = 10.2 #https://onlinelibrary.wiley.com/doi/10.1111/j.1095-8649.2012.03446.x
TimeStep = 0.5 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/MaxAge
MaxLen = 400
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityVec = NA # selectivity vector
RetenParams = c(NA, NA) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
# Linf = 664
Linf = 484
# vbK = 0.241
vbK = 0.16
CVSizeAtAge = 0.05
GrowthParams = c(Linf, vbK)
RefnceAges = NA

#* Abrolhos ####
# fit catch curve to simulated data
FishMort = 0.2 #1.5 * NatMort
SelParams = c(85, 53) # L50, L95-L50 for gear selectivity
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obsolete)
length(ObsRetCatchFreqAtLen)
length(midpt)
InitFishMort = 0.05 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 85
InitDelta = 53
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
DistnType = 1

samplesize <-  c(100,250,500,750,1000,1500,2000,3000,4000,5000,6000)

results <- as.data.frame(array(0, dim=c(length(samplesize), 3)))

for(i in 1:length(samplesize)) {
  SampleSize = samplesize[i]
  Res_sim=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
  
  ObsRetCatchFreqAtLen = Res_sim$ObsRetCatchFreqAtLen
  midpt=Res_sim$midpt
  lbnd=Res_sim$lbnd
  ubnd=Res_sim$ubnd
  
  Results_sim=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                              lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  Est <- Results_sim$ParamEst[1,1]
  uprbnd <- Results_sim$ParamEst[1,3]
  lwrbnd <- Results_sim$ParamEst[1,2]
  Results <- data.frame(Est,uprbnd,lwrbnd)
  
  results[i, ] <- Results
  
}

results$samplesize <- samplesize

results.WKW <- results %>% 
  rename(
    Fmort = V1,
    upperbnd = V2,
    lowerbnd = V3) %>% 
  mutate(mortality=Fmort) %>% 
  mutate(condition = 1-(mortality/(mortality+NatMort))) %>% 
  mutate(upr_condition = 1-(lowerbnd/(lowerbnd+NatMort))) %>% 
  mutate(lwr_condition = 1-(upperbnd/(upperbnd+NatMort))) 

# results1

#* Plot WKW sample sizes ####

WKW_Sample_Size_Abrolhos  <- results.WKW %>% 
  ggplot() +
  geom_rect(xmin=-400, xmax=6350, ymin=0.62, ymax=0.72, linetype="dotted", fill=NA, colour="grey50")+
  geom_point(aes(x=samplesize, y=condition), colour="grey20")+
  geom_errorbar(aes(x=samplesize, y=condition ,ymin=lwr_condition, ymax=upr_condition), colour="grey20")+
  theme_classic()+
  ylab("Condition Estimate")+
  xlab("Sample size")+
  geom_vline(xintercept = 479, linetype="dashed", colour="#117633", linewidth=1)
WKW_Sample_Size_Abrolhos 

setwd(fig_dir)
ggsave(WKW_Sample_Size_Abrolhos , filename="Sample-Size_C-auricularis_Abrolhos.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

# #* Metro ####

CVSizeAtAge = 0.05

# fit catch curve to simulated data
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obsolete)
length(ObsRetCatchFreqAtLen)
length(midpt)
FishMort = 0.2
SelParams = c(74, 35) # L50, L95-L50 for gear selectivity
InitFishMort = 0.05 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 76
InitDelta = 36
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
DistnType = 1

samplesize <-  c(100,250,500,750,1000,1500,2000,3000,4000,5000,6000,7000,8000)

results <- as.data.frame(array(0, dim=c(length(samplesize), 3)))

for(i in 1:length(samplesize)) {
  SampleSize = samplesize[i]
  Res_sim=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
  
  ObsRetCatchFreqAtLen = Res_sim$ObsRetCatchFreqAtLen
  midpt=Res_sim$midpt
  lbnd=Res_sim$lbnd
  ubnd=Res_sim$ubnd
  
  Results_sim=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                              lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  Est <- Results_sim$ParamEst[1,1]
  uprbnd <- Results_sim$ParamEst[1,3]
  lwrbnd <- Results_sim$ParamEst[1,2]
  Results <- data.frame(Est,uprbnd,lwrbnd)
  
  results[i, ] <- Results
  
}

results$samplesize <- samplesize

results.WKW <- results %>% 
  rename(
    Fmort = V1,
    upperbnd = V2,
    lowerbnd = V3) %>% 
  mutate(mortality=Fmort) %>% 
  mutate(condition = 1-(mortality/(mortality+NatMort))) %>% 
  mutate(upr_condition = 1-(lowerbnd/(lowerbnd+NatMort))) %>% 
  mutate(lwr_condition = 1-(upperbnd/(upperbnd+NatMort))) 

WKW_Sample_Size_Metro <- results.WKW %>% 
  ggplot() +
  geom_rect(xmin=-400, xmax=8350, ymin=0.165, ymax=0.265, linetype="dotted", fill=NA, colour="grey50")+
  geom_point(aes(x=samplesize, y=Fmort))+
  geom_errorbar(aes(x=samplesize, y=Fmort ,ymin=lowerbnd, ymax=upperbnd))+
  theme_classic()+
  ylab("Condition Estimate")+
  xlab("Sample size")+
  geom_vline(xintercept = 7852, colour="#43A999", linetype="dashed", linewidth=1)
WKW_Sample_Size_Metro

setwd(fig_dir)
ggsave(WKW_Sample_Size_Metro , filename="Sample-Size_C-auricularis_Metro.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )


#### SAMPLE SIZE FOR PINK SNAPPER ####

# Simulate data
MaxAge = 41
TimeStep = 0.5 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/MaxAge
MaxLen = 1130 
LenInc = 30
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityVec = NA # selectivity vector
DiscMort = 0 # proportion of fish that die due to natural mortality
FishMort = 0.2
SelParams = c(230, 52) # L50, L95-L50 for gear selectivity

CVSizeAtAge = 0.075

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 1136 #https://academic.oup.com/icesjms/article/74/1/180/2669555
vbK = 0.12 #https://academic.oup.com/icesjms/article/74/1/180/2669555
GrowthParams = c(Linf, vbK)
RefnceAges = NA


# fit catch curve to simulated data
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
length(ObsRetCatchFreqAtLen)
length(midpt)
InitFishMort = 0.05 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 216
InitDelta = 20
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
DistnType = 1

samplesize <-  c(100,250,500,750,1000,1500,2000,3000,4000,5000,6000)


results <- as.data.frame(array(0, dim=c(length(samplesize), 3)))

for(i in 1:length(samplesize)) {
  SampleSize = samplesize[i]
  Res_sim=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
  
  ObsRetCatchFreqAtLen = Res_sim$ObsRetCatchFreqAtLen
  midpt=Res_sim$midpt
  lbnd=Res_sim$lbnd
  ubnd=Res_sim$ubnd
  
  Results_sim=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                              lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  Est <- Results_sim$ParamEst[1,1]
  uprbnd <- Results_sim$ParamEst[1,3]
  lwrbnd <- Results_sim$ParamEst[1,2]
  Results <- data.frame(Est,uprbnd,lwrbnd)
  
  results[i, ] <- Results
  
}

results$samplesize <- samplesize

results.snapper <- results %>% 
  rename(
    Fmort = V1,
    upperbnd = V2,
    lowerbnd = V3) %>% 
  mutate(mortality=Fmort) %>% 
  mutate(condition = 1-(mortality/(mortality+NatMort))) %>% 
  mutate(upr_condition = 1-(lowerbnd/(lowerbnd+NatMort))) %>% 
  mutate(lwr_condition = 1-(upperbnd/(upperbnd+NatMort))) 

#* Plot pink snapper sample sizes ####

PinkSnapper_Sample_Size  <- results.snapper %>% 
  ggplot() +
  geom_rect(xmin=-275, xmax=6350, ymin=0.285, ymax=0.385, linetype="dotted", fill=NA, colour="grey50")+
  geom_point(aes(x=samplesize, y=condition), colour="grey20")+
  geom_errorbar(aes(x=samplesize, y=condition ,ymin=lwr_condition, ymax=upr_condition), colour="grey20")+
  theme_classic()+
  ylab("Condition Estimate")+
  xlab("Sample size")+
  geom_vline(xintercept = 872, linetype="dashed", colour="#332387", linewidth=1)
PinkSnapper_Sample_Size 

setwd(fig_dir)
ggsave(PinkSnapper_Sample_Size, filename="Sample-Size_C-auratus.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )


#### SAMPLE SIZE FOR SOUTHERN MAORI WRASSE ####

# Simulate data
MaxAge = 11.8
TimeStep = 0.5 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/MaxAge
MaxLen = 400 
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified

SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityVec = NA # selectivity vector
DiscMort = 0 # proportion of fish that die due to natural mortality

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 328 # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1095-8649.2012.03446.x
vbK = 0.28 # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1095-8649.2012.03446.x
GrowthParams = c(Linf, vbK)
RefnceAges = NA

#* Perth Metro ####
CVSizeAtAge = 0.05
FishMort = 0.2
SelParams = c(206, 122) # L50, L95-L50 for gear selectivity

# fit catch curve to simulated data
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
length(ObsRetCatchFreqAtLen)
length(midpt)
InitFishMort = 0.05 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 190
InitDelta = 110
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
DistnType = 1

samplesize <-  c(100,250,500,750,1000,1500,2000,3000,4000,5000,6000)

results <- as.data.frame(array(0, dim=c(length(samplesize), 3)))

for(i in 1:length(samplesize)) {
  SampleSize = samplesize[i]
  Res_sim=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
  
  ObsRetCatchFreqAtLen = Res_sim$ObsRetCatchFreqAtLen
  midpt=Res_sim$midpt
  lbnd=Res_sim$lbnd
  ubnd=Res_sim$ubnd
  
  Results_sim=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                              lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  Est <- Results_sim$ParamEst[1,1]
  uprbnd <- Results_sim$ParamEst[1,3]
  lwrbnd <- Results_sim$ParamEst[1,2]
  Results <- data.frame(Est,uprbnd,lwrbnd)
  
  results[i, ] <- Results
  
}

results$samplesize <- samplesize

results.maoriwrasse <- results %>% 
  rename(
    Fmort = V1,
    upperbnd = V2,
    lowerbnd = V3)%>% 
  mutate(mortality=Fmort) %>% 
  mutate(condition = 1-(mortality/(mortality+NatMort))) %>% 
  mutate(upr_condition = 1-(lowerbnd/(lowerbnd+NatMort))) %>% 
  mutate(lwr_condition = 1-(upperbnd/(upperbnd+NatMort))) 


#* Plot maori wrasse sample sizes ####

MaoriWrasse_Sample_Size  <- results.maoriwrasse %>% 
  ggplot() +
  geom_rect(xmin=-275, xmax=6350, ymin=0.6, ymax=0.7, linetype="dotted", fill=NA, colour="grey50")+
  geom_point(aes(x=samplesize, y=condition), colour="grey20")+
  geom_errorbar(aes(x=samplesize, y=condition ,ymin=lwr_condition, ymax=upr_condition), colour="grey20")+
  theme_classic()+
  ylab("Condition Estimate")+
  xlab("Sample size")+
  geom_vline(xintercept = 459, linetype="dashed", colour="#872155", linewidth=1)
MaoriWrasse_Sample_Size 

setwd(fig_dir)
ggsave(MaoriWrasse_Sample_Size, filename="Sample-Size_O-lineolatus_Metro.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#* Perth Capes ####
CVSizeAtAge = 0.05
TimeStep = 1/12
SelParams = c(109, 51) 
FishMort=0.2

# fit catch curve to simulated data
ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
length(ObsRetCatchFreqAtLen)
length(midpt)
InitFishMort = 0.05 # specify starting parameters
InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
InitL50 = 109
InitDelta = 51
params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
DistnType = 1

samplesize <-  c(100,250,500,750,1000,1500,2000,3000,4000,5000,6000)


results <- as.data.frame(array(0, dim=c(length(samplesize), 3)))

for(i in 1:length(samplesize)) {
  SampleSize = samplesize[i]
  Res_sim=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
  
  ObsRetCatchFreqAtLen = Res_sim$ObsRetCatchFreqAtLen
  midpt=Res_sim$midpt
  lbnd=Res_sim$lbnd
  ubnd=Res_sim$ubnd
  
  Results_sim=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                              lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  Est <- Results_sim$ParamEst[1,1]
  uprbnd <- Results_sim$ParamEst[1,3]
  lwrbnd <- Results_sim$ParamEst[1,2]
  Results <- data.frame(Est,uprbnd,lwrbnd)
  
  results[i, ] <- Results
  
}

results$samplesize <- samplesize

results.maoriwrasse <- results %>% 
  rename(
    Fmort = V1,
    upperbnd = V2,
    lowerbnd = V3)%>% 
  mutate(mortality=Fmort) %>% 
  mutate(condition = 1-(mortality/(mortality+NatMort))) %>% 
  mutate(upr_condition = 1-(lowerbnd/(lowerbnd+NatMort))) %>% 
  mutate(lwr_condition = 1-(upperbnd/(upperbnd+NatMort))) 


#* Plot maori wrasse sample sizes ####

MaoriWrasse_Sample_Size  <- results.maoriwrasse %>% 
  ggplot() +
  geom_rect(xmin=-275, xmax=6350, ymin=0.63, ymax=0.73, linetype="dotted", fill=NA, colour="grey50")+
  geom_point(aes(x=samplesize, y=condition), colour="grey20")+
  geom_errorbar(aes(x=samplesize, y=condition ,ymin=lwr_condition, ymax=upr_condition), colour="grey20")+
  theme_classic()+
  ylab("Condition Estimate")+
  xlab("Sample size")+
  geom_vline(xintercept = 5464, linetype="dashed", colour="#CB6778", linewidth=1)
MaoriWrasse_Sample_Size 

setwd(fig_dir)
ggsave(MaoriWrasse_Sample_Size, filename="Sample-Size_O-lineolatus_Capes.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )





