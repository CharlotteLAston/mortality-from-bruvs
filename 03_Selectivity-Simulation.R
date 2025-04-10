# ******************************************************
# Run a simulation of BRUV data to see how well we can
# retrieve the parameters 
# ******************************************************

# library(devtools)
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)

library(L3Assess)
library(tidyverse)
library(dplyr)
library(ggplot2)

rm(list=ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_data", sep="/")
sim_dir <- paste(working.dir, "Simulations", sep="/")


# *********************************************************
# specify biology and fishery parameters 
# *********************************************************
MaxAge = 10.2
MaxObsAge = 20 # max model age
TimeStep = 1/4 # model timestep (1/12 for Pink Snapper, 1/4 for Breaksea Cod, 1/4 for SMW)
NatMort =  exp(1.551-1.066*log(10.2)) #0.12 
MaxLen = 400
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
Linf = 378
vbK = 0.25
CVSizeAtAge = 0.03
GrowthParams = c(Linf, vbK)
RefnceAges = NA
FishMort = 1.5*NatMort #0.67*NatMort #1.5*NatMort
PropReleased = NA


#RetScen = 2 # 1 = low retention scenario, 2 = High retention scenario

NSIM = 200
SampleSize = 1000

param_res <- array(0, dim=c(NSIM, 10)) #Trial, Sel, F_line, L50_line, Delta_line, F_bruv, L50_bruv, Delta_bruv, Convergence - line, convergence - BRUV
Ret_scen_res <- NULL

Start <- Sys.time()
set.seed(1313)
for(RetScen in 1:3){
 
   param_res <- array(0, dim=c(NSIM, 10))
  
  for(SIM in 1:NSIM) {
    
    # **************************
    # Specify sel scenario
    # **************************
    
    if(RetScen == 1){
      SelParams = c(100, 50) # Doesn't really matter because retention is the thing that the fishing mortality applies to, we're assuming no PRM
      Ret_95 = Linf*0.5
      Ret_Delta = 50
      RetenParams = c(Ret_95-Ret_Delta, Ret_Delta)
    } else if (RetScen == 2){
      SelParams = c(100, 50) # Doesn't really matter because retention is the thing that the fishing mortality applies to, we're assuming no PRM
      Ret_95 = Linf*0.75
      Ret_Delta = 50
      RetenParams = c(Ret_95-Ret_Delta, Ret_Delta) # Issue is that even with early retention they've already stopped growing - so you won't be able to estimate F pretty much ever
    } else if (RetScen==3){
      SelParams = c(100, 50) # Doesn't really matter because retention is the thing that the fishing mortality applies to, we're assuming no PRM
      Ret_95 = Linf*0.9
      Ret_Delta = 50
      RetenParams = c(Ret_95-Ret_Delta, Ret_Delta) # Issue is that even with early retention they've already stopped growing - so you won't be able to estimate F pretty much ever
    } 
    
    # **************************
    # Simulate line data
    # **************************
    
    MLL=NA
    SampleSize = 1000 # sample size is for retained fish
    
    Res_line=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType, SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
    lendat_line1 = Res_line$ObsRetCatchFreqAtLen
    lendat_line2 = Res_line$ObsDiscCatchFreqAtLen
    
    # **************************
    # Simulate bruv data
    # **************************
    Res_bruv=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType, SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
    
    lendat_bruv1 = Res_bruv$ObsRetCatchFreqAtLen # Retained by line fishing
    lendat_bruv2 = Res_bruv$ObsDiscCatchFreqAtLen # Fish that are selected but then discarded, what we would see on the stereo-BRUVs
    
    # **************************
    # Fit to simulated line data
    # **************************

    # set random initial parameter values
    if(FishMort<0.3){
      InitFishMort = 0.2
    } else {
      InitFishMort = 0.6
    }
    InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
    
    if(RetScen==1){
      InitL50 = rnorm(1, 130, 5)
    } else if (RetScen==2){
      InitL50 = rnorm(1, 230, 5)
    } else if (RetScen==3){
      InitL50 = rnorm(1, 300, 5)
    } 
  
    InitDelta = rnorm(1, 50, 5) 
    
    params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
    
    ObsRetCatchFreqAtLen  = lendat_line1
    DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
    ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
    RetenAtLen = NA
    
    FittedRes_line=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
    # FittedRes_line$ParamEst
    
    # ******************
    # Fit to BRUV data
    # ******************
    
    L50_Ret = FittedRes_line$ParamEst[2,1] # Here saying that "retention" is the fish that get caught by the fishing gear so where mortality would be applied
    L95_Ret =  FittedRes_line$ParamEst[3,1] + L50_Ret
    RetenAtLen = 1 / (1+exp(-log(19)*(midpt-L50_Ret)/(L95_Ret-L50_Ret)))
    
    ObsRetCatchFreqAtLen = lendat_bruv1 + lendat_bruv2
    ObsDiscCatchFreqAtLen = NA
    
    InitL50 = 100 # selectivity of the stereo-BRUVs
    InitDelta = 50 # selectivity of the stereo-BRUVs
    params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
    
    FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                              lbnd, ubnd, midpt, SelectivityAtLen, RetenAtLen, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)     

    
    if (!is.na(FittedRes$convergence)) {
      param_res[SIM,1] = SIM
      param_res[SIM,2] = RetScen
      param_res[SIM,3] = FittedRes_line$ParamEst[1,1]
      param_res[SIM,4] = FittedRes_line$ParamEst[2,1]
      param_res[SIM,5] = FittedRes_line$ParamEst[3,1]
      param_res[SIM,6] = FittedRes$ParamEst[1,1]
      param_res[SIM,7] = FittedRes$ParamEst[2,1]
      param_res[SIM,8] = FittedRes$ParamEst[3,1]
      param_res[SIM,9] = FittedRes_line$convergence 
      param_res[SIM,10] = FittedRes$convergence 
    }
    
  }
   param_res <- param_res %>% 
     as.data.frame() %>% 
     rename(Trial = "V1",
            Sel_Scen = "V2",
            F_line = "V3",
            L50_line = "V4",
            Delta_line = "V5",
            F_bruv = "V6",
            L50_bruv = "V7",
            Delta_bruv = "V8",
            convergence_line = "V9",
            convergence_bruv = "V10")
  
  Ret_scen_res <- rbind(Ret_scen_res, param_res)
}
End <- Sys.time()
Time.taken <- End - Start

Full_res_low <- Ret_scen_res %>%
  mutate(F_scen = "Low")

Full_res_high <- Ret_scen_res %>%
  mutate(F_scen = "high")

Full_res <- rbind(Full_res_low, Full_res_high)

setwd(sim_dir)
saveRDS(Full_res, file="sim_sel_WKW.rds")





