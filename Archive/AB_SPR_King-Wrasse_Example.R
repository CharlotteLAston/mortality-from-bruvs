## King Wrasse ##

# Example 2: hermaphroditic species
InitRecruit <- 1 # Initial recruitment
MaxModelAge <- 20 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
Linf <- c(333, 399) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.29, 0.14) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(0.05, -3.76) # years - von Bertalanffy growth model parameters - Females, males
EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
ln_lenwt_a <- -11.82 # for log-log relationship
lenwt_b <- 3.07 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
WLrel_Type <- 2 # 1=power, 2=log-log relationship
EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 1 # Ratio of females to males at age zero
FinalSex_A50 <- 4.7 # Logistic sex change relationship parameters (inflection point)
FinalSex_A95 <- 6.9 # Logistic sex change relationship parameters (95% of max probability)
EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
mat_A50 <- c(3.8, 3.8) # females, males - Logistic length (mm) at maturity relationship parameters
mat_A95 <- c(6.7, 6.7) # females, males - Logistic length (mm) at maturity relationship parameters
EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
Gear_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
Gear_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
FemGearSelAtAge=rep(1,MaxModelAge+1)
MalGearSelAtAge = FemGearSelAtAge
EstGearSelAtAge <- data.frame(FemGearSelAtAge=FemGearSelAtAge, MalGearSelAtAge=MalGearSelAtAge) # gear selectivity at age (from age 0), inputted as values in data frame

Land_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
Land_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
FemSelLandAtAge=rep(1,MaxModelAge+1)
FemSelLandAtAge[1:6]=0
MalSelLandAtAge = FemSelLandAtAge
EstLandSelAtAge <- data.frame(FemSelLandAtAge=FemSelLandAtAge, MalSelLandAtAge=MalSelLandAtAge) # gear selectivity at age (from age 0), inputted as values in data frame


ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
# FemRetProbAtAge=rep(1,MaxModelAge+1)
# FemRetProbAtAge[1:6]=0
# MalRetProbAtAge = FemRetProbAtAge
EstRetenAtAge <- data.frame(FemRetProbAtAge=NA, MalRetProbAtAge=NA) # fish retention at age (from age 0), inputted as values in data frame
DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
NatMort = 0.3 # natural mortality  (year-1)
FMort <- 0.73 # estimate of fishing mortality, e.g. from catch curve analysis
Res=CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                             ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                             EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
                             EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
                             ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)

Res$Fem_SPR
Res$Mal_SPR
Res$Eq_CombSexRelSpBiom

Current_F <- FMort # estimate of fishing mortality, e.g. from catch curve analysis
PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel_land at age, 4=mal mat/sel_land at age,
# 5=fem sel_land/sel_disc/prob reten/gear sel, 6=fem sel_land/sel_disc/prob reten/gear sel, 7=fem F at age,
# 8=mal F at age, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=fem biom at age,
# 13=catch at age, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
# 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
# 22=plot male depletion vs eq. recruitment (plots 1-22 plotted if !is.na(EggFertParam))
PlotPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                         lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                         ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam,
                         mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge,
                         Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95,
                         EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, Current_F)

#### Logistic catch curve ####
set.seed(123)
MinAge = 1
MaxAge = 40
Ages = MinAge:MaxAge
NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
FMort = 0.1
ZMort = FMort + NatMort
SelA50 = 6
SelA95 = 8
SampleSize = 1000 # required sample size. For 2 sex model, same sample size generated for each sex.
Res=SimAgeFreqData_EqMod(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
ObsAgeFreq = unlist(as.vector(Res$CatchSample)) # input as vector for single sex model
# Specify catch curve model and required inputs for that model
NatMort = 0.104 # Logistic selectivity
MinFreq = NA # # Logistic selectivity
Init_FMort = 0.2
Init_SelA50 = 5
Init_SelDelta = 2
params = log(c(Init_FMort, Init_SelA50, Init_SelDelta))
res=GetLogisticCatchCurveResults(params, NatMort, Ages, ObsAgeFreq)
res$ParamEst

CatchCurveModel = 3 # Logistic selectivity
RecAssump = NA # Logistic selectivity
NatMort = 0.104 # Logistic selectivity
MinFreq = NA # # Logistic selectivity
Init_FMort = 0.2
Init_SelA50 = 5
Init_SelDelta = 2
params = log(c(Init_FMort, Init_SelA50, Init_SelDelta))
res=GetLogisticCatchCurveResults(params, NatMort, Ages, ObsAgeFreq)
res$ParamEst
PlotAgeBasedCatchCurveResults_NormalSpace(RecAssump, SpecRecAge=NA, MinFreq, MinAge, MaxAge, NatMort,
                                          ObsAgeFreq, CatchCurveModel, MainLabel=NA,
                                          xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
                                          ymax=NA, yint=NA, PlotCLs=T)
