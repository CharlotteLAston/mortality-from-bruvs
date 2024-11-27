###################################################

# Setting up initial population to establish an 
# equilibrium level of mortality and recruitment

###################################################
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(stringr)
library(forcats)
library(RColorBrewer)
library(geosphere)
library(abind)

rm(list = ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
pop_dir <- paste(working.dir, "Population_Files", sep="/")

#### PARAMETERS FOR THE FISH POPULATION #####

## Set timestep 
step = 1/12 # We're doing a monthly timestep here

# Beverton-Holt Parameters
h = 0.76 # Steepness 0.76
R0 = 1 # Initial recruitment

# Proportion female
prop.fem=0.5

## L. Nebulosus ##
# Von bertalanffy
Linf = 664
k = 0.241
t0 = -0.375

# Length-weight relationship
WLa = 0.000028
WLb = 2.8761

# Natural mortality
M = 0.146

# Selectivity
A50 = 3
A95 = 4

# Longevity
MaxAge = 30

# Maturity
M50 = 3.62
M95 = 5.97

## L. sebae ##
# Von bertalanffy
Linf = 623.9
k = 0.2
t0 = -1.533

# Length-weight relationship
WLa = 1.406e-08
WLb = 3.051

# Natural mortality
M = 0.104

# Selectivity
A50 = 3.762
A95 = 6.241

# Longevity
MaxAge = 32

# Maturity
M50 = 3.43
M95 = 4.37

## L. miniatus ##
# Von bertalanffy
Linf = 45
k = 0.166
t0 = -1.75

# Length-weight relationship
WLa = 8.04e-06 
WLb = 3.051

# Natural mortality
M = 0.23

# Selectivity
A50 = 2.5
A95 = 4.5 

# Longevity
MaxAge = 20

# Maturity
M50 = 1.3
M95 = 1.85

#### SET UP LIFE HISTORY VALUES FOR L. NEBULOSUS ####

Life.History <- as.data.frame(array(0, dim=c((MaxAge*12)+1, 3))) %>% 
  rename(Age = "V1") %>% 
  rename(Length = "V2") %>% 
  rename(Weight = "V3")

for (r in 1:(MaxAge*12)){                                    # This sets the ages to be every month, but in decimals
  Life.History[r+1, 1] <- Life.History[r, 1] + step 
}

Life.History <- Life.History %>% 
  mutate(Length = Linf*(1-exp(-k*(Age-t0)))) %>% 
  mutate(Weight = (WLa*(Length^WLb)/1000)) # Divide by 1000 to get in kg 

#### SET UP SURVIVAL AND MATURITY FOR UNFISHED POPULATION ####

Unfished.Pop.SetUp <- as.data.frame(array(1, dim=c((MaxAge*12)+1, 3))) %>% 
  rename(Unfish.Surv = "V1") %>% 
  rename(Unfish.Mat = "V2") %>% 
  rename(Unfish.Bio = "V3")

## Unfished Population Survival in Each Age Group
Unfished.Pop.SetUp[1, 1] <- R0

for (r in 2:((MaxAge*12)+1)){                                    # This calculates the survival in the next age based on the previous age
  Unfished.Pop.SetUp[r, 1] <- Unfished.Pop.SetUp[r-1, 1]*exp(-M/12) 
}

## Unfished Population Proportion Mature per Recruit in Each Age Group

Unfished.Pop.SetUp <- Unfished.Pop.SetUp %>% 
  mutate(Age = Life.History$Age) %>% 
  mutate(Unfish.Mat = 1/(1+(exp(-log(19)*((Age-M50)/(M95-M50))))))

## Unfished Mature Biomass per Recruit in Each Age Group

Unfished.Pop.SetUp <- Unfished.Pop.SetUp %>% 
  mutate(Unfish.Bio = Unfish.Mat * Unfish.Surv * Life.History$Weight)

## Calculate the Total Mature Unfished Biomass per Recruit

Unfish.Fem.SSB <- Unfished.Pop.SetUp %>%  # This is the total mature female spawning biomass per recruit for the unfished population
  slice(which(row_number() %% 12 == 1)) %>% 
  summarise(., sum=sum(Unfish.Bio))

#### CALCULATE BEVERTON-HOLT PARAMETERS FOR UNFISHED POPULATION ####

alpha <- (Unfish.Fem.SSB/R0)*((1-h)/(4*h))

beta <- ((h-0.2)/(0.8*h*R0))

#### SET UP SURVIVAL AND MATURITY FOR FISHED POPULATION ####

# Fishing parameters
eq.init.fish =  0.2*M

Fished.Pop.SetUp <- as.data.frame(array(1, dim=c((MaxAge*12)+1, 5))) %>% 
  rename(Selectivity = "V1") %>% 
  rename(Fishing.Mort = "V2") %>% 
  rename(TotMort = "V3") %>% 
  rename(Fish.Surv = "V4") %>% 
  rename(Fish.Mat = "V5")

## Selectivity for each age group - This is for the equilibrium population

Fished.Pop.SetUp <- Fished.Pop.SetUp %>% 
  mutate(Age = Life.History$Age) %>% 
  mutate(Length = Life.History$Length) %>% 
  mutate(Selectivity = 1/(1+(exp(-log(19)*((Age-A50)/(A95-A50))))))

## Next calculate fishing mortality 

Fished.Pop.SetUp <- Fished.Pop.SetUp %>% 
  mutate(Fishing.Mort = Selectivity * eq.init.fish)

## Calculate total mortality

Fished.Pop.SetUp <- Fished.Pop.SetUp %>% 
  mutate(TotMort = Fishing.Mort + M)

## Calculate Fished Survival
Fished.Pop.SetUp[1, 4] <- R0

for (r in 2:((MaxAge*12)+1)){                                    # This calculates the survival in the next age based on the previous age
  Fished.Pop.SetUp[r, 4] <- Fished.Pop.SetUp[r-1, 4]*exp(-Fished.Pop.SetUp[r-1,3]*step) # Divide by the time step here
}

## Calculate Mature Fished Biomass per Recruit

Fished.Pop.SetUp <- Fished.Pop.SetUp %>% 
  mutate(Fish.Mat = Unfished.Pop.SetUp$Unfish.Mat) %>% 
  mutate(Fish.Bio = Fish.Mat * Fish.Surv * Life.History$Weight) %>% 
  mutate(Fish.Died = 1-Fishing.Mort)

## Calculate the Total Mature Fished Biomass per Recruit

Fish.Fem.SSB <- Fished.Pop.SetUp %>%  # This is the total mature female spawning biomass per recruit for the unfished population
  slice(which(row_number() %% 12 == 1)) %>% 
  summarise(., sum=sum(Fish.Bio))

##### CALCULATING SPAWNER PER RECRUIT AND EQUILIRBIRUM RECRUITMENT ####

SPR <- Fish.Fem.SSB/Unfish.Fem.SSB

Equil.Rec <- (Fish.Fem.SSB-alpha)/(beta*Fish.Fem.SSB)

#### CREATING AN INITIAL EQUILIBRIUM FISHED POPULATION FOR THE MODEL ####

## Our initial population uses all the numbers from the previous step of this code 
## The maturity at age, weight at age, fished per recruit survival
## the proportion of females at age
## We also need to decide on an initial level of recruitment (in thousands)

init.recruit <- 1000 # in thousands - normally 5000 for big model, 5 for small model

## Calculate initial fished recruitment

Init.Fish.Rec <- (Fish.Fem.SSB-alpha)/(Fish.Fem.SSB*beta)*init.recruit 

## Calculate initial unfished spawning biomass

Init.Unfish.Fem.SpawnBio <- Init.Fish.Rec * Unfish.Fem.SSB

## Calculate initial fished spawning biomass

Init.Fem.SpawnBio <- Init.Fish.Rec * Fish.Fem.SSB

## Calculate new Beverton-Holt Parameters

alpha <- (Init.Unfish.Fem.SpawnBio/Init.Fish.Rec)*((1-h)/(4*h))

beta <- ((h-0.2)/(0.8*h*Init.Fish.Rec))

## Calculate number of female recruits in the next time step

N.Female.Rec <- (Init.Fem.SpawnBio/(alpha+beta*Init.Fem.SpawnBio))*prop.fem

## Calculate survival into the next time step to create a full age structured population

Starting.Pop.For.Model <- as.data.frame(array(1, dim=c((MaxAge*12)+1, 1)))

Starting.Pop.For.Model <- Starting.Pop.For.Model %>% 
  mutate(Age = Life.History$Age) 

Starting.Pop.For.Model[1,1] <- N.Female.Rec

Starting.Pop.For.Model <- Starting.Pop.For.Model %>% 
  mutate(Caught = 0)

for (r in 2:((MaxAge*12)+1)){                                    # This calculates the survival in the next age based on the previous age using total mortality from our fished population
  Starting.Pop.For.Model[r, 1] <- Starting.Pop.For.Model[r-1, 1]*exp(-Fished.Pop.SetUp[r-1,3]*step) # Divide by the time step here
  
  Starting.Pop.For.Model[r, 3] <- Starting.Pop.For.Model[r-1, 1]*(1-(exp(-Fished.Pop.SetUp[r-1,2]*step))) # Divide by the time step here
}

Starting.Pop.For.Model <- Starting.Pop.For.Model %>% 
  rename(Number = "V1") %>% 
  mutate(Length = Life.History$Length)

fished <- Starting.Pop.For.Model %>% 
  mutate(Pop.Sample = Number * Fished.Pop.SetUp$Selectivity) %>% 
  group_by(Age = ceiling(row_number()/12)) %>%
  summarise(Total.Sample = sum(Pop.Sample)) %>% 
  group_by(Age) %>% 
  summarise(Total = sum(Total.Sample)) %>% 
  mutate(Proportion = Total/sum(Total)) %>% 
  mutate(Percentage = Proportion *100)

histograms <- ggplot()+
  geom_col(data=fished, aes(y=Percentage, x=Age))

setwd(pop_dir)

saveRDS(Starting.Pop.For.Model, "L-miniatus_0.2F_Age_Structure")



#### OLD SET UP FOR ALL SPECIES #####
## Von Bertalanffy Parameters
# Linf = 945 #https://fishbase.mnhn.fr/popdyn/PopGrowthList.php?ID=1858&GenusName=Lethrinus&SpeciesName=miniatus&fc=328
# Linf = 623.9 # L sebae
# Linf = 664 # L nebulosus https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
# Linf = 666.4 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
# Linf = 590 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

# k = 0.166 # https://fishbase.mnhn.fr/popdyn/PopGrowthList.php?ID=1858&GenusName=Lethrinus&SpeciesName=miniatus&fc=328
# k = 0.2 # L sebae
# k = 0.241 # L nebulosus https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
# k = 0.221 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
# k = 0.26 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

# t0 = -1.75  # Used Wolfram Alpha and the equation from here https://www.fishbase.se/manual/key%20facts.htm with log10
# t0 = - 1.533 # L sebae
# t0 = -0.375 #  L Nebulosus https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
# t0 = -1.835 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
# t0 = -0.46 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

## Weight-Length Relationship
# WLa = 8.04e-06 #L miniatus https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=758b94902904411d6e1bb601ca365a718c208214
# WLa = 1.406e-08 # L sebae
# WLa = 0.000028 # L nebulosus This is for mm and g
# WLa = 0.932e-8 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
# Wla = 2.438e-5 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

# WLb = 3.127 # L miniatus https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=758b94902904411d6e1bb601ca365a718c208214
# WLb = 3.051 # L Sebae
# WLb = 2.8761 # L nebulosus
# WLb = 3.0924 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
# WLb = 2.9501 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

## Natural Mortality
# M = 0.23 # L Miniatus https://fish.gov.au/Archived-Reports/2012/reports/Documents/Leigh_et_al_2006_StockAssessment_RTE.pdf
# M = 0.104 # L sebae
# M = 0.146  #L nebulosus Marriot et al. 2011
# M = 0.192 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
# M = 0.164 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

# Selectivity
# A50 = 2.5 # L miniatus https://fish.gov.au/Archived-Reports/2012/reports/Documents/Leigh_et_al_2006_StockAssessment_RTE.pdf
# A50 = 3.762 # L sebae
# A50 = 3 # L nebulosus  https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
# A50 = 5.5 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
# A50 = 2.3 # P multidens Estimated in https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

# A95 =  mean(4.5, 4,4)# L miniatus https://fish.gov.au/Archived-Reports/2012/reports/Documents/Leigh_et_al_2006_StockAssessment_RTE.pdf
# A95 = 6.241 # L sebae
# A95 = 4 # L nebulosus  https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
# A95 = 8 # E multinotatus fully selected and estimated from https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
# A95 = 4.5 # Fully selected in the fishery https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

# MaxAge = 30
# MaxAge = 24 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
# MaxAge = 28 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

# Maturity
# In WA LM is probably around 371 and this paper (https://www.sciencedirect.com/science/article/abs/pii/S0165783615001885?via%3Dihub)says that
# you can estimate L95 as L50*1.15
# M50 = 1.3 # L miniatus
# M50 = 3.43 # L sebae
# M50 = 3.62 # L nebulosus From Marriott et al. 2011
# M50 = 2 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
# M50 = 4.8 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf


# M95 = 1.85 # L miniatus https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=758b94902904411d6e1bb601ca365a718c208214
# M95 = 4.37 # L sebae
# M95 = 5.97 # L nebulosus From Marriott et al. 2011
# M95 = 2.2849 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
# M95 = 11.95 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf




