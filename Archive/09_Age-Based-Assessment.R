###################################################
# Fitting age based catch curves to
# DPIRD length data for all species
###################################################
library(tidyr)
library(dplyr)
library(stringr)
library(sf)
library(ggplot2)
library(forcats)
library(ggtext)
library(tidyverse)
library(sf)
library(terra)
library(ggplot2)
library(patchwork)
library(ggnewscale)
library(scales)
library(ggspatial)
library(tidyterra)
library(L3Assess)

rm(list = ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_data", sep="/")

setwd(working.dir)
source("X_My-Functions.R")
a4.width=160

#### Read in data ####
setwd(data_dir)

Age_dpird <- read.csv("Length and age data for Charlotte.csv") %>% 
  filter(CaptureMeth == "LINE") %>% 
  filter(!is.na(Age))


#### Pink Snapper #####
set.seed(123)
MaxAge = 40
NatMort = 0.12 #4.22/MaxAge Probably more around 0.12 
# Durelle method? Dule and Froese
MaxLen = 1300
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 1136
vbK = 0.12
CVSizeAtAge = 0.075
GrowthParams = c(Linf, vbK)
RefnceAges = NA
DiscMort = 0

#* Format data ####
Age_PS <- Age_dpird %>% 
  filter(Species %in% "Pink Snapper")

MinAge = 1
MaxAge = 41

LenInterval = 1
LenCats <- seq(from=0, to=MaxAge, by=LenInterval)
LenCats

HistData <- hist(Age_PS$Age, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData
ObsAgeFreq = as.vector(HistData$counts)
Ages = MinAge:MaxAge


res=GetChapmanRobsonMortalityResults(RecAssump=1, SpecRecAge=NA, MinAge, MaxAge, ObsAgeFreq)
res$EstZMort

# Logistic catch curve 

CatchCurveModel = 3 # Logistic selectivity
RecAssump = NA # Logistic selectivity
NatMort = 0.104 # Logistic selectivity
MinFreq = NA # # Logistic selectivity
Init_FMort = 0.2
Init_SelA50 = 3
Init_SelDelta = 2
params = log(c(Init_FMort, Init_SelA50, Init_SelDelta))
res=GetLogisticCatchCurveResults(params, NatMort, Ages, ObsAgeFreq)

res$ParamEst

PlotAgeBasedCatchCurveResults_NormalSpace(RecAssump, SpecRecAge=NA, MinFreq, MinAge, MaxAge, NatMort,
                                          ObsAgeFreq, CatchCurveModel, MainLabel=NA,
                                          xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
                                          ymax=NA, yint=NA, PlotCLs=T)

#### Read in data ####
setwd(data_dir)

Age_dpird <- read.csv("Length and age data for Charlotte.csv") %>% 
  filter(CaptureMeth == "LINE") %>% 
  filter(!is.na(Age))


#### King Wrasse Jurien ####
set.seed(123)
MaxAge = 10.2
NatMort = 0.3 #4.22/MaxAge
FishMort = 0.2
MaxLen = 400
LenInc = 20


# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 378
vbK = 0.25
GrowthParams = c(Linf, vbK)
RefnceAges = NA
DiscMort = 0


#* Format data ####
Age_KW <- Age_dpird %>% 
  filter(Species %in% "King Wrasse") %>% 
  filter(Location %in% "Jurien")

MinAge = 1
MaxAge = 10

LenInterval = 1
LenCats <- seq(from=0, to=MaxAge, by=LenInterval)
LenCats

HistData <- hist(Age_KW$Age, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData
ObsAgeFreq = as.vector(HistData$counts)
Ages = MinAge:MaxAge

res=GetChapmanRobsonMortalityResults(RecAssump=0, SpecRecAge=NA, MinAge, MaxAge, ObsAgeFreq)
res$EstZMort

# Logistic catch curve
CatchCurveModel = 3 # Logistic selectivity
RecAssump = NA # Logistic selectivity
NatMort = 0.3 # Logistic selectivity
MinFreq = NA # # Logistic selectivity
Init_FMort = 0.2
Init_SelA50 = 5
Init_SelDelta = 2
params = log(c(Init_FMort, Init_SelA50, Init_SelDelta))
res=GetLogisticCatchCurveResults(params, NatMort, Ages, ObsAgeFreq)

res$ParamEst3

PlotAgeBasedCatchCurveResults_NormalSpace(RecAssump, SpecRecAge=NA, MinFreq, MinAge, MaxAge, NatMort,
                                          ObsAgeFreq, CatchCurveModel, MainLabel=NA,
                                          xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
                                          ymax=NA, yint=NA, PlotCLs=T)



#### Maori Wrasse Metro ####
set.seed(123)
MaxAge = 11.8
NatMort = 0.36 #4.22/MaxAge
FishMort = 0.2
MaxLen = 470
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 328
vbK = 0.28
CVSizeAtAge = 0.05
GrowthParams = c(Linf, vbK)
RefnceAges = NA
DiscMort = 0


#* Format data ####
Age_MW <- Age_dpird %>% 
  filter(Species %in% "Maori Wrasse") %>% 
  filter(Region %in% "Metro")

MinAge = 1
MaxAge = 11

LenInterval = 1
LenCats <- seq(from=0, to=MaxAge, by=LenInterval)
LenCats

HistData <- hist(Age_KW$Age, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData
ObsAgeFreq = as.vector(HistData$counts)
Ages = MinAge:MaxAge


res=GetChapmanRobsonMortalityResults(RecAssump=0, SpecRecAge=NA, MinAge, MaxAge, ObsAgeFreq)
res$EstZMort


# Logistic catch curve
PlotAgeBasedCatchCurveResults_NormalSpace(RecAssump, SpecRecAge=NA, MinFreq, MinAge, MaxAge, NatMort,
                                          ObsAgeFreq, CatchCurveModel, MainLabel=NA,
                                          xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
                                          ymax=NA, yint=NA, PlotCLs=T)

#### Brownspot Wrasse ####
#* Metro ####
Age_BW <- Age_dpird %>% 
  filter(Species %in% "Brownspot Wrasse") %>% 
  filter(Region %in% "Metro")

MinAge = 1
MaxAge = 12
NatMort = 4.22/12

LenInterval = 1
LenCats <- seq(from=0, to=MaxAge, by=LenInterval)
LenCats

HistData <- hist(Age_BW$Age, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData
ObsAgeFreq = as.vector(HistData$counts)
Ages = MinAge:MaxAge

res=GetChapmanRobsonMortalityResults(RecAssump=0, SpecRecAge=NA, MinAge, MaxAge, ObsAgeFreq)
res$EstZMort


# Logistic catch curve
CatchCurveModel = 3 # Logistic selectivity
RecAssump = NA # Logistic selectivity
MinFreq = NA # # Logistic selectivity
Init_FMort = 0.4
Init_SelA50 = 5
Init_SelDelta = 2
params = log(c(Init_FMort, Init_SelA50, Init_SelDelta))
res=GetLogisticCatchCurveResults(params, NatMort, Ages, ObsAgeFreq)

res$ParamEst

PlotAgeBasedCatchCurveResults_NormalSpace(RecAssump, SpecRecAge=NA, MinFreq, MinAge, MaxAge, NatMort,
                                          ObsAgeFreq, CatchCurveModel, MainLabel=NA,
                                          xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
                                          ymax=NA, yint=NA, PlotCLs=T)

#### Brownspot Wrasse ####
#* Jurien ####
Age_BW <- Age_dpird %>% 
  filter(Species %in% "Brownspot Wrasse") %>% 
  filter(Region %in% "Jurien")

MinAge = 1
MaxAge = 12
NatMort = 4.22/12

LenInterval = 1
LenCats <- seq(from=0, to=MaxAge, by=LenInterval)
LenCats

HistData <- hist(Age_BW$Age, breaks=LenCats, right=FALSE, col="light blue", main="", xlab="20 mm Length category", ylab="Frequency",las=1)
HistData
ObsAgeFreq = as.vector(HistData$counts)
Ages = MinAge:MaxAge

res=GetChapmanRobsonMortalityResults(RecAssump=0, SpecRecAge=NA, MinAge, MaxAge, ObsAgeFreq)
res$EstZMort


# Logistic catch curve
CatchCurveModel = 3 # Logistic selectivity
RecAssump = NA # Logistic selectivity
MinFreq = NA # # Logistic selectivity
Init_FMort = 0.1
Init_SelA50 = 5
Init_SelDelta = 2
params = log(c(Init_FMort, Init_SelA50, Init_SelDelta))
res=GetLogisticCatchCurveResults(params, NatMort, Ages, ObsAgeFreq)

res$ParamEst

PlotAgeBasedCatchCurveResults_NormalSpace(RecAssump, SpecRecAge=NA, MinFreq, MinAge, MaxAge, NatMort,
                                          ObsAgeFreq, CatchCurveModel, MainLabel=NA,
                                          xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
                                          ymax=NA, yint=NA, PlotCLs=T)



#### Read in data ####
setwd(data_dir)

dat <- read.csv("Age-based_Catch-Curve.csv")

#### Make data frame ####

dat <- dat %>% 
  mutate(Species = ifelse(Species %in% "King Wrasse", "Western King\nWrasse", Species)) %>% 
  mutate(Species = ifelse(Species %in% "Southern Maori Wrasse", "Southern Maori\nWrasse", Species)) %>% 
  mutate(Method  = ifelse(Method %in% "Pauly", "Length Converted\nCatch Curve", Method)) %>% 
  mutate(Method  = ifelse(Method %in% "L3", "Length Based\nCatch Curve", Method)) %>% 
  mutate(Method  = ifelse(Method %in% "Chapman and Robson", "Chapman and Robson\nMortality Estimator", Method)) %>% 
  mutate(species.location = paste0(Species, sep="_", Location)) %>% 
  mutate(species.location.method = paste0(species.location, sep="_", Method)) %>% 
  mutate(species.location.method.data = paste0(species.location.method, sep="_", Data)) %>% 
  mutate(Species=fct_relevel(Species, "Pink Snapper","Western King\nWrasse","Southern Maori\nWrasse")) %>% 
  mutate(species.location=fct_relevel(species.location, "Pink Snapper_Metro/SW","Western King\nWrasse_Jurien","Southern Maori\nWrasse_Metro")) %>%
  mutate(Method = fct_relevel(as.factor(Method), "Linear Catch Curve", "Chapman and Robson\nMortality Estimator", "Length Converted\nCatch Curve", "Length Based\nCatch Curve")) %>% 
  mutate(species.location.method.data = fct_relevel(species.location.method.data, "Pink Snapper_Metro/SW_Linear Catch Curve_Age", "Pink Snapper_Metro/SW_Chapman and Robson\nMortality Estimator_Age" ,
                                                    "Pink Snapper_Metro/SW_Length Converted\nCatch Curve_Length", "Pink Snapper_Metro/SW_Length Based\nCatch Curve_Length", "Western King\nWrasse_Jurien_Linear Catch Curve_Age",
                                                    "Western King\nWrasse_Jurien_Chapman and Robson\nMortality Estimator_Age", "Western King\nWrasse_Jurien_Length Converted\nCatch Curve_Length", 
                                                    "Western King\nWrasse_Jurien_Length Based\nCatch Curve_Length", "Southern Maori\nWrasse_Metro_Linear Catch Curve_Age", "Southern Maori\nWrasse_Metro_Chapman and Robson\nMortality Estimator_Age",
                                                    "Southern Maori\nWrasse_Metro_Length Converted\nCatch Curve_Length", "Southern Maori\nWrasse_Metro_Length Based\nCatch Curve_Length"))
 # mutate(Shape = c(16,16,17,17,16,16,17,17,16,16,17,17,16,16,17,17,16,16,17,17))

#colours <- c("#88CBED", "#A9439A", "#332387", "#117633", "#43A999", "#872155", "#CB6778")
# species.labels <- c("*L. nebulosus*<br>(Ningaloo)", "*L. nebulosus*&#42;<br>(Ningaloo)","*E. armatus*<br>(Metro<br>and SW)", "*C. auratus*<br>(Metro<br>and SW)", "*C. auratus*&#42;<br>(Metro<br>and SW)",
#                     "*C. auricularis*<br>(Metro)", "*O. lineolatus*<br>(Metro)", "*C. auricularis*<br>(Abrolhos)", "*O. lineolatus*<br>(SW)")

facet_labels <- c("Pink Snapper_Metro/SW"="Pink Snapper\n(Metro/South-west)", 
                  "Western King\nWrasse_Jurien"="Western King Wrasse\n(Jurien Bay)", "Southern Maori\nWrasse_Metro"="Southern Maori Wrasse\n(Metro)")
#### Make Plot ####

age_cc_plot <- dat %>% 
  ggplot(.)+
  geom_linerange(aes(x=species.location.method.data, ymax=Upper_95_CI, ymin=Lower_95_CI, colour=Method))+
  geom_point(aes(x=species.location.method.data, y=F_value, fill=Method, colour=Method, shape=Data), size=2)+
  scale_fill_manual(values = c("#66CCEE","#BBCC33", "#36753B", "#302383"))+
  scale_colour_manual(values = c("#66CCEE","#BBCC33", "#36753B",  "#302383"))+
  scale_shape_manual(values = c("Age"=16, "Length"=17), name="Data type", labels=c("Age", "Length"))+
  facet_grid(~species.location, scales="free", switch = "x", labeller = as_labeller(facet_labels))+
  theme_classic()+
  #scale_x_discrete(labels=species.labels)+
  ylim(0,1)+
  ylab("Fishing Mortality")+
  xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(color = NA),
        strip.placement = "outside")+
  guides(fill=guide_legend(
    keyheight=1,
    default.unit="cm",
    order=1),
    colour=guide_legend(order=1)
  )
age_cc_plot

setwd(fig_dir)
ggsave(age_cc_plot, filename="Age_assessment_plot.png", height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300 )






