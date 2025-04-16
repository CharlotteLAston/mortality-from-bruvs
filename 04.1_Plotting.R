###################################################
# Creating a plot of all the estimated F across
# all species and data types 
###################################################
library(L3Assess)
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
library(grid)
library(gridExtra)

rm(list = ls())

working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_data", sep="/")
sim_dir <- paste(working.dir, "Simulations", sep="/")

a4.width=160

#### Read in data ####
setwd(data_dir)

dat <- read.csv("catch_curve_results.csv")

#### Make data frame ####

dat <- dat %>% 
  mutate(Species = ifelse(Species %in% "King Wrasse", "Western King\nWrasse", Species)) %>% 
  mutate(Species = ifelse(Species %in% "Maori Wrasse", "Southern Maori\nWrasse", Species)) %>% 
  mutate(species.location = paste0(Species, sep="_", Location)) %>% 
  mutate(Method  = ifelse(Method %in% "Pauly", "Length Converted\nCatch Curve", Method)) %>% 
  mutate(Method  = ifelse(Method %in% "L3", "Length Based\nCatch Curve", Method)) %>% 
  mutate(species.location.method = paste0(species.location, sep="_", Method)) %>% 
  mutate(species.location.method.data = paste0(species.location.method, sep="_", Data)) %>% 
  mutate(Species=fct_relevel(Species, "Pink Snapper","Breaksea Cod","Western King\nWrasse","Southern Maori\nWrasse")) %>% 
  mutate(species.location=fct_relevel(species.location, "Pink Snapper_Metro/SW","Breaksea Cod_Metro","Western King\nWrasse_Metro","Western King\nWrasse_Jurien", "Southern Maori\nWrasse_Metro")) %>% 
  mutate(Shape = c(16,16,17,17,16,16,17,17,16,16,17,17,16,16,17,17,16,16,17,17))

#colours <- c("#88CBED", "#A9439A", "#332387", "#117633", "#43A999", "#872155", "#CB6778")
# species.labels <- c("*L. nebulosus*<br>(Ningaloo)", "*L. nebulosus*&#42;<br>(Ningaloo)","*E. armatus*<br>(Metro<br>and SW)", "*C. auratus*<br>(Metro<br>and SW)", "*C. auratus*&#42;<br>(Metro<br>and SW)",
#                     "*C. auricularis*<br>(Metro)", "*O. lineolatus*<br>(Metro)", "*C. auricularis*<br>(Abrolhos)", "*O. lineolatus*<br>(SW)")

facet_labels <- c("Pink Snapper_Metro/SW"="Pink Snapper\n(Metro/South-west)", "Breaksea Cod_Metro"="Breaksea Cod\n(Metro)",
                  "Western King\nWrasse_Metro"="Western King Wrasse\n(Metro)", "Western King\nWrasse_Jurien"="Western King Wrasse\n(Jurien Bay)", "Southern Maori\nWrasse_Metro"="Southern Maori Wrasse\n(Metro)")
#### Make Plot ####

condition_plot <- dat %>% 
  ggplot(.)+
  geom_point(aes(x=species.location.method.data, y=F_value, fill=Method, colour=Method, shape=Data), size=2)+
  #scale_fill_manual(values=c("#88CBED", "#88CBED", "#A9439A", "#332387", "#332387", "#117633","#43A999", "#872155", "#CB6778"), guide="none")+
  geom_linerange(aes(x=species.location.method.data, ymax=Upper_95_CI, ymin=Lower_95_CI, colour=Method))+
  #scale_linetype_manual(values=c("L3" = "solid", "Pauly" = "longdash"), name = "Catch-curve\nmethod")+
  scale_shape_manual(values = c("BRUV"=16, "DPIRD"=17), name="Data type", labels=c("Stereo-BRUV", "Line fishing"))+
  scale_fill_manual(values = c( "#302383", "#36753B"))+
  scale_colour_manual(values = c("#302383", "#36753B"))+
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
    default.unit="cm")
  )
  #theme(axis.text.x=ggtext::element_markdown())
condition_plot

setwd(fig_dir)
ggsave(condition_plot, filename="Condition_plot.png", height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300 )

#### Age vs. Length scatter plots ####
setwd(data_dir)

PS_colour <- "#80B1D4"

Linf = 1136
vbK = 0.12
tzero = -0.42

PS_dpird <- read.csv("Length and age data for Charlotte.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL, Age) %>% 
  filter(Species %in% "Pink Snapper") %>% 
  mutate(vonB = Linf*(1-exp(-vbK*(Age-tzero)))) %>% 
  mutate(Shape = ifelse(TL<410, 1, 19))

Age.MLL <- tzero-log(1-410/Linf)/vbK
  
PS_scatter <- PS_dpird %>% 
  ggplot()+
  geom_point(aes(x=Age, y=TL), shape=PS_dpird$Shape, colour=PS_colour, fill=PS_colour)+
  geom_line(aes(x=Age, y=vonB), colour="grey20")+
  ylab("Length (mm)")+
  geom_hline(yintercept=410, linetype="dashed", colour="grey20")+
  theme_classic()
PS_scatter

# ********************************************************************************************
BC_colour <- "#BBCC33"

Linf = 561
vbK = 0.14
tzero = -0.24

BC_dpird <- read.csv("Length and age data for Charlotte.csv") %>% 
  #dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "Breaksea Cod") %>% 
  filter(Region %in% "Metro") %>% 
  # mutate(vonB = Linf*(1-exp(-vbK*(Age-tzero)))) %>% 
  mutate(Age = tzero-log(1 - TL / Linf) / vbK) %>% 
  mutate(Shape = ifelse(TL<480, 1, 19))

# Age.MLL <- 

BC_scatter <- BC_dpird %>% 
  ggplot()+
  geom_point(aes(x=Age, y=TL), shape=BC_dpird$Shape, colour=BC_colour, fill=BC_colour)+
  #geom_line(aes(x=Age, y=vonB), colour="grey20")+
  ylab("Length (mm)")+
  geom_hline(yintercept=480, linetype="dashed", colour="grey20")+
  theme_classic()
BC_scatter



# ********************************************************************************************
WKW_colour <- "#BBCC33"

Linf = 378
vbK = 0.25
tzero = 0.10

WKW_dpird <- read.csv("Length and age data for Charlotte.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL, Age) %>% 
  filter(Species %in% "King Wrasse") %>% 
  mutate(vonB = Linf*(1-exp(-vbK*(Age-tzero)))) %>% 
  mutate(Shape = ifelse(TL<300, 1, 19))

# Age.MLL <- 

WKW_scatter <- WKW_dpird %>% 
  ggplot()+
  geom_point(aes(x=Age, y=TL), shape=WKW_dpird$Shape, colour=WKW_colour, fill=WKW_colour)+
  geom_line(aes(x=Age, y=vonB), colour="grey20")+
  ylab("Length (mm)")+
  geom_hline(yintercept=300, linetype="dashed", colour="grey20")+
  theme_classic()
WKW_scatter

# ********************************************************************************************
MW_colour <- "#FEB461"

Linf = 328
vbK = 0.28
tzero = 0.19

MW_dpird <- read.csv("Length and age data for Charlotte.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL, Age) %>% 
  filter(Species %in% "Maori Wrasse") %>% 
  mutate(vonB = Linf*(1-exp(-vbK*(Age-tzero)))) %>% 
  mutate(Shape = ifelse(TL<300, 1, 19))

# Age.MLL <- 

MW_scatter <- MW_dpird %>% 
  ggplot()+
  geom_point(aes(x=Age, y=TL), shape=MW_dpird$Shape, colour=MW_colour, fill=MW_colour)+
  geom_line(aes(x=Age, y=vonB), colour="grey20")+
  ylab("Length (mm)")+
  geom_hline(yintercept=300, linetype="dashed", colour="grey20")+
  theme_classic()
MW_scatter

#### Simulated vs actual line and stereo-BRUV data ####
line_colour <- "#BBCC33"
bruv_colour <- "#66CCEE"

#* Pink Snapper ###
SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
MaxAge = 40 # max model age
TimeStep = 1/12 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 0.12
FishMort = 1.5*NatMort
MaxLen = 1300
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(200, 50) # L50, L95-L50 for gear selectivity
#RetenParams = c(400, 50) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality
lbnd <- seq(0,MaxLen - LenInc, LenInc)
ubnd <- lbnd + LenInc
midpt <- lbnd + (LenInc/2)
nLenCl <- length(midpt)

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 1136
vbK = 0.12
CVSizeAtAge = 0.03
GrowthParams = c(Linf, vbK)
RefnceAges = NA
# generate line data
RetenParams = c((Linf*0.75)-50, 50) # L50, L95-L50 for retention #802
Res=Res_line=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, 
                                        SelectivityType, SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# plot length-frequency data (retained fish - line)
lendat_line = Res$ObsRetCatchFreqAtLen
lendat_line2 = Res$ObsDiscCatchFreqAtLen
# plot(Res$midpt, lendat_line, "o", ylim=c(0,150))

# generate BRUV data (with different realisation, i.e. starting with different random number)
MLL=NA
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# compare length-frequency data - BRUV
lendat_bruv = Res$ObsRetCatchFreqAtLen # Retained by line fishing
lendat_bruv2 = Res$ObsDiscCatchFreqAtLen # Discarded fish - difference between retained and selectivity?
lendat_bruv3 = lendat_bruv + lendat_bruv2
# lines(Res$midpt, lendat_bruv3, "o", col="red") 

setwd(data_dir)

PS_dpird <- read.csv("Length and age data for Charlotte.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL, Age) %>% 
  filter(Species %in% "Pink Snapper") %>% 
  mutate(LengthClass = trunc(TL/LenInc)+1,
         LenMidpt = midpt[LengthClass]) %>% 
  group_by(LenMidpt) %>% 
  summarise(Frequency = n()) %>% 
  complete(LenMidpt = midpt) %>% 
  mutate(Frequency = replace_na(Frequency, 0)) %>% 
  mutate(Freq.std = Frequency/max(Frequency)) 

PS_bruv <- readRDS("australian-synthesis_complete_length_chrysophrys_auratus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Marmion|TwoRocks|Two.Rocks|Capes|Geographe|SwC|south-west|Ngari|Warnbro|Rottnest"))) %>% 
  dplyr::filter(str_detect(campaign,"2006|2007|2008|2009|2010|2011")) %>% 
  dplyr::filter(!str_detect(campaign, "Perth.Canyon")) %>% 
  mutate(LengthClass = trunc(length/LenInc)+1,
         LenMidpt = midpt[LengthClass]) %>% 
  group_by(LenMidpt) %>% 
  summarise(Frequency = n()) %>% 
  complete(LenMidpt = midpt) %>% 
  mutate(Frequency = replace_na(Frequency, 0)) %>% 
  mutate(Freq.std = Frequency/max(Frequency)) 

sim_line_early_PS <- as.data.frame(lendat_line) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_line/max(lendat_line))

sim_bruv_early_PS <- as.data.frame(lendat_bruv3) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_bruv3/max(lendat_bruv3))

sim_line_late_PS <- as.data.frame(lendat_line) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_line/max(lendat_line))

sim_bruv_late_PS <- as.data.frame(lendat_bruv3) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_bruv3/max(lendat_bruv3))

comp_line_PS_high <- ggplot()+
  geom_line(sim_line_early_PS, mapping=aes(x=midpt, y=Freq.std), linetype="dotted", colour=line_colour, linewidth=1.2)+
  geom_line(sim_line_late_PS, mapping=aes(x=midpt, y=Freq.std), linetype="dashed", colour=line_colour, linewidth=0.75)+
  geom_line(PS_dpird, mapping=aes(x=LenMidpt, y=Freq.std), colour="grey50", linewidth=0.75)+
  xlab("Length (mm)")+
  ylab("Normalised frequency")+
  scale_x_continuous(breaks = seq(0, MaxLen, by=200))+
  theme_classic()
comp_line_PS_high

comp_bruv_PS_high <- ggplot()+
  geom_line(sim_bruv_early_PS, mapping=aes(x=midpt, y=Freq.std), linetype="dotted", colour=bruv_colour, linewidth=1.2)+
  geom_line(sim_bruv_late_PS, mapping=aes(x=midpt, y=Freq.std), linetype="dashed", colour=bruv_colour, linewidth=0.75)+
  geom_line(PS_bruv, mapping=aes(x=LenMidpt, y=Freq.std), colour="grey50", linewidth=0.75)+
  xlab("Length (mm)")+
  ylab("Normalised frequency")+
  scale_x_continuous(breaks = seq(0, MaxLen, by=200))+
  theme_classic()
comp_bruv_PS_high

setwd(fig_dir)
ggsave(comp_line, filename="Sim_comp_line_PS_highF.png", height=a4.width, width = a4.width, units  ="mm", dpi = 300 )
ggsave(comp_bruv, filename="Sim_comp_BRUV_PS_highF.png", height=a4.width, width = a4.width, units  ="mm", dpi = 300 )

#* Breaksea Cod ###
BC_colour <- "#CC6677"
SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
MaxAge = 20 # max model age
TimeStep = 1/4 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = exp(1.551-1.066*log(MaxAge))
FishMort = 1.5*NatMort
MaxLen = 560
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(100, 50) # L50, L95-L50 for gear selectivity

DiscMort = 0 # proportion of fish that die due to natural mortality
lbnd <- seq(0,MaxLen - LenInc, LenInc)
ubnd <- lbnd + LenInc
midpt <- lbnd + (LenInc/2)
nLenCl <- length(midpt)

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 561
vbK = 0.14
CVSizeAtAge = 0.03
GrowthParams = c(Linf, vbK)
RefnceAges = NA
# generate line data
RetenParams = c((Linf*0.75)-50, 50) # L50, L95-L50 for retention #230
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# plot length-frequency data (retained fish - line)
lendat_line = Res$ObsRetCatchFreqAtLen
lendat_line2 = Res$ObsDiscCatchFreqAtLen
# plot(Res$midpt, lendat_line, "o", ylim=c(0,150))

# generate BRUV data (with different realisation, i.e. starting with different random number)
MLL=NA
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# compare length-frequency data - BRUV
lendat_bruv = Res$ObsRetCatchFreqAtLen # Retained by line fishing
lendat_bruv2 = Res$ObsDiscCatchFreqAtLen # Discarded fish - difference between retained and selectivity?
lendat_bruv3 = lendat_bruv + lendat_bruv2
# lines(Res$midpt, lendat_bruv3, "o", col="red") 

setwd(data_dir)

BC_dpird <- read.csv("Length data for Charlotte.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "Breaksea Cod") %>% 
  filter(Region %in% "Metro") %>% 
  mutate(LengthClass = trunc(TL/LenInc)+1,
         LenMidpt = midpt[LengthClass]) %>% 
  group_by(LenMidpt) %>% 
  summarise(Frequency = n()) %>% 
  complete(LenMidpt = midpt) %>% 
  mutate(Frequency = replace_na(Frequency, 0)) %>% 
  mutate(Freq.std = Frequency/max(Frequency)) 

BC_bruv <- readRDS("australian-synthesis_complete_length_epinephelides-armatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Marmion|TwoRocks|Two.Rocks|Warnbro|Rottnest"))) %>% #Capes|Geographe|SwC|south-west|Ngari|
  dplyr::filter(str_detect(campaign,"2006|2007|2008|2009|2010|2011")) %>% 
  dplyr::filter(!str_detect(campaign, "Perth.Canyon")) %>% 
  mutate(LengthClass = trunc(length/LenInc)+1,
         LenMidpt = midpt[LengthClass]) %>% 
  group_by(LenMidpt) %>% 
  summarise(Frequency = n()) %>% 
  complete(LenMidpt = midpt) %>% 
  mutate(Frequency = replace_na(Frequency, 0)) %>% 
  mutate(Freq.std = Frequency/max(Frequency)) 

sim_line_early_BC <- as.data.frame(lendat_line) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_line/max(lendat_line))

sim_bruv_early_BC <- as.data.frame(lendat_bruv3) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_bruv3/max(lendat_bruv3))

sim_line_late_BC <- as.data.frame(lendat_line) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_line/max(lendat_line))

sim_bruv_late_BC <- as.data.frame(lendat_bruv3) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_bruv3/max(lendat_bruv3))

comp_line_BC_high <- ggplot()+
  geom_line(sim_line_early_BC, mapping=aes(x=midpt, y=Freq.std), linetype="dotted", colour=line_colour, linewidth=1.2)+
  geom_line(sim_line_late_BC, mapping=aes(x=midpt, y=Freq.std), linetype="dashed", colour=line_colour, linewidth=0.75)+
  geom_line(BC_dpird, mapping=aes(x=LenMidpt, y=Freq.std), colour="grey50", linewidth=0.75)+
  xlab("Length (mm)")+
  ylab("Normalised frequency")+
  scale_x_continuous(breaks = seq(0, MaxLen, by=100))+
  theme_classic()
comp_line_BC_high

comp_bruv_BC_high <- ggplot()+
  geom_line(sim_bruv_early_BC, mapping=aes(x=midpt, y=Freq.std), linetype="dotted", colour=bruv_colour, linewidth=1.2)+
  geom_line(sim_bruv_late_BC, mapping=aes(x=midpt, y=Freq.std), linetype="dashed", colour=bruv_colour, linewidth=0.75)+
  geom_line(BC_bruv, mapping=aes(x=LenMidpt, y=Freq.std), colour="grey50", linewidth=0.75)+
  xlab("Length (mm)")+
  ylab("Normalised frequency")+
  scale_x_continuous(breaks = seq(0, MaxLen, by=100))+
  theme_classic()
comp_bruv_BC_high

setwd(fig_dir)
ggsave(comp_line, filename="Sim_comp_line_BC_highF.png", height=a4.width, width = a4.width, units  ="mm", dpi = 300 )
ggsave(comp_bruv, filename="Sim_comp_BRUV_BC_highF.png", height=a4.width, width = a4.width, units  ="mm", dpi = 300 )

#### Put Pink Snapper and Breaksea Plots together ####
a4.height=290

setwd(fig_dir)
x.label <- textGrob("Length (mm)", gp=gpar(fontsize=12))
y.label <- textGrob("Normalised frequency", gp=gpar(fontsize=12), rot=90)


PSxBC.sim.comp <- grid.arrange(arrangeGrob(comp_line_PS_high+xlab(NULL)+ylab(NULL)+
                                             ggtitle("High mortality")+
                                             theme(plot.title = element_text(size=9, hjust=0.5))+
                                             geom_text(aes(label="a)", x=25, y=1)),
                                           comp_line_PS_low+xlab(NULL)+ylab(NULL)+
                                             ggtitle("Low mortality")+
                                             theme(plot.title = element_text(size=9, hjust=0.5))+
                                             geom_text(aes(label="b)", x=25, y=1)),
                                           comp_line_BC_high+xlab(NULL)+ylab(NULL)+
                                             ggtitle("High mortality")+
                                             theme(plot.title = element_text(size=9, hjust=0.5))+
                                             geom_text(aes(label="c)", x=25, y=1)),
                                           comp_line_BC_low+xlab(NULL)+ylab(NULL)+
                                             ggtitle("Low mortality")+
                                             theme(plot.title = element_text(size=9, hjust=0.5))+
                                             geom_text(aes(label="d)", x=25, y=1)),
                                           comp_bruv_PS_high+xlab(NULL)+ylab(NULL)+
                                             geom_text(aes(label="e)", x=25, y=1)),
                                           comp_bruv_PS_low+xlab(NULL)+ylab(NULL)+
                                             geom_text(aes(label="f)", x=25, y=1)),
                                           comp_bruv_BC_high+xlab(NULL)+ylab(NULL)+
                                             geom_text(aes(label="g)", x=25, y=1)),
                                           comp_bruv_BC_low+xlab(NULL)+ylab(NULL)+
                                             geom_text(aes(label="h)", x=25, y=1)),
                                           left=y.label,
                                           bottom=x.label,
                                           nrow=2, ncol=4))

ggsave(PSxBC.sim.comp, filename="PSxBC.sim.comp.png",height = a4.width*1.1, width = a4.height, units  ="mm", dpi = 300 )


#* Western King Wrasse ###
WKW_colour <- "#44AA99"
SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
MaxAge = 10.5 # max model age
TimeStep = 1/12 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = exp(1.551-1.066*log(10.5))
FishMort = 0.67*NatMort
MaxLen = 500
LenInc = 10
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(100, 50) # L50, L95-L50 for gear selectivity

DiscMort = 0 # proportion of fish that die due to natural mortality
lbnd <- seq(0,MaxLen - LenInc, LenInc)
ubnd <- lbnd + LenInc
midpt <- lbnd + (LenInc/2)
nLenCl <- length(midpt)

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 378
vbK = 0.25
CVSizeAtAge = 0.03
GrowthParams = c(Linf, vbK)
RefnceAges = NA
# generate line data
RetenParams = c((Linf*0.75)-50, 50) # L50, L95-L50 for retention
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# plot length-frequency data (retained fish - line)
lendat_line = Res$ObsRetCatchFreqAtLen
lendat_line2 = Res$ObsDiscCatchFreqAtLen
# plot(Res$midpt, lendat_line, "o", ylim=c(0,150))

# generate BRUV data (with different realisation, i.e. starting with different random number)
MLL=NA
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# compare length-frequency data - BRUV
lendat_bruv = Res$ObsRetCatchFreqAtLen # Retained by line fishing
lendat_bruv2 = Res$ObsDiscCatchFreqAtLen # Discarded fish - difference between retained and selectivity?
lendat_bruv3 = lendat_bruv + lendat_bruv2
# lines(Res$midpt, lendat_bruv3, "o", col="red") 

setwd(data_dir)

WKW_dpird <- read.csv("Length data for Charlotte UPDATED.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "King Wrasse") %>% 
  filter(Region %in% "Jurien") %>% 
  mutate(LengthClass = trunc(TL/LenInc)+1,
         LenMidpt = midpt[LengthClass]) %>% 
  group_by(LenMidpt) %>% 
  summarise(Frequency = n()) %>% 
  complete(LenMidpt = midpt) %>% 
  mutate(Frequency = replace_na(Frequency, 0)) %>% 
  mutate(Freq.std = Frequency/max(Frequency)) 

WKW_bruv <- readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) %>%
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009"))  %>% 
  mutate(LengthClass = trunc(length/LenInc)+1,
         LenMidpt = midpt[LengthClass]) %>% 
  group_by(LenMidpt) %>% 
  summarise(Frequency = n()) %>% 
  complete(LenMidpt = midpt) %>% 
  mutate(Frequency = replace_na(Frequency, 0)) %>% 
  mutate(Freq.std = Frequency/max(Frequency)) 

sim_line_early_WKW <- as.data.frame(lendat_line) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_line/max(lendat_line))

sim_bruv_early_WKW <- as.data.frame(lendat_bruv3) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_bruv3/max(lendat_bruv3))

sim_line_late_WKW <- as.data.frame(lendat_line) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_line/max(lendat_line))

sim_bruv_late_WKW <- as.data.frame(lendat_bruv3) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_bruv3/max(lendat_bruv3))

comp_line_WKW_low <- ggplot()+
  geom_line(sim_line_early_WKW, mapping=aes(x=midpt, y=Freq.std), linetype="dotted", colour=line_colour, linewidth=1.2)+
  geom_line(sim_line_late_WKW, mapping=aes(x=midpt, y=Freq.std), linetype="dashed", colour=line_colour, linewidth=0.75)+
  geom_line(WKW_dpird, mapping=aes(x=LenMidpt, y=Freq.std), colour="grey50", linewidth=0.75)+
  xlab("Length (mm)")+
  ylab("Normalised frequency")+
  scale_x_continuous(breaks = seq(0, MaxLen, by=100))+
  theme_classic()
comp_line_WKW_low

comp_bruv_WKW_low <- ggplot()+
  geom_line(sim_bruv_early_WKW, mapping=aes(x=midpt, y=Freq.std), linetype="dotted", colour=bruv_colour, linewidth=1.2)+
  geom_line(sim_bruv_late_WKW, mapping=aes(x=midpt, y=Freq.std), linetype="dashed", colour=bruv_colour, linewidth=0.75)+
  geom_line(WKW_bruv, mapping=aes(x=LenMidpt, y=Freq.std), colour="grey50", linewidth=0.75)+
  xlab("Length (mm)")+
  ylab("Normalised frequency")+
  scale_x_continuous(breaks = seq(0, MaxLen, by=100))+
  theme_classic()
comp_bruv_WKW_low

setwd(fig_dir)
ggsave(comp_line, filename="Sim_comp_line_WKW_highF.png", height=a4.width, width = a4.width, units  ="mm", dpi = 300 )
ggsave(comp_bruv, filename="Sim_comp_BRUV_WKW_highF.png", height=a4.width, width = a4.width, units  ="mm", dpi = 300 )

#* Southern Maori Wrasse ###
SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
MaxAge = 12 # max model age
TimeStep = 1/6 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = exp(1.551-1.066*log(11.8))
FishMort = 0.67*NatMort
MaxLen = 470
LenInc = 10
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(100, 50) # L50, L95-L50 for gear selectivity

DiscMort = 0 # proportion of fish that die due to natural mortality
lbnd <- seq(0,MaxLen - LenInc, LenInc)
ubnd <- lbnd + LenInc
midpt <- lbnd + (LenInc/2)
nLenCl <- length(midpt)

# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 328
vbK = 0.28
CVSizeAtAge = 0.03
GrowthParams = c(Linf, vbK)
RefnceAges = NA
# generate line data
RetenParams = c((Linf*0.75)-50, 50) # L50, L95-L50 for retention
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# plot length-frequency data (retained fish - line)
lendat_line = Res$ObsRetCatchFreqAtLen
lendat_line2 = Res$ObsDiscCatchFreqAtLen
# plot(Res$midpt, lendat_line, "o", ylim=c(0,150))

# generate BRUV data (with different realisation, i.e. starting with different random number)
MLL=NA
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                               SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

# compare length-frequency data - BRUV
lendat_bruv = Res$ObsRetCatchFreqAtLen # Retained by line fishing
lendat_bruv2 = Res$ObsDiscCatchFreqAtLen # Discarded fish - difference between retained and selectivity?
lendat_bruv3 = lendat_bruv + lendat_bruv2
# lines(Res$midpt, lendat_bruv3, "o", col="red") 

setwd(data_dir)

MW_dpird <- read.csv("Length data for Charlotte UPDATED.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL) %>% 
  filter(Species %in% "Maori Wrasse") %>% 
  filter(Region %in% "Metro") %>% 
  mutate(LengthClass = trunc(TL/LenInc)+1,
         LenMidpt = midpt[LengthClass]) %>% 
  group_by(LenMidpt) %>% 
  summarise(Frequency = n()) %>% 
  complete(LenMidpt = midpt) %>% 
  mutate(Frequency = replace_na(Frequency, 0)) %>% 
  mutate(Freq.std = Frequency/max(Frequency)) 

MW_bruv <- readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Rottnest|TwoRocks|Marmion|Warnbro|Two.Rocks"))) %>% 
  dplyr::filter(str_detect(campaign,"2006|2007|2008|2009|2010|2011")) %>% 
  mutate(LengthClass = trunc(length/LenInc)+1,
         LenMidpt = midpt[LengthClass]) %>% 
  group_by(LenMidpt) %>% 
  summarise(Frequency = n()) %>% 
  complete(LenMidpt = midpt) %>% 
  mutate(Frequency = replace_na(Frequency, 0)) %>% 
  mutate(Freq.std = Frequency/max(Frequency)) 

sim_line_early_SMW <- as.data.frame(lendat_line) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_line/max(lendat_line))

sim_bruv_early_SMW <- as.data.frame(lendat_bruv3) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_bruv3/max(lendat_bruv3))

sim_line_late_SMW <- as.data.frame(lendat_line) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_line/max(lendat_line))

sim_bruv_late_SMW <- as.data.frame(lendat_bruv3) %>% 
  mutate(midpt = midpt) %>% 
  mutate(Freq.std = lendat_bruv3/max(lendat_bruv3))

comp_line_low_SMW <- ggplot()+
  geom_line(sim_line_early_SMW, mapping=aes(x=midpt, y=Freq.std), linetype="dotted", colour=line_colour, linewidth=1.2)+
  geom_line(sim_line_late_SMW, mapping=aes(x=midpt, y=Freq.std), linetype="dashed", colour=line_colour, linewidth=1.2)+
  geom_line(MW_dpird, mapping=aes(x=LenMidpt, y=Freq.std), colour="grey50", linewidth=0.75)+
  xlab("Length (mm)")+
  ylab("Normalised frequency")+
  scale_x_continuous(breaks = seq(0, MaxLen, by=100))+
  theme_classic()
comp_line_low_SMW

comp_bruv_low_SMW <- ggplot()+
  geom_line(sim_bruv_early_SMW, mapping=aes(x=midpt, y=Freq.std), linetype="dotted", colour=bruv_colour, linewidth=1.2)+
  geom_line(sim_bruv_late_SMW, mapping=aes(x=midpt, y=Freq.std), linetype="dashed", colour=bruv_colour, linewidth=0.75)+
  geom_line(MW_bruv, mapping=aes(x=LenMidpt, y=Freq.std), colour="grey50", linewidth=0.75)+
  xlab("Length (mm)")+
  ylab("Normalised frequency")+
  scale_x_continuous(breaks = seq(0, MaxLen, by=100))+
  theme_classic()
comp_bruv_low_SMW

setwd(fig_dir)
ggsave(comp_line, filename="Sim_comp_line_MW_highF.png", height=a4.width, width = a4.width, units  ="mm", dpi = 300 )
ggsave(comp_bruv, filename="Sim_comp_BRUV_MW_highF.png", height=a4.width, width = a4.width, units  ="mm", dpi = 300 )


#### Put Western King Wrasse and Southern Maori Wrasse Plots together ####
a4.height=290

setwd(fig_dir)
x.label <- textGrob("Length (mm)", gp=gpar(fontsize=12))
y.label <- textGrob("Normalised frequency", gp=gpar(fontsize=12), rot=90)


WKWxSMW.sim.comp <- grid.arrange(arrangeGrob(comp_line_WKW_high+xlab(NULL)+ylab(NULL)+
                                             ggtitle("High mortality")+
                                             theme(plot.title = element_text(size=9, hjust=0.5))+
                                               geom_text(aes(label = "a)", x=25, y=1)),
                                           comp_line_WKW_low+xlab(NULL)+ylab(NULL)+
                                             ggtitle("Low mortality")+
                                             theme(plot.title = element_text(size=9, hjust=0.5))+
                                             geom_text(aes(label = "b)", x=25, y=1)),
                                           comp_line_high_SMW+xlab(NULL)+ylab(NULL)+
                                             ggtitle("High mortality")+
                                             theme(plot.title = element_text(size=9, hjust=0.5))+
                                             geom_text(aes(label = "c)", x=25, y=1)),
                                           comp_line_low_SMW+xlab(NULL)+ylab(NULL)+
                                             ggtitle("Low mortality")+
                                             theme(plot.title = element_text(size=9, hjust=0.5))+
                                             geom_text(aes(label = "d)", x=25, y=1)),
                                           comp_bruv_WKW_high+xlab(NULL)+ylab(NULL)+
                                             geom_text(aes(label = "e)", x=25, y=1)),
                                           comp_bruv_WKW_low+xlab(NULL)+ylab(NULL)+
                                             geom_text(aes(label = "f)", x=25, y=1)),
                                           comp_bruv_high_SMW+xlab(NULL)+ylab(NULL)+
                                             geom_text(aes(label = "g)", x=25, y=1)),
                                           comp_bruv_low_SMW +xlab(NULL)+ylab(NULL)+
                                             geom_text(aes(label = "h)", x=25, y=1)),
                                           left=y.label,
                                           bottom=x.label,
                                           nrow=2, ncol=4))

ggsave(WKWxSMW.sim.comp, filename="WKWxSMW.sim.comp.png",height = a4.width*1.1, width = a4.height, units  ="mm", dpi = 300 )


#### Additional plots for the paper ####
## NEED TO RUN THE PLOTS IN THEIR ORIGINAL SCRIPTS 

# Comparisons of selectivity

setwd(fig_dir)
x.label <- textGrob("Length (mm)", gp=gpar(fontsize=12))
y.label <- textGrob("Probability", gp=gpar(fontsize=12), rot=90)

a4.height=290

selectivity_plots <- grid.arrange(arrangeGrob(selectivity_plot_auratus + geom_text(aes(label = "a)", x=25, y=1)),
                                              selectivity_plot_armatus + geom_text(aes(label = "b)", x=25, y=1)),
                                              selectivity_plot_lineolatus + geom_text(aes(label = "c)", x=25, y=1)),
                                              left=y.label,
                                              bottom=x.label,
                                              nrow=2, ncol=2))

ggsave(selectivity_plots, filename="Selectivity_comp.png",height = a4.height*0.5, width = a4.width, units  ="mm", dpi = 300 )

y.label <- textGrob(expression(paste(italic("F") ~ (year^{-1}))), gp=gpar(fontsize=12), rot=90)

mortlaity_plots <- grid.arrange(arrangeGrob(mortality.at.length_auratus + xlab(NULL)+ylab(NULL) + ylim(0,1) 
                                            + scale_x_continuous(breaks=seq(0,1300, 200))+ geom_text(aes(label = "a)", x=25, y=1)),
                                            mortality.at.length_armatus + xlab(NULL)+ylab(NULL)+ ylim(0,1)
                                            + scale_x_continuous(breaks=seq(0,600, 100)) + geom_text(aes(label = "b)", x=25, y=1)),
                                            mortality.at.length_lineolatus + xlab(NULL)+ylab(NULL)+ ylim(0,1)
                                            +scale_x_continuous(breaks=seq(0,500, 100)) + geom_text(aes(label = "c)", x=25, y=1)),
                                              left=y.label,
                                              bottom=x.label,
                                              nrow=2, ncol=2))

ggsave(mortlaity_plots, filename="Mortality_at_length_comp.png",height = a4.height*0.5, width = a4.width, units  ="mm", dpi = 300 )
