###################################################

# Calculating proportions of different length
# classes and making plots 

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
library(rcartocolor)

rm(list = ls())

Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend.background = element_rect(fill="white"),
    legend.background = element_blank(),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = c(0.2, 0.8),
    text=element_text(size=10),
    strip.text.y = element_text(size = 10,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=12),
    axis.title.y=element_text(vjust=0.6, angle=90, size=12),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
pop_dir <- paste(working.dir, "Population_Files", sep="/")

#### L NEBULOSUS ####
#* Von Bertalanffy Parameters ####
Linf = 664 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
k = 0.241 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
t0 = -0.375 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr

# Maturity
M50 = 3.62 # From Marriott et al. 2011
M95 = 5.97 # From Marriott et al. 2011

#* Read in data ####
setwd(pop_dir)
F.at.0 <- readRDS("L-nebulosus_0F_Age_Structure") %>% 
  mutate(Mort.Group = "No fishing mortality")
F.at.02M <- readRDS("L-nebulosus_0.2F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 0.2M")
F.at.M23 <- readRDS("L-nebulosus_2.3F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 2/3M")
F.at.M <- readRDS("L-nebulosus_MF_Age_Structure") %>% 
  mutate(Mort.Group = "F at M")
F.at.15M <- readRDS("L-nebulosus_1.5F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 1.5M")

Full.Lengths <- rbind(F.at.0, F.at.02M ,F.at.M23,F.at.M, F.at.15M) %>%
  mutate(Mort.Group = fct_relevel(Mort.Group, c("No fishing mortality", "F at 0.2M" ,"F at 2/3M", "F at M",  "F at 1.5M")))

#* Work out the proportions in the different length groups ####
LM <- 350 #Linf*(1-exp(-k*(M50-t0))) # Age at 50% Mature is 3.62
Linf.LM <- Linf - ((Linf-LM)/2)
# Linf.LM <- Linf - ((Linf-LM)*(2/3))


Length.Props <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  mutate(Length.Group = ifelse(Length < LM, "Small",
                               ifelse(Length>Linf.LM, "Biggest", "Big"))) %>%  
  filter(!Length.Group %in% "Small") %>% 
  group_by(Mort.Group, Length.Group) %>% 
  summarise(Total = sum(Number)) %>% 
  mutate(Proportion = Total/sum(Total)) %>% 
  ungroup() 

G50 <- Linf*(1-exp(-k*(2.88811-t0))) # Worked out the age at half max growth rate using Wolfram Alpha

Growth.Props <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  mutate(Length.Group = ifelse(Length < LM, "Small",
                               ifelse(Length>Linf.LM, "Biggest", "Big"))) %>% 
  filter(!Length.Group %in% "Small") %>% 
  group_by(Mort.Group, Length.Group) %>% 
  summarise(Total = sum(Number)) %>% 
  mutate(Proportion = Total/sum(Total)) %>% 
  ungroup() 

length.plot.nebulosus <- Full.Lengths %>% 
  ggplot()+
  geom_line(data=Full.Lengths, aes(x=Length, y=Number, group=Mort.Group, colour=Mort.Group), linewidth=1)+
  scale_colour_carto_d(name = "Fishing Mortality", palette = "Safe", direction=1)+
  ylab("Abundance")+
  geom_vline(xintercept = 410.4687, linetype = "dashed")+
  geom_vline(xintercept = 494.9792, linetype = "dotted")+
  geom_vline(xintercept = 537.2344, linetype= "dotdash")+
  guides(color = guide_legend(reverse=FALSE))+
  theme_classic()+
  Theme1
length.plot.nebulosus 

#### L SEBAE ####
#* Von Bertalanffy Parameters ####
Linf = 623.9 # L sebae
k = 0.2
t0 = - 1.533

# Maturity
M50 = 3.43 # L sebae
M95 = 4.37 # L sebae

#* Read in data ####
setwd(pop_dir)
F.at.0 <- readRDS("L-sebae_0F_Age_Structure") %>% 
  mutate(Mort.Group = "No fishing mortality")
F.at.M02 <- readRDS("L-sebae_0.2F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 0.2M")
F.at.M23 <- readRDS("L-sebae_2.3F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 2/3M")
F.at.M08 <- readRDS("L-sebae_0.8F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 0.8M")
F.at.M <- readRDS("L-sebae_MF_Age_Structure") %>% 
  mutate(Mort.Group = "F at M")
F.at.15M <- readRDS("L-sebae_1.5F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 1.5M")

Full.Lengths <- rbind(F.at.0, F.at.M02 ,F.at.M23, F.at.M08, F.at.M, F.at.15M) %>% 
  mutate(Mort.Group = fct_relevel(Mort.Group, levels=c("No fishing mortality","F at 0.2M", "F at 2/3M", "F at M","F at 0.8M", "F at 1.5M")))

#* Work out the proportions in the different length groups ####
LM <- 480 #Linf*(1-exp(-k*(M50-t0))) # Age at 50% Mature is 3.62
Linf.LM <- Linf - ((Linf-LM)/2)
Linf.LM <- Linf - ((Linf-LM)*(2/3))


Length.Props <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  mutate(Length.Group = ifelse(Length < LM, "Small",
                               ifelse(Length>Linf.LM, "Biggest", "Big"))) %>%  
  filter(!Length.Group %in% "Small") %>% 
  group_by(Mort.Group, Length.Group) %>% 
  summarise(Total = sum(Number)) %>% 
  mutate(Proportion = Total/sum(Total)) %>% 
  ungroup() 

G50 <- Linf*(1-exp(-k*(3.465736-t0))) # Worked out the age at half max growth rate using Wolfram Alpha

Growth.Props <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  mutate(Length.Group = ifelse(Length < LM, "Small",
                               ifelse(Length>Linf.LM, "Biggest", "Big"))) %>% 
  filter(!Length.Group %in% "Small") %>% 
  group_by(Mort.Group, Length.Group) %>% 
  summarise(Total = sum(Number)) %>% 
  mutate(Proportion = Total/sum(Total)) %>% 
  ungroup() 

length.plot.sebae <- Full.Lengths %>% 
  ggplot()+
  geom_line(data=Full.Lengths, aes(x=Length, y=Number, group=Mort.Group, colour=Mort.Group), linewidth=1)+
  scale_colour_carto_d(name = "Fishing Mortality", palette = "Safe", direction=1)+
  ylab("Abundance")+
  geom_vline(xintercept = LM, linetype = "dashed")+
  geom_vline(xintercept = 469.7502, linetype = "dotted")+
  geom_vline(xintercept = 508.2876, linetype = "dotdash")+
  guides(color = guide_legend(reverse=FALSE))+
  theme_classic()+
  Theme1
length.plot.sebae


#### L MINIATUS ####
#* Von Bertalanffy Parameters ####
Linf = 945 # L miniatus
k = 0.166
t0 = -1.75

# Maturity
M50 = 1.3 # L miniatus
M95 = 1.85 # L miniatus

#* Read in data ####
setwd(pop_dir)
F.at.0 <- readRDS("L-miniatus_0F_Age_Structure") %>% 
  mutate(Mort.Group = "No fishing mortality")
F.at.08M <- readRDS("L-miniatus_0.8F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 0.8M")
F.at.M23 <- readRDS("L-miniatus_2.3F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 2/3M")
F.at.M <- readRDS("L-miniatus_MF_Age_Structure") %>% 
  mutate(Mort.Group = "F at M")
F.at.15M <- readRDS("L-miniatus_1.5F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 1.5M")

Full.Lengths <- rbind(F.at.0, F.at.M23, F.at.08M, F.at.M, F.at.15M) %>% 
  mutate(Mort.Group = fct_relevel(Mort.Group, levels=c("No fishing mortality", "F at M", "F at 2/3M", "F at 0.8M", "F at 1.5M")))

#* Work out the proportions in the different length groups ####
LM <- Linf*(1-exp(-k*(M50-t0))) # Age at 50% Mature is 1.2
# Linf.LM <- Linf - ((Linf-LM)/2)
Linf.LM <- Linf - ((Linf-LM)*(2/3))


Length.Props <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  mutate(Length.Group = ifelse(Length < LM, "Small",
                               ifelse(Length>Linf.LM, "Biggest", "Big"))) %>%  
  filter(!Length.Group %in% "Small") %>% 
  group_by(Mort.Group, Length.Group) %>% 
  summarise(Total = sum(Number)) %>% 
  mutate(Proportion = Total/sum(Total)) %>% 
  ungroup() 

G50 <- Linf*(1-exp(-k*(1.92541-t0))) # Worked out the age at half max growth rate using Wolfram Alpha

Growth.Props <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  mutate(Length.Group = ifelse(Length < LM, "Small",
                               ifelse(Length>G50, "Biggest", "Big"))) %>% 
  filter(!Length.Group %in% "Small") %>% 
  group_by(Mort.Group, Length.Group) %>% 
  summarise(Total = sum(Number)) %>% 
  mutate(Proportion = Total/sum(Total)) %>% 
  ungroup() 

length.plot.miniatus <- Full.Lengths %>% 
  ggplot()+
  geom_line(data=Full.Lengths, aes(x=Length, y=Number, group=Mort.Group, colour=Mort.Group), linewidth=1)+
  scale_colour_carto_d(name = "Fishing Mortality", palette = "Safe", direction=1)+
  ylab("Abundance")+
  geom_vline(xintercept = LM, linetype = "dashed")+
  geom_vline(xintercept = 565.2854, linetype = "dotted")+
  geom_vline(xintercept = 660.2141, linetype = "dotdash")+
  guides(color = guide_legend(reverse=FALSE))+
  theme_classic()+
  Theme1
length.plot.miniatus

#### E MULTINOTATUS ####
#* Von Bertalanffy Parameters ####
Linf = 666.4 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
k = 0.221 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
t0 = -1.835 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

# Maturity
M50 = 2 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
M95 = 2.2849 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

#* Read in data ####
setwd(pop_dir)
F.at.0 <- readRDS("E-multinotatus_0F_Age_Structure") %>% 
  mutate(Mort.Group = "No fishing mortality")
F.at.02M <- readRDS("E-multinotatus_0.2F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 0.2M")
F.at.M23 <- readRDS("E-multinotatus_2.3F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 2/3M")
F.at.M <- readRDS("E-multinotatus_MF_Age_Structure") %>% 
  mutate(Mort.Group = "F at M")
F.at.15M <- readRDS("E-multinotatus_1.5F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 1.5M")

Full.Lengths <- rbind(F.at.0, F.at.02M, F.at.M23, F.at.M, F.at.15M) %>%
  mutate(Mort.Group = fct_relevel(Mort.Group, levels=c("No fishing mortality", "F at M", "F at 2/3M", "F at 0.8M", "F at 1.5M")))

#* Work out the proportions in the different length groups ####
LM <- Linf*(1-exp(-k*(M50-t0))) # Age at 50% Mature is 3.62
# Linf.LM <- Linf - ((Linf-LM)/2)
Linf.LM <- Linf - ((Linf-LM)*(2/3))


Length.Props <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  mutate(Length.Group = ifelse(Length < LM, "Small",
                               ifelse(Length>Linf.LM, "Biggest", "Big"))) %>%  
  filter(!Length.Group %in% "Small") %>% 
  group_by(Mort.Group, Length.Group) %>% 
  summarise(Total = sum(Number)) %>% 
  mutate(Proportion = Total/sum(Total)) %>% 
  ungroup() 

G50 <- Linf*(1-exp(-k*(2.88811-t0))) # Worked out the age at half max growth rate using Wolfram Alpha

Growth.Props <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  mutate(Length.Group = ifelse(Length < LM, "Small",
                               ifelse(Length>Linf.LM, "Biggest", "Big"))) %>% 
  filter(!Length.Group %in% "Small") %>% 
  group_by(Mort.Group, Length.Group) %>% 
  summarise(Total = sum(Number)) %>% 
  mutate(Proportion = Total/sum(Total)) %>% 
  ungroup() 

#### P MULTIDENS ####
#* Von Bertalanffy Parameters ####
Linf = 590 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
k = 0.26 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
t0 = -0.46 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

# Maturity
M50 = 4.8 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
M95 = 11.95 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

#* Read in data ####
setwd(pop_dir)
F.at.0 <- readRDS("P-multidens_0F_Age_Structure") %>% 
  mutate(Mort.Group = "No fishing mortality")
F.at.02M <- readRDS("P-multidens_0.2F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 0.2M")
F.at.M23 <- readRDS("P-multidens_2.3F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 2/3M")
F.at.M <- readRDS("P-multidens_MF_Age_Structure") %>% 
  mutate(Mort.Group = "F at M")
F.at.15M <- readRDS("P-multidens_1.5F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 1.5M")

Full.Lengths <- rbind(F.at.0, F.at.02M, F.at.M23, F.at.M, F.at.15M) %>%
  mutate(Mort.Group = fct_relevel(Mort.Group, levels=c("No fishing mortality", "F at M", "F at 2/3M", "F at 0.8M", "F at 1.5M")))

#* Work out the proportions in the different length groups ####
LM <- Linf*(1-exp(-k*(M50-t0))) # Age at 50% Mature is 3.62
# Linf.LM <- Linf - ((Linf-LM)/2)
Linf.LM <- Linf - ((Linf-LM)*(2/3))


Length.Props <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  mutate(Length.Group = ifelse(Length < LM, "Small",
                               ifelse(Length>Linf.LM, "Biggest", "Big"))) %>%  
  filter(!Length.Group %in% "Small") %>% 
  group_by(Mort.Group, Length.Group) %>% 
  summarise(Total = sum(Number)) %>% 
  mutate(Proportion = Total/sum(Total)) %>% 
  ungroup() 
