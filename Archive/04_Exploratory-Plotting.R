###################################################

# Creating different plots to represent the 
# populations

###################################################
library(tidyverse)
library(plyr)
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


#### L NEBULOSUS ####
#* Von Bertalanffy Parameters ####
Linf = 664 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
k = 0.241 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
t0 = -0.375 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr

# Maturity
M50 = 3.62 # From Marriott et al. 2011
M95 = 5.97 # From Marriott et al. 2011

# Selectivity
A50 = 3
A95 = 4

# Length-weight relationship
WLa = 0.000028
WLb = 2.8761

# Age
MaxAge = 30

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

Full.Lengths <- rbind(F.at.0, F.at.02M ,F.at.M23, F.at.M, F.at.15M) %>%
  mutate(Mort.Group = fct_relevel(Mort.Group, c("No fishing mortality", "F at 0.2M" ,"F at 2/3M", "F at M", "F at 1.5M")))

Only.Lengths <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup()

Selectivity <- as.data.frame(array(1, dim=c((MaxAge*12)+1, 2))) %>% 
  rename(Age = "V1") %>% 
  rename(selectivity = "V2") %>% 
  mutate(Age = F.at.0$Age) %>% 
  mutate(selectivity = 1/(1+(exp(-log(19)*((Age-A50)/(A95-A50))))))

# Full.Lengths <- Full.Lengths %>% 
#   mutate(Pop.Sample = Number * Selectivity$selectivity)

Mort.Levels <-  c("No fishing mortality", "F at 0.2M" ,"F at 2/3M", "F at M", "F at 1.5M")

Pop.Dist <- list()

for(M in 1:5){
  
  length <- Only.Lengths %>% 
    filter(Mort.Group %in% Mort.Levels[M])
  
  fished <- Full.Lengths %>% 
    filter(Mort.Group %in% Mort.Levels[M]) %>% 
    mutate(Pop.Sample = Number * Selectivity$selectivity) %>% 
    group_by(Age = base::ceiling(row_number()/12)) %>%
    dplyr::summarise(Total.Sample = sum(Pop.Sample)) %>% 
    mutate(length = Linf*(1-exp(-k*(Age-t0)))) %>% 
    mutate(length.group = round(length/5)*5) %>% 
    ungroup() %>%
    filter(length > 480) %>% 
    group_by(length) %>%
    dplyr::summarise(Total = sum(Total.Sample)) %>%
    mutate(Proportion = Total/sum(Total)) %>%
    mutate(Percentage = Proportion *100) %>%
    ungroup() %>% 
    mutate(Mort.Group = Mort.Levels[M])
  
  Pop.Dist[[M]] <- fished
}


F.at.0 <- ggplot()+
  geom_col(data=Pop.Dist[[1]], aes(y=Percentage, x=Age))

F.at.0.2M <- ggplot()+
  geom_col(data=Pop.Dist[[2]], aes(y=Percentage, x=Age))

F.at.2.3M <- ggplot()+
  geom_col(data=Pop.Dist[[3]], aes(y=Percentage, x=Age))

F.at.M <- ggplot()+
  geom_col(data=Pop.Dist[[4]], aes(y=Percentage, x=Age))

F.at.1.5M <- ggplot()+
  geom_col(data=Pop.Dist[[5]], aes(y=Percentage, x=Age))

## Adding in real data
setwd(data_dir)
Shallow2015 <- read.csv("2015-08_Ningaloo_stereoBRUVs_length.csv") %>% 
  filter(!campaignid %in% "2015-08_Ningaloo.deep.sanctuaries_stereoBRUVs") %>% 
  mutate(scientific = paste0(genus, sep=" ", species)) %>% 
  dplyr::select(sample, length, status, campaignid)
Deep2019 <- read.csv("ningaloo.complete.length.csv") %>% 
  mutate(scientific = paste0(genus, sep=" ", species)) %>% 
  filter(scientific %in% "Lethrinus nebulosus") %>% 
  filter(!is.na(length)) 
Deep2019Zones <- read.csv("final.data.csv") %>% 
  dplyr::select(sample, status)

Deep2019 <- Deep2019 %>% 
  left_join(., Deep2019Zones, by="sample") %>% 
  dplyr::select(sample, length, status.y, campaignid) %>% 
  rename(status="status.y") %>% 
  mutate(status = ifelse(status %in% "NT", "No-take", "Fished"))

# Deep2021 <- readRDS("Ningaloo_PtCloates_BOSS-BRUV_dat_length.rds") %>% 
#   filter(scientific %in% "Lethrinidae Lethrinus nebulosus") %>% 
#   filter(!is.na(length)) %>% 
#   dplyr::select(sample, length, status, campaignid) %>% 
#   mutate(status = ifelse(status %in% "No-Take", "No-take", status))


Full.Lengths.Real.NTZ <- rbind(Shallow2015, Deep2019) %>% 
  mutate(length.group = round(length/10)*10) %>% 
  filter(status %in% c("No-take")) %>% 
  group_by(length.group) %>% 
  dplyr::summarise(n.group = n()) %>%
  ungroup() %>% 
  mutate(Total = sum(n.group)) %>% 
  group_by(length.group) %>%
  # dplyr::summarise(Total = sum(n.group)) %>%
  mutate(Proportion = n.group/Total) %>%
  mutate(Percentage = Proportion *100)
  

combined.plot.nebulosus.ntz <- ggplot()+
  geom_area(data=Pop.Dist[[1]], aes(x=length, y=Percentage), fill="darkseagreen1", colour="darkseagreen1", alpha=0.5)+ # F = 0
  geom_area(data=Pop.Dist[[2]], aes(x=length, y=Percentage), fill="skyblue2", colour="skyblue2", alpha=0.5)+ # F = 0.2M
  geom_area(data=Pop.Dist[[3]], aes(x=length, y=Percentage), fill="skyblue2", colour="skyblue2", alpha=0.5)+ # F = 2/3M
  geom_area(data=Pop.Dist[[4]], aes(x=length, y=Percentage), fill="moccasin", colour="moccasin", alpha=0.5)+ # F = M
  geom_area(data=Pop.Dist[[5]], aes(x=length, y=Percentage), fill="salmon", colour="salmon", alpha=0.4)+ # F = 1.5M
  geom_col(data=Full.Lengths.Real.NTZ, aes(x=length.group, y=Percentage), fill="grey30", colour="grey20", width=10)+
  ylab("Percentage Frequency")+
  xlab("Length (mm)")+
  theme_classic()
combined.plot.nebulosus.ntz 

Full.Lengths.Real.F <- rbind(Shallow2015, Deep2019) %>% 
  mutate(length.group = round(length/10)*10) %>% 
  filter(status %in% c("Fished")) %>% 
  group_by(length.group) %>% 
  dplyr::summarise(n.group = n()) %>%
  ungroup() %>% 
  mutate(Total = sum(n.group)) %>% 
  group_by(length.group) %>%
  # dplyr::summarise(Total = sum(n.group)) %>%
  mutate(Proportion = n.group/Total) %>%
  mutate(Percentage = Proportion *100)


combined.plot.nebulosus.f <- ggplot()+
  geom_area(data=Pop.Dist[[1]], aes(x=length, y=Percentage), fill="seagreen4", colour="seagreen4", alpha=0.5)+ # F = 0
  geom_area(data=Pop.Dist[[2]], aes(x=length, y=Percentage), fill="darkseagreen1", colour="darkseagreen1", alpha=0.5)+ # F = 0.2M
  geom_area(data=Pop.Dist[[3]], aes(x=length, y=Percentage), fill="skyblue2", colour="skyblue2", alpha=0.5)+ # F = 2/3M
  geom_area(data=Pop.Dist[[4]], aes(x=length, y=Percentage), fill="moccasin", colour="moccasin", alpha=0.5)+ # F = M
  geom_area(data=Pop.Dist[[5]], aes(x=length, y=Percentage), fill="salmon", colour="salmon", alpha=0.4)+ # F = 1.5M
  geom_col(data=Full.Lengths.Real.F, aes(x=length.group, y=Percentage), fill="grey30", colour="grey20", width=10)+
  ylab("Percentage Frequency")+
  xlab("Length (mm)")+
  theme_classic()
combined.plot.nebulosus.f

#* Condition plot for L nebulosus ####
# Create values for condition from spreadsheet 
fish.condition <- as.data.frame(array(1, dim=c(6, 3))) %>% 
  rename(Species = "V1") %>% 
  rename(Percent = "V2") %>% 
  rename(Zone = "V3") %>% 
  mutate(Species = c("L. nebulosus", "L. nebulosus", "L. sebae","L. sebae", "L. miniatus", "L. miniatus")) %>% 
  mutate(Percent = c(24.1, 19.8, 71.6, 63.4, 3.06, NA)) %>% 
  mutate(Zone = c("NPZ", "MUZ", "NPZ", "MUZ", "NPZ", "MUZ")) %>% 
  mutate(Zone = as.factor(Zone)) %>% 
  mutate(Zone = fct_relevel(Zone, "NPZ", "HPZ", "MUZ", "SPZ"))

plot.df <-   as.data.frame(array(0, dim=c(4, 2))) %>% 
  rename(Percentage = "V1") %>% 
  rename(Zone = "V2") %>% 
  mutate(Zone = c("NPZ", "HPZ", "MUZ", "SPZ")) %>% 
  mutate(Zone = as.factor(Zone)) %>% 
  mutate(Zone = fct_relevel(Zone, "NPZ", "HPZ", "MUZ", "SPZ"))

nebs.condition <- fish.condition %>% 
  filter(Species %in% "L. nebulosus")

nebulosus.condition <- 
  ggplot()+
  geom_point(data=plot.df, aes(y=Percentage, x=Zone, group=Zone), size=0, colour=NA)+
  geom_point(data = nebs.condition, aes(y=Percent, x=Zone, group=Zone), colour=NA) + 
  annotate("rect", xmin = 0, xmax = 4.5, ymin = 0, ymax = 33,fill = "salmon", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 4.5, ymin = 33, ymax = 42,fill = "moccasin", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 4.5, ymin = 42, ymax = 48, fill = "skyblue2", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 48, ymax = 54,fill = "skyblue2", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 54, ymax = 100, fill = "darkseagreen1", alpha=0.5)+
  annotate("rect", xmin = 1.5, xmax = 4.5, ymin = 48, ymax = 100,fill = "darkseagreen1", alpha=0.5)+
  geom_point(data = nebs.condition, aes(y=Percent, x=Zone, group=Zone), size=4, colour="grey30" ) + 
  geom_vline(xintercept = 1.5, linetype="dashed")+
  geom_vline(xintercept = 2.5, linetype="dashed")+
  geom_vline(xintercept = 3.5, linetype="dashed")+
  ylim(0,100)+
  ylab("% of large size class")+
  theme_classic()
nebulosus.condition

#* L. nebulosus biomass ####

Full.Lengths.Mat <- Full.Lengths %>% 
  mutate(Mature = 1/(1+(exp(-log(19)*((Age-M50)/(M95-M50)))))) %>% 
  mutate(No.Mat = Mature*Number) %>% 
  mutate(Weight = (WLa*(Length^WLb)/1000)) %>% 
  mutate(Mat.Bio = No.Mat*Weight)

Prop.Mat <- Full.Lengths.Mat %>% 
  group_by(Mort.Group) %>% 
  #slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  group_by(Mort.Group) %>% 
  summarise(Total = sum(Mat.Bio)) %>% 
  mutate(Proportion = Total/Total[Mort.Group=="No fishing mortality"]) %>% 
  ungroup() 

# Working out proportion of mature from real data

Full.Mature.Real <- rbind(Shallow2015, Deep2019) %>% 
  mutate(length.group = round(length/10)*10) %>% 
  mutate(Weight = (WLa*(length^WLb)/1000)) %>% 
  mutate(Age = ((-log(1 - length/Linf))/k + t0)) %>% 
  mutate(Age = ifelse(is.nan(Age), 30, Age)) %>% 
  mutate(Prop.Mat = 1/(1+(exp(-log(19)*((Age-M50)/(M95-M50)))))) %>% 
  mutate(Mat.Bio = Prop.Mat*Weight) %>% 
  mutate(Year = ifelse(campaignid %in% c("2015-08_Ningaloo.lagoon.sanctuaries_stereoBRUVs", "2015-08_Ningaloo.shallow.sanctuaries_stereoBRUVs"), 2015, 2019)) %>% 
  group_by(Year) %>% 
  summarise(Prop.Bio = (sum(Mat.Bio)/56701.955)*100)   # 0.002688516


plot.df <-   as.data.frame(array(0, dim=c(5, 2))) %>% 
  rename(Percentage = "V1") %>% 
  rename(Year = "V2") %>% 
  mutate(Year = c(2015, 2016, 2017, 2018, 2019))

nebulosus.biomass <- 
  ggplot()+
  geom_point(data=plot.df, aes(y=Percentage, x=Year, group=Year), size=0, colour=NA)+
  #geom_point(data = Full.Mature.Real, aes(y=Prop.Bio, x=Year, group=Year), colour=NA) + 
  annotate("rect", xmin = 2014.5, xmax = 2019.5, ymin = 0, ymax = 17,fill = "salmon", alpha=0.5)+ #1.5M
  annotate("rect", xmin = 2014.5, xmax = 2019.5, ymin = 17, ymax = 28,fill = "moccasin", alpha=0.5)+ #2/3 M
  annotate("rect", xmin = 2014.5, xmax = 2019.5, ymin = 28, ymax = 42, fill = "skyblue2", alpha=0.5)+ #M = F
  annotate("rect", xmin = 2014.5, xmax = 2019.5, ymin = 42, ymax = 75,fill = "skyblue2", alpha=0.5)+ #0.2M
  annotate("rect", xmin = 2014.5, xmax = 2019.5, ymin = 75, ymax = 100, fill = "darkseagreen1", alpha=0.5)+
  geom_point(data = Full.Mature.Real, aes(y=Prop.Bio, x=Year),  size=4, colour="grey30") +
  # geom_vline(xintercept = 1.5, linetype="dashed")+
  # geom_vline(xintercept = 2.5, linetype="dashed")+
  # geom_vline(xintercept = 3.5, linetype="dashed")+
  ylim(0,100)+
  ylab("SB/Ksb")+
  theme_classic()
nebulosus.biomass



#### L SEBAE ####
#* Von Bertalanffy Parameters ####
Linf = 623.9 # L sebae
k = 0.2
t0 = - 1.533

# Maturity
M50 = 3.43 # L sebae
M95 = 4.37 # L sebae

# Selectivity
A50 = 3.762
A95 = 6.241

# Weight length
WLa = 1.406e-08
WLb = 3.051

# Longevity
MaxAge = 32

#* Read in data ####
setwd(pop_dir)
F.at.0 <- readRDS("L-sebae_0F_Age_Structure") %>% 
  mutate(Mort.Group = "No fishing mortality")
F.at.02M <- readRDS("L-sebae_0.2F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 0.2M")
F.at.M23 <- readRDS("L-sebae_2.3F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 2/3M")
F.at.M <- readRDS("L-sebae_MF_Age_Structure") %>% 
  mutate(Mort.Group = "F at M")
F.at.15M <- readRDS("L-sebae_1.5F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 1.5M")


Full.Lengths <- rbind(F.at.0, F.at.02M ,F.at.M23, F.at.M, F.at.15M) %>%
  dplyr::mutate(Mort.Group = fct_relevel(Mort.Group, c("No fishing mortality", "F at 0.2M" ,"F at 2/3M", "F at M", "F at 1.5M")))

Only.Lengths <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup()

Selectivity <- as.data.frame(array(1, dim=c((MaxAge*12)+1, 2))) %>% 
  dplyr::rename(Age = "V1") %>% 
  dplyr::rename(selectivity = "V2") %>% 
  mutate(Age = F.at.0$Age) %>% 
  mutate(selectivity = 1/(1+(exp(-log(19)*((Age-A50)/(A95-A50))))))

# Full.Lengths <- Full.Lengths %>%
#   mutate(Pop.Sample = Number * Selectivity$selectivity)

Mort.Levels <-  c("No fishing mortality", "F at 0.2M" ,"F at 2/3M", "F at M", "F at 1.5M")

Pop.Dist <- list()

for(M in 1:5){
  
  length <- Only.Lengths %>% 
    dplyr::filter(Mort.Group %in% Mort.Levels[M])
  
  fished <- Full.Lengths %>% 
    dplyr::filter(Mort.Group %in% Mort.Levels[M]) %>% 
    mutate(Pop.Sample = Number * Selectivity$selectivity) %>% 
    dplyr::group_by(Age = base::ceiling(row_number()/12)) %>%
    dplyr::summarise(Total.Sample = sum(Pop.Sample)) %>% 
    mutate(length = Linf*(1-exp(-k*(Age-t0)))) %>% 
    mutate(length.group = round(length/5)*5) %>% 
    ungroup() %>%
    filter(length>480) %>% 
    group_by(length) %>%
    dplyr::summarise(Total = sum(Total.Sample)) %>%
    mutate(Proportion = Total/sum(Total)) %>%
    mutate(Percentage = Proportion *100) %>%
    ungroup() %>% 
    mutate(Mort.Group = Mort.Levels[M])
  
  Pop.Dist[[M]] <- fished
}


F.at.0 <- ggplot()+
  geom_col(data=Pop.Dist[[1]], aes(y=Percentage, x=Age))

F.at.0.2M <- ggplot()+
  geom_col(data=Pop.Dist[[2]], aes(y=Percentage, x=Age))

F.at.2.3M <- ggplot()+
  geom_col(data=Pop.Dist[[3]], aes(y=Percentage, x=Age))

F.at.M <- ggplot()+
  geom_col(data=Pop.Dist[[4]], aes(y=Percentage, x=Age))

F.at.1.5M <- ggplot()+
  geom_col(data=Pop.Dist[[5]], aes(y=Percentage, x=Age))

## Adding in real data
setwd(data_dir)
Pilbara <- readRDS("Pilbara_data") %>% 
  #filter(Length < 2000) %>% 
  filter(scientific %in% "Lutjanus sebae") %>% 
  dplyr::select(Sample, Length, Status)

Full.Lengths.Real.ntz <- Pilbara %>% 
  filter(Status %in% "No-take") %>% 
  mutate(length.group = round(Length/15)*15) %>% 
  group_by(length.group) %>% 
  dplyr::summarise(n.group = n()) %>%
  ungroup() %>% 
  mutate(Total = sum(n.group)) %>% 
  group_by(length.group) %>%
  # dplyr::summarise(Total = sum(n.group)) %>%
  mutate(Proportion = n.group/Total) %>%
  mutate(Percentage = Proportion *100)


combined.plot.sebae <- ggplot()+
  geom_area(data=Pop.Dist[[1]], aes(x=length, y=Percentage), fill="darkseagreen1", colour="darkseagreen1", alpha=0.5)+ # F = 0
  geom_area(data=Pop.Dist[[2]], aes(x=length, y=Percentage), fill="skyblue2", colour="skyblue2", alpha=0.5)+ # F = 0.2M
  geom_area(data=Pop.Dist[[3]], aes(x=length, y=Percentage), fill="skyblue2", colour="skyblue2", alpha=0.5)+ # F = 2/3M
  geom_area(data=Pop.Dist[[4]], aes(x=length, y=Percentage), fill="moccasin", colour="moccasin", alpha=0.5)+ # F = M
  geom_area(data=Pop.Dist[[5]], aes(x=length, y=Percentage), fill="salmon", colour="salmon", alpha=0.4)+ # F = 1.5M
  geom_col(data=Full.Lengths.Real.ntz, aes(x=length.group, y=Percentage), fill="grey30", colour="grey20", width=10)+
  ylab("Percentage Frequency")+
  xlab("Length (mm)")+
  theme_classic()
combined.plot.sebae

Full.Lengths.Real.f <- Pilbara %>% 
  filter(Status %in% "Fished") %>% 
  mutate(length.group = round(Length/15)*15) %>% 
  group_by(length.group) %>% 
  dplyr::summarise(n.group = n()) %>%
  ungroup() %>% 
  mutate(Total = sum(n.group)) %>% 
  group_by(length.group) %>%
  # dplyr::summarise(Total = sum(n.group)) %>%
  mutate(Proportion = n.group/Total) %>%
  mutate(Percentage = Proportion *100)


combined.plot.sebae.f <- ggplot()+
  geom_area(data=Pop.Dist[[1]], aes(x=length, y=Percentage), fill="seagreen4", colour="seagreen4", alpha=0.5)+ # F = 0
  geom_area(data=Pop.Dist[[2]], aes(x=length, y=Percentage), fill="darkseagreen1", colour="darkseagreen1", alpha=0.5)+ # F = 0.2M
  geom_area(data=Pop.Dist[[3]], aes(x=length, y=Percentage), fill="skyblue2", colour="skyblue2", alpha=0.5)+ # F = 2/3M
  geom_area(data=Pop.Dist[[4]], aes(x=length, y=Percentage), fill="moccasin", colour="moccasin", alpha=0.5)+ # F = M
  geom_area(data=Pop.Dist[[5]], aes(x=length, y=Percentage), fill="salmon", colour="salmon", alpha=0.4)+ # F = 1.5M
  geom_col(data=Full.Lengths.Real.f, aes(x=length.group, y=Percentage), fill="grey30", colour="grey20", width=10)+
  ylab("Percentage Frequency")+
  xlab("Length (mm)")+
  theme_classic()
combined.plot.sebae.f

#* L sebae condition plot ####
seb.condition <- fish.condition %>%
  filter(Species %in% "L. sebae")

sebae.condition <- 
  ggplot()+
  geom_point(data=plot.df, aes(y=Percentage, x=Zone, group=Zone), size=0, colour=NA)+
  geom_point(data = seb.condition, aes(y=Percent, x=Zone, group=Zone), colour=NA) + 
  annotate("rect", xmin = 0, xmax = 4.5, ymin = 0, ymax = 33,fill = "salmon", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 4.5, ymin = 33, ymax = 42,fill = "moccasin", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 4.5, ymin = 42, ymax = 48, fill = "skyblue2", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 48, ymax = 54,fill = "skyblue2", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 54, ymax = 100, fill = "darkseagreen1", alpha=0.5)+
  annotate("rect", xmin = 1.5, xmax = 4.5, ymin = 48, ymax = 100,fill = "darkseagreen1", alpha=0.5)+
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 100,fill = "grey90")+
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 0, ymax = 100, fill = "grey90")+
  geom_point(data = seb.condition, aes(y=Percent, x=Zone, group=Zone), size=4, colour="grey30" ) + 
  geom_vline(xintercept = 1.5, linetype="dashed")+
  geom_vline(xintercept = 2.5, linetype="dashed")+
  geom_vline(xintercept = 3.5, linetype="dashed")+
  ylim(0,100)+
  ylab("% of large size class")+
  theme_classic()
sebae.condition

#* L sebae Biomass ####
Full.Lengths.Mat <- Full.Lengths %>% 
  mutate(Mature = 1/(1+(exp(-log(19)*((Age-M50)/(M95-M50)))))) %>% 
  mutate(No.Mat = Mature*Number) %>% 
  mutate(Weight = (WLa*(Length^WLb)/1000)) %>% 
  mutate(Mat.Bio = No.Mat*Weight)

Prop.Mat <- Full.Lengths.Mat %>% 
  group_by(Mort.Group) %>% 
  #slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  group_by(Mort.Group) %>% 
  summarise(Total = sum(Mat.Bio)) %>% 
  mutate(Proportion = Total/Total[Mort.Group=="No fishing mortality"]) %>% 
  ungroup() 


# Working out proportion of mature from real data
setwd(data_dir)
Pilbara <- readRDS("Pilbara_data") %>% 
  filter(Length < 2000) %>% 
  filter(scientific %in% "Lutjanus sebae") %>% 
  dplyr::select(Sample, Length, Status)


Full.Mature.Real <- Pilbara %>% 
  mutate(length.group = round(Length/10)*10) %>% 
  mutate(Weight = (WLa*(Length^WLb)/1000)) %>% 
  mutate(Age = ((-log(1 - Length/Linf))/k + t0)) %>% 
  mutate(Age = ifelse(is.nan(Age), MaxAge, Age)) %>% 
  mutate(Prop.Mat = 1/(1+(exp(-log(19)*((Age-M50)/(M95-M50)))))) %>% 
  mutate(Mat.Bio = Prop.Mat*Weight) %>% 
  mutate(Year = 2010) %>% 
  group_by(Year) %>% 
  summarise(Prop.Bio = (sum(Mat.Bio)/122.84299)*100)   

plot.df <-   as.data.frame(array(0, dim=c(5, 2))) %>% 
  rename(Percentage = "V1") %>% 
  rename(Year = "V2") %>% 
  mutate(Year = c(2008, 2009, 2010, 2011, 2012))

sebae.biomass <- 
  ggplot()+
  geom_point(data=plot.df, aes(y=Percentage, x=Year, group=Year), size=0, colour=NA)+
  #geom_point(data = Full.Mature.Real, aes(y=Prop.Bio, x=Year, group=Year), colour=NA) + 
  annotate("rect", xmin = 2008, xmax = 2012.5, ymin = 0, ymax = 24,fill = "salmon", alpha=0.5)+ #1.5M
  annotate("rect", xmin = 2008, xmax = 2012.5, ymin = 24, ymax = 36,fill = "moccasin", alpha=0.5)+ #2/3 M
  annotate("rect", xmin = 2008, xmax = 2012.5, ymin = 36, ymax = 49, fill = "skyblue2", alpha=0.5)+ #M = F
  annotate("rect", xmin = 2008, xmax = 2012.5, ymin = 49, ymax = 79,fill = "skyblue2", alpha=0.5)+ #0.2M
  annotate("rect", xmin = 2008, xmax = 2012.5, ymin = 79, ymax = 100, fill = "darkseagreen1", alpha=0.5)+
  geom_point(data = Full.Mature.Real, aes(y=Prop.Bio, x=Year),  size=4, colour="grey30") +
  # geom_vline(xintercept = 1.5, linetype="dashed")+
  # geom_vline(xintercept = 2.5, linetype="dashed")+
  # geom_vline(xintercept = 3.5, linetype="dashed")+
  ylim(0,100)+
  ylab("SB/Ksb")+
  theme_classic()
sebae.biomass


#### L MINIATUS ####
#* Von Bertalanffy Parameters ####
Linf = 945 # L miniatus
k = 0.166
t0 = -1.75

# Maturity
M50 = 1.3 # L miniatus
M95 = 1.85 # L miniatus

# Selectivity
A50 = 2.5
A95 = 4.5 

# Length-weight relationship
WLa = 8.04e-06 
WLb = 3.051

# Longevity
MaxAge = 20

#* Read in data ####
setwd(pop_dir)
F.at.0 <- readRDS("L-miniatus_0F_Age_Structure") %>% 
  mutate(Mort.Group = "No fishing mortality")
F.at.02M <- readRDS("L-miniatus_0.2F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 0.2M")
F.at.M23 <- readRDS("L-miniatus_2.3F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 2/3M")
F.at.M <- readRDS("L-miniatus_MF_Age_Structure") %>% 
  mutate(Mort.Group = "F at M")
F.at.15M <- readRDS("L-miniatus_1.5F_Age_Structure") %>% 
  mutate(Mort.Group = "F at 1.5M")

Full.Lengths <- rbind(F.at.0, F.at.02M , F.at.M23, F.at.M, F.at.15M) %>% 
  mutate(Mort.Group = fct_relevel(Mort.Group, c("No fishing mortality", "F at 0.2M", "F at M", "F at 2/3M", "F at 1.5M")))

Only.Lengths <- Full.Lengths %>% 
  group_by(Mort.Group) %>% 
  slice(which(row_number() %% 12 == 1)) %>% 
  ungroup()

Selectivity <- as.data.frame(array(1, dim=c((MaxAge*12)+1, 2))) %>% 
  dplyr::rename(Age = "V1") %>% 
  dplyr::rename(selectivity = "V2") %>% 
  mutate(Age = F.at.0$Age) %>% 
  mutate(selectivity = 1/(1+(exp(-log(19)*((Age-A50)/(A95-A50))))))

# Full.Lengths <- Full.Lengths %>%
#   mutate(Pop.Sample = Number * Selectivity$selectivity)

Mort.Levels <-  c("No fishing mortality", "F at 0.2M" ,"F at 2/3M", "F at M", "F at 1.5M")

Pop.Dist <- list()

for(M in 1:5){
  
  length <- Only.Lengths %>% 
    dplyr::filter(Mort.Group %in% Mort.Levels[M])
  
  fished <- Full.Lengths %>% 
    dplyr::filter(Mort.Group %in% Mort.Levels[M]) %>% 
    mutate(Pop.Sample = Number * Selectivity$selectivity) %>% 
    dplyr::group_by(Age = base::ceiling(row_number()/12)) %>%
    dplyr::summarise(Total.Sample = sum(Pop.Sample)) %>% 
    mutate(length = Linf*(1-exp(-k*(Age-t0)))) %>% 
    mutate(length.group = (round(length/5)*5)*10) %>% 
    ungroup() %>%
    filter(length>580) %>% 
    group_by(length) %>%
    dplyr::summarise(Total = sum(Total.Sample)) %>%
    mutate(Proportion = Total/sum(Total)) %>%
    mutate(Percentage = Proportion *100) %>%
    ungroup()
  
  Pop.Dist[[M]] <- fished
}

## Adding in real data
setwd(data_dir)

Ab2021 <- readRDS("2021-05_Abrolhos_BOSS-BRUV_dat_length.rds") %>% 
  filter(scientific %in% "Lethrinidae Lethrinus miniatus") %>% 
  filter(!is.na(length)) %>% 
  dplyr::select(sample, length, status)

Full.Lengths.Real <- Ab2021 %>% 
  mutate(length.group = round(length/15)*15) %>% 
  group_by(length.group) %>% 
  dplyr::summarise(n.group = n()) %>%
  ungroup() %>% 
  mutate(Total = sum(n.group)) %>% 
  group_by(length.group) %>%
  # dplyr::summarise(Total = sum(n.group)) %>%
  mutate(Proportion = n.group/Total) %>%
  mutate(Percentage = Proportion *100)


combined.plot.miniatus <- ggplot()+
  geom_area(data=Pop.Dist[[1]], aes(x=length, y=Percentage), fill="darkseagreen1", colour="darkseagreen1", alpha=0.5)+ # F = 0
  geom_area(data=Pop.Dist[[2]], aes(x=length, y=Percentage), fill="skyblue2", colour="skyblue2", alpha=0.5)+ # F = 0.2M
  geom_area(data=Pop.Dist[[3]], aes(x=length, y=Percentage), fill="skyblue2", colour="skyblue2", alpha=0.5)+ # F = 2/3M
  geom_area(data=Pop.Dist[[4]], aes(x=length, y=Percentage), fill="moccasin", colour="moccasin", alpha=0.5)+ # F = M
  geom_area(data=Pop.Dist[[5]], aes(x=length, y=Percentage), fill="salmon", colour="salmon", alpha=0.4)+ # F = 1.5M
  geom_col(data=Full.Lengths.Real, aes(x=length.group, y=Percentage), fill="palegreen4", colour="grey20", width=15)+
  ylab("Percentage Frequency")+
  xlab("Length (mm)")+
  theme_classic()
combined.plot.miniatus

#* L. miniatus condition plot ####

min.condition <- fish.condition %>% 
  filter(Species %in% "L. miniatus")

miniatus.condition <- 
  ggplot()+
  geom_point(data=plot.df, aes(y=Percentage, x=Zone, group=Zone), colour=NA)+
  geom_point(data = min.condition, aes(y=Percent, x=Zone, group=Zone), colour=NA) + 
  annotate("rect", xmin = 0, xmax = 4.5, ymin = 0, ymax = 33,fill = "salmon", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 4.5, ymin = 33, ymax = 42,fill = "moccasin", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 4.5, ymin = 42, ymax = 48, fill = "skyblue2", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 48, ymax = 54,fill = "skyblue2", alpha=0.5)+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 54, ymax = 100, fill = "darkseagreen1", alpha=0.5)+
  annotate("rect", xmin = 1.5, xmax = 4.5, ymin = 48, ymax = 100,fill = "darkseagreen1", alpha=0.5)+
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 100,fill = "grey90")+
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 0, ymax = 100, fill = "grey90")+
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = 0, ymax = 100, fill = "grey90")+
  geom_point(data = min.condition, aes(y=Percent, x=Zone, group=Zone), size=4, colour="grey30" ) + 
  geom_vline(xintercept = 1.5, linetype="dashed")+
  geom_vline(xintercept = 2.5, linetype="dashed")+
  geom_vline(xintercept = 3.5, linetype="dashed")+
  ylim(0,100)+
  ylab("% of large size class")+
  theme_classic()
miniatus.condition

#* L. miniatus Biomass ####

Full.Lengths.Mat <- Full.Lengths %>% 
  mutate(Length = Length*10) %>% 
  mutate(Mature = 1/(1+(exp(-log(19)*((Age-M50)/(M95-M50)))))) %>% 
  mutate(No.Mat = Mature*Number) %>% 
  mutate(Weight = (WLa*(Length^WLb)/1000)) %>% 
  mutate(Mat.Bio = No.Mat*Weight)

Prop.Mat <- Full.Lengths.Mat %>% 
  group_by(Mort.Group) %>% 
  #slice(which(row_number() %% 12 == 1)) %>% 
  ungroup() %>% 
  group_by(Mort.Group) %>% 
  summarise(Total = sum(Mat.Bio)) %>% 
  mutate(Proportion = Total/Total[Mort.Group=="No fishing mortality"]) %>% 
  ungroup() 


# Working out proportion of mature from real data
setwd(data_dir)

Ab2021 <- readRDS("2021-05_Abrolhos_BOSS-BRUV_dat_length.rds") %>% 
  filter(scientific %in% "Lethrinidae Lethrinus miniatus") %>% 
  filter(!is.na(length)) %>% 
  dplyr::select(sample, length, status)

Full.Mature.Real <- Ab2021 %>% 
  mutate(length.group = round(length/10)*10) %>% 
  mutate(Weight = (WLa*(length^WLb)/1000)) %>% 
  mutate(Age = ((-log(1 - length/Linf))/k + t0)) %>% 
  mutate(Age = ifelse(is.nan(Age), MaxAge, Age)) %>% 
  mutate(Prop.Mat = 1/(1+(exp(-log(19)*((Age-M50)/(M95-M50)))))) %>% 
  mutate(Mat.Bio = Prop.Mat*Weight) %>% 
  mutate(Year = 2021) %>% 
  group_by(Year) %>% 
  summarise(Prop.Bio = (sum(Mat.Bio)/5972.672)*100)   

plot.df <-   as.data.frame(array(0, dim=c(5, 2))) %>% 
  rename(Percentage = "V1") %>% 
  rename(Year = "V2") %>% 
  mutate(Year = c(2019, 2020, 2021, 2022, 2023))

miniatus.biomass <- 
  ggplot()+
  geom_point(data=plot.df, aes(y=Percentage, x=Year, group=Year), size=0, colour=NA)+
  #geom_point(data = Full.Mature.Real, aes(y=Prop.Bio, x=Year, group=Year), colour=NA) + 
  annotate("rect", xmin = 2019, xmax = 2023, ymin = 0, ymax = 24,fill = "salmon", alpha=0.5)+ #1.5M
  annotate("rect", xmin = 2019, xmax = 2023, ymin = 24, ymax = 36,fill = "moccasin", alpha=0.5)+ #2/3 M
  annotate("rect", xmin = 2019, xmax = 2023, ymin = 36, ymax = 49, fill = "skyblue2", alpha=0.5)+ #M = F
  annotate("rect", xmin = 2019, xmax = 2023, ymin = 49, ymax = 79,fill = "skyblue2", alpha=0.5)+ #0.2M
  annotate("rect", xmin = 2019, xmax = 2023, ymin = 79, ymax = 100, fill = "darkseagreen1", alpha=0.5)+
  geom_point(data = Full.Mature.Real, aes(y=Prop.Bio, x=Year),  size=4, colour="grey30") +
  # geom_vline(xintercept = 1.5, linetype="dashed")+
  # geom_vline(xintercept = 2.5, linetype="dashed")+
  # geom_vline(xintercept = 3.5, linetype="dashed")+
  ylim(0,100)+
  ylab("SB/Ksb")+
  theme_classic()
miniatus.biomass



