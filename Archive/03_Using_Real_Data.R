###################################################

# Calculating proportions of different length
# classes from real data

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

#### L NEBULOSUS ####
#* Von Bertalanffy Parameters ####
Linf = 664 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
k = 0.241 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
t0 = -0.375 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr

# Maturity
M50 = 3.62 # From Marriott et al. 2011
M95 = 5.97 # From Marriott et al. 2011

#* Read in Data ####
setwd(data_dir)
# Shallow2015 <- read.csv("2015-08_Ningaloo_stereoBRUVs_length.csv") %>% 
#   mutate(scientific = paste0(genus, sep=" ", species)) %>% 
#   dplyr::select(sample, length, status, campaignid)
Deep2019 <- read.csv("ningaloo.complete.length.csv") %>% 
  mutate(scientific = paste0(genus, sep=" ", species)) %>% 
  filter(scientific %in% "Lethrinus nebulosus") %>% 
  filter(!is.na(length)) %>% 
  dplyr::select(campaignid, sample, scientific, length)


Deep2021 <- readRDS("Ningaloo_PtCloates_BOSS-BRUV_dat_length.rds") %>% 
  filter(scientific %in% "Lethrinidae Lethrinus nebulosus") %>% 
  filter(!is.na(length)) %>% 
  dplyr::select(sample, length, campaignid, scientific) 

DBCA_data <- readRDS("DBCA_lengths.RDS") %>% 
  filter(scientific %in% "Lethrinus nebulosus") %>% 
  filter(!campaignid %in% "2021-05_Jurien.Bay.MP.Monitoring_stereoBRUVs") %>% 
  dplyr::select(sample, length, campaignid, scientific) 

Shallow_Data <- read.csv("2015-08_Ningaloo_stereoBRUVs_length.csv")%>% 
  mutate(scientific = paste0(genus, sep=" ", species)) %>% 
  dplyr::select(sample, length, campaignid, scientific) 

Full.Lengths <- rbind(DBCA_data, Deep2019, Deep2021)

setwd(data_dir)
write.csv(Full.Lengths, "Nebulosus_lengths.csv")

histograms <- ggplot()+
  geom_histogram(data=Full.Lengths, aes(x=length))+
  facet_grid(~status)

#* Work out the proportions in the different length groups ####
LM <- 350 #Linf*(1-exp(-k*(M50-t0))) # Age at 50% Mature is 3.62
Linf.LM <- Linf - ((Linf-LM)/2)
Linf.LM <- Linf - ((Linf-LM)*(2/3))


Length.Props <- Full.Lengths %>% 
  mutate(Length.Group = ifelse(length < LM, "Small",
                               ifelse(length>Linf.LM, "Biggest", "Big"))) %>%  
  filter(!Length.Group %in% "Small") %>% 
  group_by(status, Length.Group) %>% 
  summarise(n.group = n()) %>% 
  mutate(Proportion = n.group/sum(n.group)) %>% 
  ungroup() 

#### L MINIATUS ####
#* Von Bertalanffy Parameters ####
Linf = 945 # L miniatus
k = 0.166
t0 = -1.75

# Maturity
M50 = 1.3 # L miniatus
M95 = 1.85 # L miniatus

#* Read in Data ####
setwd(data_dir)
Shallow2015 <- read.csv("2015-08_Ningaloo_stereoBRUVs_length.csv") %>%
  mutate(scientific = paste0(genus, sep=" ", species)) %>%
  dplyr::select(sample, length, status, campaignid)
Deep2019 <- read.csv("ningaloo.complete.length.csv") %>% 
  mutate(scientific = paste0(genus, sep=" ", species)) %>% 
  filter(scientific %in% "Lethrinus miniatus") %>% 
  filter(!is.na(length)) 
Deep2019Zones <- read.csv("final.data.csv") %>% 
  dplyr::select(sample, status)

Deep2019 <- Deep2019 %>% 
  left_join(., Deep2019Zones, by="sample") %>% 
  dplyr::select(sample, length, status.y, campaignid) %>% 
  rename(status="status.y") %>% 
  mutate(status = ifelse(status %in% "NT", "No-take", "Fished"))

Deep2021 <- readRDS("Ningaloo_PtCloates_BOSS-BRUV_dat_length.rds") %>% 
  filter(scientific %in% "Lethrinidae Lethrinus miniatus") %>% 
  filter(!is.na(length)) %>% 
  dplyr::select(sample, length, status, campaignid) %>% 
  mutate(status = ifelse(status %in% "No-Take", "No-take", status))

Full.Lengths <- rbind(Deep2019, Deep2021)

#* Work out the proportions in the different length groups ####
LM <- Linf*(1-exp(-k*(M50-t0))) # Age at 50% Mature is ~1.3
# Linf.LM <- Linf - ((Linf-LM)/2)
Linf.LM <- Linf - ((Linf-LM)*(2/3))

Length.Props <- Full.Lengths %>% 
  mutate(Length.Group = ifelse(length < LM, "Small",
                               ifelse(length>Linf.LM, "Biggest", "Big"))) %>%  
  filter(!Length.Group %in% "Small") %>% 
  group_by(status, Length.Group) %>% 
  summarise(n.group = n()) %>% 
  mutate(Proportion = n.group/sum(n.group)) %>% 
  ungroup() 

#### L SEBAE ####
#* Von Bertalanffy Parameters ####
Linf = 623.9 # L sebae
k = 0.2
t0 = - 1.533

# Maturity
M50 = 3.43 # L sebae
M95 = 4.37 # L sebae

#* Read in Data ####
setwd(data_dir)
Shallow2015 <- read.csv("2015-08_Ningaloo_stereoBRUVs_length.csv") %>%
  mutate(scientific = paste0(genus, sep=" ", species)) %>%
  dplyr::select(sample, length, status, campaignid)
Deep2019 <- read.csv("ningaloo.complete.length.csv") %>% 
  mutate(scientific = paste0(genus, sep=" ", species)) %>% 
  filter(scientific %in% "Lutjanus sebae") %>% 
  filter(!is.na(length)) 
Deep2019Zones <- read.csv("final.data.csv") %>% 
  dplyr::select(sample, status)

Deep2019 <- Deep2019 %>% 
  left_join(., Deep2019Zones, by="sample") %>% 
  dplyr::select(sample, length, status.y, campaignid) %>% 
  rename(status="status.y") %>% 
  mutate(status = ifelse(status %in% "NT", "No-take", "Fished"))

Deep2021 <- read.csv("Ningaloo_PtCloates_BOSS.checked.length.csv") %>% 
  mutate(scientific = paste0(genus, sep=" ", species)) %>% 
  filter(scientific %in% "Lutjanus sebae") %>% 
  filter(!is.na(length)) %>% 
  dplyr::select(sample, length, campaignid) 

Deep2021Zones <- readRDS("Ningaloo_PtCloates_BOSS-BRUV_dat_cti.rds") %>% 
  dplyr::select(status, sample, campaignid)

Deep2021 <- Deep2021 %>% 
  left_join(., Deep2021Zones, by=c("campaignid", "sample")) %>% 
  dplyr::select(sample, length, status, campaignid) %>% 
  mutate(status = ifelse(status %in% "No-Take", "No-take", status))

Full.Lengths <- rbind(Deep2019, Deep2021)

#* Work out the proportions in the different length groups ####
LM <- Linf*(1-exp(-k*(M50-t0))) # Age at 50% Mature is 3.62
# Linf.LM <- Linf - ((Linf-LM)/2)
Linf.LM <- Linf - ((Linf-LM)*(2/3))

Length.Props <- Full.Lengths %>% 
  mutate(Length.Group = ifelse(length < LM, "Small",
                               ifelse(length>Linf.LM, "Biggest", "Big"))) %>%  
  filter(!Length.Group %in% "Small") %>% 
  group_by(status, Length.Group) %>% 
  summarise(n.group = n()) %>% 
  mutate(Proportion = n.group/sum(n.group)) %>% 
  ungroup() 

#### ABROLHOS ####
Linf = 945 # L miniatus
k = 0.166
t0 = -1.75

# Maturity
M50 = 1.3 # L miniatus
M95 = 1.85 # L miniatus


#* Read in Data ####
setwd(data_dir)
Ab2021 <- readRDS("2021-05_Abrolhos_BOSS-BRUV_dat_length.rds") %>% 
  filter(scientific %in% "Lethrinidae Lethrinus miniatus") %>% 
  filter(!is.na(length)) %>% 
  dplyr::select(sample, length, status)

Full.Lengths <- Ab2021

histograms <- ggplot()+
  geom_histogram(data=Full.Lengths, aes(x=length))+
  facet_grid(~status)


#* Work out the proportions in the different length groups ####
LM <- Linf*(1-exp(-k*(M50-t0))) # Age at 50% Mature is 3.62
# Linf.LM <- Linf - ((Linf-LM)/2)
Linf.LM <- Linf - ((Linf-LM)*(2/3))


Length.Props <- Full.Lengths %>% 
  mutate(Length.Group = ifelse(length < LM, "Small",
                               ifelse(length>Linf.LM, "Biggest", "Big"))) %>%  
  filter(!Length.Group %in% "Small") %>% 
  group_by(status, Length.Group) %>% 
  summarise(n.group = n()) %>% 
  mutate(Proportion = n.group/sum(n.group)) %>% 
  ungroup() 

#### PILBARA ####
#* L sebae ####
Linf =  524.77 # L sebae https://pdf.sciencedirectassets.com/272396/1-s2.0-S0272771400X01096/1-s2.0-S027277140190887X/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEL%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJIMEYCIQD5TfFWbdlZqdgqRd8enp%2FQDMCr6KF8EEKxPiihysqWtgIhAJ9K%2F5Vdpe62iy6%2FTLi486IY%2FS5RIVOq8ZfUvDeFt7xbKrwFCLj%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQBRoMMDU5MDAzNTQ2ODY1Igy4Dy9u0Y7T1G%2FWAHcqkAUcbJVlJffuuphzTI13cHcot40N1A4wyBkFpXNbBu6rI79FmtEM7FPTFEIwZxL01KH%2BL3csPAxokkWxNyVkhaZWdTtVk%2F8RZPQMp5lKi%2B5yWFJjfWNHqki0%2FfWwV5G5SQwVgl8k6zD0Gk1%2FimS0kVeJgVW%2BTDff%2FJN%2B%2F%2FUIDabbcV0CGzjL%2Bd74hvegDI57S701z1FuJYfVgk2qaNmaMRgUgRsPRPNi9Cp5E7ya4G%2FcFJUnTzKs%2BciTdhvoQGCj7ZR9V3xAyF3TbWYUs%2Fk4kKtZ9Raas5P85%2F6gAvZijjD02W59d4HsrfaGTa1%2F3DBA5ITD%2BWdOoOzCPXIgt39a68YeTv87f9OW37GJAjSJ3R2Ed3XeOazhNJag7NbKbYQrO8a8N9UC%2FEzWGegv2P5HrD81ycI4s1%2Fm2VNJ5rKL5TJb96NWJFXhmzRLeKDTdMGcgWeo4%2B43amxt5rqA1Q62KsmX1aD6Xz46FTWUjjjfJulKW4kIXwhxA5hi0I591StMCSjMNxAwHPZzomj1%2Bu7%2F1TZBegse%2BTqZ0MuxCu4WQYNzpoS5a0KAqmWFm4JDGv70ysu2D80x6cKF%2FH9Rj76fdsho%2F8kGk2KaGXKQPI8qC64TCuHYr63Kn9CGnl%2FNGPn5J10XHz89%2FkXhbPpLzM%2Bi5ILXkQuVz0MNvp3YnhCLCxYQX91MovDxLnIGPKahhAtzD2EuUxEL8qh9jVE61V28%2BQInHoCYv%2BJOd7%2BX%2BKeDW%2B%2BOiKnCXRuGlYKmkViTAAcyBgJzoELu4SeDAmXWZ9mGlLX9JvtXp4XNT%2F5ssB8hhLwfwtNWEkBqLI36VccnjpfEa6w%2B7rHJqkdG%2B6oWEpk61SiFJ7ry8eLVLMvx3jfdIa9wkDCBwaWvBjqwAcsyslAnrZylMASgommY4HoU6u0iSWOI4O7RyJYoWUT3OHvYMMJypNFb9dRCyo%2BbHKFwa7d%2Fv4PCzkqCr927uMyS46PhUZW5hJsetzDEhEpP0P81lAlWjXbMLse%2FH3AR%2FS%2BLRWFBxHoIAFJ0fl84tWazeC9F26ta5pZ6vlTsYOYeSr%2FIko8Yh%2F5lbOvKbNfKweiD99d6XxTgWT9%2B3qidXvUaInflrtnowvo22cjLLTdP&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240307T065546Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY6JFP3BHD%2F20240307%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=820b394fc34f86763ba8f00526fe71fead5a1c61b6e6f1ea10bcaac5a411872d&hash=ff25a4a06d359e65d24dbef00a8f60fa5734291cc1930e5c54c3e63b6dc9932b&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S027277140190887X&tid=spdf-ab407c3e-20ce-4b5d-855d-4e24ba3805b0&sid=828ef6fd661c544e512a70b724c203cd7a45gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=04125e5552590e055706&rr=8608ae8a58d55aa8&cc=au
k =  0.2330
t0 = - 1.533

# Maturity
M50 = 3.43 # L sebae
M95 = 4.37 # L sebae

#* P multidens ###
Linf = 590 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
k = 0.26 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
t0 = -0.46 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

# Maturity
M50 = 4.8 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
M95 = 11.95 # P multidens https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

#* E multinotatus ###
Linf = 666.4 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
k = 0.221 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
t0 = -1.835 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf

# Maturity
M50 = 2 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf
M95 = 2.2849 # E multinotatus https://www.fish.wa.gov.au/Documents/resource_assessment/resource_assessment_report_013.pdf


#* Read in Data ####
setwd(data_dir)
Pilbara <- readRDS("Pilbara_data") %>% 
  #filter(Length < 2000) %>% 
  filter(scientific %in% "Lutjanus sebae") %>% 
  dplyr::select(Sample, Length, Status)

Full.Lengths <- Pilbara %>% 
  glimpse()

histograms <- ggplot()+
  geom_histogram(data=Full.Lengths, aes(x=Length))+
  facet_grid(~Status)

#* Work out the proportions in the different length groups ####
LM <- 480 #Linf*(1-exp(-k*(M50-t0))) # Age at 50% Mature is 3.62
Linf.LM <- Linf - ((Linf-LM)/2)
Linf.LM <- Linf - ((Linf-LM)*(2/3))


Length.Props <- Full.Lengths %>% 
  mutate(Length.Group = ifelse(Length < LM, "Small",
                               ifelse(Length>Linf.LM, "Biggest", "Big"))) %>%  
  filter(!Length.Group %in% "Small") %>% 
  group_by(Status, Length.Group) %>% 
  summarise(n.group = n()) %>% 
  mutate(Proportion = n.group/sum(n.group)) %>% 
  ungroup() 
