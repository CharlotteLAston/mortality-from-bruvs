###################################################

# Setting up real data to go into the next step

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

#### READ IN DATA ####
setwd(data_dir)

lengths <- read.csv("2010-09_Pilbara_stereoBRUV_Length (1).csv")
metadata <- read.csv("2010-09_Pilbara_stereoBRUV_Metadata (3).csv")

#### FORMAT DATA ####

full.data <- lengths %>% 
  left_join(., metadata, by="Sample") %>% 
  mutate(scientific = paste0(Genus, sep=" ", Species)) %>% 
  filter(scientific %in% c("Lutjanus sebae", "Epinephelus multinotatus", "Pristipomoides multidens"))

setwd(data_dir)
saveRDS(full.data, "Pilbara_data")

#### Western King Wrasse #####

dat <- readRDS("australian-synthesis_complete_length.RDS") %>% 
  filter(length <= 400) %>% 
  mutate(ID = row_number())

mps <- st_read("Rottnest_Sanctuaries.shp") %>% 
  st_transform(gdacrs)

plot(mps$geometry)

sp.dat <- st_as_sf(dat, coords=c("longitude", "latitude")) 
gdacrs <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
st_crs(sp.dat) <- gdacrs
plot(sp.dat$geometry, add=T)

NTZ_points <- st_intersection(sp.dat, mps) %>% 
  mutate(status = "NTZ") %>% 
  dplyr::select(ID, status) %>% 
  st_drop_geometry()

dat <- dat %>% 
  left_join(., NTZ_points, by="ID") %>% 
  mutate(status.x = ifelse(!is.na(status.y), status.y, status.x)) %>% 
  mutate(status.x = ifelse(is.na(status.x), "Fished", status.x)) %>%  # anything that didn't intersect is in the fished zone
  dplyr::select(!status.y)

dat <- dat %>% 
  rename(status = "status.x")

setwd(data_dir)
saveRDS(dat, "Western-King-Wrasse_Status")



