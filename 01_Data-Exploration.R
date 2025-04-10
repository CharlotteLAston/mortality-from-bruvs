###################################################
# Having a look at DPIRD data and plotting 
# preliminary plots to see how it lines up with
# BRUV data
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
library(patchwork)
library(ggnewscale)
library(scales)
library(ggspatial)
library(tidyterra)
library(grid)
library(gridExtra)
library(gtable)

rm(list = ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_data", sep="/")

a4.width <- 160

#* Read in broad spatial data needed for plots ####
setwd(sp_dir)

e <- ext(110.8, 127, -36, -20)

aus <- st_read("cstauscd_r.mif", crs=4283) %>% 
  filter(FEAT_CODE %in% c("mainland", "island")) 
ausc <- st_crop(aus, e)

aumpa <- st_read("AustraliaNetworkMarineParks.shp") 
mpa <- st_crop(aumpa,e) %>% 
  dplyr::mutate(Area_KM2 = as.character(Area_KM2)) %>%
  dplyr::filter(!ResName %in% c("Perth Canyon"),
                !Area_KM2 %in% c("2107.97750059", "14448.0532117",
                                 "22493.1162347", "14231.1456468"))


mpa$ZoneName <- factor(mpa$ZoneName, levels = c("Recreational Use Zone",
                                                "Multiple Use Zone", 
                                                "Special Purpose Zone",
                                                "National Park Zone"))
npz <- mpa[mpa$ZoneName %in% "National Park Zone", ] 

wampa <- st_read("WA_MPA_2020.shp")
st_crs(wampa) <- st_crs(aumpa) 
wampa <- st_make_valid(wampa)
wampa <- st_crop(wampa, e)  

# simplify zone names
wampa$waname <- gsub("( \\().+(\\))", "", wampa$ZONE_TYPE)
wampa$waname <- gsub(" [1-4]", "", wampa$waname)
wampa$waname[wampa$NAME == "Ngari Capes"]     <- "General Use"
wampa$waname <- dplyr::recode(wampa$waname, 
                              "General Use" = "General Use Zone",
                              # "MMA" = "Marine Management Area",
                              # "Recreation Area" = "Recreation Zone",
                              # "Conservation Area" = "Sanctuary Zone",
                              "Special Purpose Zone (Shore Based Activities)" = 
                                "Special Purpose Zone\n(Shore Based Activities)")

unique(wampa$waname)
sanc <- wampa %>%
  dplyr::filter(waname %in% "Sanctuary Zone")

bathy <- rast("bath_250_good.tif") %>%
  crop(e)

#* Read in DPIRD data ####
setwd(data_dir)
dpird <- read.csv("Length data for Charlotte with brownspot.csv") %>% 
  dplyr::select(Species, Year, Region, Lat, Lon, TL)

dpird_sp <- dpird %>% 
  filter(!is.na(Lat))%>% 
  mutate(sample_type = ifelse(Species %in% c("Pink Snapper", "Breaksea Cod"), "Fishery dependent", "Fishery independent"))


#* Read in BRUV data ####
setwd(data_dir)

PS_data <- readRDS("australian-synthesis_complete_length_chrysophrys_auratus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Marmion|TwoRocks|Two.Rocks|Capes|Geographe|SwC|south-west|Ngari|Warnbro|Rottnest"))) %>% 
  dplyr::filter(str_detect(campaign,"2006|2007|2008|2009|2010|2011")) %>% 
  dplyr::filter(!str_detect(campaign, "Perth.Canyon")) %>% 
  dplyr::select(sample, campaign, latitude, longitude)%>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)

# WKW_metro_data <- readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
#   dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Two.Rocks|TwoRocks|Warnbro|Marmion|Rottnest"))) %>%
#   dplyr::filter(str_detect(campaign, "2006|2007|2008|2009")) %>% 
#   dplyr::select(sample, campaign, latitude, longitude) %>% 
#   rename(x = longitude, y = latitude) %>% 
#   distinct(sample, .keep_all=TRUE)

WKW_jurien_data <- readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) %>%
  # dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Abrolhos_stereo-BRUVs"))) %>%
  # filter(length<=MaxLen) %>% 
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009")) %>% 
  dplyr::select(sample, campaign, latitude, longitude) %>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)

# WKW_capes_data <- readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
#   dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Capes"))) %>%
#   # dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Abrolhos_stereo-BRUVs"))) %>%
#   # filter(length<=MaxLen) %>% 
#   dplyr::filter(str_detect(campaign, "2006|2007|2008|2009|2010|2011")) %>% 
#   dplyr::select(sample, campaign, latitude, longitude) %>% 
#   rename(x = longitude, y = latitude) %>% 
#   distinct(sample, .keep_all=TRUE)

MW_metro_data <-readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Rottnest|TwoRocks|Marmion|Warnbro|Two.Rocks"))) %>%
  dplyr::filter(str_detect(campaign,"2006|2007|2008|2009")) %>% 
  dplyr::select(sample, campaign, latitude, longitude) %>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)

# MW_jurien_data <- readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
#   dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) %>% 
#   dplyr::filter(str_detect(campaign,"2006|2007|2008|2009")) %>% 
#   dplyr::select(sample, campaign, latitude, longitude) %>% 
#   rename(x = longitude, y = latitude) %>% 
#   distinct(sample, .keep_all=TRUE)

# MW_capes_data <- readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
#   dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Capes"))) %>% 
#   dplyr::filter(str_detect(campaign, "2006|2007|2008|2009|2010|2011")) %>% 
#   dplyr::select(sample, campaign, latitude, longitude) %>% 
#   rename(x = longitude, y = latitude) %>% 
#   distinct(sample, .keep_all=TRUE)

BS_Data <- readRDS("australian-synthesis_complete_length_epinephelides-armatus.RDS") %>% 
  # filter(length<MaxLen) %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Marmion|TwoRocks|Two.Rocks|Capes|Geographe|SwC|south-west|Ngari|Warnbro|Rottnest"))) %>% 
  dplyr::filter(str_detect(campaign,"2006|2007|2008|2009|2010|2011")) %>% 
  dplyr::select(sample, campaign, latitude, longitude) %>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)

BW_Data <- readRDS("australian-synthesis_complete_length_notolabrus-parilus.RDS") %>%
  # filter(length<MaxLen) %>%
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"2006|2007|2008|2009"))) %>%
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"JMarmion|TwoRocks|Two.Rocks|Warnbro|Rottnest"))) %>%
  dplyr::select(sample, campaign, latitude, longitude) %>%
  rename(x = longitude, y = latitude) %>%
  distinct(sample, .keep_all=TRUE)

all_data <- rbind(PS_data, WKW_jurien_data, MW_metro_data, BS_Data, BW_Data) 

#* plot heat map of all samples ####

p1 <- ggplot() +
  geom_spatraster(data = bathy) +
  scale_fill_gradient(high = "#b9cded", low = "#09357d",  name = "Depth (m)") +
  new_scale_fill() +
  stat_density2d(data = all_data, aes(x = x, y = y, fill = after_stat(level),
                                           alpha = after_stat(level)),
                 geom = "polygon", binwidth = 0.05, show.legend = T)+
  scale_fill_gradient2(low = "#f5ece1", mid = "#f2b46d", high = "#f78a05", name="Area covered by\nstereo-BRUV sampling") +
  scale_alpha(range = c(0.5, 0.8), guide = "none") +   
  geom_point(data = dpird_sp, aes(x=Lon, y=Lat, shape=sample_type), size=1, fill=NA)+
  scale_shape_manual(values=c(4,19), name="Line fishing\nsampling type")+
  geom_sf(data = ausc, fill = "seashell2", colour = "black") +
  coord_sf(xlim = c(114, 118), #124.9889
           ylim = c(-36, -30)) + #-35.5, 21.9
  theme_minimal()+
  theme(legend.position = c(1.4,0.3),
        legend.justification = "right",
        legend.key = element_rect(fill=NA, colour=NA))+
  annotation_scale(location="tr")
p1

inset.map <- ggplot() +
  geom_sf(data = aus, fill = "seashell2", colour = "grey90", size = 0.075) +
  #geom_sf(data = aumpa, alpha = 5/6, colour = "grey85", size = 0.02) +
  coord_sf(xlim = c(108, 125), ylim = c(-37, -13)) +
  annotate(geom = "text", x=c(110.2), y=c(-29.4), label = c("Indian\nOcean"), size=3)+
  annotate(geom = "text", x=c(120.8), y=c(-25.94), label = c("Western\nAustralia"), size=3.5)+
  # annotate("rect", xmin = 113, xmax = 114.35, ymin = -22.8, ymax = -21.5,   # Change here
  #          colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  # annotate("rect", xmin = 113, xmax = 115.061, ymin = -29.25, ymax = -28,
  #          colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  annotate("rect", xmin=114, xmax=118, ymin=-36, ymax=-30,
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2)+
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "grey70"))+
  ylab(NULL)+
  xlab(NULL)
inset.map

sample.plot <- p1 + inset_element(inset.map, left = -1.29, right = 3.8, top = 1, bottom = 0.6)  
sample.plot

setwd(fig_dir)
ggsave(sample.plot, filename="Sample_plot.png", height = a4.width*1.5, width = a4.width*2.5, units  ="mm", dpi = 300 )


#### Spatial plots of sample overlap ####
#* Pink Snapper

snapper <- ggplot() +
  geom_spatraster(data = bathy) +
  scale_fill_gradient(high = "#b9cded", low = "#09357d",  name = "Depth") +
  new_scale_fill() +
  stat_density2d(data = PS_data, aes(x = x, y = y, fill = after_stat(level),
                                           alpha = after_stat(level)),
                 geom = "polygon", binwidth = 0.05, show.legend = F) +
  scale_fill_gradient2(low = "#f5ece1", mid = "#f2b46d", high = "#f78a05") +
  scale_alpha(range = c(0.5, 0.8), guide = "none") +   
  geom_sf(data = ausc, fill = "seashell2", colour = "black") +
  geom_point(data = dpird_sp %>% filter(Species %in% "Pink Snapper"), aes(x=Lon, y=Lat))+
  coord_sf(xlim = c(112, 118), #124.9889
           ylim = c(-36, -30.5)) + #-35.5, 21.9
  theme_minimal()+
  annotation_scale(location="tr")
snapper

setwd(fig_dir)
ggsave(snapper, filename="snapper_plot.png", height = a4.width*1, width = a4.width*2, units  ="mm", dpi = 300 )

#* Breaksea

breaksea <- ggplot() +
  geom_spatraster(data = bathy) +
  scale_fill_gradient(high = "#b9cded", low = "#09357d",  name = "Depth") +
  new_scale_fill() +
  stat_density2d(data = BS_Data, aes(x = x, y = y, fill = after_stat(level),
                                     alpha = after_stat(level)),
                 geom = "polygon", binwidth = 0.05, show.legend = F) +
  scale_fill_gradient2(low = "#f5ece1", mid = "#f2b46d", high = "#f78a05") +
  scale_alpha(range = c(0.5, 0.8), guide = "none") +   
  geom_sf(data = ausc, fill = "seashell2", colour = "black") +
  geom_point(data = dpird_sp %>% filter(Species %in% "Breaksea Cod"), aes(x=Lon, y=Lat))+
  coord_sf(xlim = c(112, 118), #124.9889
           ylim = c(-36, -30.5)) + #-35.5, 21.9
  theme_minimal()+
  annotation_scale(location="tr")
breaksea

setwd(fig_dir)
ggsave(breaksea, filename="breaksea_plot.png", height = a4.width*1, width = a4.width*2, units  ="mm", dpi = 300 )

#* Western King Wrasse
temp_wkw <- rbind(WKW_jurien_data) 


KingWrasse <- ggplot() +
  geom_spatraster(data = bathy) +
  scale_fill_gradient(high = "#b9cded", low = "#09357d",  name = "Depth") +
  new_scale_fill() +
  stat_density2d(data = temp_wkw, aes(x = x, y = y, fill = after_stat(level),
                                     alpha = after_stat(level)),
                 geom = "polygon", binwidth = 0.05, show.legend = F) +
  scale_fill_gradient2(low = "#f5ece1", mid = "#f2b46d", high = "#f78a05") +
  scale_alpha(range = c(0.5, 0.8), guide = "none") +   
  geom_sf(data = ausc, fill = "seashell2", colour = "black") +
  geom_point(data = dpird_sp %>% filter(Species %in% "King Wrasse" & Region %in% "Jurien"), aes(x=Lon, y=Lat))+
  coord_sf(xlim = c(113, 117), #124.9889
           ylim = c(-29.5, -33)) + #-35.5, 21.9
  theme_minimal()+
  annotation_scale(location="tr")
KingWrasse

setwd(fig_dir)
ggsave(KingWrasse, filename="KingWrasse_plot.png", height = a4.width*1, width = a4.width*2, units  ="mm", dpi = 300 )

#* Maori Wrasse
temp_mw <- rbind(MW_metro_data) 

MaoriWrasse <- ggplot() +
  geom_spatraster(data = bathy) +
  scale_fill_gradient(high = "#b9cded", low = "#09357d",  name = "Depth") +
  new_scale_fill() +
  stat_density2d(data = temp_mw, aes(x = x, y = y, fill = after_stat(level),
                                      alpha = after_stat(level)),
                 geom = "polygon", binwidth = 0.05, show.legend = F) +
  scale_fill_gradient2(low = "#f5ece1", mid = "#f2b46d", high = "#f78a05") +
  scale_alpha(range = c(0.5, 0.8), guide = "none") +   
  geom_sf(data = ausc, fill = "seashell2", colour = "black") +
  geom_point(data = dpird_sp %>% filter(Species %in% "Maori Wrasse"&Region %in% "Metro"), aes(x=Lon, y=Lat))+
  coord_sf(xlim = c(114, 117), #124.9889
           ylim = c(-32.5, -30.5)) + #-35.5, 21.9
  theme_minimal()+
  annotation_scale(location="tr")
MaoriWrasse

setwd(fig_dir)
ggsave(MaoriWrasse, filename="MaoriWrasse_plot.png", height = a4.width*1, width = a4.width*2, units  ="mm", dpi = 300 )

#* Brownspot Wrasse
temp_mw <- rbind(BW_Data) 

BrownspotWrasse <- ggplot() +
  geom_spatraster(data = bathy) +
  scale_fill_gradient(high = "#b9cded", low = "#09357d",  name = "Depth") +
  new_scale_fill() +
  stat_density2d(data = BW_Data, aes(x = x, y = y, fill = after_stat(level),
                                     alpha = after_stat(level)),
                 geom = "polygon", binwidth = 0.05, show.legend = F) +
  scale_fill_gradient2(low = "#f5ece1", mid = "#f2b46d", high = "#f78a05") +
  scale_alpha(range = c(0.5, 0.8), guide = "none") +   
  geom_sf(data = ausc, fill = "seashell2", colour = "black") +
  geom_point(data = dpird_sp %>% filter(Species %in% "Brownspot Wrasse"&Region %in% "Metro"), aes(x=Lon, y=Lat))+
  coord_sf(xlim = c(113, 118), #124.9889
           ylim = c(-32.5, -30.5)) + #-35.5, 21.9
  theme_minimal()+
  annotation_scale(location="tr")
BrownspotWrasse

setwd(fig_dir)
ggsave(BrownspotWrasse, filename="BrownspottedWrasse_plot.png", height = a4.width*1, width = a4.width*2, units  ="mm", dpi = 300 )


#### Plots of length distributions ####
#* Pink snapper ####
setwd(data_dir)
PS_length <- readRDS("australian-synthesis_complete_length_chrysophrys_auratus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Marmion|TwoRocks|Two.Rocks|Warnbro|Rottnest|Capes|Geographe|SwC|south-west|Ngari"))) %>% #Capes|Geographe|SwC|south-west|Ngari|
  dplyr::filter(str_detect(campaign,"2006|2007|2008|2009|2010|2011")) %>% 
  dplyr::filter(!str_detect(campaign, "Perth.Canyon")) %>% 
  dplyr::select(species, length) %>% 
  mutate(Method = "BRUV")

PS_dpird <- dpird %>% 
  filter(Species %in% "Pink Snapper") %>% 
  filter(Region %in% c("Metro", " South West")) %>% 
  select(Species, TL) %>% 
  rename(length = "TL",
         species = "Species") %>% 
  mutate(Method = "DPIRD")
  
snapper_lengths_metro <- rbind(PS_length, PS_dpird)

snapper_metro_plot <- ggplot()+
  geom_histogram(data=snapper_lengths_metro, aes(x=length, fill=Method), position="identity", alpha=0.5)+
  scale_fill_manual(values=c("#66CCEE", "#BBCC33"), labels=c("Stereo-BRUV", "Line fishing"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()



#* Breaksea Cod  ####
setwd(data_dir)

BS_length <-  readRDS("australian-synthesis_complete_length_epinephelides-armatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Marmion|TwoRocks|Two.Rocks|Warnbro|Rottnest"))) %>% #Capes|Geographe|SwC|south-west|Ngari|
  dplyr::filter(str_detect(campaign,"2006|2007|2008|2009|2010|2011")) %>% 
  dplyr::select(species, length) %>% 
  mutate(Method = "BRUV")

BS_dpird <- dpird %>% 
  filter(Species %in% "Breaksea Cod") %>% 
  filter(Region %in% "Metro") %>% 
  select(Species, TL) %>% 
  rename(length = "TL",
         species = "Species") %>% 
  mutate(Method = "DPIRD")

Breaksea_lengths_metro <- rbind(BS_length, BS_dpird)

Breaksea_metro_plot <- ggplot()+
  geom_histogram(data=Breaksea_lengths_metro, aes(x=length, fill=Method), position="identity", alpha=0.5)+
  scale_fill_manual(values=c("#66CCEE", "#BBCC33"), labels=c("Stereo-BRUV", "Line fishing"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()
Breaksea_metro_plot
#_____________________________________________________________________________________________________________

BS_length <- readRDS("australian-synthesis_complete_length_epinephelides-armatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"2006|2007|2008|2009|2010|2011"))) %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Capes|Geographe|SwC|south-west|Ngari"))) %>% 
  dplyr::filter(!str_detect(campaign, "Perth.Canyon")) %>% 
  dplyr::select(species, length) %>% 
  mutate(Method = "BRUV")

BS_dpird <- dpird %>% 
  filter(Species %in% "Breaksea Cod") %>% 
  filter(Region %in% "South-West") %>% 
  select(Species, TL) %>% 
  rename(length = "TL",
         species = "Species") %>% 
  mutate(Method = "DPIRD")

Breaksea_lengths_SW <- rbind(BS_length, BS_dpird)

Breaksea_SW_plot <- ggplot()+
  geom_histogram(data=Breaksea_lengths_SW, aes(x=length, fill=Method))+
  theme_classic()
Breaksea_SW_plot


#* Western King Wrasse ####
setwd(data_dir)

WKW_length <-  readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Two.Rocks|TwoRocks|Warnbro|Marmion|Rottnest"))) %>%
  dplyr::filter(str_detect(campaign, "2005|2006|2007|2008")) %>% 
  dplyr::select(species, length) %>% 
  mutate(Method = "BRUV")

WKW_dpird <- dpird %>% 
  filter(Species %in% "King Wrasse") %>% 
  filter(Region %in% "Metro") %>% 
  select(Species, TL) %>% 
  rename(length = "TL",
         species = "Species") %>% 
  mutate(Method = "DPIRD")

KingWrasse_lengths_metro <- rbind(WKW_length, WKW_dpird)

KingWrasse_metro_plot <- ggplot()+
  geom_histogram(data=KingWrasse_lengths_metro, aes(x=length, fill=Method), position="identity", alpha=0.5)+
  scale_fill_manual(values=c("#66CCEE", "#BBCC33"), labels=c("Stereo-BRUV", "Line fishing"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()
KingWrasse_metro_plot
#_____________________________________________________________________________________________________________

WKW_length <-  readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) %>% 
  dplyr::filter(str_detect(campaign,"2005|2006|2007|2008")) %>% 
  dplyr::select(species, length) %>% 
  mutate(Method = "BRUV")

WKW_dpird <- dpird %>% 
  filter(Species %in% "King Wrasse") %>% 
  filter(Region %in% "Jurien") %>% 
  select(Species, TL) %>% 
  rename(length = "TL",
         species = "Species") %>% 
  mutate(Method = "DPIRD")

KingWrasse_lengths_Jur <- rbind(WKW_length, WKW_dpird)

KingWrasse_Jur_plot <- ggplot()+
  geom_histogram(data=KingWrasse_lengths_Jur, aes(x=length, fill=Method), position="identity", alpha=0.5)+
  scale_fill_manual(values=c("#66CCEE", "#BBCC33"), labels=c("Stereo-BRUV", "Line fishing"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()
KingWrasse_Jur_plot

#* Maori Wrasse ####
setwd(data_dir)

MW_length <-readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
  dplyr::filter(str_detect(campaign, "2005|2006|2007|2008")) %>%
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Rottnest|TwoRocks|Marmion|Warnbro|Two.Rocks"))) %>% 
  dplyr::select(species, length) %>% 
  mutate(Method = "BRUV")

MW_dpird <- dpird %>% 
  filter(Species %in% "Maori Wrasse") %>% 
  filter(Region %in% "Metro") %>% 
  select(Species, TL) %>% 
  rename(length = "TL",
         species = "Species") %>% 
  mutate(Method = "DPIRD")

MaoriWrasse_lengths_metro <- rbind(MW_length, MW_dpird)

MaoriWrasse_metro_plot <- ggplot()+
  geom_histogram(data=MaoriWrasse_lengths_metro, aes(x=length, fill=Method), position="identity", alpha=0.5)+
  scale_fill_manual(values=c("#66CCEE", "#BBCC33"), labels=c("Stereo-BRUV", "Line fishing"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()
MaoriWrasse_metro_plot

#* Brownspot Wrasse ####
setwd(data_dir)

BW_length <-readRDS("australian-synthesis_complete_length_notolabrus-parilus.RDS") %>%
  dplyr::filter(str_detect(campaign, "2006|2006|2007|2009")) %>%
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Rottnest|TwoRocks|Marmion|Warnbro|Two.Rocks"))) %>%
  dplyr::select(species, length) %>%
  mutate(Method = "BRUV")
range(BW_length$length)
estimate_mode(BW_length$length)

BW_dpird <- dpird %>%
  filter(Species %in% "Brownspot Wrasse") %>%
  filter(Region %in% "Metro") %>%
  select(Species, TL) %>%
  rename(length = "TL",
         species = "Species") %>%
  mutate(Method = "DPIRD")
range(BW_dpird$length)
estimate_mode(BW_dpird$length)

BrownspotWrasse_lengths_metro <- rbind(BW_length, BW_dpird)

BrownspotWrasse_metro_plot <- ggplot()+
  geom_histogram(data=BrownspotWrasse_lengths_metro, aes(x=length, fill=Method), position="identity",alpha=0.5)+
  scale_fill_manual(values=c("#66CCEE", "#BBCC33"), labels=c("Stereo-BRUV", "Line fishing"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()
BrownspotWrasse_metro_plot

setwd(data_dir)

BW_length <-readRDS("australian-synthesis_complete_length_notolabrus-parilus.RDS") %>%
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009")) %>%
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) %>%
  dplyr::select(species, length) %>%
  mutate(Method = "BRUV")

BW_dpird <- dpird %>%
  filter(Species %in% "Brownspot Wrasse") %>%
  filter(Region %in% "Jurien") %>%
  select(Species, TL) %>%
  rename(length = "TL",
         species = "Species") %>%
  mutate(Method = "DPIRD")

BrownspotWrasse_lengths_jurien <- rbind(BW_length, BW_dpird)

BrownspotWrasse_jurien_plot <- ggplot()+
  geom_histogram(data=BrownspotWrasse_lengths_jurien, aes(x=length, fill=Method), position="identity",alpha=0.5)+
  scale_fill_manual(values=c("#66CCEE", "#BBCC33"), labels=c("Stereo-BRUV", "Line fishing"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()
BrownspotWrasse_jurien_plot

#### Put plots together ####
setwd(fig_dir)
x.label <- textGrob("Length (mm)", gp=gpar(fontsize=10))
y.label <- textGrob("Count", gp=gpar(fontsize=10), rot=90)
legend <- gtable_filter(ggplotGrob(snapper_metro_plot), "guide-box")


Length_Dist_Plot <- grid.arrange(arrangeGrob(snapper_metro_plot + theme(legend.position="none"),
                                          Breaksea_metro_plot + theme(legend.position="none"),
                                          KingWrasse_Jur_plot + theme(legend.position="none"),
                                          MaoriWrasse_metro_plot + theme(legend.position="none"),
                                          BrownspotWrasse_metro_plot + theme(legend.position="none"),
                                          left=y.label,
                                          bottom=x.label,
                                          right=legend,
                                          ncol=2))
ggsave(Length_Dist_Plot, filename="Length_Dist_Plot.png",height = a4.width*1, width = a4.width*1.1, units  ="mm", dpi = 300 )




