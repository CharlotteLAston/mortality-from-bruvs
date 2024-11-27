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

working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_data", sep="/")

a4.width=160

#### Make data frame ####

dat <- data.frame(species=as.factor(c("L. nebulosus","L. nebulosus" ,"E. armatus", "C. auratus", "C. auratus","C. auricularis", "C. auricularis", "O. lineolatus", "O. lineolatus")),
                  location = c("Ningaloo", "Ningaloo*","Perth Metro and South West", "Perth Metro and South West", "Perth Metro and South West*","Abrolhos", "Perth Metro", "Perth Metro", "South West"),
                  mortality = c(0.27, 0.227, 0.1, 0.23, 0.17, 0.13, 0.32, 0.12, 0.088),
                  #nat.mort = c(0.146, 0.3, 0.102, 0.41, 0.41, 0.36, 0.36),
                  upp_bnd = c(0.35, 0.351, 0.034, 0.26, 0.2, 0.23, 0.294, 0.24, 0.115),
                  lwr_bnd = c(0.2, 0.148, 0.24, 0.2, 0.12, 0.071, 0.337, 0.06, 0.067),
                  Shape = c(16,17,16,16,17,16,16,16,16)) %>% 
  mutate(species.location = paste0(species, sep="_", location)) %>% 
  mutate(species.location=fct_relevel(species.location, "L. nebulosus_Ningaloo","L. nebulosus_Ningaloo*","E. armatus_Perth Metro and South West","C. auratus_Perth Metro and South West",
                                      "C. auratus_Perth Metro and South West*","C. auricularis_Perth Metro", "O. lineolatus_Perth Metro", "C. auricularis_Abrolhos","O. lineolatus_South West")) 

#colours <- c("#88CBED", "#A9439A", "#332387", "#117633", "#43A999", "#872155", "#CB6778")
species.labels <- c("*L. nebulosus*<br>(Ningaloo)", "*L. nebulosus*&#42;<br>(Ningaloo)","*E. armatus*<br>(Metro<br>and SW)", "*C. auratus*<br>(Metro<br>and SW)", "*C. auratus*&#42;<br>(Metro<br>and SW)",
                    "*C. auricularis*<br>(Metro)", "*O. lineolatus*<br>(Metro)", "*C. auricularis*<br>(Abrolhos)", "*O. lineolatus*<br>(SW)")
#### Make Plot ####

condition_plot <- dat %>% 
  ggplot(.)+
  # geom_rect(xmin=0, xmax=8.5,ymin=0, ymax=0.4, colour=NA, fill="#e06666ff", alpha=0.15)+
  # geom_rect(xmin=0, xmax=8.5,ymin=0.4, ymax=0.5, colour=NA, fill="#ffe599ff", alpha=0.15)+
  # geom_rect(xmin=0, xmax=8.5,ymin=0.5, ymax=0.6, colour=NA, fill="#a4c2f4ff", alpha=0.15)+
  # geom_rect(xmin=0, xmax=8.5,ymin=0.6, ymax=0.9, colour=NA, fill="#6aa84fff", alpha=0.15)+
  # geom_rect(xmin=0, xmax=8.5,ymin=0.9, ymax=1, colour=NA, fill="#539653ff", alpha=0.15)+
  # geom_rect(xmin=0, xmax=7.5,ymin=0, ymax=1, colour=NA, fill="grey100", alpha=0.05)+
  geom_point(aes(x=species.location, y=mortality, fill=species.location, colour=species.location), shape=dat$Shape, size=2)+
  scale_fill_manual(values=c("#88CBED", "#88CBED", "#A9439A", "#332387", "#332387", "#117633","#43A999", "#872155", "#CB6778"), guide="none")+
  geom_linerange(aes(x=species.location, ymax=upp_bnd, ymin=lwr_bnd, colour=species.location))+
  scale_colour_manual(values=c("#88CBED", "#88CBED", "#A9439A", "#332387", "#332387", "#117633", "#43A999", "#872155", "#CB6778"), guide="none")+
  # geom_vline(xintercept = 7.5, colour="white", linewidth=1.5)+
  # geom_vline(xintercept = 1.5, colour="white", linewidth=1.5)+
  theme_classic()+
  scale_x_discrete(labels=species.labels)+
  ylim(0,1)+
  ylab("Fishing Mortality")+
  xlab(NULL)+
  theme(axis.text.x=ggtext::element_markdown())
condition_plot

setwd(fig_dir)
ggsave(condition_plot, filename="Condition_plot.png", height = a4.width*1, width = a4.width*1.2, units  ="mm", dpi = 300 )


#### Plotting Sample Map ####
#* Load in data ####
setwd(data_dir)
# SE_data <- readRDS("australian-synthesis_complete_length_lethrinus_nebulosus.RDS") %>% 
#   dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"2006|2007|2008|2009|2010|2011|2012|2013|2014|2015|2016"))) %>% 
#   dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Ningaloo|Exmouth")))  %>% 
#   dplyr::select(sample, campaign, latitude, longitude) %>% 
#   rename(x = longitude, y = latitude) %>% 
#   filter(x <115.1) %>% 
#   distinct(sample, .keep_all=TRUE)

PS_data <- readRDS("australian-synthesis_complete_length_chrysophrys_auratus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"2006|2007|2008|2009|2010|2011"))) %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Marmion|TwoRocks|Two.Rocks|Capes|Geographe|SwC|south-west|Ngari|Warnbro|Rottnest"))) %>% 
  dplyr::filter(!str_detect(campaign, "Perth.Canyon")) %>% 
  dplyr::select(sample, campaign, latitude, longitude)%>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)

WKW_metro_data <- readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Two.Rocks|TwoRocks|Warnbro|Marmion|Rottnest"))) %>%
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009")) %>% 
  dplyr::select(sample, campaign, latitude, longitude) %>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)
  
WKW_jurien_data <- readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) %>%
  # dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Abrolhos_stereo-BRUVs"))) %>%
  # filter(length<=MaxLen) %>% 
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009")) %>% 
  dplyr::select(sample, campaign, latitude, longitude) %>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)

WKW_capes_data <- readRDS("australian-synthesis_complete_length_coris_auricularis.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Capes"))) %>%
  # dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Abrolhos_stereo-BRUVs"))) %>%
  # filter(length<=MaxLen) %>% 
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009|2010|2011")) %>% 
  dplyr::select(sample, campaign, latitude, longitude) %>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)

MW_metro_data <-readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"2006|2007|2008|2009|2010|2011"))) %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Rottnest|TwoRocks|Marmion|Warnbro|Two.Rocks"))) %>%
  dplyr::select(sample, campaign, latitude, longitude) %>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)

MW_jurien_data <- readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"2006|2007|2008|2009|2010|2011"))) %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Jurien"))) %>% 
  dplyr::select(sample, campaign, latitude, longitude) %>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)

MW_capes_data <- readRDS("australian-synthesis_complete_length_opthalmolepis_lineolatus.RDS") %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"Capes"))) %>% 
  dplyr::filter(str_detect(campaign, "2006|2007|2008|2009|2010|2011")) %>% 
  dplyr::select(sample, campaign, latitude, longitude) %>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)

BS_Data <- readRDS("australian-synthesis_complete_length_epinephelides-armatus.RDS") %>% 
  # filter(length<MaxLen) %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"2006|2007|2008|2009|2010|2011"))) %>% 
  dplyr::filter_all(.vars_predicate = any_vars(str_detect(.,"JBMP|Jurien|Marmion|TwoRocks|Two.Rocks|Capes|Geographe|SwC|south-west|Ngari|Warnbro|Rottnest"))) %>% 
  dplyr::select(sample, campaign, latitude, longitude) %>% 
  rename(x = longitude, y = latitude) %>% 
  distinct(sample, .keep_all=TRUE)

#* Load spatial data ####
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

#* Sort out metadata for Ningaloo ####
# setwd(data_dir)
# metadata_2015 <- metadata_2015 <- read.csv("2015-08_Ningaloo_stereoBRUVs_length.csv")
# metadata_2019 <- read.csv("ningaloo.complete.length.csv") %>% 
#   mutate(sample = as.character(sample)) %>% 
#   filter(genus %in% ("Lethrinus")) %>% 
#   filter(species %in% c("nebulosus"))
# dbca_metadata <- readRDS("DBCA_lengths.RDS") %>% 
#   mutate(sample = as.character(sample)) %>% 
#   filter(scientific %in% "Lethrinus nebulosus")
# metadata_2022 <- readRDS("Ningaloo_PtCloates_BOSS-BRUV_dat_length.rds")
# 
# SE_2019 <- SE_data %>% 
#   filter(campaignid %in% c("2019-08_Ningaloo.MP.Monitoring_stereoBRUVs", "2019-08_Ningaloo_stereo-BRUVs")) %>% 
#   left_join(., metadata_2019, by=c("sample")) %>% 
#   left_join(., dbca_metadata, by=c("sample")) %>% 
#   distinct(sample,.keep_all=TRUE) %>% 
#   mutate(latitude.x=ifelse(is.na(latitude.x), latitude.y, latitude.x),
#          longitude.x = ifelse(is.na(longitude.x), longitude.y, longitude.x)) %>% 
#   dplyr::select(sample, campaignid.x, latitude.x, longitude.x) %>% 
#   rename(campaign = "campaignid.x",
#          latitude = "latitude.x",
#          longitude = "longitude.x")
# 
# SE_2015 <- SE_data %>% 
#   filter(campaignid %in% c("2015-08_Ningaloo.deep.sanctuaries_stereoBRUVs")) %>% 
#   left_join(., metadata_2015, by=c("sample")) %>% 
#   distinct(sample,.keep_all=TRUE) %>% 
#   dplyr::select(sample, campaignid.x, latitude, longitude) %>% 
#   rename(campaign = "campaignid.x")
# 
# SE_2022 <- SE_data %>% 
#   filter(campaignid %in% c("2022-05_PtCloates_stereo-BRUVS", "2021-08_Pt-Cloates_stereo-BRUVs", "2021-08_PtCloates_BOSS")) %>% 
#   left_join(., metadata_2022, by=c("sample")) %>% 
#   distinct(sample,.keep_all=TRUE) %>% 
#   dplyr::select(sample, campaignid.x, latitude, longitude) %>% 
#   rename(campaign = "campaignid.x")
# 
# SE_data <- rbind(SE_2015, SE_2019, SE_2022)

#* Put data all together ####

# Full_samples <- rbind(SE_data, BC_data, PS_data, WKW_data, MW_data) %>% 
#   rename(x = longitude, y = latitude) %>% 
#   distinct(sample, .keep_all=TRUE)

#* plot heat map of samples ####

p1 <- ggplot() +
  geom_spatraster(data = bathy) +
  scale_fill_gradient(high = "#b9cded", low = "#09357d",  name = "Depth") +
  new_scale_fill() +
  stat_density2d(data = MW_capes_data, aes(x = x, y = y, fill = after_stat(level),
                                 alpha = after_stat(level)),
                 geom = "polygon", binwidth = 0.05, show.legend = F) + # , breaks = c(0, 0.05, 0.1, 0.99, 1)
  # scale_fill_viridis_c(option = "magma") +
  # scale_fill_gradient2(low = "#03051AFF", mid = "#F4875EFF", high = "#EC4B3E") +
  scale_fill_gradient2(low = "#f5ece1", mid = "#f2b46d", high = "#f78a05") +
  scale_alpha(range = c(0.5, 0.8), guide = "none") +   
  #geom_sf(data = pred_polys, fill = NA, colour = "grey90", linewidth = 0.5) +
  geom_sf(data = ausc, fill = "seashell2", colour = "black") +
  # geom_sf(data = ausc, fill = "gray80", colour = "black") +
  # geom_sf(data = npz, fill = NA, colour = "#7bbc63", linewidth = 0.2) +                          # Add national park zones
  # geom_sf(data = sanc, fill = NA, colour = "#bfd054", linewidth = 0.2) +
  # geom_point(data = hab, 
  #            aes(x = x, y = y),
  #            size = 0.25, alpha = 1) +
  # labs(x = "", y = "", title = "Benthic") +
  coord_sf(xlim = c(112, 118), #124.9889
           ylim = c(-36, -33)) + #-35.5, 21.9
  theme_minimal()+
  annotation_scale(location="tr")
p1

inset.map <- ggplot() +
  geom_sf(data = aus, fill = "seashell2", colour = "grey90", size = 0.075) +
  #geom_sf(data = aumpa, alpha = 5/6, colour = "grey85", size = 0.02) +
  coord_sf(xlim = c(108, 125), ylim = c(-37, -13)) +
  annotate(geom = "text", x=c(110.2), y=c(-29.4), label = c("Indian\nOcean"), size=3)+
  annotate(geom = "text", x=c(120.8), y=c(-25.94), label = c("Western\nAustralia"), size=3.5)+
  annotate("rect", xmin = 113, xmax = 114.35, ymin = -22.8, ymax = -21.5,   # Change here
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  annotate("rect", xmin = 113, xmax = 115.061, ymin = -29.25, ymax = -28,
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  annotate("rect", xmin=114.5, xmax=115.55, ymin=-34.5, ymax=-33.25,
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2)+
  annotate("rect", xmin=115, xmax=116, ymin=-33, ymax=-31,
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2)+
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "grey70"))+
  ylab(NULL)+
  xlab(NULL)
inset.map

sample.plot <- p1 + inset_element(inset.map, left = -1.29, right = 4.125, top = 1, bottom = 0.6)  
sample.plot

setwd(fig_dir)
ggsave(sample.plot, filename=paste0("Deployment_heatmap.png"),height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )




