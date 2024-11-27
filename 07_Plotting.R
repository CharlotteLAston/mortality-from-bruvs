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

rm(list = ls())

working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_data", sep="/")

a4.width=160

#### Read in data ####
setwd(data_dir)

dat <- read.csv("catch_curve_results.csv")

#### Make data frame ####

dat <- dat %>% 
  mutate(Species = ifelse(Species %in% "King Wrasse", "Western King\nWrasse", Species)) %>% 
  mutate(species.location = paste0(Species, sep="_", Location)) %>% 
  mutate(species.location.method = paste0(species.location, sep="_", Method)) %>% 
  mutate(species.location.method.data = paste0(species.location.method, sep="_", Data)) %>% 
  mutate(Species=fct_relevel(Species, "Pink Snapper","Breaksea Cod","Western King\nWrasse","Maori Wrasse")) %>% 
  mutate(species.location=fct_relevel(species.location, "Pink Snapper_Metro/SW","Breaksea Cod_Metro","Western King\nWrasse_Metro","Western King\nWrasse_Jurien", "Maori Wrasse_Metro")) %>% 
  mutate(Shape = c(16,16,17,17,16,16,17,17,16,16,17,17,16,16,17,17,16,16,17,17))

#colours <- c("#88CBED", "#A9439A", "#332387", "#117633", "#43A999", "#872155", "#CB6778")
# species.labels <- c("*L. nebulosus*<br>(Ningaloo)", "*L. nebulosus*&#42;<br>(Ningaloo)","*E. armatus*<br>(Metro<br>and SW)", "*C. auratus*<br>(Metro<br>and SW)", "*C. auratus*&#42;<br>(Metro<br>and SW)",
#                     "*C. auricularis*<br>(Metro)", "*O. lineolatus*<br>(Metro)", "*C. auricularis*<br>(Abrolhos)", "*O. lineolatus*<br>(SW)")

facet_labels <- c("Pink Snapper_Metro/SW"="Pink Snapper\n(Metro/South-west)", "Breaksea Cod_Metro"="Breaksea Cod\n(Metro)",
                  "Western King\nWrasse_Metro"="Western King Wrasse\n(Metro)", "Western King\nWrasse_Jurien"="Western King Wrasse\n(Jurien Bay)", "Maori Wrasse_Metro"="Maori Wrasse\n(Metro)")
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
        strip.placement = "outside")
  #theme(axis.text.x=ggtext::element_markdown())
condition_plot

setwd(fig_dir)
ggsave(condition_plot, filename="Condition_plot.png", height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300 )
