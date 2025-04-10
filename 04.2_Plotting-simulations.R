# ******************************************************
# Run a simulation of BRUV data to see how well we can
# retrieve the parameters 
# ******************************************************


library(L3Assess)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(gtable)
library(ggtext)
library(grid)

rm(list=ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_data", sep="/")
sim_dir <- paste(working.dir, "Simulations", sep="/")

line_colour <- "#BBCC33"
bruv_colour <- "#66CCEE"

a4.width=160

#### Read in data #####
setwd(sim_dir)

# WKW_sim <- readRDS("sim_sel_WKW.rds")
PS_sim <- readRDS("sim_sel_PS.rds")
BC_sim <- readRDS("sim_sel_BC.rds")
MW_sim <- readRDS("sim_sel_SMW.rds")
WKW_sim <- readRDS("sim_sel_WKW.rds")

#### Pink Snapper ####
# Pink Snapper will fall down because of their slower growth 
# Especially at higher Fs
# They grow so slowly so if you start fishing too many of hte big ones you basically wipe out that
# section of the length distribution and then can't fit anything/becomes highly uncertain 

NatMort = 0.12
Linf = 1136
FMort_true_low = 0.67*NatMort
FMort_true_high = 1.5*NatMort

Ret_delta = 50
Sel_true = 200

Ret_0.25_true = (Linf*0.25)-Ret_delta 
Ret_0.5_true = (Linf*0.5)-Ret_delta 
Ret_0.75_true = (Linf*0.75)-Ret_delta 
Ret_0.9_true = (Linf*0.9)-Ret_delta 

Fmort_plot <- c(FMort_true_low, FMort_true_low, FMort_true_high, FMort_true_high)

#Fmort_for_plot <- data.frame(Sel_Scen, Fmort_plot)


PS_sim <- PS_sim %>% 
  mutate(rel_err_L50_line = ifelse(Sel_Scen==1,(L50_line-Ret_0.25_true)/Ret_0.25_true, 
                               ifelse(Sel_Scen==2,  (L50_line-Ret_0.5_true)/Ret_0.5_true, 
                                      ifelse(Sel_Scen==3,(L50_line-Ret_0.75_true)/Ret_0.75_true,
                                             (L50_line-Ret_0.9_true)/Ret_0.9_true))),
         rel_err_L50_bruv = (L50_bruv-Sel_true)/Sel_true,
         rel_err_delta_line = (Delta_line-Ret_delta)/Ret_delta,
         rel_err_delta_bruv = (Delta_bruv-Ret_delta)/Ret_delta,
         rel_err_F_line = ifelse(F_scen %in% "Low", (F_line-FMort_true_low)/FMort_true_low,
                                 (F_line-FMort_true_high)/FMort_true_high),
         rel_err_F_bruv = ifelse(F_scen %in% "Low", (F_bruv-FMort_true_low)/FMort_true_low,
                                 (F_bruv-FMort_true_high)/FMort_true_high)) %>% 
  mutate(abs_err_L50_line = abs(ifelse(Sel_Scen==1,(L50_line-Ret_0.25_true)/Ret_0.25_true, 
                               ifelse(Sel_Scen==2,  (L50_line-Ret_0.5_true)/Ret_0.5_true, 
                                      ifelse(Sel_Scen==3,(L50_line-Ret_0.75_true)/Ret_0.75_true,
                                             (L50_line-Ret_0.9_true)/Ret_0.9_true)))),
         abs_err_L50_bruv = abs((L50_bruv-Sel_true)/Sel_true),
         abs_err_delta_line = abs((Delta_line-Ret_delta)/Ret_delta),
         abs_err_delta_bruv = abs((Delta_bruv-Ret_delta)/Ret_delta),
         abs_err_F_line = ifelse(F_scen %in% "Low", abs((F_line-FMort_true_low)/FMort_true_low),
                                 abs((F_line-FMort_true_high)/FMort_true_high)),
         abs_err_F_bruv = ifelse(F_scen %in% "Low", abs((F_bruv-FMort_true_low)/FMort_true_low),
                                 abs((F_bruv-FMort_true_high)/FMort_true_high))) %>% 
  filter(Sel_Scen>1) %>% 
  mutate(Sel_Scen = as.factor(Sel_Scen)) %>% 
  filter(convergence_line==0|convergence_bruv==0)

single_values_PS <- PS_sim %>% 
  group_by(F_scen, Sel_Scen) %>% 
  summarise(bruvL50MRE = median(rel_err_L50_bruv),
            bruvDeltaMRE = median(rel_err_delta_bruv),
            bruvFMRE = median(rel_err_F_bruv),
            lineL50MRE = median(rel_err_L50_line),
            lineDeltaMRE = median(rel_err_delta_line),
            lineFMRE = median(rel_err_F_line),
            bruvL50MARE = median(abs_err_L50_bruv),
            bruvDeltaMARE = median(abs_err_delta_bruv),
            bruvFMARE = median(abs_err_F_bruv),
            lineL50MARE = median(abs_err_L50_line),
            lineDeltaMARE = median(abs_err_delta_line),
            lineFMARE = median(abs_err_F_line))

#*Create violin plots ####

F_scen_names <- as_labeller(c("high"="High mortality", "Low"="Low mortality"))

bruv_F_MRE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=rel_err_F_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=bruvFMRE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), limits=c(-1,1))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
       # strip.text = element_text(size=5))
bruv_F_MRE_PS

line_F_MRE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=rel_err_F_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=lineFMRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), limits=c(NA,1))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic() +
  theme(strip.background = element_blank())
  #       strip.text = element_text(size=5))
line_F_MRE_PS

bruv_F_MARE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=abs_err_F_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=bruvFMARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(NA,1))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_F_MARE_PS

line_F_MARE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=abs_err_F_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=lineFMARE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(NA,1))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_F_MARE_PS  

# setwd(fig_dir)
# ggsave(line_F_MRE_PS, filename="MRE_line_PS.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(bruv_F_MRE_PS, filename="MRE_bruv_PS.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(line_F_MARE_PS, filename="MARE_line_PS.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(bruv_F_MARE_PS, filename="MARE_bruv_PS.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )


#### Breaksea Cod ####

NatMort = exp(1.551-1.066*log(20))
Linf = 560
FMort_true_low = 0.67*NatMort
FMort_true_high = 1.5*NatMort


Ret_delta = 50
Sel_true = 100

Ret_0.25_true = (Linf*0.25)-Ret_delta 
Ret_0.5_true = (Linf*0.5)-Ret_delta 
Ret_0.75_true = (Linf*0.75)-Ret_delta 
Ret_0.9_true = (Linf*0.9)-Ret_delta 


Sel_Scen <- c(1,2,3,4,1,2,3,4)
Fmort_plot <- c(FMort_true_low, FMort_true_low, FMort_true_high, FMort_true_high)

### Need to fix the selectivity scenarios
BC_sim<- BC_sim %>% 
  mutate(rel_err_L50_line = ifelse(Sel_Scen==1,(L50_line-Ret_0.25_true)/Ret_0.25_true, 
                                   ifelse(Sel_Scen==2,  (L50_line-Ret_0.5_true)/Ret_0.5_true, 
                                          ifelse(Sel_Scen==3,(L50_line-Ret_0.75_true)/Ret_0.75_true,
                                                 (L50_line-Ret_0.9_true)/Ret_0.9_true))),
         rel_err_L50_bruv = (L50_bruv-Sel_true)/Sel_true,
         rel_err_delta_line = (Delta_line-Ret_delta)/Ret_delta,
         rel_err_delta_bruv = (Delta_bruv-Ret_delta)/Ret_delta,
         rel_err_F_line = ifelse(F_scen %in% "Low", (F_line-FMort_true_low)/FMort_true_low,
                                 (F_line-FMort_true_high)/FMort_true_high),
         rel_err_F_bruv = ifelse(F_scen %in% "Low", (F_bruv-FMort_true_low)/FMort_true_low,
                                 (F_bruv-FMort_true_high)/FMort_true_high)) %>% 
  mutate(abs_err_L50_line = abs(ifelse(Sel_Scen==1,(L50_line-Ret_0.25_true)/Ret_0.25_true, 
                                       ifelse(Sel_Scen==2,  (L50_line-Ret_0.5_true)/Ret_0.5_true, 
                                              ifelse(Sel_Scen==3,(L50_line-Ret_0.75_true)/Ret_0.75_true,
                                                     (L50_line-Ret_0.9_true)/Ret_0.9_true)))),
         abs_err_L50_bruv = abs((L50_bruv-Sel_true)/Sel_true),
         abs_err_delta_line = abs((Delta_line-Ret_delta)/Ret_delta),
         abs_err_delta_bruv = abs((Delta_bruv-Ret_delta)/Ret_delta),
         abs_err_F_line = ifelse(F_scen %in% "Low", abs((F_line-FMort_true_low)/FMort_true_low),
                                 abs((F_line-FMort_true_high)/FMort_true_high)),
         abs_err_F_bruv = ifelse(F_scen %in% "Low", abs((F_bruv-FMort_true_low)/FMort_true_low),
                                 abs((F_bruv-FMort_true_high)/FMort_true_high))) %>% 
  filter(Sel_Scen>1) %>% 
  mutate(Sel_Scen = as.factor(Sel_Scen)) %>% 
  filter(convergence_line==0|convergence_bruv==0)

single_values_BC <- BC_sim %>% 
  group_by(F_scen, Sel_Scen) %>% 
  summarise(bruvL50MRE = median(rel_err_L50_bruv),
            bruvDeltaMRE = median(rel_err_delta_bruv),
            bruvFMRE = median(rel_err_F_bruv),
            lineL50MRE = median(rel_err_L50_line),
            lineDeltaMRE = median(rel_err_delta_line),
            lineFMRE = median(rel_err_F_line),
            bruvL50MARE = median(abs_err_L50_bruv),
            bruvDeltaMARE = median(abs_err_delta_bruv),
            bruvFMARE = median(abs_err_F_bruv),
            lineL50MARE = median(abs_err_L50_line),
            lineDeltaMARE = median(abs_err_delta_line),
            lineFMARE = median(abs_err_F_line))


#*Create violin plots ####
  
F_scen_names <- as_labeller(c("high"="High mortality", "Low"="Low mortality"))

bruv_F_MRE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=rel_err_F_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=bruvFMRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.25", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  geom_hline(yintercept=0, colour="grey20", linetype="dashed")+
  scale_y_continuous(breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), limits=c(-1,1))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_F_MRE_BC

line_F_MRE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=rel_err_F_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=lineFMRE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  geom_hline(yintercept=0, colour="grey20", linetype="dashed")+
  scale_y_continuous(breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), limits=c(-1,1))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_F_MRE_BC

bruv_F_MARE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=abs_err_F_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=bruvFMARE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(NA,1))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_F_MARE_BC

line_F_MARE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=abs_err_F_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=lineFMARE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(NA,1))+
  theme_classic()+
  theme(strip.background = element_blank())
line_F_MARE_BC  

# setwd(fig_dir)
# ggsave(line_F_MRE_BC, filename="MRE_line_BC.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(bruv_F_MRE_BC, filename="MRE_bruv_BC.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(line_F_MARE_BC, filename="MARE_line_BC.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(bruv_F_MARE_BC, filename="MARE_bruv_BC.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )

#### Southern Maori Wrasse ####
NatMort = exp(1.551-1.066*log(11.8))
Linf = 328
FMort_true_low = 0.67*NatMort
FMort_true_high = 1.5*NatMort

Ret_delta = 50
Sel_true = 100

Ret_0.25_true = (Linf*0.25)-Ret_delta 
Ret_0.5_true = (Linf*0.5)-Ret_delta 
Ret_0.75_true = (Linf*0.75)-Ret_delta 
Ret_0.9_true = (Linf*0.9)-Ret_delta 


Sel_Scen <- c(1,2,1,2)
Fmort_plot <- c(FMort_true_low, FMort_true_low, FMort_true_high, FMort_true_high)

MW_sim <- MW_sim %>% 
  mutate(rel_err_L50_line = ifelse(Sel_Scen==1,(L50_line-Ret_0.25_true)/Ret_0.25_true, 
                                   ifelse(Sel_Scen==2,  (L50_line-Ret_0.5_true)/Ret_0.5_true, 
                                          ifelse(Sel_Scen==3,(L50_line-Ret_0.75_true)/Ret_0.75_true,
                                                 (L50_line-Ret_0.9_true)/Ret_0.9_true))),
         rel_err_L50_bruv = (L50_bruv-Sel_true)/Sel_true,
         rel_err_delta_line = (Delta_line-Ret_delta)/Ret_delta,
         rel_err_delta_bruv = (Delta_bruv-Ret_delta)/Ret_delta,
         rel_err_F_line = ifelse(F_scen %in% "Low", (F_line-FMort_true_low)/FMort_true_low,
                                 (F_line-FMort_true_high)/FMort_true_high),
         rel_err_F_bruv = ifelse(F_scen %in% "Low", (F_bruv-FMort_true_low)/FMort_true_low,
                                 (F_bruv-FMort_true_high)/FMort_true_high)) %>% 
  mutate(abs_err_L50_line = abs(ifelse(Sel_Scen==1,(L50_line-Ret_0.25_true)/Ret_0.25_true, 
                                       ifelse(Sel_Scen==2,  (L50_line-Ret_0.5_true)/Ret_0.5_true, 
                                              ifelse(Sel_Scen==3,(L50_line-Ret_0.75_true)/Ret_0.75_true,
                                                     (L50_line-Ret_0.9_true)/Ret_0.9_true)))),
         abs_err_L50_bruv = abs((L50_bruv-Sel_true)/Sel_true),
         abs_err_delta_line = abs((Delta_line-Ret_delta)/Ret_delta),
         abs_err_delta_bruv = abs((Delta_bruv-Ret_delta)/Ret_delta),
         abs_err_F_line = ifelse(F_scen %in% "Low", abs((F_line-FMort_true_low)/FMort_true_low),
                                 abs((F_line-FMort_true_high)/FMort_true_high)),
         abs_err_F_bruv = ifelse(F_scen %in% "Low", abs((F_bruv-FMort_true_low)/FMort_true_low),
                                 abs((F_bruv-FMort_true_high)/FMort_true_high))) %>% 
  filter(Sel_Scen>1) %>% 
  mutate(Sel_Scen = as.factor(Sel_Scen)) %>% 
  filter(convergence_line==0|convergence_bruv==0)

single_values_MW <- MW_sim %>% 
  group_by(F_scen, Sel_Scen) %>% 
  summarise(bruvL50MRE = median(rel_err_L50_bruv),
            bruvDeltaMRE = median(rel_err_delta_bruv),
            bruvFMRE = median(rel_err_F_bruv),
            lineL50MRE = median(rel_err_L50_line),
            lineDeltaMRE = median(rel_err_delta_line),
            lineFMRE = median(rel_err_F_line),
            bruvL50MARE = median(abs_err_L50_bruv),
            bruvDeltaMARE = median(abs_err_delta_bruv),
            bruvFMARE = median(abs_err_F_bruv),
            lineL50MARE = median(abs_err_L50_line),
            lineDeltaMARE = median(abs_err_delta_line),
            lineFMARE = median(abs_err_F_line))


#*Create violin plots ####

F_scen_names <- as_labeller(c("high"="High mortality", "Low"="Low mortality"))

bruv_F_MRE_SMW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=rel_err_F_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=bruvFMRE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), limits=c(-1,1))+
  geom_hline(yintercept=0, colour="grey20", linetype="dashed")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_F_MRE_SMW

line_F_MRE_SMW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=rel_err_F_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=lineFMRE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), limits=c(-1,1))+
  geom_hline(yintercept=0, colour="grey20", linetype="dashed")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_F_MRE_SMW

bruv_F_MARE_SMW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=abs_err_F_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=bruvFMARE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(NA,1))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_F_MARE_SMW

line_F_MARE_SMW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=abs_err_F_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=lineFMARE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(NA,1))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_F_MARE_SMW  

# setwd(fig_dir)
# ggsave(line_F_MRE_SMW, filename="MRE_line_MW.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(bruv_F_MRE_SMW, filename="MRE_bruv_MW.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(line_F_MARE_SMW, filename="MARE_line_MW.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(bruv_F_MARE_SMW, filename="MARE_bruv_MW.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )

#### Western King Wrasse ####
NatMort = exp(1.551-1.066*log(10.2))
Linf = 378
FMort_true_low = 0.67*NatMort
FMort_true_high = 1.5*NatMort

Ret_delta = 50
Sel_true = 100

Ret_0.5_true = (Linf*0.5)-Ret_delta 
Ret_0.75_true = (Linf*0.75)-Ret_delta 
Ret_0.9_true = (Linf*0.9)-Ret_delta 


Sel_Scen <- c(1,2,1,2)
Fmort_plot <- c(FMort_true_low, FMort_true_low, FMort_true_high, FMort_true_high)

WKW_sim <- WKW_sim %>% 
  mutate(rel_err_L50_line = ifelse(Sel_Scen==1,(L50_line-Ret_0.5_true)/Ret_0.5_true, 
                                    ifelse(Sel_Scen==2,(L50_line-Ret_0.75_true)/Ret_0.75_true,
                                           (L50_line-Ret_0.9_true)/Ret_0.9_true)),
         rel_err_L50_bruv = (L50_bruv-Sel_true)/Sel_true,
         rel_err_delta_line = (Delta_line-Ret_delta)/Ret_delta,
         rel_err_delta_bruv = (Delta_bruv-Ret_delta)/Ret_delta,
         rel_err_F_line = ifelse(F_scen %in% "Low", (F_line-FMort_true_low)/FMort_true_low,
                                 (F_line-FMort_true_high)/FMort_true_high),
         rel_err_F_bruv = ifelse(F_scen %in% "Low", (F_bruv-FMort_true_low)/FMort_true_low,
                                 (F_bruv-FMort_true_high)/FMort_true_high)) %>% 
  mutate(abs_err_L50_line = abs(ifelse(Sel_Scen==1,(L50_line-Ret_0.5_true)/Ret_0.5_true, 
                                       ifelse(Sel_Scen==2,(L50_line-Ret_0.75_true)/Ret_0.75_true,
                                              (L50_line-Ret_0.9_true)/Ret_0.9_true))),
         abs_err_L50_bruv = abs((L50_bruv-Sel_true)/Sel_true),
         abs_err_delta_line = abs((Delta_line-Ret_delta)/Ret_delta),
         abs_err_delta_bruv = abs((Delta_bruv-Ret_delta)/Ret_delta),
         abs_err_F_line = ifelse(F_scen %in% "Low", abs((F_line-FMort_true_low)/FMort_true_low),
                                 abs((F_line-FMort_true_high)/FMort_true_high)),
         abs_err_F_bruv = ifelse(F_scen %in% "Low", abs((F_bruv-FMort_true_low)/FMort_true_low),
                                 abs((F_bruv-FMort_true_high)/FMort_true_high))) %>% 
  mutate(Sel_Scen = as.factor(Sel_Scen)) %>% 
  filter(convergence_line==0|convergence_bruv==0)

single_values_WKW <- WKW_sim %>% 
  group_by(F_scen, Sel_Scen) %>% 
  summarise(bruvL50MRE = median(rel_err_L50_bruv),
            bruvDeltaMRE = median(rel_err_delta_bruv),
            bruvFMRE = median(rel_err_F_bruv),
            lineL50MRE = median(rel_err_L50_line),
            lineDeltaMRE = median(rel_err_delta_line),
            lineFMRE = median(rel_err_F_line),
            bruvL50MARE = median(abs_err_L50_bruv),
            bruvDeltaMARE = median(abs_err_delta_bruv),
            bruvFMARE = median(abs_err_F_bruv),
            lineL50MARE = median(abs_err_L50_line),
            lineDeltaMARE = median(abs_err_delta_line),
            lineFMARE = median(abs_err_F_line))


#*Create violin plots ####

F_scen_names <- as_labeller(c("high"="High mortality", "Low"="Low mortality"))

bruv_F_MRE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=rel_err_F_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=bruvFMRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), limits=c(-1,1))+
  geom_hline(yintercept=0, colour="grey20", linetype="dashed")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_F_MRE_WKW

line_F_MRE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=rel_err_F_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=lineFMRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), limits=c(-1,1))+
  geom_hline(yintercept=0, colour="grey20", linetype="dashed")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_F_MRE_WKW

bruv_F_MARE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=abs_err_F_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=bruvFMARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(NA,1))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_F_MARE_WKW

line_F_MARE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=abs_err_F_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=lineFMARE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(NA,1))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_F_MARE_WKW  

# setwd(fig_dir)
# ggsave(line_F_MRE_WKW, filename="MRE_line_WKW.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(bruv_F_MRE_WKW, filename="MRE_bruv_WKW.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(line_F_MARE_WKW, filename="MARE_line_WKW.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )
# ggsave(bruv_F_MARE_WKW, filename="MARE_bruv_WKW.png", height=a4.width*0.5, width = a4.width*0.5, units  ="mm", dpi = 300 )


#### Save MRE and MARE Plots ####
setwd(fig_dir)
x.label <- textGrob(expression("Proportion of L"*infinity*" = to L"[95]), gp=gpar(fontsize=12))
y.label <- textGrob("Median relative error", gp=gpar(fontsize=12), rot=90)

a4.height = 250
PS_line <- data.frame(F_scen= c("high", "Low"),
                      label = c("a)", ""))
BC_line <- data.frame(F_scen= c("high", "Low"),
                      label = c("b)", ""))
WKW_line <- data.frame(F_scen= c("high", "Low"),
                      label = c("c)", ""))
SMW_line <- data.frame(F_scen= c("high", "Low"),
                      label = c("d)", ""))
PS_bruv <- data.frame(F_scen= c("high", "Low"),
                      label = c("e)", ""))
BC_bruv <- data.frame(F_scen= c("high", "Low"),
                      label = c("f)", ""))
WKW_bruv <- data.frame(F_scen= c("high", "Low"),
                      label = c("g)", ""))
SMW_bruv <- data.frame(F_scen= c("high", "Low"),
                      label = c("h)", ""))


MRE.F <- grid.arrange(arrangeGrob(line_F_MRE_PS + geom_text(PS_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                 line_F_MRE_BC + geom_text(BC_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                 line_F_MRE_WKW + geom_text(WKW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                 line_F_MRE_SMW + geom_text(SMW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                 bruv_F_MRE_PS + geom_text(PS_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                 bruv_F_MRE_BC + geom_text(BC_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                 bruv_F_MRE_WKW + geom_text(WKW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                 bruv_F_MRE_SMW + geom_text(SMW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                 left=y.label,
                                 bottom=x.label,
                                 ncol=4, nrow=2))

ggsave(MRE.F, filename="All_MRE.png",height = a4.width, width = a4.height, units  ="mm", dpi = 300 )

y.label <- textGrob("Median absolute relative error", gp=gpar(fontsize=12), rot=90)

MARE.F <- grid.arrange(arrangeGrob(line_F_MARE_PS + geom_text(PS_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  line_F_MARE_BC + geom_text(BC_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  line_F_MARE_WKW + geom_text(WKW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  line_F_MARE_SMW + geom_text(SMW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_F_MARE_PS + geom_text(PS_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_F_MARE_BC + geom_text(BC_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_F_MARE_WKW + geom_text(WKW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_F_MARE_SMW + geom_text(SMW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  left=y.label,
                                  bottom=x.label,
                                  ncol=4, nrow=2))

ggsave(MARE.F, filename="All_MARE.png",height = a4.width, width = a4.height, units  ="mm", dpi = 300 )


#### Pink Snapper Selectivity #####
#* L50
F_scen_names <- as_labeller(c("high"="High mortality", "Low"="Low mortality"))

bruv_L50_MRE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=rel_err_L50_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=bruvL50MRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.1,-0.05,0,0.05, 0.1), limits=c(-0.1,0.1))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
# strip.text = element_text(size=5))
bruv_L50_MRE_PS

line_L50_MRE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=rel_err_L50_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=lineL50MRE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.1,-0.05,0,0.05, 0.1), limits=c(-0.1,0.1))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic() +
  theme(strip.background = element_blank())
#       strip.text = element_text(size=5))
line_L50_MRE_PS

bruv_L50_MARE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=abs_err_L50_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=bruvL50MARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2), limits=c(NA,0.2))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_L50_MARE_PS

line_L50_MARE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=abs_err_L50_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=lineL50MARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2), limits=c(NA,0.2))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_L50_MARE_PS  


#* Delta
bruv_delta_MRE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=rel_err_delta_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=bruvDeltaMRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.6,-0.4,-0.2,0,0.2, 0.4, 0.6), limits=c(-0.6,0.6))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
# strip.text = element_text(size=5))
bruv_delta_MRE_PS

line_delta_MRE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=rel_err_delta_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=lineDeltaMRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.6,-0.4,-0.2,0,0.2, 0.4, 0.6), limits=c(-0.6,0.6))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic() +
  theme(strip.background = element_blank())
#       strip.text = element_text(size=5))
line_delta_MRE_PS

bruv_delta_MARE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=abs_err_delta_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=bruvDeltaMARE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8), limits=c(NA,0.8))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_delta_MARE_PS

line_delta_MARE_PS <- ggplot()+
  geom_violin(dat=PS_sim, aes(x=Sel_Scen, y=abs_err_delta_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_PS, aes(x=Sel_Scen, y=lineDeltaMARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8), limits=c(NA,0.8))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_delta_MARE_PS  



#### Breaksea Cod Selectivity #####
#* L50
F_scen_names <- as_labeller(c("high"="High mortality", "Low"="Low mortality"))

bruv_L50_MRE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=rel_err_L50_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=bruvL50MRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.1,-0.05,0,0.05, 0.1), limits=c(-0.1,0.1))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_L50_MRE_BC

line_L50_MRE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=rel_err_L50_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=lineL50MRE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.1,-0.05,0,0.05, 0.1), limits=c(-0.1,0.1))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic() +
  theme(strip.background = element_blank())
line_L50_MRE_BC

bruv_L50_MARE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=abs_err_L50_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=bruvL50MARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2), limits=c(NA,0.2))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_L50_MARE_BC

line_L50_MARE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=abs_err_L50_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=lineL50MARE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2), limits=c(NA,0.2))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_L50_MARE_BC  


#* Delta
bruv_delta_MRE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=rel_err_delta_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=bruvDeltaMRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.6,-0.4,-0.2,0,0.2, 0.4, 0.6), limits=c(-0.6,0.6))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_delta_MRE_BC

line_delta_MRE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=rel_err_delta_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=lineDeltaMRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.6,-0.4,-0.2,0,0.2, 0.4, 0.6), limits=c(-0.6,0.6))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic() +
  theme(strip.background = element_blank())
line_delta_MRE_BC

bruv_delta_MARE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=abs_err_delta_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=bruvDeltaMARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8), limits=c(NA,0.8))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_delta_MARE_BC

line_delta_MARE_BC <- ggplot()+
  geom_violin(dat=BC_sim, aes(x=Sel_Scen, y=abs_err_delta_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_BC, aes(x=Sel_Scen, y=lineDeltaMARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8), limits=c(NA,0.8))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_delta_MARE_BC  

#### Southern Maori Wrasse Selectivity #####
#* L50
F_scen_names <- as_labeller(c("high"="High mortality", "Low"="Low mortality"))

bruv_L50_MRE_MW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=rel_err_L50_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=bruvL50MRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.1,-0.05,0,0.05, 0.1), limits=c(-0.1,0.1))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_L50_MRE_MW

line_L50_MRE_MW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=rel_err_L50_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=lineL50MRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.1,-0.05,0,0.05, 0.1), limits=c(-0.1,0.1))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic() +
  theme(strip.background = element_blank())
line_L50_MRE_MW

bruv_L50_MARE_MW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=abs_err_L50_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=bruvL50MARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2), limits=c(NA,0.2))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_L50_MARE_MW

line_L50_MARE_MW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=abs_err_L50_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=lineL50MARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2), limits=c(NA,0.2))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_L50_MARE_MW  


#* Delta
bruv_delta_MRE_MW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=rel_err_delta_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=bruvDeltaMRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.6,-0.4,-0.2,0,0.2, 0.4, 0.6), limits=c(-0.6,0.6))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_delta_MRE_MW

line_delta_MRE_MW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=rel_err_delta_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=lineDeltaMRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.6,-0.4,-0.2,0,0.2, 0.4, 0.6), limits=c(-0.6,0.6))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic() +
  theme(strip.background = element_blank())
line_delta_MRE_MW

bruv_delta_MARE_MW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=abs_err_delta_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=bruvDeltaMARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8), limits=c(NA,0.8))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_delta_MARE_MW

line_delta_MARE_MW <- ggplot()+
  geom_violin(dat=MW_sim, aes(x=Sel_Scen, y=abs_err_delta_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_MW, aes(x=Sel_Scen, y=lineDeltaMARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8), limits=c(NA,0.8))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_delta_MARE_MW  

#### Western King Wrasse Selectivity #####
#* L50
F_scen_names <- as_labeller(c("high"="High mortality", "Low"="Low mortality"))

bruv_L50_MRE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=rel_err_L50_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=bruvL50MRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.1,-0.05,0,0.05, 0.1), limits=c(-0.1,0.1))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_L50_MRE_WKW

line_L50_MRE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=rel_err_L50_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=lineL50MRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.1,-0.05,0,0.05, 0.1), limits=c(-0.1,0.1))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic() +
  theme(strip.background = element_blank())
line_L50_MRE_WKW

bruv_L50_MARE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=abs_err_L50_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=bruvL50MARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2), limits=c(NA,0.2))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_L50_MARE_WKW

line_L50_MARE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=abs_err_L50_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=lineL50MARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2), limits=c(NA,0.2))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_L50_MARE_WKW  


#* Delta
bruv_delta_MRE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=rel_err_delta_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=bruvDeltaMRE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.6,-0.4,-0.2,0,0.2, 0.4, 0.6), limits=c(-0.6,0.6))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_delta_MRE_WKW

line_delta_MRE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=rel_err_delta_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=lineDeltaMRE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9",  "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(-0.6,-0.4,-0.2,0,0.2, 0.4, 0.6), limits=c(-0.6,0.6))+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20")+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic() +
  theme(strip.background = element_blank())
line_delta_MRE_WKW

bruv_delta_MARE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=abs_err_delta_bruv), fill=bruv_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=bruvDeltaMARE), size=0.5)+
  scale_x_discrete(labels = c( "0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8), limits=c(NA,0.8))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
bruv_delta_MARE_WKW

line_delta_MARE_WKW <- ggplot()+
  geom_violin(dat=WKW_sim, aes(x=Sel_Scen, y=abs_err_delta_line), fill=line_colour, linewidth=0.25)+
  geom_point(dat=single_values_WKW, aes(x=Sel_Scen, y=lineDeltaMARE), size=0.5)+
  scale_x_discrete(labels = c("0.5", "0.75", "0.9", "0.5", "0.75", "0.9"))+
  ylab(NULL)+
  xlab(NULL)+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8), limits=c(NA,0.8))+
  facet_wrap(~F_scen, labeller=F_scen_names)+
  theme_classic()+
  theme(strip.background = element_blank())
line_delta_MARE_WKW  


#### Save MRE and MARE Plots ####
setwd(fig_dir)
x.label <- textGrob(expression("Proportion of L"*infinity*" = to L"[95]), gp=gpar(fontsize=12))
y.label <- textGrob("Median relative error", gp=gpar(fontsize=12), rot=90)

a4.height = 250


MRE.L50 <- grid.arrange(arrangeGrob(line_L50_MRE_PS + geom_text(PS_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  line_L50_MRE_BC + geom_text(BC_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  line_L50_MRE_WKW + geom_text(WKW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  line_L50_MRE_MW + geom_text(SMW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_L50_MRE_PS + geom_text(PS_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_L50_MRE_BC + geom_text(BC_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_L50_MRE_WKW + geom_text(WKW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_L50_MRE_MW + geom_text(SMW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  left=y.label,
                                  bottom=x.label,
                                  ncol=4, nrow=2))

ggsave(MRE.L50, filename="All_MRE_L50.png",height = a4.width, width = a4.height, units  ="mm", dpi = 300 )

y.label <- textGrob("Median absolute relative error", gp=gpar(fontsize=12), rot=90)

MARE.L50 <-grid.arrange(arrangeGrob(line_L50_MARE_PS + geom_text(PS_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                    line_L50_MARE_BC + geom_text(BC_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                    line_L50_MARE_WKW + geom_text(WKW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                    line_L50_MARE_MW + geom_text(SMW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                    bruv_L50_MARE_PS + geom_text(PS_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                    bruv_L50_MARE_BC + geom_text(BC_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                    bruv_L50_MARE_WKW + geom_text(WKW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                    bruv_L50_MARE_MW + geom_text(SMW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                    left=y.label,
                                    bottom=x.label,
                                    ncol=4, nrow=2))


ggsave(MARE.L50, filename="All_MARE_L50.png",height = a4.width, width = a4.height, units  ="mm", dpi = 300 )


## Delta
y.label <- textGrob("Median relative error", gp=gpar(fontsize=12), rot=90)

a4.height = 250


MRE.delta <- grid.arrange(arrangeGrob(line_delta_MRE_PS + geom_text(PS_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  line_delta_MRE_BC + geom_text(BC_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  line_delta_MRE_WKW + geom_text(WKW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  line_delta_MRE_MW + geom_text(SMW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_delta_MRE_PS + geom_text(PS_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_delta_MRE_BC + geom_text(BC_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_delta_MRE_WKW + geom_text(WKW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  bruv_delta_MRE_MW + geom_text(SMW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                  left=y.label,
                                  bottom=x.label,
                                  ncol=4, nrow=2))

ggsave(MRE.delta, filename="All_MRE_delta.png",height = a4.width, width = a4.height, units  ="mm", dpi = 300 )

y.label <- textGrob("Median absolute relative error", gp=gpar(fontsize=12), rot=90)

MARE.delta <- grid.arrange(arrangeGrob(line_delta_MARE_PS + geom_text(PS_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                      line_delta_MARE_BC + geom_text(BC_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                      line_delta_MARE_WKW + geom_text(WKW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                      line_delta_MARE_MW + geom_text(SMW_line, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                      bruv_delta_MARE_PS + geom_text(PS_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                      bruv_delta_MARE_BC + geom_text(BC_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                      bruv_delta_MARE_WKW + geom_text(WKW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                      bruv_delta_MARE_MW + geom_text(SMW_bruv, mapping=aes(label=label, x=-Inf, y=Inf, hjust=-0.5, vjust=1.2)),
                                      left=y.label,
                                      bottom=x.label,
                                      ncol=4, nrow=2))

ggsave(MARE.delta, filename="All_MARE_delta.png",height = a4.width, width = a4.height, units  ="mm", dpi = 300 )




