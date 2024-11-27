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
library(FSSgam)
library(mgcv)
library(MuMIn)
library(car)
library(doBy)
library(doSNOW)
library(dfadjust)

working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_data", sep="/")

dat <- data.frame(species=as.factor(c("L. nebulosus","E. armatus ", "C. auratus", "C. auricularis", "C. auricularis", "O. lineolatus", "O. lineolatus")),
                  mortality = c(0.27, 0.1, 0.23, 0.13, 0.32, 0.12, 0.088),
                  nat.mort = c(0.146, 0.3, 0.102, 0.41, 0.41, 0.36, 0.36),
                  upp_bnd = c(0.35, 0.034, 0.26, 0.23, 0.294, 0.24, 0.115),
                  lwr_bnd = c(0.2, 0.24, 0.2, 0.071, 0.337, 0.06, 0.067)) %>% 
  mutate(condition = 1-(mortality/(mortality+nat.mort))) %>% 
  mutate(condition_upp = 1-(lwr_bnd/(lwr_bnd+nat.mort))) %>% 
  mutate(condition_lwr = 1-(upp_bnd/(upp_bnd+nat.mort))) %>% 
  mutate(variance = (condition_lwr+condition_upp)/2) %>% 
  # group_by(species) %>% 
  # mutate(average = mean(condition)) %>% 
  # ungroup() %>% 
  distinct(species, .keep_all=T) %>% 
  filter(!species %in% "L. nebulosus") %>%
  mutate(max.length = c(500,1130,397,400),
         trophic.level = c(3.9,3.3,3.4,3.5)) %>% 
  as.data.frame()
  

Length_model <- lm(condition~max.length, data=dat, weights = 1/variance)
summary(Length_model)
SEadj <- dfadjustSE(Length_model)

SEint <- SEadj$coefficients[1,4]
SElength <- SEadj$coefficients[2,4]

coefs <- Length_model$coefficients

mean_length <- coefs[2]
intercept <- coefs[1]

lower_length <- coefs[2] - SElength
upper_length <- coefs[2] + SElength

lower_inter <- coefs[1] - SEint
upper_int <- coefs[1] + SEint

Trophic_model <- lm(condition~trophic.level, data=dat, weights = 1/variance)
summary(Trophic_model)

Null_model <- lm(condition~1, data=dat, weights = 1/variance)
summary(Null_model)


dat.new <- dat %>% 
  distinct(species, .keep_all=T)

testdata <- expand.grid(max.length=mean(dat.new$max.length)) %>%
  distinct()%>%
  glimpse()


mean_length_pred <- intercept+(testdata$max.length*mean_length)
upper_length_pred <- upper_int+(testdata$max.length*upper_length)
lower_length_pread <- lower_inter+(testdata$max.length*lower_length)




