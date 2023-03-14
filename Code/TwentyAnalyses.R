# --------------------------------------------------------------------------------------- #
# - FILE NAME:   TwentyAnalyses.R         
# - DATE:        07/01/2023
# - DESCRIPTION: Supplementary analyses to test the sensitivity of our results
#                when increasing the number of 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(ggplot2)
library(tidybayes)
library(tidyverse)
library(dplyr)
library(rstan)
library(rstanarm)
library(brms)
library(bayesplot)
library(bayestestR)
library(ggdist)
library(cowplot)
library(wesanderson)
library(tidyr)
library(ggalt)
library(ggmap)
library(ggforce)
library(ggridges)
library(patchwork)
library(ggalluvial)
library(ggimage)
library(ggpubr)
library(stringr)
library(forcats)
library(igraph)
library(MetBrewer)

# Set default ggplot theme

theme_set(theme_minimal()+
            theme(axis.title.x = element_text(size=12,
                                              margin = margin(t = 10, r = 0, b = 0, l = 0)), 
                  axis.title.y = element_text(size=12,
                                              margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.line.x = element_line(color="black", linewidth =  0.5),
                  axis.line.y = element_line(color="black", linewidth = 0.5),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(color="black", size = 12),
                  axis.text.y = element_text(color="black", size = 12),
                  strip.text.x = element_text(size = 12),
                  axis.ticks = element_line(color="black"),
                  plot.title = element_text(hjust = 0.5)))

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 

load(paste0(DataPath, "/FiveData.RData"))

# Subset the 10 data points -----------------------------------------------------

pop_data10 <- pop_data %>% 
  mutate(Year=as.numeric(Year)) %>% 
  group_by(ID) %>% 
  filter(length(unique(Year))>=10)

pops_data10 <- pops_data %>% 
  filter(ID%in%pop_data10$ID)

# MODELLING #########################################################

# Set modelling parameters

iter = 8000
thin = 0.0005*iter
warmup = 0.1*iter

# Set priors

priors <- c(prior(normal(0, 1), class = b),
            prior(exponential(1), class = sigma))

# Threats effects on mean population trend --------

# Combined stressors

m2 <- brm(mu | se(sigma, sigma = TRUE)  ~ threats-1 + (1|SpeciesName)+ (1|Realm), 
          iter = iter, thin = thin, warmup = warmup,
          control = list(adapt_delta = .975, max_treedepth = 20),
          data = pops_data10, family = gaussian, prior = priors,
          cores = 18)

# System 

ms2 <- brm(mu | se(sigma, sigma = TRUE) ~ threats:System-1 + (1|SpeciesName)+ (1|Realm), 
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = pops_data10, family = gaussian, prior = priors,
           cores = 18)

# Taxon 

mt2 <- brm(mu | se(sigma, sigma = TRUE) ~ threats:Taxon-1 + (1|SpeciesName)+ (1|Realm),
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = pops_data10, family = gaussian, prior = priors,
           cores = 18)

# Subset only the interactive ones

pops_data10 <- pops_data10 %>% 
  filter(none == "None" | n.threat>1)

# Individual models ############################################################

## System ---------------------------------------------------------------------- 

# Marine 

mmu_mar2 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                 pollution + habitatl+climatechange + invasive+ 
                 exploitation+disease-1 +(1|SpeciesName),  
               iter = iter, thin = thin, warmup = warmup,
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = pops_data10 %>% filter(System=="Marine"), 
               family = gaussian, prior = priors,
               cores = 18)

# Freshwater 

mmu_fr2 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                pollution + habitatl+climatechange + invasive+ 
                exploitation+disease-1 +(1|SpeciesName),  
              iter = iter, thin = thin, warmup = warmup,
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = pops_data10 %>% filter(System=="Freshwater"), 
              family = gaussian, prior = priors,
              cores = 18)

# Terrestrial 

mmu_ter2 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                 pollution + habitatl+climatechange + invasive+ 
                 exploitation+disease-1 +(1|SpeciesName),  
               iter = iter, thin = thin, warmup = warmup,
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = pops_data10 %>% filter(System=="Terrestrial"), 
               family = gaussian, prior = priors,
               cores = 18)

## Taxa -------------------------------------------------------------------------

# Amphibians 

mmu_am2 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                habitatl+ exploitation +
                disease-1 +(1|SpeciesName),  
              iter = iter, thin = thin, warmup = warmup,
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = pops_data10 %>% filter(Taxon=="Amphibians"), 
              family = gaussian, 
              prior = priors,
              cores = 18)

# Birds 

mmu_bi2 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution + 
                habitatl+climatechange + invasive + exploitation +
                disease -1 +(1|SpeciesName),  
              iter = iter, thin = thin, warmup = warmup,
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = pops_data10 %>% filter(Taxon=="Birds"), 
              family = gaussian, 
              prior = priors,
              cores = 18)

# Fishes 

mmu_f2 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution +  
               habitatl+climatechange + invasive -1 +(1|SpeciesName),  
             iter = iter, thin = thin, warmup = warmup,
             control = list(adapt_delta = .975, max_treedepth = 20),
             data = pops_data10 %>% filter(Taxon=="Fish"), 
             family = gaussian, 
             prior = priors,
             cores = 18)

# Mammals 

mmu_ma2 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution + 
                habitatl+climatechange + invasive + exploitation +
                disease -1 +(1|SpeciesName),  
              iter = iter, thin = thin, warmup = warmup,
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = pops_data10 %>% filter(Taxon=="Mammals"), 
              family = gaussian, 
              prior = priors,
              cores = 18)

# Reptiles 

mmu_rep2 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution + 
                 habitatl+climatechange + invasive + exploitation +
                 disease -1 +(1|SpeciesName),  
               iter = iter, thin = thin, warmup = warmup,
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = pops_data10 %>% filter(Taxon=="Reptiles"), 
               family = gaussian, 
               prior = priors,
               cores = 18)


# Save the models --------------------------------------------------------------

setwd(ResultPath)
save(file = "TenResults.RData", 
     m2, ms2, mt2,
     mmu_am2, mmu_bi2, mmu_f2, 
     mmu_fr2, mmu_ma2, mmu_mar2,
     mmu_rep2, mmu_ter2)  

# Subset the 20 years data -----------------------------------------------------

pop_data20 <- pop_data %>% 
  mutate(Year=as.numeric(Year)) %>% 
  group_by(ID) %>% 
  filter(length(unique(Year))>=20)

pops_data20 <- pops_data %>% 
  filter(ID%in%pop_data20$ID)

# Threats effects on mean population trend --------

# Combined stressors

m3 <- brm(mu | se(sigma, sigma = TRUE)  ~ threats-1 + (1|SpeciesName)+ (1|Realm), 
          iter = iter, thin = thin, warmup = warmup,
          control = list(adapt_delta = .975, max_treedepth = 20),
          data = pops_data20, family = gaussian, prior = priors,
          cores = 18)

# System 

ms3 <- brm(mu | se(sigma, sigma = TRUE) ~ threats:System-1 + (1|SpeciesName)+ (1|Realm), 
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = pops_data20, family = gaussian, prior = priors,
           cores = 18)

# Taxon 

mt3 <- brm(mu | se(sigma, sigma = TRUE) ~ threats:Taxon-1 + (1|SpeciesName)+ (1|Realm),
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = pops_data20, family = gaussian, prior = priors,
           cores = 18)

# Subset only the interactive ones

pops_data20 <- pops_data20 %>% 
  filter(none == "None" | n.threat>1)

# Individual models ############################################################

## System ---------------------------------------------------------------------- 

# Marine 

mmu_mar3 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                  pollution + habitatl+climatechange + invasive+ 
                  exploitation+disease-1 +(1|SpeciesName),  
                iter = iter, thin = thin, warmup = warmup,
                control = list(adapt_delta = .975, max_treedepth = 20),
                data = pops_data10 %>% filter(System=="Marine"), 
                family = gaussian, prior = priors,
                cores = 18)

# Freshwater 

mmu_fr3 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                 pollution + habitatl+climatechange + invasive+ 
                 exploitation+disease-1 +(1|SpeciesName),  
               iter = iter, thin = thin, warmup = warmup,
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = pops_data10 %>% filter(System=="Freshwater"), 
               family = gaussian, prior = priors,
               cores = 18)

# Terrestrial 

mmu_ter3 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                  pollution + habitatl+climatechange + invasive+ 
                  exploitation+disease-1 +(1|SpeciesName),  
                iter = iter, thin = thin, warmup = warmup,
                control = list(adapt_delta = .975, max_treedepth = 20),
                data = pops_data10 %>% filter(System=="Terrestrial"), 
                family = gaussian, prior = priors,
                cores = 18)

## Taxa -------------------------------------------------------------------------

# Amphibians 

mmu_am3 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                 habitatl+ exploitation +
                 disease-1 +(1|SpeciesName),  
               iter = iter, thin = thin, warmup = warmup,
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = pops_data10 %>% filter(Taxon=="Amphibians"), 
               family = gaussian, 
               prior = priors,
               cores = 18)

# Birds 

mmu_bi3 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution + 
                 habitatl+climatechange + invasive + exploitation +
                 disease -1 +(1|SpeciesName),  
               iter = iter, thin = thin, warmup = warmup,
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = pops_data10 %>% filter(Taxon=="Birds"), 
               family = gaussian, 
               prior = priors,
               cores = 18)

# Fishes 

mmu_f3 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution +  
                habitatl+climatechange + invasive -1 +(1|SpeciesName),  
              iter = iter, thin = thin, warmup = warmup,
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = pops_data10 %>% filter(Taxon=="Fish"), 
              family = gaussian, 
              prior = priors,
              cores = 18)

# Mammals 

mmu_ma3 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution + 
                 habitatl+climatechange + invasive + exploitation +
                 disease -1 +(1|SpeciesName),  
               iter = iter, thin = thin, warmup = warmup,
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = pops_data10 %>% filter(Taxon=="Mammals"), 
               family = gaussian, 
               prior = priors,
               cores = 18)

# Reptiles 

mmu_rep3 <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution + 
                  habitatl+climatechange + invasive + exploitation +
                  disease -1 +(1|SpeciesName),  
                iter = iter, thin = thin, warmup = warmup,
                control = list(adapt_delta = .975, max_treedepth = 20),
                data = pops_data10 %>% filter(Taxon=="Reptiles"), 
                family = gaussian, 
                prior = priors,
                cores = 18)

# Save the models --------------------------------------------------------------

setwd(ResultPath)
save(file = "TwentyResults.RData", 
     m3, ms3, mt3,
     mmu_am3, mmu_bi3, mmu_f3, 
     mmu_fr3, mmu_ma3, mmu_mar3,
     mmu_rep3, mmu_ter3) 

