# --------------------------------------------------------------------------------------- #
# - FILE NAME:   TenAnalyses.R         
# - DATE:        07/01/2021
# - DESCRIPTION: Supplementary analyses. 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(tidyverse)
library(brms)
library(data.table)
library(zoo)
library(dplyr)
library(tidybayes)
library(car)

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 

load(paste0(DataPath, "/FiveData.RData"))

# MODELLING #########################################################

# Set modelling parameters

iter = 8000
thin = 0.0005*iter
warmup = 0.1*iter

# Set priors

priors <- c(prior(normal(0, 1), class = b),
            prior(exponential(1), class = sigma))

# Transform to a data frame

pops_data <- as.data.frame(pops_data)

# Threats effects on mean population trend --------

# Combined stressors

m1 <- brm(mu | se(sigma, sigma = TRUE)  ~ threats-1 + (1|SpeciesName)+ (1|Realm), 
         iter = iter, thin = thin, warmup = warmup,
         control = list(adapt_delta = .975, max_treedepth = 20),
         data = pops_data, family = gaussian, prior = priors,
         cores = 18)

# System 

ms1 <- brm(mu | se(sigma, sigma = TRUE) ~ threats:System-1 + (1|SpeciesName)+ (1|Realm), 
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = pops_data, family = gaussian, prior = priors,
           cores = 18)

# Taxon 

mt1 <- brm(mu | se(sigma, sigma = TRUE) ~ threats:Taxon-1 + (1|SpeciesName)+ (1|Realm),
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
            data = pops_data, family = gaussian, prior = priors,
            cores = 18)

# Effect of protection

mp1 <- brm(mu | se(sigma, sigma = TRUE) ~ threats:Protected_status-1 + (1|SpeciesName)+ (1|Realm),
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = pops_data, family = gaussian, prior = priors,
           cores = 18)

# Effect of management

mm1 <- brm(mu | se(sigma, sigma = TRUE) ~ threats:as.factor(Managed)-1 + (1|SpeciesName)+ (1|Realm),
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = pops_data, family = gaussian, prior = priors,
           cores = 18)


# Subset only the interactive ones

pops_data <- pops_data %>% 
  filter(none == "None" | n.threat>1)

# Combined stressors ###########################################################

# System -----------------------------------------------------------------------

mmus <- brm(mu | se(sigma, sigma = TRUE)  ~ none:System + pollution:System + 
              habitatl:System+climatechange:System+ invasive:System+ 
              exploitation:System+disease:System-1 +(1|SpeciesName),  
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = pops_data, family = gaussian, prior = priors,
           cores = 18)

# Taxa -------------------------------------------------------------------------

mmut <- brm(mu | se(sigma, sigma = TRUE)  ~ none:Taxon + pollution:Taxon + 
              habitatl:Taxon+climatechange:Taxon+ invasive:Taxon+ exploitation:Taxon+
              disease:Taxon-1 +(1|SpeciesName),  
            iter = iter, thin = thin, warmup = warmup,
            control = list(adapt_delta = .975, max_treedepth = 20),
            data = pops_data, family = gaussian, prior = priors,
            cores = 18)

# Individual models ############################################################

## System ---------------------------------------------------------------------- 

# Marine 

mmu_mar <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                 pollution + habitatl+climatechange + invasive+ 
                 exploitation+disease-1 +(1|SpeciesName),  
               iter = iter, thin = thin, warmup = warmup,
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = pops_data %>% filter(System=="Marine"), 
               family = gaussian, prior = priors,
               cores = 18)

# Freshwater 

mmu_fr <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                pollution + habitatl+climatechange + invasive+ 
                exploitation+disease-1 +(1|SpeciesName),  
              iter = iter, thin = thin, warmup = warmup,
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = pops_data %>% filter(System=="Freshwater"), 
              family = gaussian, prior = priors,
              cores = 18)

# Terrestrial 

mmu_ter <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                 pollution + habitatl+climatechange + invasive+ 
                 exploitation+disease-1 +(1|SpeciesName),  
               iter = iter, thin = thin, warmup = warmup,
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = pops_data %>% filter(System=="Terrestrial"), 
               family = gaussian, prior = priors,
               cores = 18)

## Taxa -------------------------------------------------------------------------

# Amphibians 

mmu_am <- brm(mu | se(sigma, sigma = TRUE)  ~ none + 
                habitatl+ exploitation +
                disease-1 +(1|SpeciesName),  
              iter = iter, thin = thin, warmup = warmup,
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = pops_data %>% filter(Taxon=="Amphibians"), 
              family = gaussian, 
              prior = priors,
              cores = 18)

# Birds 

mmu_bi <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution + 
                habitatl+climatechange + invasive + exploitation +
                disease -1 +(1|SpeciesName),  
              iter = iter, thin = thin, warmup = warmup,
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = pops_data %>% filter(Taxon=="Birds"), 
              family = gaussian, 
              prior = priors,
              cores = 18)

# Fishes 

mmu_f <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution +  
               habitatl+climatechange + invasive -1 +(1|SpeciesName),  
             iter = iter, thin = thin, warmup = warmup,
             control = list(adapt_delta = .975, max_treedepth = 20),
             data = pops_data %>% filter(Taxon=="Fish"), 
             family = gaussian, 
             prior = priors,
             cores = 18)

# Mammals 

mmu_ma <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution + 
                habitatl+climatechange + invasive + exploitation +
                disease -1 +(1|SpeciesName),  
              iter = iter, thin = thin, warmup = warmup,
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = pops_data %>% filter(Taxon=="Mammals"), 
              family = gaussian, 
              prior = priors,
              cores = 18)

# Reptiles 

mmu_rep <- brm(mu | se(sigma, sigma = TRUE)  ~ none + pollution + 
                 habitatl+climatechange + invasive + exploitation +
                 disease -1 +(1|SpeciesName),  
               iter = iter, thin = thin, warmup = warmup,
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = pops_data %>% filter(Taxon=="Reptiles"), 
               family = gaussian, 
               prior = priors,
               cores = 18)

# Save the models --------------------------------------------------------------

setwd(ResultPath)
save(file = "Results.RData", 
     m1, ms1, mt1,
     mp1, mm1, 
     mmus, mmut,
     mmu_am, mmu_bi, mmu_f, 
     mmu_fr, mmu_ma, mmu_mar,
     mmu_rep,mmu_ter)  
