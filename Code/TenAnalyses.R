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

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 

load(paste0(DataPath, "/TenYearsData.RData"))

# MODELLING #########################################################

# Set modelling parameters

iter = 8000
thin = 0.0005*iter
warmup = 0.1*iter

# Set priors

priors <- c(prior(normal(0, 1), class = b),
            prior(exponential(1), class = sigma))

# Transform to a data frame

pops_data10 <- as.data.frame(pops_data10)

# Threats effects on mean population trend --------

# Combined stressors

m1 <- brm(mu | se(tau, sigma = TRUE)  ~ threats-1 + (1|SpeciesName)+ (1|Realm), 
         iter = iter, thin = thin, warmup = warmup,
         control = list(adapt_delta = .975, max_treedepth = 20),
          data = pops_data10, family = gaussian, prior = priors,
          cores = 18)

# System 

ms1 <- brm(mu | se(tau, sigma = TRUE) ~ threats:System-1 + (1|SpeciesName)+ (1|Realm), 
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = pops_data10, family = gaussian, prior = priors,
           cores = 18)

# Taxon 

mt1 <- brm(mu | se(tau, sigma = TRUE) ~ threats:Taxon-1 + (1|SpeciesName)+ (1|Realm),
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
            data = pops_data10, family = gaussian, prior = priors,
            cores = 18)

# Presence absence model 

load(paste0(DataPath, "/LivingPlanetData.RData"))

# Keep only the threat data

dd_long <- dd_long %>% distinct(ID, .keep_all=T)

# Match the two databases including the presence absence of the disturbance 

pops_data10 <- pops_data10 %>% 
  left_join(dd_long[,c("ID", "pollution", "habitatl","climatechange", "habitatd", 
                       "invasive", "exploitation", "disease")])
pops_data10 <- pops_data10 %>% 
  mutate(none=if_else(threats=="None", "None", "NA"),
         none=gsub("NA", NA, none))

# Replace NAs for "no"

pops_data10 <- pops_data10 %>% replace_na(list(pollution = "No",
                                               habitatl = "No",
                                               climatechange = "No",
                                               habitatd = "No", 
                                               invasive="No", 
                                               exploitation="No",
                                               disease="No",
                                               none="No"))
# Correct habitat loss 

pops_data10 <- pops_data10 %>% mutate(habitatl=ifelse(habitatd!="No", 
                                                      "Habitat loss", habitatl))

# Relevel the factor

pops_data10 <- pops_data10 %>% 
  mutate(pollution=fct_relevel(pollution, "No"),
         habitatl=fct_relevel(habitatl, "No"),
         climatechange=fct_relevel(climatechange, "No"),
         habitatd=fct_relevel(habitatd, "No"),
         invasive=fct_relevel(invasive, "No"),
         exploitation=fct_relevel(exploitation, "No"),
         disease=fct_relevel(disease, "No"),
         none=fct_relevel(none, "No"))

pops_data10 <- pops_data10 %>% 
  mutate(habitatl=fct_relevel(habitatl, "No"),
         invasive=fct_relevel(invasive, "No"))

# Subset only the interactive ones

pops_data10 <- pops_data10 %>% 
  filter(none == "None" | n.threat>1)

# Combined stressors ###########################################################

# System -----------------------------------------------------------------------

mmus <- brm(mu | se(tau, sigma = TRUE)  ~ none:System + pollution:System + 
              habitatl:System+climatechange:System+ invasive:System+ 
              exploitation:System+disease:System-1 +(1|SpeciesName),  
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = pops_data10, family = gaussian, prior = priors,
           cores = 18)

# Taxa -------------------------------------------------------------------------

mmut <- brm(mu | se(tau, sigma = TRUE)  ~ none:Taxon + pollution:Taxon + 
              habitatl:Taxon+climatechange:Taxon+ invasive:Taxon+ exploitation:Taxon+
              disease:Taxon-1 +(1|SpeciesName),  
            iter = iter, thin = thin, warmup = warmup,
            control = list(adapt_delta = .975, max_treedepth = 20),
            data = pops_data10, family = gaussian, prior = priors,
            cores = 18)


# Save the models --------------------------------------------------------------

setwd(ResultPath)
save(file = "TenResults.RData", 
     m1, ms1, mt1,
     mp1, mm1, 
     mmus, mmut)  
