# --------------------------------------------------------------------------------------- #
# - FILE NAME:   AnalysisData.R         
# - DATE:        07/01/2021
# - DESCRIPTION: Prepare data for the analyses
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(tidyverse)
library(data.table)
library(zoo)
library(dplyr)
library(naniar)

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 

load(paste0(DataPath, "/SSResults.RData"))
load(paste0(DataPath, "/LivingPlanetData2.RData"))

# Classify taxonomic groups into major ones

pops_data <- pops_data %>% 
  mutate(Taxon= ifelse(Class=="Holocephali"|Class=="Elasmobranchii" | 
                         Class=="Myxini"|Class=="Cephalaspidomorphi"|
                         Class=="Actinopterygii"|Class=="Sarcopterygii",
                       "Fish", 
                       ifelse(Class=="Aves", "Birds",
                              ifelse(Class=="Mammalia", 
                                     "Mammals",
                                     ifelse(Class=="Amphibia", 
                                            "Amphibians",
                                            ifelse(Class=="Reptilia",
                                                   "Reptiles",Class))))))

# Keep only the threat data

dd_long <- dd_long %>% distinct(ID, .keep_all=T)

# Remove no threats

pops_data <- pops_data %>% 
  mutate(threats=replace_na(threats, "None")) %>% 
  filter(!is.na(threats))

pop_data <- pop_data %>% 
  mutate(threats=replace_na(threats, "None")) %>% 
  filter(!is.na(threats))

# Match the two databases including the presence absence of the disturbance 

pops_data <- pops_data %>% 
  left_join(dd_long[,c("ID", "pollution", "habitatl","climatechange", "habitatd", 
                       "invasive", "exploitation", "disease")])
pops_data <- pops_data %>% 
  mutate(none=if_else(threats=="None", "None", "NA"),
         none=gsub("NA", NA, none))

# Replace NAs for "no"

pops_data <- pops_data %>% replace_na(list(pollution = "No",
                                           habitatl = "No",
                                           climatechange = "No",
                                           habitatd = "No", 
                                           invasive="No", 
                                           exploitation="No",
                                           disease="No",
                                           none="No"))
# Correct habitat loss 

pops_data <- pops_data %>% mutate(habitatl=ifelse(habitatd!="No", 
                                                  "Habitat loss", habitatl))

# Relevel the factor

pops_data <- pops_data %>% 
  mutate(pollution=fct_relevel(pollution, "No"),
         habitatl=fct_relevel(habitatl, "No"),
         climatechange=fct_relevel(climatechange, "No"),
         habitatd=fct_relevel(habitatd, "No"),
         invasive=fct_relevel(invasive, "No"),
         exploitation=fct_relevel(exploitation, "No"),
         disease=fct_relevel(disease, "No"),
         none=fct_relevel(none, "No"))

pops_data <- pops_data %>% 
  mutate(habitatl=fct_relevel(habitatl, "No"),
         invasive=fct_relevel(invasive, "No"))

# Save the data 

setwd(DataPath)
save(pops_data, pop_data,
     file = "FiveData.RData")
