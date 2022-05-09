# --------------------------------------------------------------------------------------- #
# - FILE NAME:   MultiEffects.R         
# - DATE:        15/03/2020
# - DESCRIPTION: Use results from the models to estimate wheter multiple
#               stressors are causing additive, synergistic or antagonistic
#               effects.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(ggplot2) #Now I load the latest version of ggplot2
library(tidybayes)
library(tidyverse)
library(dplyr)
library(brms)
library(bayesplot)
library(tidyr)
library(stringr)
library(cowplot)
library(ggforce)

# Set default ggplot theme

theme_set(theme_minimal()+
            theme(axis.title.x = element_text(size=15, margin = margin(t = 10, r = 0, b = 0, l = 0)), 
                  axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.line.x = element_line(color="black", size = 0.5),
                  axis.line.y = element_line(color="black", size = 0.5),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(color="black", size = 12),
                  axis.text.y = element_text(color="black", size = 12),
                  strip.text.x = element_text(size = 12),
                  axis.ticks = element_line(color="black"),
                  plot.margin = unit(c(0, 0, 0, 0), "cm")))


# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 

load(paste0(ResultPath, "/TenResults.RData"))
load(paste0(DataPath, "/TenYearsData.RData"))

# Change name invasive species 

pops_data10 <- pops_data10 %>% 
  mutate(threats = gsub("Invasive spp/genes", "Invasive", threats),)

# Estimate the multiple stressors effects in population trend ##################

# First we estimate the effect of each individual stressor
# We will set an arbitrary threshold of 95% cofindence

dat <- m1 %>%
  gather_draws(`b_.*`, regex = TRUE) %>% 
  mutate(.variable = gsub("b_threats", "", .variable), 
         .variable = gsub("InvasivesppDgenes", "Invasive", .variable),
         .variable = gsub("HabitatdegradationDchange", "Habitat degradation", .variable),
         .variable = gsub("Habitatloss", "Habitat loss", .variable),
         .variable = gsub("Climatechange", "Climate change", .variable),
         number = str_count(.variable,"[A-Z]")) %>% 
  filter(number==1) %>% 
  group_by(.variable) %>% 
  median_qi(.value, .width = c(.95))

dat

# Now we do the same but for more than one stressor
# we will filter according to the available ones from above

dat2 <- m1 %>%
  gather_draws(`b_.*`, regex = TRUE) %>% 
  mutate(.variable = gsub("b_threats", "", .variable), 
         .variable = gsub("InvasivesppDgenes", "Invasive", .variable),
         .variable = gsub("HabitatdegradationDchange", "Habitat degradation", .variable),
         .variable = gsub("Habitatloss", "Habitat loss", .variable),
         .variable = gsub("Climatechange", "Climate change", .variable),
         number = str_count(.variable,"[A-Z]")) %>% 
  filter(number>1) %>% 
  group_by(.variable) %>% 
  median_qi(.value, .width = c(.95))

dat2

#create a data frame

data_int <- expand.grid(var1=unique(dat$.variable), 
                        var2=unique(dat$.variable),
                        var3=NA,
                        value= NA,
                        lower=NA,
                        upper=NA,
                        .width=NA,
                        n.stress=2)

data_int <- data_int %>% mutate(.variable= paste0(var1, ":", 
                                                  var2))

#For loop for two interactions

count <- 0
for(i in 1:(length(dat$.variable))){
     for(j in 1:(length(dat$.variable))){
      count <- count+1
      data_int$value[count] <- dat$.value[i]+dat$.value[j]
      data_int$lower[count] <- dat$.lower[i]+dat$.lower[j]
      data_int$upper[count] <- dat$.upper[i]+dat$.upper[j]
      data_int$.width[count] <- dat$.width[i]
    }
  }

# Now we calculate the interaction among three stressors
# this is to set the expectations of how the additive effects should be

#create a data frame
data_int3 <- expand.grid(var1=unique(dat$.variable), 
                        var2=unique(dat$.variable),
                        var3=unique(dat$.variable),
                        value= NA,
                        lower=NA,
                        upper=NA,
                        .width=NA,
                        n.stress=3)

data_int3 <- data_int3 %>% mutate(.variable= paste0(var1, ":",
                                                    var2, ":", var3))

#For loop for three interactions
count <- 0
for(i in 1:(length(dat$.variable))){
  for(j in 1:(length(dat$.variable))){
    for(k in 1:(length(dat$.variable))){
    count <- count+1
    data_int3$value[count] <- dat$.value[i]+dat$.value[j]+dat$.value[k]
    data_int3$lower[count] <- dat$.lower[i]+dat$.lower[j]+dat$.lower[k]
    data_int3$upper[count] <- dat$.upper[i]+dat$.upper[j]+dat$.upper[k]
    data_int3$.width[count] <- dat$.width[i]
    }
    }
  }

# Join with the previous one 

data_int <- rbind(data_int, data_int3)

# Now we need to match the stressors with the current list of interactions
# Notice that for some the deduction of the low CI migh be higher than the 
# deduction from the upper CI, so we will correct for that

data_int2 <- data_int %>% 
  left_join(pops_data10 %>% dplyr::select(ID, SpeciesName, System, Taxon, threats, mu, uCI, lCI), 
            by=c(".variable"="threats")) %>%
  drop_na(ID) %>% 
  mutate(effect=#ifelse(uCI>0 & lCI<0,"Additive", 
           ifelse(upper<lCI, "Antagonistic", 
                  ifelse(uCI<lower,
                         "Synergy", "Additive")))

# Proportion of interaction types by System ####################################

# First, I am going to create separate each stressor interaction, according 
# to the presence of each one of them. 
# Then I am going to classify the interactions according to the values of the
# confidence intervals. 

data_prop <- data_int2 %>% 
  mutate(interaction.type=factor(effect, 
                                  levels = c("Additive", 
                                             "Synergy", "Antagonistic")),
         interaction.type=fct_relevel(interaction.type,"Antagonistic"),
    climate=str_extract(.variable, "Climate change"),
         exploitation=str_extract(.variable, "Exploitation"),
         habitatd=str_extract(.variable, "Habitat degradation"),
         habitatl=str_extract(.variable, "Habitat loss"),
         invasive=str_extract(.variable, "Invasive"),
         disease=str_extract(.variable, "Disease"),
         pollution=str_extract(.variable, "Pollution"),
    n.stress= str_count(.variable,"[A-Z]"))

# Now we are going to estimate the proportion of each interaction type 
# according to each group

climate_change <- data_prop %>% 
  drop_na(climate) %>%
  group_by(climate, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>% 
  rename(variable=climate)
  
habitat_loss <- data_prop %>% 
  drop_na(habitatl) %>%
  group_by(habitatl,System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=habitatl)

exploitation <- data_prop %>% 
  drop_na(exploitation) %>%
  group_by(exploitation, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=exploitation)

invasive <- data_prop %>% 
  drop_na(invasive) %>%
  group_by(invasive, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=invasive)

disease <- data_prop %>% 
  drop_na(disease) %>%
  group_by(disease, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=disease)

pollution <- data_prop %>% 
  drop_na(pollution) %>%
  group_by(pollution, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=pollution)

# Now we join them 

data_freq <- rbind(climate_change,
                   habitat_loss, exploitation, 
                   disease, invasive, 
                   pollution)

# Proportion by taxon ----------------------------------------------------------

climate_change <- data_prop %>% 
  drop_na(climate) %>%
  group_by(climate, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>% 
  rename(variable=climate)

habitat_loss <- data_prop %>% 
  drop_na(habitatl) %>%
  group_by(habitatl,Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=habitatl)

exploitation <- data_prop %>% 
  drop_na(exploitation) %>%
  group_by(exploitation, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=exploitation)

invasive <- data_prop %>% 
  drop_na(invasive) %>%
  group_by(invasive, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=invasive)

disease <- data_prop %>% 
  drop_na(disease) %>%
  group_by(disease, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=disease)

pollution <- data_prop %>% 
  drop_na(pollution) %>%
  group_by(pollution, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=pollution)

# Now we join them 

data_freq_taxon <- rbind(climate_change,
                         habitat_loss, exploitation, 
                         disease, invasive, 
                         pollution)

# Save 

setwd(ResultPath)
save(data_int2, data_prop, 
     data_freq, data_freq_taxon,
     file="PopMultEffects.RData")

