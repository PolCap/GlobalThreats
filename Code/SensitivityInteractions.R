# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SensitivityInteractions.R         
# - DATE:        07/01/2023
# - DESCRIPTION: Test the sensitivity of the models to recalculate the interactions
#                based on the mean rather than the CI.
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
                  axis.line.x = element_line(color="black", linewidth = 0.5),
                  axis.line.y = element_line(color="black", linewidth =  0.5),
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

load(paste0(ResultPath, "/Results.RData"))
load(paste0(DataPath, "/FiveData.RData"))

# Change name invasive species 

pops_data<- pops_data%>% 
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
  left_join(pops_data%>% dplyr::select(ID, SpeciesName, System, Taxon, threats, mu, uCI, lCI), 
            by=c(".variable"="threats")) %>%
  drop_na(ID) %>% 
  mutate(effect=#ifelse(uCI>0 & lCI<0,"Additive", 
           ifelse(value<lCI, "Antagonistic", 
                  ifelse(uCI<value,
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
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>% 
  rename(variable=climate)

habitat_loss <- data_prop %>% 
  drop_na(habitatl) %>%
  group_by(habitatl,System, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=habitatl)

exploitation <- data_prop %>% 
  drop_na(exploitation) %>%
  group_by(exploitation, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=exploitation)

invasive <- data_prop %>% 
  drop_na(invasive) %>%
  group_by(invasive, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=invasive)

disease <- data_prop %>% 
  drop_na(disease) %>%
  group_by(disease, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=disease)

pollution <- data_prop %>% 
  drop_na(pollution) %>%
  group_by(pollution, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
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
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>% 
  rename(variable=climate)

habitat_loss <- data_prop %>% 
  drop_na(habitatl) %>%
  group_by(habitatl,Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=habitatl)

exploitation <- data_prop %>% 
  drop_na(exploitation) %>%
  group_by(exploitation, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=exploitation)

invasive <- data_prop %>% 
  drop_na(invasive) %>%
  group_by(invasive, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=invasive)

disease <- data_prop %>% 
  drop_na(disease) %>%
  group_by(disease, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=disease)

pollution <- data_prop %>% 
  drop_na(pollution) %>%
  group_by(pollution, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
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
     file="PopMeanEffects.RData")

# Multiplicative null model ####################################################

# First we estimate the effect of each individual stressor

dat <- m1 %>%
  gather_draws(`b_.*`, regex = TRUE) %>% 
  mutate(.variable = gsub("b_threats", "", .variable), 
         .variable = gsub("InvasivesppDgenes", "Invasive", .variable),
         .variable = gsub("HabitatdegradationDchange", "Habitat degradation", 
                          .variable),
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

# Create a data frame

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
    data_int$value[count] <- dat$.value[i]+dat$.value[j]-(dat$.value[i]*dat$.value[j])
    data_int$lower[count] <- dat$.lower[i]+dat$.lower[j]-(dat$.lower[i]*dat$.lower[j])
    data_int$upper[count] <- dat$.upper[i]+dat$.upper[j]-(dat$.upper[i]*dat$.upper[j])
    data_int$.width[count] <- dat$.width[i]
  }
}

# Now we calculate the interaction among three stressors
# this is to set the expectations of how the additive effects should be

# create a data frame

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
      data_int3$value[count] <- dat$.value[i]+dat$.value[j]+dat$.value[k]-(dat$.value[i]*dat$.value[j]*dat$.value[k])
      data_int3$lower[count] <- dat$.lower[i]+dat$.lower[j]+dat$.lower[k]-(dat$.lower[i]*dat$.lower[j]*dat$.lower[k])
      data_int3$upper[count] <- dat$.upper[i]+dat$.upper[j]+dat$.upper[k]-(dat$.upper[i]*dat$.upper[j]*dat$.upper[k])
      data_int3$.width[count] <- dat$.width[i]
    }
  }
}

# Join with the previous one 

data_int <- rbind(data_int, data_int3)

# Using mu from the state-space estimation -------------------------------------

data_int_m <- data_int %>% 
  left_join(pops_data %>% dplyr::select(ID, SpeciesName, System, Taxon, threats, mu, uCI, lCI), 
            by=c(".variable"="threats")) %>%
  drop_na(ID) %>% 
  mutate(effect=ifelse(upper<lCI, "Antagonistic", 
                       ifelse(uCI<lower,
                              "Synergy", "Additive")))

# Proportion of interaction types by System ####################################

# First, I am going to create separate each stressor interaction, according 
# to the presence of each one of them. 
# Then I am going to classify the interactions according to the values of the
# confidence intervals. 

data_prop_m <- data_int_m %>% 
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

climate_change <- data_prop_m %>% 
  drop_na(climate) %>%
  group_by(climate, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>% 
  rename(variable=climate)

habitat_loss <- data_prop_m %>% 
  drop_na(habitatl) %>%
  group_by(habitatl,System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=habitatl)

exploitation <- data_prop_m %>% 
  drop_na(exploitation) %>%
  group_by(exploitation, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=exploitation)

invasive <- data_prop_m %>% 
  drop_na(invasive) %>%
  group_by(invasive, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=invasive)

disease <- data_prop_m %>% 
  drop_na(disease) %>%
  group_by(disease, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=disease)

pollution <- data_prop_m %>% 
  drop_na(pollution) %>%
  group_by(pollution, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=pollution)

# Now we join them 

data_freq_m <- rbind(climate_change,
                     habitat_loss, exploitation, 
                     disease, invasive, 
                     pollution)

# Proportion by taxon ----------------------------------------------------------

climate_change <- data_prop_m %>% 
  drop_na(climate) %>%
  group_by(climate, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>% 
  rename(variable=climate)

habitat_loss <- data_prop_m %>% 
  drop_na(habitatl) %>%
  group_by(habitatl,Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=habitatl)

exploitation <- data_prop_m %>% 
  drop_na(exploitation) %>%
  group_by(exploitation, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=exploitation)

invasive <- data_prop_m %>% 
  drop_na(invasive) %>%
  group_by(invasive, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=invasive)

disease <- data_prop_m %>% 
  drop_na(disease) %>%
  group_by(disease, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=disease)

pollution <- data_prop_m %>% 
  drop_na(pollution) %>%
  group_by(pollution, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=pollution)

# Now we join them 

data_freq_m_taxon <- rbind(climate_change,
                           habitat_loss, exploitation, 
                           disease, invasive, 
                           pollution)

# Dominant null model ##########################################################

# Create a data frame

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
    data_int$value[count] <- max(dat$.value[i],dat$.value[j])
    data_int$lower[count] <- max(dat$.lower[i],dat$.lower[j])
    data_int$upper[count] <- max(dat$.upper[i],dat$.upper[j])
    data_int$.width[count] <- dat$.width[i]
  }
}

# Now we calculate the interaction among three stressors
# this is to set the expectations of how the additive effects should be

# create a data frame

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
      data_int3$value[count] <- max(dat$.value[i],dat$.value[j],dat$.value[k])
      data_int3$lower[count] <- max(dat$.lower[i],dat$.lower[j],dat$.lower[k])
      data_int3$upper[count] <- max(dat$.upper[i],dat$.upper[j],dat$.upper[k])
      data_int3$.width[count] <- dat$.width[i]
    }
  }
}

# Join with the previous one 

data_int <- rbind(data_int, data_int3)

# Compare with the upper and lower bounds of the  -------------------------------------

data_int_d <- data_int %>% 
  left_join(pops_data %>% dplyr::select(ID, SpeciesName, System, Taxon, threats, mu, uCI, lCI), 
            by=c(".variable"="threats")) %>%
  drop_na(ID) %>% 
  mutate(effect=ifelse(upper<lCI, "Antagonistic", 
                       ifelse(uCI<lower,
                              "Synergy", "Additive")))

# Proportion of interaction types by System ####################################

# First, I am going to create separate each stressor interaction, according 
# to the presence of each one of them. 
# Then I am going to classify the interactions according to the values of the
# confidence intervals. 

data_prop_d <- data_int_d %>% 
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

climate_change <- data_prop_d %>% 
  drop_na(climate) %>%
  group_by(climate, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>% 
  rename(variable=climate)

habitat_loss <- data_prop_d %>% 
  drop_na(habitatl) %>%
  group_by(habitatl,System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=habitatl)

exploitation <- data_prop_d %>% 
  drop_na(exploitation) %>%
  group_by(exploitation, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=exploitation)

invasive <- data_prop_d %>% 
  drop_na(invasive) %>%
  group_by(invasive, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=invasive)

disease <- data_prop_d %>% 
  drop_na(disease) %>%
  group_by(disease, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=disease)

pollution <- data_prop_d %>% 
  drop_na(pollution) %>%
  group_by(pollution, System, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=pollution)

# Now we join them 

data_freq_d <- rbind(climate_change,
                     habitat_loss, exploitation, 
                     disease, invasive, 
                     pollution)

# Proportion by taxon ----------------------------------------------------------

climate_change <- data_prop_d %>% 
  drop_na(climate) %>%
  group_by(climate, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>% 
  rename(variable=climate)

habitat_loss <- data_prop_d %>% 
  drop_na(habitatl) %>%
  group_by(habitatl,Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=habitatl)

exploitation <- data_prop_d %>% 
  drop_na(exploitation) %>%
  group_by(exploitation, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=exploitation)

invasive <- data_prop_d %>% 
  drop_na(invasive) %>%
  group_by(invasive, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=invasive)

disease <- data_prop_d %>% 
  drop_na(disease) %>%
  group_by(disease, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=disease)

pollution <- data_prop_d %>% 
  drop_na(pollution) %>%
  group_by(pollution, Taxon, interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  rename(variable=pollution)

# Now we join them 

data_freq_d_taxon <- rbind(climate_change,
                           habitat_loss, exploitation, 
                           disease, invasive, 
                           pollution)



# Fig. S5: Compare null models across systems ----------------------------------

load(paste0(ResultPath,"/PopMultEffects.RData"))

# Color palette for taxons

taxon_pal <- wes_palette("Cavalcanti1", n = 5)

# System palette

system_pal <- c("#A1D6E2","#336B87", "#CC954E")

# Panel a: Simple additive -----------------------------------------------------

data <- data_freq %>% mutate(System=as.factor(System),
                             freq=freq, 
                             interaction.type=fct_relevel(interaction.type,
                                                          c("Synergy", 
                                                            "Antagonistic", 
                                                            "Additive"))) 

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nObsType <- nlevels(as.factor(data$interaction.type))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$System)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$System <- rep(levels(data$System), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(System, variable)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label

label_data <- data %>% group_by(id, variable) %>% summarize(tot=sum(freq))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(System) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start 
grid_data <- grid_data[-1,]

# Make the plot

(gs11a <- ggplot(data) +
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=freq, fill=interaction.type), 
             stat="identity", alpha=0.8) +
    # Add a personalised grid
    geom_segment(data=grid_data, aes(x = end-0.5, y = 100, xend = start-0.5, yend = 100),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 80, xend = start-0.5, yend = 80),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 60, xend = start-0.5, yend = 60),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 40, xend = start-0.5, yend = 40),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 20, xend = start-0.5, yend = 20),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    # Colour of the fill
    scale_fill_manual(name = "",
                      values = c("#B32315",
                                          "#1E63B3",
                                          "#694364"))+
                                            # Add text showing the freq of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$id),5), y = c(20,40, 60, 80,100), 
             label = c("20", "40", "60", "80","100"), 
             color="grey", size=5, angle=0, fontface="bold", hjust=1.5) +
    ylim(-100,220) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    # Add labels on top of each bar
    geom_text(data=label_data, aes(x=id, y=tot+5, label=variable, hjust=hjust), 
              color="black", fontface="bold",alpha=0.6, size=5, 
              angle= label_data$angle, inherit.aes = FALSE ) +
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), 
                 colour = system_pal, alpha=0.8, size=2, inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -22, label=System), 
              hjust=c(1,0.5,0), colour = system_pal, alpha=0.8, size=4, 
              fontface="bold", inherit.aes = FALSE))

# Panel b: Multiplicative ------------------------------------------------------

data <- data_freq_m %>% mutate(System=as.factor(System),
                               freq=freq, 
                               interaction.type=fct_relevel(interaction.type,
                                                            c("Synergy", 
                                                              "Antagonistic", 
                                                              "Additive"))) 

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nObsType <- nlevels(as.factor(data$interaction.type))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$System)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$System <- rep(levels(data$System), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(System, variable)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label

label_data <- data %>% group_by(id, variable) %>% summarize(tot=sum(freq))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(System) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start 
grid_data <- grid_data[-1,]

# Make the plot

(gs11b <- ggplot(data) +
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=freq, fill=interaction.type), 
             stat="identity", alpha=0.8) +
    # Add a personalised grid
    geom_segment(data=grid_data, aes(x = end-0.5, y = 100, xend = start-0.5, yend = 100),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 80, xend = start-0.5, yend = 80),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 60, xend = start-0.5, yend = 60),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 40, xend = start-0.5, yend = 40),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 20, xend = start-0.5, yend = 20),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    # Colour of the fill
    scale_fill_manual(name = "",
                      values = c("#B32315",
                                          "#1E63B3",
                                          "#694364"))+
                                            # Add text showing the freq of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$id),5), y = c(20,40, 60, 80,100), 
             label = c("20", "40", "60", "80","100"), 
             color="grey", size=5, angle=0, fontface="bold", hjust=1.5) +
    ylim(-100,220) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    # Add labels on top of each bar
    geom_text(data=label_data, aes(x=id, y=tot+5, label=variable, hjust=hjust), 
              color="black", fontface="bold",alpha=0.6, size=5, 
              angle= label_data$angle, inherit.aes = FALSE ) +
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), 
                 colour = system_pal, alpha=0.8, size=2, inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -22, label=System), 
              hjust=c(1,0.5,0), colour = system_pal, alpha=0.8, size=4, 
              fontface="bold", inherit.aes = FALSE))

# Panel c: Dominace ------------------------------------------------------

data <- data_freq_d %>% mutate(System=as.factor(System),
                               freq=freq, 
                               interaction.type=fct_relevel(interaction.type,
                                                            c("Synergy", 
                                                              "Antagonistic", 
                                                              "Additive"))) 

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nObsType <- nlevels(as.factor(data$interaction.type))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$System)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$System <- rep(levels(data$System), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(System, variable)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label

label_data <- data %>% group_by(id, variable) %>% summarize(tot=sum(freq))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(System) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start 
grid_data <- grid_data[-1,]

# Make the plot

(gs11c <- ggplot(data) +
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=freq, fill=interaction.type), 
             stat="identity", alpha=0.8) +
    # Add a personalised grid
    geom_segment(data=grid_data, aes(x = end-0.5, y = 100, xend = start-0.5, yend = 100),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 80, xend = start-0.5, yend = 80),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 60, xend = start-0.5, yend = 60),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 40, xend = start-0.5, yend = 40),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 20, xend = start-0.5, yend = 20),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    # Colour of the fill
    scale_fill_manual(name = "",
                      values = c("#B32315",
                                          "#1E63B3",
                                          "#694364"))+
                                            # Add text showing the freq of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$id),5), y = c(20,40, 60, 80,100), 
             label = c("20", "40", "60", "80","100"), 
             color="grey", size=5, angle=0, fontface="bold", hjust=1.5) +
    ylim(-100,220) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    # Add labels on top of each bar
    geom_text(data=label_data, aes(x=id, y=tot+5, label=variable, hjust=hjust), 
              color="black", fontface="bold",alpha=0.6, size=5, 
              angle= label_data$angle, inherit.aes = FALSE ) +
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), 
                 colour = system_pal, alpha=0.8, size=2, inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -22, label=System), 
              hjust=c(1,0.5,0), colour = system_pal, alpha=0.8, size=4, 
              fontface="bold", inherit.aes = FALSE))

# Now by taxa ------------------------------------------------------------------

# Panel d: Simple additive -----------------------------------------------------

data <- data_freq_taxon %>% 
  mutate(Taxon=as.factor(Taxon),
         freq=freq,
         interaction.type=fct_relevel(interaction.type,
                                      c("Synergy", 
                                        "Antagonistic", 
                                        "Additive"))) 

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nObsType <- nlevels(as.factor(data$interaction.type))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$Taxon)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$Taxon <- rep(levels(data$Taxon), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(Taxon, variable)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label

label_data <- data %>% group_by(id, variable) %>% summarize(tot=sum(freq))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(Taxon) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start 
grid_data <- grid_data[-1,]

# Make the plot

(gs11d <- ggplot(data) +
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=freq, fill=interaction.type), 
             stat="identity", alpha=0.8) +
    # Add a personalised grid
    geom_segment(data=grid_data, aes(x = end-0.5, y = 100, xend = start-0.5, yend = 100),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 80, xend = start-0.5, yend = 80),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 60, xend = start-0.5, yend = 60),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 40, xend = start-0.5, yend = 40),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 20, xend = start-0.5, yend = 20),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    # Colour of the fill
    scale_fill_manual(name = "",values = c("#B32315",
                                                    "#1E63B3",
                                                    "#694364"))+
                                                      # Add text showing the freq of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$id),5), y = c(20,40, 60, 80,100), 
             label = c("20", "40", "60", "80","100"), 
             color="grey", size=5, angle=0, fontface="bold", hjust=1) +
    ylim(-100,220) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    # Add labels on top of each bar
    geom_text(data=label_data, aes(x=id, y=tot+5, label=variable, hjust=hjust), 
              color="black", fontface="bold", alpha=0.6, size=5, 
              angle= label_data$angle, inherit.aes = FALSE ) +
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), 
                 colour = taxon_pal, alpha=0.8, size=2 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -22, label=Taxon), 
              hjust=c(0.9,0.8,0.5,0,0), 
              colour = taxon_pal, alpha=0.8, size=4, 
              fontface="bold",inherit.aes = FALSE))

# Panel e: Multiplicative ------------------------------------------------------

data <- data_freq_m_taxon %>% 
  mutate(Taxon=as.factor(Taxon),
         freq=freq,
         interaction.type=fct_relevel(interaction.type,
                                      c("Synergy", 
                                        "Antagonistic", 
                                        "Additive"))) 

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nObsType <- nlevels(as.factor(data$interaction.type))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$Taxon)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$Taxon <- rep(levels(data$Taxon), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(Taxon, variable)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label

label_data <- data %>% group_by(id, variable) %>% summarize(tot=sum(freq))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(Taxon) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start 
grid_data <- grid_data[-1,]

# Make the plot

(gs11e <- ggplot(data) +
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=freq, fill=interaction.type), 
             stat="identity", alpha=0.8) +
    # Add a personalised grid
    geom_segment(data=grid_data, aes(x = end-0.5, y = 100, xend = start-0.5, yend = 100),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 80, xend = start-0.5, yend = 80),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 60, xend = start-0.5, yend = 60),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 40, xend = start-0.5, yend = 40),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 20, xend = start-0.5, yend = 20),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    # Colour of the fill
    scale_fill_manual(name = "",values = c("#B32315",
                                                    "#1E63B3",
                                                    "#694364"))+
                                                      # Add text showing the freq of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$id),5), y = c(20,40, 60, 80,100), 
             label = c("20", "40", "60", "80","100"), 
             color="grey", size=5, angle=0, fontface="bold", hjust=1) +
    ylim(-100,220) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    # Add labels on top of each bar
    geom_text(data=label_data, aes(x=id, y=tot+5, label=variable, hjust=hjust), 
              color="black", fontface="bold", alpha=0.6, size=5, 
              angle= label_data$angle, inherit.aes = FALSE ) +
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), 
                 colour = taxon_pal, alpha=0.8, size=2 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -22, label=Taxon), 
              hjust=c(0.9,0.8,0.5,0,0), 
              colour = taxon_pal, alpha=0.8, size=4, 
              fontface="bold",inherit.aes = FALSE))

# Panel f: Dominance -----------------------------------------------------------

data <- data_freq_d_taxon %>% 
  mutate(Taxon=as.factor(Taxon),
         freq=freq,
         interaction.type=fct_relevel(interaction.type,
                                      c("Synergy", 
                                        "Antagonistic", 
                                        "Additive"))) 

# Set a number of 'empty bar' to add at the end of each group

empty_bar <- 2
nObsType <- nlevels(as.factor(data$interaction.type))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$Taxon)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$Taxon <- rep(levels(data$Taxon), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(Taxon, variable)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label

label_data <- data %>% group_by(id, variable) %>% summarize(tot=sum(freq))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines

base_data <- data %>% 
  group_by(Taxon) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)

grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start 
grid_data <- grid_data[-1,]

# Make the plot

(gs11f <- ggplot(data) +
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=freq, fill=interaction.type), 
             stat="identity", alpha=0.8) +
    # Add a personalised grid
    geom_segment(data=grid_data, aes(x = end-0.5, y = 100, xend = start-0.5, yend = 100),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 80, xend = start-0.5, yend = 80),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 60, xend = start-0.5, yend = 60),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 40, xend = start-0.5, yend = 40),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 20, xend = start-0.5, yend = 20),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    # Colour of the fill
    scale_fill_manual(name = "",values = c("#B32315",
                                                    "#1E63B3",
                                                    "#694364"))+
                                                      # Add text showing the freq of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$id),5), y = c(20,40, 60, 80,100), 
             label = c("20", "40", "60", "80","100"), 
             color="grey", size=5, angle=0, fontface="bold", hjust=1) +
    ylim(-100,220) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    # Add labels on top of each bar
    geom_text(data=label_data, aes(x=id, y=tot+5, label=variable, hjust=hjust), 
              color="black", fontface="bold", alpha=0.6, size=5, 
              angle= label_data$angle, inherit.aes = FALSE ) +
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), 
                 colour = taxon_pal, alpha=0.8, size=2 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -22, label=Taxon), 
              hjust=c(0.9,0.8,0.5,0,0), 
              colour = taxon_pal, alpha=0.8, size=4, 
              fontface="bold",inherit.aes = FALSE))


# Combine them -----------------------------------------------------------------

#Create first row 

(row1 <- plot_grid(gs11a+ggtitle("Simple additive"), 
                   gs11b+ggtitle("Multiplicative"), 
                   gs11c+ggtitle("Dominance"),
                   gs11d+ggtitle("Simple additive"), 
                   gs11e+ggtitle("Multiplicative"), 
                   gs11f+ggtitle("Dominance"), 
                   nrow=2, labels="auto"))

# Get the legend

legend <- get_legend(gs11a+theme(legend.position = "bottom",
                                 legend.text = element_text(size=14)))
# Plot it 

(figS11<- plot_grid(row1, legend, 
                    nrow = 2, rel_heights = c(1,0.05)))


# Save it

ggsave("Figure S11.pdf", figS11, 
       width = 16, height = 14,
       path = ResultPath)

