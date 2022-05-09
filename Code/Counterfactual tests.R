# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Counterfactual tests.R         
# - DATE:        15/09/2020
# - DESCRIPTION: Code to test scenarios where the threats have been removed.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(ggplot2)
library(tidybayes)
library(tidyverse)
library(dplyr)
library(brms)
library(bayesplot)
library(bayestestR)
library(ggdist)
library(cowplot)
library(tidyr)
library(patchwork)
library(MetBrewer)
library(foreach)
library(doSNOW)
library(parallel)
library(naniar)

# Set default ggplot theme

theme_set(theme_minimal()+
            theme(axis.title.x = element_text(size=14,
                                              margin = margin(t = 10, r = 0, b = 0, l = 0)), 
                  axis.title.y = element_text(size=14,
                                              margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.line.x = element_line(color="black", size = 0.5),
                  axis.line.y = element_line(color="black", size = 0.5),
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

load(paste0(ResultPath, "/TenResults.RData"))
load(paste0(DataPath, "/TenYearsData.RData"))

# Change invasive species

pops_data10 <- pops_data10 %>% 
  mutate(threats=gsub("Invasive spp/genes", "Invasive", threats))

# Palette for the threats 

threat_palette<-c(met.brewer(name="Hokusai1", n=6, type="continuous"), "grey60")

# Conterfactual tests ----------------------------------------------------------

# Obtain the predicted value for each population 

pop_perd <- fitted(m1, scale = "response")

# Edit the data to have the necessary variables and join with data values

pop_perd <- pop_perd %>%
  as.data.frame() %>% 
  bind_cols(pops_data10[c("threats", "ID")]) %>% # Add the threat affecting the population 
  mutate(lambda=exp(Estimate), # Backtransform to lambda
         number = str_count(threats,"[A-Z]")) # Get the number of threats


# Obtain the mean coefficient for each threat 

mus_mean <- m1 %>% 
  gather_draws(`b_.*`, regex = TRUE) %>%
  mean_qi(.width = .95) %>%
  mutate(.variable = gsub("b_threats", "", .variable), 
         .variable = gsub("InvasivesppDgenes", "Invasive", .variable),
         .variable = gsub("Habitatloss", "Habitat loss", .variable),
         .variable = gsub("Climatechange", "Climate change", .variable),# Create an initial fake population
         number = str_count(.variable,"[A-Z]")) 

# Substract the mean coefficient from the predicted values
# We will parallelise the code, so this will take multiple steps

# Define the  threats

threats <- unique(pops_data10$threats[which(pops_data10$n.threat>=1)]) %>% 
  as_tibble() %>% 
  mutate(number=str_count(value,"[A-Z]"))

# Define the number of cores to use

clus  <- makeCluster(detectCores() - 1)

# register the cluster with doParallel

registerDoSNOW(clus)

# load in the options using .options.snow
mu_scenarios <- foreach(i = threats$value, 
                        .combine = 'rbind',
                        .packages = c("tidyverse", "naniar")) %dopar% {
                          
                          x <- threats$number[threats$value==i]
                          
                          if(x==1){
                            new <- pop_perd %>% 
                              filter(grepl(i,threats)) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==i])) %>% 
                              rbind(pop_perd %>% 
                                      filter(!grepl(i,threats))) %>% 
                              mutate(scenario=paste("No", i))
                          } else if(x==2){
                            x2 <- unlist(strsplit(i, 
                                                  ":", 
                                                  fixed=T))
                            
                            new1 <- pop_perd %>% 
                              filter(grepl(x2[1],threats), 
                                     grepl(x2[2],threats), number!=1) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==i])) 
                            
                            new2 <- pop_perd %>% 
                              filter(grepl(x2[1],threats), 
                                     number!=2, !grepl(i,threats),
                                     !ID%in%new1$ID) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==x2[1]]))
                            
                            new <- pop_perd %>% 
                              filter(grepl(x2[2],threats), number!=2,
                                     !grepl(i,threats), 
                                     !ID%in%new1$ID,
                                     !ID%in%new2$ID) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==x2[2]])) %>% 
                              rbind(new1) %>% 
                              rbind(new2)  
                              
                            new <- new %>% 
                              rbind(pop_perd %>% 
                                      filter(!ID%in%new$ID)) %>% 
                              mutate(scenario=paste("No", i))
                          }else if(x==3){
                            
                            # Separate the threats 
                            x3 <- unlist(strsplit(i, 
                                                  ":", 
                                                  fixed=T))
                            
                            # Create new names
                            
                            two1 <- paste0(x3[1], ":", x3[2]) 
                            two2 <- paste0(x3[2], ":", x3[3]) 
                            two3 <- paste0(x3[1], ":", x3[3])
                            
                            # Create a data frame with the three threats
                            new1 <- pop_perd %>% 
                              filter(grepl(x3[1],threats), 
                                     grepl(x3[2],threats),
                                     grepl(x3[3],threats), 
                                     number==3) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==i])) 
                            
                            # Subset the mus with two threats
                            two_threat1 <- mus_mean %>%  filter(grepl(x3[1], .variable),
                                                                grepl(x3[2], .variable),
                                                                number==2)
                            two_threat2 <- mus_mean %>%  filter(grepl(x3[2], .variable),
                                                                grepl(x3[3], .variable),
                                                                number==2)
                            two_threat3 <- mus_mean %>%  filter(grepl(x3[1], .variable),
                                                                grepl(x3[3], .variable),
                                                                number==2)
                            
                            # subset two threats and substract
                            if(dim(two_threat1)[1]==0){
                              lam <- mus_mean$.value[mus_mean$.variable==x3[1]]+mus_mean$.value[mus_mean$.variable==x3[2]]
                              new2 <- tryCatch({pop_perd %>% 
                                filter(grepl(x3[1], threats),
                                       grepl(x3[2], threats),
                                       !ID%in%new1$ID) %>% 
                                mutate(lambda=exp(Estimate-lam))},
                                error=function(e){e <- new1$lambda=NA})
                            }else{new2 <- tryCatch({pop_perd %>% 
                              filter(grepl(x3[1], threats),
                                      grepl(x3[2], threats),
                                     !ID%in%new1$ID) %>% 
                              mutate(lambda=exp(Estimate-two_threat1$.value))},
                            error=function(e){e <- new1$lambda=NA})}
                            if(dim(two_threat2)[1]==0){
                              lam <- mus_mean$.value[mus_mean$.variable==x3[2]]+mus_mean$.value[mus_mean$.variable==x3[3]]
                              new3 <- tryCatch({pop_perd %>% 
                                filter(grepl(x3[2], threats),
                                       grepl(x3[3], threats),
                                       !ID%in%new1$ID) %>% 
                                mutate(lambda=exp(Estimate-lam))},
                          error=function(e){e <- new1$lambda=NA})
                            }else{
                              new3 <- tryCatch({pop_perd %>% 
                                filter(grepl(x3[2], threats),
                                       grepl(x3[3], threats),
                                       !ID%in%new1$ID) %>% 
                                mutate(lambda=exp(Estimate-two_threat2$.value))},
                                error=function(e){e <- new1$lambda=NA})  
                            }
                            if(dim(two_threat3)[1]==0){ 
                              lam <- mus_mean$.value[mus_mean$.variable==x3[1]]+mus_mean$.value[mus_mean$.variable==x3[3]]
                            new4 <- tryCatch({pop_perd %>% 
                                filter(grepl(x3[1], threats),
                                       grepl(x3[3], threats),
                                       !ID%in%new1$ID) %>% 
                                mutate(lambda=exp(Estimate-lam))},
                                error=function(e){e <- new1$lambda=NA})
                            }else{
                              new4 <- tryCatch({pop_perd %>% 
                                  filter(grepl(x3[1], threats),
                                         grepl(x3[3], threats),
                                         !ID%in%new1$ID) %>% 
                                  mutate(lambda=exp(Estimate-two_threat3$.value))},
                                  error=function(e){e <- new1$lambda=NA})
                            }
                            new1 %>% replace_with_na()
                            # Remove single effects
                            
                              new5 <- pop_perd %>% 
                                filter(grepl(x3[1], threats), 
                                       !ID%in%new1$ID, 
                                       !ID%in%new2$ID,
                                       !ID%in%new4$ID) %>% 
                                mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==x3[1]]))
                            
                            
                            new6 <- pop_perd %>% 
                              filter(grepl(x3[2], threats),
                                     !ID%in%new1$ID, 
                                     !ID%in%new2$ID,
                                     !ID%in%new3$ID) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==x3[2]]))
                            
                            # Final subset
                            
                            new <- pop_perd %>% 
                              filter(grepl(x3[3], threats),
                                       !ID%in%new1$ID, 
                                       !ID%in%new3$ID,
                                       !ID%in%new4$ID) %>% 
                                mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==x3[3]])) %>% 
                                rbind(new1) %>% 
                                rbind(new2) %>% 
                                rbind(new3) %>% 
                                rbind(new4) %>% 
                                rbind(new5) %>% 
                                rbind(new6)
                            
                              
                            
                            new <- new %>% 
                              rbind(pop_perd %>% 
                                      filter(!ID%in%new$ID)) %>% 
                              mutate(scenario=paste("No", i)) %>% 
                              drop_na(lambda)
                          }
                          
                        }

# Stop the parallelisation

stopCluster(clus)

# Add the control scenario 

scenarios_data <- 
  pop_perd %>% mutate(scenario="No Management") %>% 
  bind_rows(mu_scenarios) %>% 
  filter(threats!="None") %>% 
  mutate(prop=((log(lambda)-Estimate)/abs(Estimate))*100,
         prop=round(prop,2),
         number_scen= str_count(scenario,"[A-Z]")) %>% 
  group_by(scenario) %>% 
  mutate(pos=median(lambda), 
         pos2=median(prop)) %>% 
  ungroup()  

# Add the control scenario 

scenarios_median <- scenarios_data %>% 
  group_by(scenario) %>% 
  summarise(mcontrol=median(Estimate),
            mcounter=median(log(lambda))) %>% 
  mutate(prop=((mcounter-mcontrol)/abs(mcontrol))*100,
         prop=round(prop,2))

# Figure 4: Counterfactual single threats -----------------------------------------------

(figure4 <- scenarios_data %>% 
   filter(scenario!="No Management", 
          number_scen==2) %>%
   ggplot(aes(x=reorder(scenario, pos), y=log(lambda), 
              fill=scenario, 
              colour=scenario)) + 
   stat_halfeye(adjust = .5, 
                width = .6, 
                .width = 0,
                alpha=.8,
                justification = -.3, 
                point_colour = NA) + 
   geom_boxplot(width = .2, 
                outlier.shape = NA, 
                alpha=.8,
                colour="black") +
   theme(legend.position = "none")+
   scale_fill_manual(values=c(threat_palette))+
   scale_colour_manual(values=threat_palette)+
   gghalves::geom_half_point(show.legend = FALSE,
                             side = "r",
                             range_scale = .2,
                             shape=21, alpha=0.5,
                             colour="grey25") +
   geom_hline(yintercept = 0)+
   geom_hline(yintercept = median(pop_perd$Estimate[pop_perd$threats!="None"]), 
              linetype="dashed",size=.8)+
   coord_cartesian(xlim = c(1.2, NA), clip = "off") +
   labs(y=expression(paste("Population trend (", mu,")")),
        x="Counterfactual scenarios")+
   # Add annotation of the horizontal line
   geom_curve(aes(x = 7, 
                  y = median(pop_perd$Estimate[pop_perd$threats!="None"]),
                  xend = 7.3, 
                  yend = -0.015),
              arrow = arrow(length = unit(0.07, "inch")), size = 0.4, 
              curvature = -0.3, colour="black")+
   annotate("text", x = 7.3, y = -0.03, 
            label = "Median trend of \n non-managed scenarios")+
   geom_curve(aes(x = 7, 
                  y = 0,
                  xend = 7.3, 
                  yend = 0.01),
              arrow = arrow(length = unit(0.07, "inch")), size = 0.4, 
              curvature = 0.3, colour="black")+
   annotate("text", x = 7.3, y = 0.015, 
            label = "Population trend = 0")+
   theme(plot.margin=unit(c(0,2,0,0), "cm")))

# Save it

ggsave("Figure 4.pdf", figure4, path=ResultPath,
       width = 11, height = 8)

# Figure S6: Counterfactual plots -----------------------------------------------


# Panel a

(gs6a <- scenarios_data %>% 
    filter(scenario!="No Management", 
           number_scen==3) %>%
    ggplot(aes(x=reorder(scenario, pos), y=log(lambda), 
               fill=scenario, 
               colour=scenario)) + 
    stat_halfeye(adjust = .5, 
                 width = .6, 
                 .width = 0,
                 alpha=.8,
                 justification = -.3, 
                 point_colour = NA) + 
    geom_boxplot(width = .15, 
                 outlier.shape = NA, 
                 alpha=.8,
                 colour="black") +
    theme(legend.position = "none")+
    scale_fill_manual(values=met.brewer(name="Tiepolo", n=14, type="continuous")) +
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = median(pop_perd$Estimate[pop_perd$threats!="None"]), 
               linetype="dashed",size=.8)+
      labs(y=expression(paste("Population trend (", mu,")")),
         x="Counterfactual scenarios")+
    theme(plot.margin=unit(c(0,0,0,0), "cm")))

# Panel b

(gs6b <- scenarios_data %>% 
    filter(scenario!="No management", 
           number_scen==4) %>%
    ggplot(aes(x=reorder(scenario, pos), y=log(lambda), 
               fill=scenario, 
               colour=scenario)) + 
    stat_halfeye(adjust = .5, 
                 width = .6, 
                 .width = 0,
                 alpha=.8,
                 justification = -.3, 
                 point_colour = NA) + 
    geom_boxplot(width = .15, 
                 outlier.shape = NA, 
                 alpha=.8,
                 colour="black") +
    theme(legend.position = "none")+
    scale_fill_manual(values=met.brewer(name="Nattier", n=18, type="continuous")) +
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = median(pop_perd$Estimate[pop_perd$threats!="None"]), 
               linetype="dashed",size=.8)+
    labs(y=expression(paste("Population trend (", mu,")")),
         x="")+
    theme(plot.margin=unit(c(0,0,0,-2), "cm")))

# Combined plot

(figureS6 <- (gs6a+coord_flip())+
    (gs6b+coord_flip())+
    plot_annotation(tag_levels = c('a')))

# Save the plot

ggsave("Figure S6.pdf",figureS6,
       path = ResultPath, width = 14, height = 12)