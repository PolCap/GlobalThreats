# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Figures.R         
# - DATE:        15/09/2020
# - DESCRIPTION: Code to produce the figures.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

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
                  axis.line.x = element_line(color="black", linewidth = 0.5),
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

load(paste0(ResultPath, "/Results.RData"))
load(paste0(DataPath, "/FiveData.RData"))

# Change name invasive species 

pops_data<- pops_data%>% 
  mutate(threats = gsub("Invasive spp/genes", "Invasive", threats),)

# Load functions for the plots 

source(paste0(CodePath,"/PlotFunctions.R"))

#Create a palette for the threats 

threat_palette<-c(met.brewer(name="Hokusai1", n=6, type="continuous"), "grey60")

# Color palette for taxons

taxon_pal <- wes_palette("Cavalcanti1", n = 5)

# System palette

system_pal <- c("#A1D6E2","#336B87", "#CC954E")

# Number of time series --------------------------------------------------------

# Number of species 

length(unique(pops_data$SpeciesName))

# Number of population time series per taxon

pops_data%>% 
  group_by(Taxon) %>% 
  summarise(n=n())

# Number of population time series per system

pops_data%>% 
  group_by(System) %>% 
  summarise(n=n())

# Min, max, and mean duration 

max(pop_data$Duration)
min(pop_data$Duration)
mean(pop_data$Duration)

# Min and max years

max(pop_data$Year)
min(pop_data$Year)

# Figure 1: Map ################################################################

# Set the world map

world <- map_data("world")

# Create map

(g2a <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    #coord_proj("+proj=wintri") +
    theme_map() + 
    geom_point(data = pops_data, 
               aes(x = Longitude, y = Latitude, 
                   fill = mu),
               alpha=0.7, shape=21, size=3) +
    scale_y_continuous(limits = c(-80, 80)) +
    scale_fill_gradient2(expression(paste("Population trend (", mu, ")")),
                          midpoint = 0, 
                          low = "#c1666b", 
                          mid = "#e4dfda",
                          high = "#4281a4",
                          guide =guide_colourbar(nbin=100,
                                                 barwidth = 9,
                                                 title.position="top"),
                          breaks=c(-0.10, -0.05, 0, 0.05, 0.10),
                         limits=c(-0.10, 0.10))+
    expand_limits(x = 0, y = 0)+
    theme(legend.position = "none",
          legend.direction = "horizontal",
          legend.title = element_text(size = 12, hjust =0.5),
          legend.text = element_text(size = 10),
          plot.margin = unit(c(0,0,0,0), units = , "cm")))

# Legend

legend1 <- get_legend(g2a+theme(legend.position = "bottom",
                               legend.title = element_text(size = 14)))

# Create distribution 1

anot <- data.frame(x=median(pops_data$Duration), 
                   y=300, 
                   label =paste("Median = ", median(pops_data$Duration)))


(gd1 <- pops_data%>% 
    ggplot() +
    geom_histogram(aes(x = Duration),
                   binwidth = 2,
                       alpha = .6,
                 position = "stack",
                 fill = "grey40", color = "white") +
    geom_vline(aes(xintercept = median(Duration)),
             linetype = "longdash", size=0.5) +
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(breaks = seq(10, 60, 10),
                       expand = c(0,0))+
    ggrepel::geom_text_repel(data = anot, aes(x=x, y=y+100,
                                      label =label),
                    nudge_x = 10, nudge_y = 0)+
    labs(x="Duration (years)", 
       y="Number of datasets")+
    theme(axis.title.x = element_text(margin = margin(t = 0)),
          axis.text.y = element_text(margin = margin(r = 0))))


# Combined map ----

(g2a<- ggdraw(g2a) +
   draw_plot(gd1, x = -0.32, y = -0.3143, scale = 0.35) +
   draw_plot(legend1, 
             x=0.4, y=-0.435, scale = 1)) 

# Panel b: frequency of the different threats by system ------------------------

# Readjust the data frame to contain the threats as individual columns

threat_freq_sys <- pop_data %>% 
  distinct(ID, .keep_all=T) %>% 
  mutate(none=ifelse(threats=="None", "none", NA)) %>% 
  dplyr::select(ID, System, none, pollution, habitatl, climatechange, 
         invasive, exploitation, disease) %>% 
  pivot_longer(3:9, names_to="threats") %>% 
  drop_na(value) %>% 
  group_by(System, threats) %>% 
  summarise(n=n()) %>% 
  mutate(freq=(n/sum(n)),
         threats=ifelse(threats=="climatechange", "Climate change",
                        ifelse(threats=="exploitation", 
                               "Exploitation",
                               ifelse(threats=="invasive", 
                                      "Invasive",ifelse(threats=="disease", 
                                                        "Disease",ifelse(threats=="habitatl",
                                                                         "Habitat loss", 
                                                                         ifelse(threats=="pollution",
                                                                                "Pollution", 
                                                                                ifelse(threats=="none", 
                                                                                       "None", threats)))))))) 
# Plot

(g2b <- threat_freq_sys %>% 
    mutate(threats=factor(threats,
                          levels=c("Habitat loss", "Exploitation", 
                                   "Climate change","Pollution", 
                                   "Invasive",
                                   "Disease",
                                   "None"))) %>% 
    ggplot(aes(fill = threats,
               y=freq, x=System)) + 
    geom_bar(stat = "identity",
             position = "fill") +
    scale_fill_manual(values = threat_palette)+
    scale_y_continuous(breaks = seq(0, 1, .2), 
                       label = scales::percent,
                       expand = c(0,0))+ 
    labs(y="Proportion of threats (%)", x="",fill="") +
    # geom_text(aes(label = paste0(round(freq*100, 0), "%")),
    #           size=5,
    #           position = position_stack(vjust = 0.5)) +
    theme(legend.position="bottom",
          legend.text = element_text(size=12),
          strip.text = element_text(hjust = 0),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size = 14), 
          plot.margin = unit(c(0.5,0,0,0), units = "cm")))


# Panel c: frequency of the different threats by taxon -------------------------

pop_data <- pop_data %>% 
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
                                                   "Reptiles", "NA"))))))

threat_freq_tax <- pop_data %>% 
  distinct(ID, .keep_all=T) %>% 
  mutate(none=ifelse(threats=="None", "none", NA)) %>% 
  dplyr::select(ID, Taxon, none, pollution, habitatl, climatechange, 
         invasive, exploitation, disease) %>% 
  pivot_longer(3:9, names_to="threats") %>% 
  drop_na(value) %>% 
  group_by(Taxon, threats) %>% 
  summarise(n=n()) %>% 
  mutate(freq=(n/sum(n)),
         threats=ifelse(threats=="climatechange", "Climate change",
                        ifelse(threats=="exploitation", 
                               "Exploitation",
                               ifelse(threats=="invasive", 
                                      "Invasive",ifelse(threats=="disease", 
                                                        "Disease",ifelse(threats=="habitatl",
                                                                         "Habitat loss", 
                                                                         ifelse(threats=="pollution",
                                                                                "Pollution", 
                                                                                ifelse(threats=="none", 
                                                                                       "None", threats)))))))) 
# Plot it  

(g2c <- threat_freq_tax %>% 
    mutate(threats=factor(threats,
                          levels=c("Habitat loss", "Exploitation", 
                                   "Climate change","Pollution", 
                                   "Invasive",
                                   "Disease",
                                   "None"))) %>% 
    ggplot(aes(fill = threats,
               y=freq, x=Taxon)) + 
    geom_bar(stat = "identity",
             position = "fill") +
    scale_fill_manual(values = threat_palette)+
    scale_y_continuous(breaks = seq(0, 1, .2), 
                       label = scales::percent,
                       expand = c(0,0))+ 
    labs(y="Proportion of threats (%)", x="",fill="") +
    # geom_text(aes(label = paste0(round(freq*100, 0), "%")),
    #           size=3,
    #           position = position_stack(vjust = 0.5)) +
    theme(legend.position="none",
          legend.text = element_text(size=12),
          strip.text = element_text(hjust = 0),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size = 14), 
          plot.margin = unit(c(0.5,0,0,0), units = "cm")))


# Left plot

(left_panel <- plot_grid(g2b+theme(legend.position = "none"),
                         g2c, 
                         labels = c("b", "c"), nrow = 2))

# Add legend

legend2 <- get_legend(g2b)

#Column2

(col2 <- plot_grid(left_panel, legend2,
                  ncol = 1,
                  rel_heights = c(1, .1)))

# Combine all plots ------------------------------------------------------------

(figure1 <- plot_grid(g2a, NULL, col2, 
                      rel_widths = c(2.1, -0.1, 1), 
                      nrow=1,
                      labels = c("a", "")))

# Save it

ggsave("Figure 1.pdf", figure1, 
       width = 18, height = 8,
       path = ResultPath)

# Figure 2: Single and multiple stressors by system and taxa ###################
# Panel a: Single stressors by system ------------------------------------------

(g3a <- ms1 %>%
   gather_draws(`b_.*`, regex = TRUE) %>% 
   median_qi(.width = .95) %>%
   mutate(.variable = gsub("b_threats", "", .variable), 
          .variable = gsub("InvasivesppDgenes", "Invasive", .variable),
          .variable = gsub("HabitatdegradationDchange", "Habitat degradation", .variable),
          .variable = gsub("Habitatloss", "Habitat loss", .variable),
          .variable = gsub("Climatechange", "Climate change", .variable),
          System = gsub(".*System", "", .variable),
          System =gsub(":", "", System),
          .variable = gsub(":System.*", "", .variable),
          number = str_count(.variable,"[A-Z]")) %>%
   filter(number==1) %>% 
   ggplot(aes(x = reorder(.variable, -.lower),
              y = .value, colour=System, group=System)) +
    geom_point(size=4,position =position_dodge(.7)) +
    geom_errorbar( aes(ymin=.lower, ymax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                   position=position_dodge(0.7))+
  geom_hline(yintercept = 0, lty = 2, size = 0.5) +
   labs(x = "",
        y = expression(paste("Population trend (", mu, ")",sep = ""))) +
   scale_y_continuous(labels = scaleFUN,limits = c(-.25,.25)) +
   scale_colour_manual("", values = system_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Table S1 ---------------------------------------------------------------------

TableS1 <- as.data.frame(describe_posterior(ms1, ci = 0.95, test="none")) %>% 
  mutate(Parameter = gsub("b_threats", "", Parameter), 
         Parameter = gsub("InvasivesppDgenes", "Invasive", Parameter),
         Parameter = gsub("HabitatdegradationDchange", "Habitat degradation", Parameter),
         Parameter = gsub("Habitatloss", "Habitat loss", Parameter),
         Parameter = gsub("Climatechange", "Climate change", Parameter),
         Parameter =gsub("\\.", " ", Parameter),
         System = gsub(".*System", "", Parameter),
         System =gsub("\\.", "", System),
         Parameter = gsub("System.*", "", Parameter),
         number = str_count(Parameter,"[A-Z]")) 

# Save it 

setwd(ResultPath)

write.csv2(TableS1, "Table S1.csv")

# Panel b: Single stressors by taxon -------------------------------------------

(g3b <- mt1 %>%
   gather_draws(`b_.*`, regex = TRUE) %>% 
   median_qi(.width = .95) %>%
   mutate(.variable = gsub("b_threats", "", .variable), 
          .variable = gsub("InvasivesppDgenes", "Invasive", .variable),
          .variable = gsub("HabitatdegradationDchange", "Habitat degradation", .variable),
          .variable = gsub("Habitatloss", "Habitat loss", .variable),
          .variable = gsub("Climatechange", "Climate change", .variable),
          Taxon = gsub(".*Taxon", "", .variable),
          Taxon =gsub(":", "", Taxon),
          .variable = gsub(":Taxon.*", "", .variable),
          number = str_count(.variable,"[A-Z]")) %>%
    filter(number==1, 
           Taxon!="Amphibians"|.variable!="Disease",
           Taxon!="Amphibians"|.variable!="Exploitation",
           Taxon!="Amphibians"|.variable!="Climate change",
           Taxon!="Amphibians"|.variable!="Pollution",
           Taxon!="Amphibians"|.variable!="Invasive",
           Taxon!="Reptiles"|.variable!="Disease",
           Taxon!="Fish"|.variable!="Disease",
           Taxon!="Reptiles"|.variable!="Pollution") %>%
   ggplot(aes(x = reorder(.variable, -.lower),
              y = .value, colour=Taxon, group=Taxon)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar( aes(ymin=.lower, ymax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                  position=position_dodge(0.7))+
   geom_hline(yintercept = 0, lty = 2, size = 0.5) +
   labs(x = "",
        y = "") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.25,.25)) +
   scale_colour_manual("", 
                       values = taxon_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Table S2 ---------------------------------------------------------------------

TableS2 <- as.data.frame(describe_posterior(mt1, ci = 0.95, test="none")) %>% 
  mutate(Parameter = gsub("b_threats", "", Parameter), 
         Parameter = gsub("InvasivesppDgenes", "Invasive", Parameter),
         Parameter = gsub("HabitatdegradationDchange", "Habitat degradation", Parameter),
         Parameter = gsub("Habitatloss", "Habitat loss", Parameter),
         Parameter = gsub("Climatechange", "Climate change", Parameter),
         Parameter =gsub("\\.", " ", Parameter),
         Taxon = gsub(".*Taxon", "", Parameter),
         Taxon =gsub("\\.", "", Taxon),
         Parameter = gsub("Taxon.*", "", Parameter),
         number = str_count(Parameter,"[A-Z]"))

# Save it 

write.csv2(TableS2, "TableS2.csv")

# Panel c: Multiple stressors by system ----------------------------------------

(g3c <- mmu_fr %>%
   gather_draws(`b_.*`, regex = TRUE) %>%
   median_qi(.width = .95) %>%
   mutate(System="Freshwater") %>% 
   rbind(mmu_mar %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(System="Marine")) %>% 
   rbind(mmu_ter %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(System="Terrestrial")) %>% 
   mutate(.variable = gsub("b_", "", .variable),
          .variable = gsub("InvasivesppDgenes", "Invasive", .variable),
          .variable = gsub("HabitatdegradationDchange", "Habitat degradation", .variable),
          .variable = gsub("Habitatloss", "Habitat loss", .variable),
          .variable = gsub("Climatechange", "Climate change", .variable),
          .variable = gsub("noneNone", "None", .variable),
          threat = gsub("pollution", "", .variable),
          threat = gsub("habitatl", "", threat),
          threat = gsub("invasive", "", threat),
          threat = gsub("climatechange", "", threat),
          threat = gsub("exploitation", "", threat),
          threat = gsub("none", "", threat),
          threat = gsub("disease", "", threat),
          threat=ifelse(grepl("None",.variable), "None",
                        ifelse(grepl("No", .variable), "No", 
                               threat)),
          threat =factor(threat, levels = c("None",
                                            "Habitat loss",
                                            "Exploitation",
                                            "Pollution",
                                            "Climate change",
                                            "Disease",
                                            "Invasive"))) %>%
   filter(threat!="No") %>%
   ggplot(aes(x = reorder(threat, -.lower),
              y = .value, group=System, colour=System)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar( aes(ymin=.lower, ymax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                  position=position_dodge(0.7))+
   geom_hline(yintercept = 0, lty = 2, size = 0.5) +
   labs(x = "",
        y = expression(paste("Population trend (", mu, ")",sep = ""))) +
   scale_y_continuous(labels = scaleFUN,limits = c(-.4,.4)) +
   scale_colour_manual("", values = system_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Table S3 ---------------------------------------------------------------------

TableS3 <- as.data.frame(describe_posterior(mmu_fr, ci = 0.95, test="none")) %>% 
  mutate(System="Freshwater") %>% 
  rbind(as.data.frame(describe_posterior(mmu_mar, ci = 0.95, test="none")) %>% 
          mutate(System="Marine")) %>% 
  rbind(as.data.frame(describe_posterior(mmu_ter, ci = 0.95, test="none")) %>% 
          mutate(System="Terrestrial")) %>% 
  mutate(Parameter = gsub("b_", "", Parameter),
         Parameter = gsub("InvasivesppDgenes", "Invasive", Parameter),
         Parameter = gsub("HabitatdegradationDchange", "Habitat degradation", Parameter),
         Parameter = gsub("Habitatloss", "Habitat loss", Parameter),
         Parameter = gsub("Climatechange", "Climate change", Parameter),
         Parameter = gsub("noneNone", "None", Parameter),
         threat = gsub("pollution", "", Parameter),
         threat = gsub("habitatl", "", threat),
         threat = gsub("invasive", "", threat),
         threat = gsub("climatechange", "", threat),
         threat = gsub("exploitation", "", threat),
         threat = gsub("none", "", threat),
         threat = gsub("disease", "", threat),
         threat=ifelse(grepl("None",Parameter), "None",
                       ifelse(grepl("No", Parameter), "No", 
                              threat))) 
# Save it 

write.csv2(TableS3, "TableS3.csv")

# Panel d: Multiple stressors by taxon -----------------------------------------

(g3d <- mmu_am %>%
   gather_draws(`b_.*`, regex = TRUE) %>%
   median_qi(.width = .95) %>%
   mutate(Taxon="Amphibians") %>% 
   rbind(mmu_bi %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(Taxon="Birds")) %>% 
   rbind(mmu_f %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(Taxon="Fish")) %>% 
   rbind(mmu_ma %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(Taxon="Mammals")) %>% 
   rbind(mmu_rep %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(Taxon="Reptiles")) %>% 
   mutate(.variable = gsub("b_", "", .variable),
          .variable = gsub("InvasivesppDgenes", "Invasive", .variable),
          .variable = gsub("HabitatdegradationDchange", "Habitat degradation", .variable),
          .variable = gsub("Habitatloss", "Habitat loss", .variable),
          .variable = gsub("Climatechange", "Climate change", .variable),
          .variable = gsub("noneNone", "None", .variable),
          threat = gsub("pollution", "", .variable),
          threat = gsub("habitatl", "", threat),
          threat = gsub("invasive", "", threat),
          threat = gsub("climatechange", "", threat),
          threat = gsub("exploitation", "", threat),
          threat = gsub("none", "", threat),
          threat = gsub("disease", "", threat),
          threat=ifelse(grepl("None",.variable), "None",
                        ifelse(grepl("No", .variable), "No", 
                               threat)),
          threat =factor(threat, levels = c("None",
                                            "Habitat loss",
                                            "Exploitation",
                                            "Pollution",
                                            "Climate change",
                                            "Disease",
                                            "Invasive"))) %>%
   filter(threat!="No",
          threat!="Disease"|Taxon!="Fish") %>%
   ggplot(aes(x = reorder(threat, -.lower),
              y = .value, group=Taxon, colour=Taxon)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar( aes(ymin=.lower, ymax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                  position=position_dodge(0.7))+
   geom_hline(yintercept = 0, lty = 2, size = 0.5) +
   labs(x = "",
        y = "") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.4,.4)) +
   scale_colour_manual("", values = taxon_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Table S4 ---------------------------------------------------------------------

TableS4 <- as.data.frame(describe_posterior(mmu_am, ci = 0.95, test="none")) %>% 
  mutate(Taxon="Amphibians") %>% 
  rbind(as.data.frame(describe_posterior(mmu_bi, ci = 0.95, test="none")) %>% 
          mutate(Taxon="Birds")) %>% 
  rbind(as.data.frame(describe_posterior(mmu_f, ci = 0.95, test="none")) %>% 
          mutate(Taxon="Fish"))  %>% 
  rbind(as.data.frame(describe_posterior(mmu_ma, ci = 0.95, test="none")) %>% 
          mutate(Taxon="Mammals")) %>% 
  rbind(as.data.frame(describe_posterior(mmu_rep, ci = 0.95, test="none")) %>% 
          mutate(Taxon="Reptiles")) %>% 
  mutate(Parameter = gsub("b_", "", Parameter),
         Parameter = gsub("InvasivesppDgenes", "Invasive", Parameter),
         Parameter = gsub("HabitatdegradationDchange", "Habitat degradation", Parameter),
         Parameter = gsub("Habitatloss", "Habitat loss", Parameter),
         Parameter = gsub("Climatechange", "Climate change", Parameter),
         Parameter = gsub("noneNone", "None", Parameter),
         threat = gsub("pollution", "", Parameter),
         threat = gsub("habitatl", "", threat),
         threat = gsub("invasive", "", threat),
         threat = gsub("climatechange", "", threat),
         threat = gsub("exploitation", "", threat),
         threat = gsub("none", "", threat),
         threat = gsub("disease", "", threat),
         threat = gsub("\\.", "", threat),  
         threat=ifelse(grepl("None",Parameter), "None",
                       ifelse(grepl("No", Parameter), "No", 
                              threat))) 

# Save it 

write.csv2(TableS4, "TableS4.csv")

# Combine ----------------------------------------------------------------------

# Create title

title <- ggdraw() + 
  draw_label("Single threats",
             fontface = 'bold',size = 18,
             x = 0,
             hjust = 0) 

# Create subtitles

subtitle1 <- ggdraw() + 
  draw_label("System",
             size = 16,
             x = 0.5,
             vjust = 0) 

subtitle2 <- ggdraw() + 
  draw_label("Taxon",
             x = 0.5,
             size = 16,
             vjust = 0) 
subtitles <- plot_grid(subtitle1, subtitle2, nrow = 1)

# Create row 1

(row1 <- plot_grid(g3a+theme(legend.position = "none",
                             axis.title.x = element_blank()),
                   g3b+theme(legend.position = "none",
                             axis.title.x = element_blank()),
                   labels = "auto"))
# Create title 2

title2 <- ggdraw() + 
  draw_label("Interacting threats",
             fontface = 'bold',
             x = 0, size = 18,
             hjust = 0) 

# Create row 2

(row2 <- plot_grid(g3c+theme(legend.position = "none",
                             axis.title.x = element_blank()),
                   g3d+theme(legend.position = "none",
                             axis.title.x = element_blank()),
                   labels = c("c","d")))

# Greate two legends 

legends <- plot_grid(get_legend(g3a+theme(legend.position = "bottom",
                                          legend.text = element_text(size = 14))),
          get_legend(g3b+theme(legend.position = "bottom",
                               legend.text = element_text(size = 14))))

# Put them together 

(figure2<- plot_grid(title, subtitles, row1, title2, row2, legends,
                    rel_heights = c(0.1,0.1,1,.1,1,.1), ncol=1))


ggsave("Figure 2.pdf",figure2,
       width = 14, height = 10, path = ResultPath)

# Figure 3: Interaction frequency ##############################################

load(paste0(ResultPath,"/PopMultEffects.RData"))

# General non-additive plot

# First classify non-additive 

data_ad <- data_int2 %>% 
  mutate(interaction.type=ifelse(effect!="Additive", 
                                 "Non-additive", "Additive")) %>% 
  group_by(interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) 

# Change scale to make it easier for the plot


data_freq<- data_freq %>% mutate(freq=freq/100)

data_freq_taxon<- data_freq_taxon %>% mutate(freq=freq/100)

# Save them 

setwd(ResultPath)
write.csv2(data_freq, "Table S5.csv")
write.csv2(data_freq_taxon, "Table S6.csv")

# Panel a: System --------------------------------------------------------------

data <- data_freq %>% mutate(System=as.factor(System),
                             freq=freq*100, 
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

(g3a <- ggplot(data) +
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

# Panel b: Taxa ----------------------------------------------------------------

data <- data_freq_taxon %>% mutate(Taxon=as.factor(Taxon),
                             freq=freq*100,
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

(g3b <- ggplot(data) +
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


# Combine ----------------------------------------------------------------------

#Create first row 

(row1 <- plot_grid(g3a, g3b, labels="auto"))

# Get the legend

legend <- get_legend(g3a+theme(legend.position = "bottom",
                               legend.text = element_text(size=14)))
# Plot it 

(figure3<- plot_grid(row1, legend, 
                     nrow = 2, rel_heights = c(1,0.05)))


# Save it

ggsave("Figure 3.pdf", figure3, 
       width = 10, height = 8,
       path = ResultPath)

