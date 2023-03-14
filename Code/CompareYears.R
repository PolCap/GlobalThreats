# --------------------------------------------------------------------------------------- #
# - FILE NAME:   CompareYears.R         
# - DATE:        07/01/2023
# - DESCRIPTION: Supplementary analyses to test the sensitivity of our results
#                when increasing the number of studied years in our time series. 
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
library(patchwork)
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

# Load the results 

load(paste0(ResultPath, "/Results.RData"))
load(paste0(ResultPath, "/TenResults.RData"))
load(paste0(ResultPath, "/TwentyResults.RData"))

# Load functions for the plots 

source(paste0(CodePath,"/PlotFunctions.R"))

#Create a palette for the threats 

threat_palette<-c(met.brewer(name="Hokusai1", n=6, type="continuous"), "grey60")

# Color palette for taxons

taxon_pal <- wes_palette("Cavalcanti1", n = 5)

# System palette

system_pal <- c("#A1D6E2","#336B87", "#CC954E")

# Figure S8: Compare results for systems  --------------------------------------
# Panel a: Single stressors by system 5 data points ----------------------------

(gs8a <- ms1 %>%
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
        y = expression(paste("Population trend (", mu, ")",sep = "")),
        title = "5 data points") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.4,.4)) +
   scale_colour_manual("", values = system_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Panel b: Single stressors by system 10 years data ----------------------------

(gs8b <- ms2 %>%
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
        y = expression(paste("Population trend (", mu, ")",sep = "")),
        title = "10 data points") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.4,.4)) +
   scale_colour_manual("", values = system_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Panel c: Single stressors by system 20 years data ----------------------------

(gs8c <- ms3 %>%
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
   filter(number==1,
          System!="Freshwater"|.variable!="Disease",
          System!="Marine"|.variable!="Disease") %>% 
   ggplot(aes(x = reorder(.variable, -.lower),
              y = .value, colour=System, group=System)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar( aes(ymin=.lower, ymax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                  position=position_dodge(0.7))+
   geom_hline(yintercept = 0, lty = 2, size = 0.5) +
   labs(x = "",
        y = "",
        title = "20 data points") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.4,.4)) +
   scale_colour_manual("", values = system_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Panel d: Multiple stressors by system 5 years -------------------------------

(gs8d <- mmu_fr %>%
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
        y = expression(paste("Population trend (", mu, ")",sep = "")),
        title="5 data points") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.18,.18)) +
   scale_colour_manual("", values = system_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Panel e: Multiple stressors by system 10 years -------------------------------

(gs8e <- mmu_fr2 %>%
   gather_draws(`b_.*`, regex = TRUE) %>%
   median_qi(.width = .95) %>%
   mutate(System="Freshwater") %>% 
   rbind(mmu_mar2 %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(System="Marine")) %>% 
   rbind(mmu_ter2 %>%
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
        y = expression(paste("Population trend (", mu, ")",sep = "")),
        title="10 data points") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.18,.18)) +
   scale_colour_manual("", values = system_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))


# Panel f: Multiple stressors by system 20 years -------------------------------

(gs8f <- mmu_fr3 %>%
   gather_draws(`b_.*`, regex = TRUE) %>%
   median_qi(.width = .95) %>%
   mutate(System="Freshwater") %>% 
   rbind(mmu_mar3 %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(System="Marine")) %>% 
   rbind(mmu_ter3 %>%
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
        y = "",
        title="20 data points") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.18,.18)) +
   scale_colour_manual("", values = system_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Combine ----------------------------------------------------------------------

(figureS8 <- (gs8a+gs8b+gs8c)/(gs8d+gs8e+gs8f)+
   plot_annotation(tag_levels = "a")+
   plot_layout(guides = 'collect'))

ggsave("Figure S8.pdf", figureS8, 
       height = 10, width = 16, 
       path=ResultPath)

# Figure S9: Compare results by Taxon  -----------------------------------------
# Panel a: Single stressors by taxon 5 points ----------------------------------

(gs9a <- mt1 %>%
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
          Taxon!="Amphibians"|.variable!="Exploitation",
          Taxon!="Amphibians"|.variable!="Invasive",
          Taxon!="Amphibians"|.variable!="Disease",
          Taxon!="Amphibians"|.variable!="Pollution",
          Taxon!="Amphibians"|.variable!="Climate change",
          Taxon!="Reptiles"|.variable!="Disease",
          Taxon!="Reptiles"|.variable!="Pollution",
          Taxon!="Fish"|.variable!="Disease") %>%
   group_by(.variable, Taxon) %>% 
   mutate(m=min(.value)) %>% 
   ungroup() %>% 
   ggplot(aes(x = reorder(.variable, -.lower),
              y = .value, colour=Taxon, group=Taxon)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar(aes(ymin=.lower, ymax=.upper), 
                 width=0, alpha=0.9, size=1.3,
                 position=position_dodge(0.7))+
   geom_hline(yintercept = 0, lty = 2, size = 0.5) +
   labs(x = "",
        y = expression(paste("Population trend (", mu, ")",sep = "")), 
        title = "5 points data") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.4,.4)) +
   scale_colour_manual("", 
                       values = taxon_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Panel b: Single stressors by taxon 10 points ----------------------------------

(gs9b <- mt2 %>%
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
          Taxon!="Amphibians"|.variable!="Exploitation",
          Taxon!="Amphibians"|.variable!="Invasive",
          Taxon!="Amphibians"|.variable!="Disease",
          Taxon!="Amphibians"|.variable!="Pollution",
          Taxon!="Amphibians"|.variable!="Climate change",
          Taxon!="Reptiles"|.variable!="Disease",
          Taxon!="Reptiles"|.variable!="Pollution",
          Taxon!="Fish"|.variable!="Disease") %>%
   group_by(.variable, Taxon) %>% 
   mutate(m=min(.value)) %>% 
   ungroup() %>% 
   ggplot(aes(x = reorder(.variable, -.lower),
              y = .value, colour=Taxon, group=Taxon)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar(aes(ymin=.lower, ymax=.upper), 
                 width=0, alpha=0.9, size=1.3,
                 position=position_dodge(0.7))+
   geom_hline(yintercept = 0, lty = 2, size = 0.5) +
   labs(x = "",
        y = expression(paste("Population trend (", mu, ")",sep = "")), 
        title = "10 points data") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.4,.4)) +
   scale_colour_manual("", 
                       values = taxon_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Panel c: Single stressors by taxon 20 years ----------------------------------

(gs9c <- mt3 %>%
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
          Taxon!="Amphibians"|.variable!="Exploitation",
          Taxon!="Amphibians"|.variable!="Invasive",
          Taxon!="Amphibians"|.variable!="Disease",
          Taxon!="Amphibians"|.variable!="Pollution",
          Taxon!="Amphibians"|.variable!="Climate change",
          Taxon!="Reptiles"|.variable!="Disease",
          Taxon!="Reptiles"|.variable!="Pollution",
          Taxon!="Fish"|.variable!="Disease",
          Taxon!="Fish"|.variable!="Climate change",
          Taxon!="Fish"|.variable!="Pollution",
          Taxon!="Birds"|.variable!="Disease",
          Taxon!="Mammals"|.variable!="Invasive",
          Taxon!="Mammals"|.variable!="Climate change",
          Taxon!="Reptiles"|.variable!="Invasive") %>%
   ggplot(aes(x = reorder(.variable, -.lower),
              y = .value, colour=Taxon, group=Taxon)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar( aes(ymin=.lower, ymax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                  position=position_dodge(0.7))+
   geom_hline(yintercept = 0, lty = 2, size = 0.5) +
   labs(x = "",
        y = "",title = "20 data points") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.4,.4)) +
   scale_colour_manual("", 
                       values = taxon_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Panel d: Multiple stressors by taxon 5 years --------------------------------

(gs9d <- mmu_am %>%
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
        y = expression(paste("Population trend (", mu, ")",sep = "")),
        title="5 data points") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.6,.6)) +
   scale_colour_manual("", values = taxon_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Panel e: Multiple stressors by taxon 20 years --------------------------------

(gs9e <- mmu_am2 %>%
   gather_draws(`b_.*`, regex = TRUE) %>%
   median_qi(.width = .95) %>%
   mutate(Taxon="Amphibians") %>% 
   rbind(mmu_bi2 %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(Taxon="Birds")) %>% 
   rbind(mmu_f2 %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(Taxon="Fish")) %>% 
   rbind(mmu_ma2 %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(Taxon="Mammals")) %>% 
   rbind(mmu_rep2 %>%
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
        y = "",
        title="10 data points") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.6,.6)) +
   scale_colour_manual("", values = taxon_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Panel f: Multiple stressors by taxon 20 years --------------------------------

(gs9f <- mmu_am2 %>%
   gather_draws(`b_.*`, regex = TRUE) %>%
   median_qi(.width = .95) %>%
   mutate(Taxon="Amphibians") %>% 
   rbind(mmu_bi2 %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(Taxon="Birds")) %>% 
   rbind(mmu_f2 %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(Taxon="Fish")) %>% 
   rbind(mmu_ma2 %>%
           gather_draws(`b_.*`, regex = TRUE) %>%
           median_qi(.width = .95) %>%
           mutate(Taxon="Mammals")) %>% 
   rbind(mmu_rep2 %>%
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
        y = "",
        title="20 data points") +
   scale_y_continuous(labels = scaleFUN,limits = c(-.6,.6)) +
   scale_colour_manual("", values = taxon_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))


# Combine ----------------------------------------------------------------------

(figureS9 <- (gs9a+gs9b+gs9c)/(gs9d+gs9e+gs9f)+
   plot_annotation(tag_levels = "a")+
   plot_layout(guides = 'collect'))

ggsave("Figure S9.pdf", figureS9, 
       height = 10, width = 16,
       path = ResultPath)

