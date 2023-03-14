# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupFigures.R         
# - DATE:        21/02/2022
# - DESCRIPTION: Supplementary figures 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(ggplot2) #Now I load the latest version of ggplot2
library(tidybayes)
library(dplyr)
library(rstan)
library(rstanarm)
library(brms)
library(bayesplot)
library(bayestestR)
library(ggdist)
library(magrittr)
library(forcats)
library(modelr)
library(emmeans)
library(cowplot)
library(rphylopic)
library(RCurl)
library(stringr)
library(wesanderson)
library(tidyr)
library(patchwork)
library(igraph)
library(MetBrewer)

# Colour palettes

taxon_pal <- wes_palette("Cavalcanti1", n = 5)
system_pal <- c("#A1D6E2","#336B87", "#CC954E")


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
                  axis.ticks = element_line(color="black")))


# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 

load(paste0(ResultPath, "/Results.RData"))
load(paste0(DataPath, "/FiveYearsData.RData"))

# Change name invasive species 

pop_data <- pop_data %>% 
  mutate(threats = gsub("Invasive spp/genes", "Invasive", threats),)

# Load functions for the plots 

source(paste0(CodePath,"/PlotFunctions.R"))

# Figure S1: Conceptual interactions ############################################

# Simulate the data 

sim_data <- data.frame(value=c(-1,   -2,  -1.5,-0.5, -2.5, -1.5,-3.5, 0.2), 
                       upper=c(-0.5, -1.5,-0.5, 0,-2,  -1,  -3,   0.7),
                       lower=c(-1.5, -2.5,-2.5,-1, -3, -2,  -4,   -0.3),
                       id=c(1:8),
                       stressor= c("A", "B", "A+B","A:B","A:B","A:B", "A:B","A:B"), 
                       interaction=c("","", "", "Additive","Additive","Additive", "Synergistic", "Antagonistic"),
                       group=c("Individual threats","Individual threats", 
                               "Null model", 
                               "Interactive effects","Interactive effects","Interactive effects",
                               "Interactive effects","Interactive effects")) %>% 
  mutate(group=factor(group, levels=c("Individual threats", "Null model", "Interactive effects"))) 

# We draw the arrows that we want to show

(figureS1 <- sim_data %>% 
    ggplot(aes(y=value, x=stressor, 
               group=id, colour=interaction))+
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = -0.5, linetype="dashed") +
    geom_hline(yintercept = -2.5, linetype="dashed") +
    geom_point(size=5,
               position = position_dodge(0.9)) +
    geom_errorbar(aes(ymax=upper, ymin=lower), 
                  width=0.15,position = position_dodge(0.9))+
    scale_colour_manual(name = "",
                        values = c("grey45",
                                   "#694364",
                                   "#1E63B3",
                                   "#B32315"))+
    facet_wrap(~group,scales = "free_x")+
    geom_text(data=sim_data %>% filter(id==6),
              aes(x = 0.75, y = value),
              size = 4, color = "#694364", lineheight = .9,
              label = "Additive")+
    geom_text(data=sim_data %>% filter(group=="Interactive effects", 
                                       interaction=="Antagonistic"),
              aes(x = 0.92, y = upper+0.2),
              size = 4, color = "#1E63B3", lineheight = .9,
              label = "Antagonistic")+
    geom_curve(data=sim_data %>% filter(group=="Interactive effects", 
                                        interaction=="Antagonistic"),
               aes(x = 0.92, y = upper, 
                   xend = 1.31, yend = value),
               arrow = arrow(length = unit(0.07, "inch")), size = 0.4,
               color = "#1E63B3", curvature = 0.3)+
    geom_text(data=sim_data %>% filter(group=="Interactive effects", 
                                       interaction=="Synergistic"),
              aes(x = 0.75, y = lower),
              size = 4, color = "#B32315", lineheight = .9,
              label = "Synergistic")+
    geom_curve(data=sim_data %>% filter(group=="Interactive effects", 
                                        interaction=="Synergistic"),
               aes(x = 0.75, y = lower+0.2, 
                   xend = 1.129, yend = value),
               arrow = arrow(length = unit(0.07, "inch")), size = 0.4,
               color = "#B32315", curvature = -0.3)+
    labs(x="Threats", y="Population trend") +
    ylim(-4, 1.2) +
    theme(panel.spacing = unit(0, units="cm"),
          legend.position = "none",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          plot.margin = unit(c(0, 1, 0, 1), units = "cm")))

# Save it

ggsave("Figure S1.pdf", figureS1,
       path = ResultPath, height = 4, width = 8)

# Figure S6: threats across systems --------------------------------------------  

# One threat -----------------------------------------------------------------

(gs1a <- ms1 %>%
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
   ggplot(aes(y = reorder(.variable, -.value),
              x = .value, colour=System, group=System)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar( aes(xmin=.lower, xmax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                  position=position_dodge(0.7))+
   geom_vline(xintercept = 0, lty = 2, size = 0.5) +
   labs(y = "",
        x = expression(paste("Population trend (", mu, ")",sep = ""))) +
   scale_y_discrete(labels = scales::wrap_format(12)) +
   scale_colour_manual("", values = system_pal)+
   scale_x_continuous(labels = scaleFUN) +
   theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm")))

# Two threats ------------------------------------------------------------------

(gs1b <- ms1 %>%
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
    filter(number==2, 
           .variable!="Exploitation:Disease"|System!="Freshwater",
           .variable!="Pollution:Disease"|System!="Terrestrial",
           .variable!="Climate change:Disease"|System!="Terrestrial",
           .variable!="Climate change:Disease"|System!="Marine",
           .variable!="Pollution:Invasive"|System!="Marine",
           .variable!="Pollution:Invasive"|System!="Freshwater",
           .variable!="Invasive:Disease"|System!="Freshwater",
           .variable!="Invasive:Disease"|System!="Marine",
           .variable!="Invasive:Climate change"|System!="Freshwater",
           .variable!="Invasive:Climate change"|System!="Marine") %>% 
   ggplot(aes(y = reorder(.variable, -.value),
              x = .value, colour=System, group=System)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar( aes(xmin=.lower, xmax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                  position=position_dodge(0.7))+
   geom_vline(xintercept = 0, lty = 2, size = 0.5) +
   labs(y = "",
        x = expression(paste("Population trend (", mu, ")",sep = ""))) +
   scale_y_discrete(labels = scales::wrap_format(12)) +
   scale_colour_manual("", values = system_pal)+
   scale_x_continuous(labels = scaleFUN) +
   theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm")))

# Three threats ----------------------------------------------------------------

(gs1c <- ms1 %>%
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
   filter(number==3, 
          .variable!="Habitat loss:Climate change:Disease"|System!="Marine",
          .variable!="Invasive:Climate change:Disease"|System!="Marine",
          .variable!="Invasive:Climate change:Disease"|System!="Terrestrial",
           .variable!="Pollution:Habitat loss:Invasive"|System!="Marine",
   #        .variable!="Pollution:Invasive:Exploitation"|System!="Freshwater",
           .variable!="Pollution:Climate change:Disease"|System!="Terrestrial",
           .variable!="Invasive:Climate change:Exploitation"|System!="Freshwater",
   .variable!="Pollution:Exploitation:Disease"|System!="Terrestrial",
   .variable!="Pollution:Exploitation:Disease"|System!="Freshwater",
           .variable!="Invasive:Climate change:Exploitation"|System!="Terrestrial",
           .variable!="Invasive:Exploitation:Disease"|System!="Freshwater",
   .variable!="Invasive:Exploitation:Disease"|System!="Terrestrial",
   #        .variable!="Pollution:Exploitation:Disease"|System!="Terrestrial",
           .variable!="Habitat loss:Invasive:Disease"|System!="Freshwater",
          .variable!="Pollution:Climate change:Exploitation"|System!="Terrestrial",
          .variable!="Pollution:Climate change:Exploitation"|System!="Freshwater",
           .variable!="Pollution:Invasive:Climate change"|System!="Terrestrial",
           .variable!="Pollution:Invasive:Climate change"|System!="Freshwater",
           .variable!="Habitat loss:Invasive:Climate change"|System!="Marine") %>% 
   #        .variable!="Pollution:Invasive:Exploitation"|System!="Terrestrial") %>% 
   ggplot(aes(y = reorder(.variable, -.value),
              x = .value, colour=System, group=System)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar( aes(xmin=.lower, xmax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                  position=position_dodge(0.7))+
   geom_vline(xintercept = 0, lty = 2, size = 0.5) +
   labs(y = "",
        x = expression(paste("Population trend (", mu, ")",sep = ""))) +
   scale_y_discrete(labels = scales::wrap_format(25)) +
   scale_colour_manual("", values = system_pal)+
   scale_x_continuous(labels = scaleFUN) +
   theme(plot.margin = unit(c(0, 0.75, 0, 0), "cm")))

# We combine them --------------------------------------------------------------

# First row

row1 <- plot_grid(gs1a+theme(legend.position = "none"),
                  gs1b+theme(legend.position = "none"),
                  gs1c+theme(legend.position = "none"),
                  nrow = 1,labels = "auto")

# get the legend

legend <- get_legend(gs1a+theme(legend.position = "bottom",
                                legend.text = element_text(size = 16)))

# Figure 

(figS6 <- plot_grid(row1, legend,
                  rel_heights = c(1,.1),
                  ncol = 1))

# Save

ggsave("Figure S6.pdf", figS6,
       height = 10, width = 17,
       path = ResultPath )

# Figure S7: threats across taxon ----------------------------------------------  

# One threat -------------------------------------------------------------------

(gs2a <- mt1 %>%
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
   ggplot(aes(y = reorder(.variable, -.value),
              x = .value, colour=Taxon, group=Taxon)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar( aes(xmin=.lower, xmax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                  position=position_dodge(0.7))+
   geom_vline(xintercept = 0, lty = 2, size = 0.5) +
   labs(x = "",
        y = expression(paste("Population trend (", mu, ")",sep = ""))) +
   scale_y_discrete(labels = scales::wrap_format(25)) +
   scale_colour_manual("", values = taxon_pal)+
   scale_x_continuous(labels = scaleFUN) +
   theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm")))

# Two threats ------------------------------------------------------------------

(gs2b <- mt1 %>%
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
          number = str_count(.variable,"[A-Z]"))%>%
   filter(number==2, 
           Taxon!="Reptiles"|.variable!="Habitat loss:Invasive",
           Taxon!="Reptiles"|.variable!="Exploitation:Disease",
           Taxon!="Reptiles"|.variable!="Climate change:Disease",
           Taxon!="Reptiles"|.variable!="Pollution:Invasive",
           Taxon!="Reptiles"|.variable!="Invasive:Disease",
           Taxon!="Reptiles"|.variable!="Habitat loss:Disease",
           Taxon!="Reptiles"|.variable!="Pollution:Disease",
           Taxon!="Amphibians"|.variable!="Exploitation:Disease",
           Taxon!="Reptiles"|.variable!="Habitat loss:Disease",
           Taxon!="Amphibians"|.variable!="Climate change:Disease",
           Taxon!="Amphibians"|.variable!="Pollution:Invasive",
           Taxon!="Amphibians"|.variable!="Habitat loss:Exploitation",
           Taxon!="Amphibians"|.variable!="Pollution:Disease",
           Taxon!="Amphibians"|.variable!="Invasive:Climate change",
           Taxon!="Amphibians"|.variable!="Invasive:Exploitation",
           Taxon!="Amphibians"|.variable!="Climate change:Exploitation",
           Taxon!="Amphibians"|.variable!="Pollution:Exploitation",
           Taxon!="Fish"|.variable!="Exploitation:Disease",
           Taxon!="Fish"|.variable!="Pollution:Invasive",
           Taxon!="Fish"|.variable!="Pollution:Disease",
           Taxon!="Fish"|.variable!="Climate change:Disease",
           Taxon!="Fish"|.variable!="Invasive:Disease",
           Taxon!="Fish"|.variable!="Habitat loss:Disease",
           Taxon!="Fish"|.variable!="Invasive:Climate change",
           Taxon!="Birds"|.variable!="Pollution:Invasive",
           Taxon!="Birds"|.variable!="Invasive:Disease",
           Taxon!="Mammals"|.variable!="Pollution:Disease",
           Taxon!="Mammals"|.variable!="Invasive:Exploitation",
           Taxon!="Mammals"|.variable!="Climate change:Disease",
           Taxon!="Mammals"|.variable!="Invasive:Climate change") %>% 
   ggplot(aes(y = reorder(.variable, -.value),
              x = .value, colour=Taxon, group=Taxon)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar( aes(xmin=.lower, xmax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                  position=position_dodge(0.7))+
   geom_vline(xintercept = 0, lty = 2, size = 0.5) +
   labs(y = "",
        x = expression(paste("Population trend (", mu, ")",sep = ""))) +
   scale_y_discrete(labels = scales::wrap_format(25)) +
   scale_colour_manual("", values = taxon_pal)+
   scale_x_continuous(labels = scaleFUN) +
   theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm")))

# Three threats ----------------------------------------------------------------

(gs2c <- mt1 %>%
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
          number = str_count(.variable,"[A-Z]"))%>%
   filter(number==3, 
           Taxon!="Reptiles"|.variable!="Invasive:Climate change:Disease",
           Taxon!="Mammals"|.variable!="Invasive:Climate change:Disease",
           Taxon!="Fish"|.variable!="Invasive:Climate change:Disease",
           Taxon!="Birds"|.variable!="Invasive:Climate change:Disease",
           Taxon!="Reptiles"|.variable!="Pollution:Climate change:Disease",
           Taxon!="Mammals"|.variable!="Pollution:Climate change:Disease",
           Taxon!="Fish"|.variable!="Pollution:Climate change:Disease",
          Taxon!="Reptiles"|.variable!="Pollution:Invasive:Exploitation",
          Taxon!="Mammals"|.variable!="Pollution:Invasive:Exploitation",
          Taxon!="Fish"|.variable!="Pollution:Invasive:Exploitation",
          Taxon!="Amphibians"|.variable!="Pollution:Invasive:Exploitation",
           Taxon!="Amphibians"|.variable!="Habitat loss:Climate change:Exploitation",
           Taxon!="Reptiles"|.variable!="Habitat loss:Invasive:Disease",
           Taxon!="Mammals"|.variable!="Habitat loss:Invasive:Disease",
           Taxon!="Fish"|.variable!="Habitat loss:Invasive:Disease",
           Taxon!="Reptiles"|.variable!="Pollution:Habitat loss:Climate change",
           Taxon!="Amphibians"|.variable!="Pollution:Habitat loss:Climate change",
           Taxon!="Fish"|.variable!="Pollution:Habitat loss:Climate change",
           Taxon!="Reptiles"|.variable!="Invasive:Exploitation:Disease",
           Taxon!="Mammals"|.variable!="Invasive:Exploitation:Disease",
           Taxon!="Fish"|.variable!="Invasive:Exploitation:Disease",
           Taxon!="Amphibians"|.variable!="Invasive:Exploitation:Disease",
           Taxon!="Reptiles"|.variable!="Habitat loss:Climate change:Disease",
           Taxon!="Mammals"|.variable!="Habitat loss:Climate change:Disease",
           Taxon!="Fish"|.variable!="Habitat loss:Climate change:Disease",
          Taxon!="Reptiles"|.variable!="Pollution:Habitat loss:Disease",
          Taxon!="Fish"|.variable!="Pollution:Habitat loss:Disease",
          Taxon!="Amphibians"|.variable!="Pollution:Habitat loss:Disease",
          Taxon!="Reptiles"|.variable!="Pollution:Climate change:Exploitation",
          Taxon!="Mammals"|.variable!="Pollution:Climate change:Exploitation",
          Taxon!="Fish"|.variable!="Pollution:Climate change:Exploitation",
          Taxon!="Amphibians"|.variable!="Pollution:Climate change:Exploitation",
           Taxon!="Reptiles"|.variable!="Pollution:Habitat loss:Invasive",
           Taxon!="Mammals"|.variable!="Pollution:Habitat loss:Invasive",
           Taxon!="Fish"|.variable!="Pollution:Habitat loss:Invasive",
           Taxon!="Amphibians"|.variable!="Pollution:Habitat loss:Invasive",
           Taxon!="Amphibians"|.variable!="Habitat loss:Invasive:Exploitation",
          Taxon!="Reptiles"|.variable!="Habitat loss:Invasive:Climate change",
          Taxon!="Amphibians"|.variable!="Habitat loss:Invasive:Climate change",
           Taxon!="Amphibians"|.variable!="Invasive:Climate change:Exploitation",
           Taxon!="Mammals"|.variable!="Invasive:Climate change:Exploitation",
           Taxon!="Reptiles"|.variable!="Invasive:Climate change:Exploitation",
           Taxon!="Fish"|.variable!="Invasive:Climate change:Exploitation",
          Taxon!="Reptiles"|.variable!="Pollution:Invasive:Climate change",
          Taxon!="Mammals"|.variable!="Pollution:Invasive:Climate change",
          Taxon!="Fish"|.variable!="Pollution:Invasive:Climate change",
          Taxon!="Amphibians"|.variable!="Pollution:Invasive:Climate change",
          Taxon!="Fish"|.variable!="Pollution:Exploitation:Disease",
          Taxon!="Amphibians"|.variable!="Pollution:Exploitation:Disease",
           Taxon!="Amphibians"|.variable!="Habitat loss:Invasive:Exploitation",
           Taxon!="Reptiles"|.variable!="Habitat loss:Exploitation:Disease",
           Taxon!="Fish"|.variable!="Habitat loss:Exploitation:Disease",
           Taxon!="Mammals"|.variable!="Habitat loss:Exploitation:Disease",
           Taxon!="Amphibians"|.variable!="Habitat loss:Exploitation:Disease",
          Taxon!="Birds"|.variable!="Pollution:Climate change:Exploitation",
          Taxon!="Amphibians"|.variable!="Invasive:Climate change:Disease") %>% 
   ggplot(aes(y = reorder(.variable, -.value),
              x = .value, colour=Taxon, group=Taxon)) +
   geom_point(size=4,position =position_dodge(.7)) +
   geom_errorbar( aes(xmin=.lower, xmax=.upper), 
                  width=0, alpha=0.9, size=1.3,
                  position=position_dodge(0.7))+
   geom_vline(xintercept = 0, lty = 2, size = 0.5) +
   labs(y = "",
        x = expression(paste("Population trend (", mu, ")",sep = ""))) +
   scale_y_discrete(labels = scales::wrap_format(25)) +
   scale_colour_manual("", values = taxon_pal)+
   scale_x_continuous(labels = scaleFUN) +
   theme(plot.margin = unit(c(0, 0.75, 0, 0), "cm")))

# We combine them --------------------------------------------------------------

# First row

row1 <- plot_grid(gs2a+theme(legend.position = "none"),
                  gs2b+theme(legend.position = "none"),
                  gs2c+theme(legend.position = "none"),
                  nrow = 1,labels = "auto")

# get the legend

legend <- get_legend(gs2a+theme(legend.position = "bottom",
                                legend.text = element_text(size = 16)))

# Figure S7

(figs7 <- plot_grid(row1, legend,
                   rel_heights = c(1,.1),
                   ncol = 1))

# Save

ggsave("Figure S7.pdf", figs7,
       height = 10, width = 16,
       path = ResultPath )

# Figure S10: Network of interactions ##########################################

load(paste0(ResultPath,"/PopMultEffects.RData"))

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

# Panel a: network of interactions ---------------------------------------------

# Now we are going to estimate the proportion of each interaction type 
# according to each group

data_non <- data_freq %>% 
  mutate(interaction.type=ifelse(interaction.type!="Additive", 
                                 "Non-additive", "Additive")) %>% 
  group_by(variable, System, interaction.type) %>% 
  summarise(number=sum(n)) %>% 
  group_by(variable,System) %>% 
  summarise(interaction.type,
            freq=(number/sum(number))) %>% 
  ungroup()

# Subset the pairs of interactions

pairs <- data_prop %>% 
  drop_na(.variable) %>%
  #filter(interaction.type!="Additive") %>% 
  group_by(.variable, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>% 
  # Now we separate the threats 
  separate(.variable, c("Primary", "Secondary", "Tertiary"), ":") %>% 
  # Create a variable with the first and second
  mutate(Primary_Secondary=paste0(Primary, ":", Secondary))

# Subset only two threats 

two <- pairs %>% filter(is.na(Tertiary), 
                        interaction.type!="Additive") %>% 
  select(Primary, Secondary, interaction.type, freq)

# Create the network

Net <- graph_from_data_frame(d= two, 
                             directed=F) 

# Generate colors based on threats

V(Net)$color <- c("#DF865D", "#BF5B4D", "#6D2F20",
                  "#A5B88D", "#E7B266", "#224B5E")

# Compute node degrees (number of interacting threats) and use that to set node size

deg <- degree(Net, mode="all")
V(Net)$size <- deg*5

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(Net)$label.color <- "black"
#V(Net)$label <- NA

# Set edge width based on weight:
E(Net)$width <- E(Net)$freq/2

# Delete the links with 0 freq

Net <- delete.edges(Net, which(E(Net)$freq==0))

# Set a new palette for the interactions

palette <- c("#B32315", "#1E63B3")

# Set the layout

l <- layout_with_kk(Net)

# Interaction types 
dev.off()

par(mar = c(0, 0, 0, 0), mgp = c(0, 0, 0), bg=NA)
plot(Net, edge.color=palette[(E(Net)$interaction.type=="Antagonistic")+1],
     vertex.label.family="Helvetica", 
     edge.curved=curve_multiple(Net), layout=l)

# Keep the plot

gs3a <- recordPlot()

# Panel d: Three way interactions ----------------------------------------------

# Three pairs top part

three_up <- pairs %>% filter(!is.na(Tertiary), 
                             interaction.type!="Additive") %>% 
  select(Primary, Secondary, interaction.type, freq)

# Three pairs bottom part

three_bot <- pairs %>% filter(!is.na(Tertiary), 
                              interaction.type!="Additive") %>% 
  select(Secondary, Tertiary, interaction.type, freq) %>% 
  rename(Primary=Secondary, 
         Secondary=Tertiary)

# Join them 

three <- rbind(three_up, three_bot)

# Create the network

Net2 <- graph_from_data_frame(d= three, directed=F) 

# Generate colors based on threats

V(Net2)$color <- c("#6D2F20", "#A5B88D", "#E7B266", "#DF865D","#BF5B4D","#224B5E")

# Compute node degrees (number of interacting threats) and use that to set node size

deg <- degree(Net2, mode="all")
V(Net2)$size <- deg*2

# The labels are currently node IDs.

V(Net2)$label.color <- "black"

# Set edge width based on weight:

E(Net2)$width <- E(Net2)$freq/2

# Delete frequency 0

Net2 <- delete.edges(Net2, which(E(Net2)$freq==0))

# Set the layout for multiple edges

l2 <- layout_with_kk(Net2)

# Plot

dev.off()
par(mar = c(0, 0, 0, 0), mgp = c(0, 0, 0), bg=NA)
plot(Net2, edge.color=palette[(E(Net2)$interaction.type=="Antagonistic")+1],
     vertex.label.family="Helvetica",
     edge.curved=curve_multiple(Net2), layout=l2)

gs3b <- recordPlot()

# Combine ----------------------------------------------------------------------

# Second row 

(figS10 <- plot_grid(gs3a, 
                    gs3b,nrow = 1, 
                    labels = "auto"))

# Save it

ggsave("Figure S10.pdf", figS10, 
       width = 10, height = 8,
       path = ResultPath)

# Table S7: two way interactions ------------------------------------------------

# Two stressors 
# Re-structure the data set for the plot 

two_pairs<- gather_set_data(pairs, 1:2)

two_pairs <- two_pairs %>% 
  mutate(x=as.factor(x),
         x = fct_inorder(x),
         freq=freq/100) %>% 
  #filter(interaction.type!="Additive") %>% 
  distinct(id, .keep_all=T) %>% 
  filter(is.na(Tertiary))%>% arrange(interaction.type)

# Save it as a supplementary table 

setwd(ResultPath)

write.csv2(two_pairs, "TableS7.csv")

# Now we plot it

# (g4c<- two_pairs %>% 
#     filter(interaction.type!="Additive") %>% 
#     ggplot(aes(axis1 = Primary, axis2 = Secondary, 
#            y = freq)) +
#   geom_alluvium(aes(fill = interaction.type),width = .55,
#                 curve_type = "cubic", alpha=0.8) +
#   geom_stratum(width = .55) +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum)), size=5) +
#   scale_fill_manual(name = "",
#                     values = c(Antagonistic="#1E63B3",
#                                Additive="#694364",
#                                Synergy="#B32315"))+
#     theme_void()+
#     guides(fill = "none")+
#     theme(plot.margin = unit(c(0,0,0,0), "cm"))) 

# (g3d <- ggplot(two_pairs, aes(x, id = id, split = y, value = freq)) +
#     geom_parallel_sets(aes(fill = interaction.type), 
#                        alpha = 0.8, axis.width = 0.6) +
#     geom_parallel_sets_axes(axis.width = 0.6, 
#                             fill="white", colour="black") +
#     geom_parallel_sets_labels(angle = 0, colour="black")+
#     scale_fill_manual(name = "",
#                       values = c(Antagonistic="#1E63B3",
#                                  Additive="#694364",
#                                  Synergy="#B32315"))+
#     theme_void() + 
#     theme(plot.margin = margin(-0.5, -3, -0.5,-3, "cm"))+
#     guides(fill = "none"))


# Re-structure the data set for the plot 

three_pairs<- gather_set_data(pairs, c(7,3))
three_pairs <- three_pairs %>% 
  mutate(x=as.factor(x),
         x = fct_inorder(x),
         freq=freq/100)  %>%
  drop_na() %>% 
  distinct(id, .keep_all=T) %>% 
  arrange(interaction.type)

# Save it as a supplementary table 

setwd(ResultPath)

write.csv2(three_pairs, "TableS8.csv")

#Now we plot it

# (g4d<- three_pairs %>% 
#     filter(interaction.type!="Additive") %>% 
#     mutate(interaction.type=factor(interaction.type, levels=c("Antagonistic", "Synergy"))) %>% 
#     ggplot(aes(axis1 = Primary_Secondary, axis2 = Tertiary, 
#                   y = freq)) +
#     geom_alluvium(aes(fill = interaction.type),
#                   width = 0.65,
#                   curve_type = "cubic", alpha=0.8) +
#     geom_stratum(width = 0.65) +
#     geom_text(stat = "stratum",
#               aes(label = after_stat(stratum)),size=5) +
#     scale_fill_manual(name = "",
#                       values = c(Antagonistic="#1E63B3",
#                                  Synergy="#B32315"))+
#     theme_void()+
#     theme(legend.position = "none",
#           plot.margin = unit(c(0,0,0,0), "cm"))) 


# Figure S12: Sensitivity interactions #########################################

load(paste0(ResultPath,"/PopMeanEffects.RData"))

# Color palette for taxons

taxon_pal <- wesanderson::wes_palette("Cavalcanti1", n = 5)

# System palette

system_pal <- c("#A1D6E2","#336B87", "#CC954E")

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

(gs12a <- ggplot(data) +
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

(gs12b <- ggplot(data) +
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

(row1 <- plot_grid(gs12a, gs12b, labels="auto"))

# Get the legend

legend <- get_legend(gs12a+theme(legend.position = "bottom",
                               legend.text = element_text(size=14)))
# Plot it 

(figS12<- plot_grid(row1, legend, 
                       nrow = 2, rel_heights = c(1,0.05)))


# Save it

ggsave("Figure S12.pdf", figS12, 
       width = 10, height = 8,
       path = ResultPath)

# Tables S9 & S10 --------------------------------------------------------------

# Save them 

setwd(ResultPath)
write.csv2(data_freq %>% mutate(freq=freq/100), "Table S9.csv")
write.csv2(data_freq_taxon%>% mutate(freq=freq/100), "Table S10.csv")
