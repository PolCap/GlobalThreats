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

load(paste0(ResultPath, "/TenResults.RData"))
load(paste0(DataPath, "/TenYearsData.RData"))

# Change name invasive species 

pops_data10 <- pops_data10 %>% 
  mutate(threats = gsub("Invasive spp/genes", "Invasive", threats),)

# Load functions for the plots 

source(paste0(CodePath,"/PlotFunctions.R"))

# Figure S1: Conceptual interactions ############################################

# Simulate the data 

sim_data <- data.frame(value=c(-1, -2, -1.5,-1.5, -3.5, 0.2), 
                       upper=c(-0.5, -1.5,-0.5, -1, -3, 0.7),
                       lower=c(-1.5, -2.5, -2.5, -2, -4, -0.3),
                       id=c(1:6),
                       stressor= c("A", "B", "A+B","A:B", "A:B","A:B"), 
                       interaction=c("","", "", "Additive", "Synergistic", "Antagonistic"),
                       group=c("Individual threats","Individual threats", 
                               "Null model", 
                               "Interactive effects", "Interactive effects","Interactive effects")) %>% 
  mutate(group=factor(group, levels=c("Individual threats", "Null model", "Interactive effects"))) 

# We draw the arrows that we want to show

(figure1 <- sim_data %>% 
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
    geom_text(data=sim_data %>% filter(group=="Interactive effects", 
                                       interaction=="Additive"),
              aes(x = 1, y = upper+0.2),
              size = 4, color = "#694364", lineheight = .9,
              label = "Additive")+
    geom_curve(data=sim_data %>% filter(group=="Interactive effects", 
                                        interaction=="Additive"),
               aes(x = 1, y = upper, xend = 0.75, 
                   yend = value+0.01),
               arrow = arrow(length = unit(0.07, "inch")), size = 0.4,
               color = "#694364", curvature = -0.3)+
    geom_text(data=sim_data %>% filter(group=="Interactive effects", 
                                       interaction=="Antagonistic"),
              aes(x = 0.92, y = upper+0.2),
              size = 4, color = "#1E63B3", lineheight = .9,
              label = "Antagonistic")+
    geom_curve(data=sim_data %>% filter(group=="Interactive effects", 
                                        interaction=="Antagonistic"),
               aes(x = 0.92, y = upper, xend = 1.24, yend = value),
               arrow = arrow(length = unit(0.07, "inch")), size = 0.4,
               color = "#1E63B3", curvature = 0.3)+
    geom_text(data=sim_data %>% filter(group=="Interactive effects", 
                                       interaction=="Synergistic"),
              aes(x = 1.30, y = upper+0.2),
              size = 4, color = "#B32315", lineheight = .9,
              label = "Synergistic")+
    geom_curve(data=sim_data %>% filter(group=="Interactive effects", 
                                        interaction=="Synergistic"),
               aes(x = 1.30, y = upper, xend = 1.06, yend = value),
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

ggsave("Figure 1.pdf", figure1,
       path = ResultPath, height = 4, width = 8)

# Figure S2: threats across systems --------------------------------------------  

# One threat -----------------------------------------------------------------

(gs2a <- ms1 %>%
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

(gs2b <- ms1 %>%
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
   filter(number==2, .variable!="Exploitation:Disease"|System!="Freshwater",
          .variable!="Climate change:Disease"|System!="Marine",
          .variable!="Climate change:Disease"|System!="Terrestrial",
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

(gs2c <- ms1 %>%
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
   filter(number==3, .variable!="Habitat loss:Climate change:Disease"|System!="Marine",
          .variable!="Invasive:Climate change:Disease"|System!="Marine",
          .variable!="Invasive:Climate change:Disease"|System!="Terrestrial",
          .variable!="Pollution:Habitat loss:Invasive"|System!="Marine",
          .variable!="Pollution:Invasive:Exploitation"|System!="Freshwater",
          .variable!="Pollution:Climate change:Disease"|System!="Terrestrial",
          .variable!="Invasive:Climate change:Exploitation"|System!="Freshwater",
          .variable!="Invasive:Climate change:Exploitation"|System!="Terrestrial",
          .variable!="Invasive:Exploitation:Disease"|System!="Freshwater",
          .variable!="Pollution:Exploitation:Disease"|System!="Terrestrial",
          .variable!="Habitat loss:Invasive:Disease"|System!="Freshwater",
          .variable!="Pollution:Climate change:Exploitation"|System!="Terrestrial",
          .variable!="Pollution:Climate change:Exploitation"|System!="Freshwater",
          .variable!="Pollution:Invasive:Climate change"|System!="Terrestrial",
          .variable!="Pollution:Invasive:Climate change"|System!="Freshwater") %>% 
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
   theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm")))

# We combine them --------------------------------------------------------------

# First row

row1 <- plot_grid(gs2a+theme(legend.position = "none"),
                  gs2b+theme(legend.position = "none"),
                  gs2c+theme(legend.position = "none"),
                  nrow = 1,labels = "auto")

# get the legend

legend <- get_legend(gs2a+theme(legend.position = "bottom",
                                legend.text = element_text(size = 16)))

# Figure 1

(figs2 <- plot_grid(row1, legend,
                  rel_heights = c(1,.1),
                  ncol = 1))

# Save

ggsave("Fig.S2.pdf", figs2,
       height = 10, width = 16,
       path = ResultPath )

# Figure S3: threats across taxon ----------------------------------------------  

# One threat -------------------------------------------------------------------

(gs3a <- mt1 %>%
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
          Taxon!="Reptiles"|.variable!="Pollution",
          Taxon!="Reptiles"|.variable!="Climate change",
          Taxon!="Reptiles"|.variable!="Disease",
          Taxon!="Amphibians"|.variable!="Pollution",
          Taxon!="Amphibians"|.variable!="Climate change",
          Taxon!="Amphibians"|.variable!="Disease",
          Taxon!="Amphibians"|.variable!="Exploitation") %>% 
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
   scale_x_continuous(labels = scaleFUN, limits = c(-.2,.2)) +
   theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm")))

# Two threats ------------------------------------------------------------------

(gs3b <- mt1 %>%
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
          Taxon!="Amphibians"|.variable!="Climate change:Disease",
          Taxon!="Amphibians"|.variable!="Pollution:Invasive",
          Taxon!="Amphibians"|.variable!="Habitat loss:Exploitation",
          Taxon!="Amphibians"|.variable!="Pollution:Disease",
          Taxon!="Amphibians"|.variable!="Invasive:Climate change",
          Taxon!="Amphibians"|.variable!="Invasive:Exploitation",
          Taxon!="Amphibians"|.variable!="Climate change:Exploitation",
          Taxon!="Amphibians"|.variable!="Pollution:Exploitation",
          Taxon!="Fish"|.variable!="Exploitation:Disease",
          Taxon!="Fish"|.variable!="Habitat loss:Invasive",
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
   scale_x_continuous(labels = scaleFUN, limits = c(-.2,.2)) +
   theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm")))

# Three threats ----------------------------------------------------------------

(gs3c <- mt1 %>%
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
          Taxon!="Reptiles"|.variable!="Habitat loss:Climate change:Disease",
          Taxon!="Mammals"|.variable!="Habitat loss:Climate change:Disease",
          Taxon!="Fish"|.variable!="Habitat loss:Climate change:Disease",
          Taxon!="Reptiles"|.variable!="Pollution:Habitat loss:Disease",
          Taxon!="Fish"|.variable!="Pollution:Habitat loss:Disease",
          Taxon!="Amphibians"|.variable!="Pollution:Habitat loss:Disease",
          Taxon!="Reptiles"|.variable!="Pollution:Climate change:Exploitation",
          Taxon!="Fish"|.variable!="Pollution:Climate change:Exploitation",
          Taxon!="Amphibians"|.variable!="Pollution:Climate change:Exploitation",
          Taxon!="Reptiles"|.variable!="Invasive:Climate change:Exploitation",
          Taxon!="Mammals"|.variable!="Invasive:Climate change:Exploitation",
          Taxon!="Fish"|.variable!="Invasive:Climate change:Exploitation",
          Taxon!="Amphibians"|.variable!="Invasive:Climate change:Exploitation",
          Taxon!="Reptiles"|.variable!="Pollution:Habitat loss:Invasive",
          Taxon!="Fish"|.variable!="Pollution:Habitat loss:Invasive",
          Taxon!="Amphibians"|.variable!="Pollution:Habitat loss:Invasive",
          Taxon!="Reptiles"|.variable!="Habitat loss:Invasive:Climate change",
          Taxon!="Amphibians"|.variable!="Habitat loss:Invasive:Climate change",
          Taxon!="Amphibians"|.variable!="Invasive:Climate change:Exploitation",
          Taxon!="Reptiles"|.variable!="Pollution:Invasive:Climate change",
          Taxon!="Mammals"|.variable!="Pollution:Invasive:Climate change",
          Taxon!="Fish"|.variable!="Pollution:Invasive:Climate change",
          Taxon!="Amphibians"|.variable!="Pollution:Invasive:Climate change",
          Taxon!="Fish"|.variable!="Pollution:Exploitation:Disease",
          Taxon!="Amphibians"|.variable!="Pollution:Exploitation:Disease",
          Taxon!="Amphibians"|.variable!="Habitat loss:Invasive:Exploitation",
          Taxon!="Reptiles"|.variable!="Habitat loss:Exploitation:Disease",
          Taxon!="Fish"|.variable!="Habitat loss:Exploitation:Disease",
          Taxon!="Amphibians"|.variable!="Habitat loss:Exploitation:Disease",
          Taxon!="Birds"|.variable!="Pollution:Climate change:Exploitation") %>% 
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
   scale_x_continuous(labels = scaleFUN, limits = c(-.2,.2)) +
   theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm")))

# We combine them --------------------------------------------------------------

# First row

row1 <- plot_grid(gs3a+theme(legend.position = "none"),
                  gs3b+theme(legend.position = "none"),
                  gs3c+theme(legend.position = "none"),
                  nrow = 1,labels = "auto")

# get the legend

legend <- get_legend(gs3a+theme(legend.position = "bottom",
                                legend.text = element_text(size = 16)))

# Figure S2

(figs3 <- plot_grid(row1, legend,
                   rel_heights = c(1,.1),
                   ncol = 1))

# Save

ggsave("Fig.S3.pdf", figs3,
       height = 10, width = 16,
       path = ResultPath )



# Figure S4: Network of interactions ###########################################

load(paste0(ResultPath,"/PopMultEffects.RData"))

# First classify non-additive 

data_ad <- data_int2 %>% 
  mutate(interaction.type=ifelse(effect!="Additive", 
                                 "Non-additive", "Additive")) %>% 
  group_by(interaction.type) %>% 
  summarise(n=n()) %>% 
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) 


# Change scale to make it easier for the plot


data_freq<- data_freq %>% mutate(freq=freq/100)

data_freq_taxon<- data_freq_taxon %>% mutate(freq=freq/100)

# Save them 

setwd(ResultPath)
write.csv(data_freq, "Table S5.csv")
write.csv(data_freq_taxon, "Table S6.csv")

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

par(mar = c(0, 0, 0, 0), mgp = c(0, 0, 0), bg=NA)
plot(Net, edge.color=palette[(E(Net)$interaction.type=="Antagonistic")+1],
     vertex.label.family="Helvetica", 
     edge.curved=curve_multiple(Net), layout=l)

# Keep the plot

gs4a <- recordPlot()

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

par(mar = c(0, 0, 0, 0), mgp = c(0, 0, 0), bg=NA)
plot(Net2, edge.color=palette[(E(Net2)$interaction.type=="Antagonistic")+1],
     vertex.label.family="Helvetica",
     edge.curved=curve_multiple(Net2), layout=l2)

gs4b <- recordPlot()

# Combine ----------------------------------------------------------------------

# Second row 

(figS4 <- plot_grid(gs4a, 
                    gs4b,nrow = 1, 
                    labels = "auto"))

# Save it

ggsave("Figure S4.pdf", figS4, 
       width = 10, height = 8,
       path = ResultPath)

# Table S7: two way interactions ------------------------------------------------

# Two stressors 
# Re-structure the data set for the plot 

two_pairs<- gather_set_data(pairs, 1:2)

two_pairs <- two_pairs %>% 
  mutate(x = fct_inorder(x)) %>% 
  #filter(interaction.type!="Additive") %>% 
  distinct(id, .keep_all=T) %>% 
  filter(is.na(Tertiary))%>% arrange(interaction.type)

# Save it as a supplementary table 

setwd(ResultPath)

write.csv(two_pairs, "TableS7.csv")

# Re-structure the data set for the plot 

three_pairs<- gather_set_data(pairs, c(7,3))
three_pairs <- three_pairs %>% 
  mutate(x = fct_inorder(x)) %>%
  drop_na() %>% 
  distinct(id, .keep_all=T) %>% 
  arrange(interaction.type)

# Save it as a supplementary table 

setwd(ResultPath)

write.csv(three_pairs, "TableS8.csv")

# Across systems 

pairs <- data_prop %>% 
  drop_na(.variable) %>%
  group_by(System, .variable, interaction.type) %>% 
  summarise(n=n()) %>% 
  tidyr::complete(System, .variable, interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100)

# Now we separate the threats 

pairs <- pairs %>% separate(.variable, c("Primary", "Secondary", "Tertiary"), ":")

# Create a variable with the first and second

pairs <- pairs %>% mutate(Primary_Secondary=paste0(Primary, ":", Secondary))

# Two stressors 
# Re-structure the data set for the plot 

two_pairs<- gather_set_data(pairs, 1:2)

two_pairs <- two_pairs %>% 
  mutate(x = fct_inorder(x)) %>% 
  distinct(id, .keep_all=T) %>% 
  filter(is.na(Tertiary))

# Save it as a supplementary table 

setwd(ResultPath)

write.csv(two_pairs, "TableS9.csv")

# Three stressors 
# Re-structure the data set for the plot 

three_pairs<- gather_set_data(pairs, c(8,3))
three_pairs <- three_pairs %>% 
  mutate(x = fct_inorder(x)) %>%
  drop_na() %>% 
  distinct(id, .keep_all=T) %>% 
  arrange(interaction.type)

# Save it as a supplementary table 

setwd(ResultPath)

write.csv(three_pairs, "TableS10.csv")