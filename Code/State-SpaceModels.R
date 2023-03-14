# --------------------------------------------------------------------------------------- #
# - FILE NAME:   State-SpaceModels.R         
# - DATE:        02/07/2020
# - DESCRIPTION: Code to estimate the population change using state-space models 
#                following Dennis et al. (2014)
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(piecewiseSEM)
library(tidyverse)
library(data.table)
library(zoo)
library(dplyr)
library(expss)
library(MASS)
library(VGAM)
library(goeveg)
library(RcppRoll)
library(patchwork)

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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Read in data

load(paste0(DataPath, "/LivingPlanetData2.RData"))

# Estimate population change

dd_long2 <- dd_long %>%
  group_by(ID) %>%  
  # Calculate population change
  mutate(popchange=log(Count+(max(Count, na.rm=T) /100)),
         #ifelse(Count<1, log(Count/(max(Count,na.rm = T))*100+1), 
         #        log(Count+1)),
         Protected_status=gsub(" .*", "", Protected_status)) %>% 
  # Remove any groupings we have created in the pipe
  ungroup()

# Add the biogeographical region

dd_long2 <- dd_long2 %>%
  mutate(T_realm=na_if(T_realm, "NULL"),
         M_realm=na_if(M_realm, "NULL"),
         FW_realm=na_if(FW_realm, "NULL"),
         Realm=coalesce(T_realm, M_realm, FW_realm))

# Find the data ids

dd_final <- dd_long2 %>%
  group_by(ID) %>%
  drop_na(popchange) %>%
  mutate(start=min(Year),
         end=max(Year),
        Duration=end-start) %>% 
  filter(Duration>9)

# Spread the data again

pops <- dd_final %>%
  mutate(year = as.factor(as.character(Year))) %>% 
  dplyr::select(ID, SpeciesName, Class, Order,
                System, Latitude, Longitude,
                Region,
                Protected_status, n.threat, 
                Primary_threat,
                Secondary_threat,
                Tertiary_threat,
                Managed,
                Realm,
                threats, Duration, 
                Year, year, popchange) %>%
  drop_na(popchange) %>%
  group_by(ID) %>%
  dplyr::select(-Year) %>%
  filter(length(unique(year))>=5) %>% 
  pivot_wider(names_from = year, values_from = popchange)

pops_data <- pops %>%
  dplyr::select(ID, SpeciesName, Class, Order,
                System, Latitude, Longitude,
                Region,
                Protected_status, n.threat, 
                Primary_threat,
                Secondary_threat,
                Tertiary_threat,
                Managed,
                Realm,
                threats, Duration) 

# Compile all time-series into a list to estimate state-space population trends##
pop <- list()
for(i in 1:nrow(pops)){
  data <- pops[i, 18:85]
  Year <- colnames(data)[is.na(data) == FALSE]
  data <- data[Year]
  N <- as.vector(t(data))
  pop[[i]] <- data.frame(pops$ID[i], Year, N)
}

# Load the functions to estimate state-space trends from Dennis et al. 2014

# source("State-space modelling.R")
source("Humbert_function.R")

# Loop through every time-series, estimating mu from state-space (REML) models ----

for(i in 1:nrow(pops)) {
  Observed.t <- as.numeric(pop[[i]]$N)
  Time.t <- as.numeric(pop[[i]]$Year)
  print(i)
  
  T.t = Time.t - Time.t[1] # Time starts at zero
  
  Y.t <-  Observed.t # The observations
  q <-  length(Y.t) - 1 # Number of time series transitions
  qp1 <-  q + 1 # q+1
  S.t <-  T.t[2:qp1] - T.t[1:q] # Time intervals
  m <-  rep(1, qp1) # Room for Kalman means
  v <-  rep(1, qp1) # Room form variances for Kalman calculations
  
  # Calculating EGOE AND EGPN estimates for use as initial values ----
  
  # The EGOE estimates
  Ybar <-  mean(Y.t)
  Tbar <-  mean(T.t)
  mu.egoe <-  sum((T.t - Tbar)*(Y.t - Ybar))/sum((T.t - Tbar)*(T.t - Tbar))
  x0.egoe <-  Ybar - mu.egoe*Tbar
  ssq.egoe <-  0
  Yhat.egoe <-  x0.egoe + mu.egoe*T.t
  tsq.egoe <-  sum((Y.t - Yhat.egoe)*(Y.t - Yhat.egoe))/(q - 1)
  
  # The EGPN estimates
  Ttr <-  sqrt(S.t)
  Ytr <-  (Y.t[2:qp1] - Y.t[1:q])/Ttr
  mu.egpn <-  sum(Ttr*Ytr)/sum(Ttr*Ttr)
  Ytrhat <-  mu.egpn*Ttr
  ssq.egpn <-  sum((Ytr - Ytrhat)*(Ytr - Ytrhat))/(q - 1)
  tsq.egpn <-  0
  x0.egpn <-  Y.t[1]
  
  # Initial values for EGSS are averages of EGOE and EGPN values 
  # Initial values for EGSS are averages of EGOE and EGPN values 
  mu0 <- (mu.egoe+mu.egpn)/2    #  For ML only 
  x00 <- x0.egoe                #  For ML only     
  ssq0 <-  ssq.egpn/2 # For ML and REML
  tsq0 <-  tsq.egoe/2 # For ML and REML
  
  # Calculate ML & REML parameter estimates
  
  # The REML estimates.
  EGSSreml <-  optim(par =  c(ssq0, tsq0),
                   negloglike.reml, NULL, method = "Nelder-Mead", yt = Y.t, tt = T.t);
  params.reml <-  c(exp(EGSSreml$par[1]), exp(EGSSreml$par[2]))
  ssq.reml <-  params.reml[1]   
  #  These are the REML estimates
  tsq.reml <-  params.reml[2]   
  vx <-  matrix(0, qp1, qp1)
  for (ti in 1:q){
    vx[(ti + 1):qp1, (ti + 1):qp1] <-  matrix(1, 1, (qp1 - ti))*T.t[ti + 1]
    }
  Sigma.mat <-  ssq.reml*vx
  Itausq <-  matrix(rep(0,(qp1*qp1)), nrow = q + 1, ncol = q + 1)
  diag(Itausq) <-  rep(tsq.reml, q + 1)
  V <-  Sigma.mat + Itausq
  D1mat <-  cbind(-diag(1/S.t), matrix(0, q, 1)) + cbind(matrix(0, q, 1), diag(1/S.t))
  V1mat <-  D1mat%*%V%*%t(D1mat)
  W.t <-  (Y.t[2:qp1] - Y.t[1:q])/S.t
  j1 <-  matrix(1, q, 1)
  V1inv <-  ginv(V1mat)
  mu.reml <-  (t(j1)%*%V1inv%*%W.t)/(t(j1)%*%V1inv%*%j1)
  j <-  matrix(1, qp1, 1)
  Vinv <-  ginv(V)
  x0.reml <-  (t(j)%*%Vinv%*%(Y.t-as.vector(mu.reml)*T.t))/(t(j)%*%Vinv%*%j)
  Var_mu.reml <-  1/(t(j1)%*%V1inv%*%j1) # Variance of mu
  mu_hi.reml <-  mu.reml + 1.96*sqrt(Var_mu.reml) # 95% CI for mu
  mu_lo.reml <-  mu.reml - 1.96*sqrt(Var_mu.reml)
  
  #  Calculate estimated population sizes for EGSS model with Kalman filter
  #  Choose REML estimates here for calculating model values
  mu <-  mu.reml  
  ssq <-  ssq.reml  
  tsq <-  tsq.reml
  x0 <-  x0.reml
  m[1] <-  x0       
  #  Initial mean of Y(t)
  v[1] <-  tsq 
  #  Initial variance of Y(t)
  for (ti in 1:q){
    # Loop to generate estimated population abundances
    # using Kalman filter (see equations 6 & 7, Dennis et al. (2006)).
    m[ti + 1] <-  mu + (m[ti] + ((v[ti] - tsq)/v[ti])*(Y.t[ti] - m[ti]))
    v[ti + 1] <-  tsq*((v[ti] - tsq)/v[ti]) + ssq + tsq
  }
  Predict.t <-  exp(m + ((v - tsq)/v)*(Y.t - m))
  pop[[i]]$Pred.N <- Predict.t
  pop[[i]]$var <- v
  
  #  The following statement calculates exp{E[X(t) | Y(t), Y(t-1),...,Y(0)]};
  #  Print the parameter estimates
  parms.egoe <-  c(mu.egoe, ssq.egoe, tsq.egoe, x0.egoe) #  Collect for printing
  parms.egpn <-  c(mu.egpn, ssq.egpn, tsq.egpn, x0.egpn)
  parms.reml <-  c(mu.reml, ssq.reml, tsq.reml, x0.reml) 
  names <-  c("mu", "ssq", "tsq", "x0")             
  types <-  c("EGOE","EGPN","EGSS-ML","EGSS-REML")
  
  # Add to dataframe
  pops_data[i,18] <- parms.reml[1]
  pops_data[i,19] <- mu_lo.reml[1,1]
  pops_data[i,20] <- mu_hi.reml[1,1]
  pops_data[i,21] <- Var_mu.reml[1,1]
  pops_data[i,22] <- parms.reml[2]
  pops_data[i,23] <- parms.reml[3]
  pops_data[i,24] <- parms.reml[4]
  pops_data[i,25] <- min(Time.t)
  pops_data[i,26] <- max(Time.t)
  
}

colnames(pops_data)[18:26]<- c("mu", "lCI", "uCI", "var", "sigma",
                               "tau", "x0", "Start", "End")


# Create data frame over time for each population

pop_pred <- do.call("rbind", pop)

# Estimate one-step-ahead residuals

pop_pred <- pop_pred %>% 
  group_by(pops.ID.i.) %>% 
  mutate(residuals=N-lag(log(Pred.N)),
         std_residuals= residuals/sd(residuals,na.rm = T)) 

# Create unique ids for species and year

dd_final$ID2 <- paste(dd_final$ID, dd_final$Year)
pop_pred$ID2 <- paste(pop_pred$pops.ID.i., pop_pred$Year)

#bind both data frames

pop_data <- left_join(dd_final, pop_pred %>% 
                        dplyr::select(ID2,Pred.N,
                                     residuals,std_residuals),
                  by = "ID2")

# General checks ---------------------------------------------------------------

# Normality

norm_pop <- pop_data %>% 
  drop_na(std_residuals) %>% 
  group_by(ID) %>% 
  summarise(norm_test=shapiro.test(std_residuals)[[2]], 
            resid_max=max(std_residuals)) %>% 
  filter(norm_test>=0.01, 
         resid_max<10) 

# Plot it 

(g1<- pop_data %>% 
    filter(ID%in%norm_pop$ID) %>% 
    ggplot() +
    geom_line(aes(x = std_residuals, 
                  y = ..scaled.., 
                  group = as.character(pops.ID.i.)), 
              stat="density", size=0.1, alpha=0.05, adjust = 5) +
    coord_cartesian(xlim = c(-5,5)) +
    labs(x = "Standardized residuals", y = "Density"))

# Autocorrelation 

autoc <- pop_data %>% 
  drop_na(std_residuals) %>% 
  group_by(ID) %>% 
  summarise(auto_f=acf(std_residuals, plot = F)[[1]], 
            lag = acf(std_residuals, plot = F)[[4]]) %>% 
  filter(lag!=0)

# Plot it 

(g2 <- autoc %>% 
    ggplot() +
    geom_point(aes(x = as.factor(lag), y = auto_f), 
               position = "jitter", alpha = 0.03) +
    geom_boxplot(aes(x = as.factor(lag), y = auto_f), 
                 alpha = 0.5, 
                 outlier.shape = NA) +
    geom_hline(aes(yintercept = 0.5), colour = "blue", linetype = "dashed") +
    geom_hline(aes(yintercept = -0.5), colour = "blue", linetype = "dashed") +
    coord_cartesian(xlim = c(0.5,19), ylim = c(-1,1)) +
    scale_x_discrete(expand = c(0,0), breaks = c(2,4,6,8,10,12,14,16,18)) +
    labs(x = "Lag", y= "ACF"))


# Plot the diagnostics 

(figS2 <- g1+g2+plot_annotation(tag_levels = "a"))

ggsave("FigureS2.pdf", figS2, 
       width = 10, height = 6, 
       path = ResultPath)

# Keep only the good models

pops_data <- pops_data %>% 
  filter(ID%in%norm_pop$ID)

pop_data <- pop_data %>% 
  filter(ID%in%norm_pop$ID)

# Save the data ----
setwd(DataPath)
save(pops_data, pop_data, file = "SSResults.RData")

