# Multiple threats drive the decline of global vertebrate populations

Pol Capdevila<sup>1</sup>*, Louise McRae<sup>2</sup>, Valentina Marconi<sup>2</sup>, Thomas F. Johnson<sup>3</sup>, Robin Freeman<sup>2</sup>, Christopher F Clements<sup>1</sup>

<sup>1</sup>School of Biological Sciences, University of Bristol, 24 Tyndall Ave, BS8 1TQ, Bristol, UK. 

<sup>2</sup>Institute of Zoology, Zoological Society of London, Regent’s Park, London NW1 4RY, UK.

#### Contact: pcapdevila.pc[at]gmail.com

---

## Abstract

_Anthropogenic threats are reshaping Earth’s biodiversity at an unprecedented pace and scale. Conservation policies aiming to halt the loss of biodiversity often rely on global rankings of threats based on their prevalence, where habitat loss and exploitation are often considered the most important. However, these global assessments rarely quantify the impacts of single or multiple threats, which could mask the true effects of the Anthropocene. Here, we quantify the impacts of threats by analysing 3,415 vertebrate populations worldwide, where information about the threats affecting these populations is known. We show that, despite being the most prevalent threats, habitat loss and exploitation are not causing the most rapid population declines. Instead, less prevalent threats, such as invasive species, diseases, and climate change have the strongest negative impacts. However, habitat loss and exploitation have the strongest impacts when acting in conjunction with other threats, even if their cumulative effects are mostly additive (equal to the sum of their individual effects). Finally, we show that mitigating the effects of climate change is as important as mitigating exploitation and/or habitat loss in slowing population declines. These results imply that the most prevalent threats from previous reports are not necessarily the most impactful at local scales and that consequences of climate change on vertebrate populations may have been underestimated._

---

## Data

- __`TenYearsData.RData`__: contains the population trends of 3,145 vertebrate time-series from the state-space models. 

---

# Code

To run the statistical analyses we used different R scripts: 

- __`TenAnalyses.R`__: code to analyse the factors influencing resistance and recovery loss. These analyses require a lot of computing power.
- __`MultiEffects.R`__: code to calculate the different proportions of the cummulative effects. 
- __`SensiMultiEffects.R`__: code to calculate the sensitivity of the cummulative effects analyses to the null model choice. 
- __`Counterfactual tests.R`__: code to simulate the counterfactual scenarios. 
- __`Figures.R`__: code to create the figures 1-3 and tables S1-S4 of the study. 
- __`SupFigures.R`__: code to create the supplementary figures. 

---

# Software

_R version 4.0.2 or greater_

To download `R`, go to https://www.r-project.org and for `RStudio`, go to https://www.rstudio.com/products/rstudio/download/ .
