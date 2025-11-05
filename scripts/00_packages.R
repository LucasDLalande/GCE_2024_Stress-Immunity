### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## Early-life glucocorticoids accelerate the senescence rate of lymphocyte count in roe deer
##
## Package loading
## 
## Lucas Lalande et al. 2024
## May 2024
##
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
rm(list=ls())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#### 1. Load packages  #### 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Add package used in the R script
packages <- c("interplot", "lemon", "texreg", "tiff", "lme4", "wesanderson", "MuMIn", "ggplot2", "sjmisc", "ggeffects", "dplyr", "lubridate", "ggpubr", "gridExtra", "egg", "grid", "OnAge", "cowplot")

# Package function

require <- function(x) { 
  if (!base::require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE) ; 
    base::require(x, character.only = TRUE)
  } 
}

# Load package

base::lapply(packages, require)

