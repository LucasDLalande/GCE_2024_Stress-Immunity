### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## Early-life glucocorticoids accelerate the senescence rate of lymphocyte count in roe deer
##
## Juvenile FGM correction
## 
## Lucas Lalande et al. 2024
## May 2024
##
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
rm(list=ls())

# DATASET PREPARATION ----
{
  # Data cleaning and keeping only individuals aged 2 or more years of age (see Cheynel et al., 2017 Scientific Reports)
  {
    data <- read.csv("dataset/dataset_anonym/00_dataset-immunosenescence-FGM_anonym.csv", h=T, sep=";", dec=".")
    data <- data[!is.na(data$ageannee),]
    data <- data[!is.na(data$masse),]
    data$pop <- as.factor(data$pop)
    data$sexe <- factor(data$sexe)
    data$numind <- as.factor(data$numind)
    data$sexe <- as.factor(data$sexe)
    data$idkit <- as.factor(data$idkit)
    data$id_bioch <- as.factor(data$id_bioch)
    data$annee <- as.factor(data$annee)
    data$date <- dmy(data$date)
    
    # FGM of Priority = 3 --> NA
    # FGMs with priority = 3 are samples that had been poorly stored so are no reliable samples
    data <- data %>%
      mutate((replace(FGMngg, FGM_Priority == "3", NA)))
    data <- data[,-41]
    names(data)[41] <- "FGMngg"
  }
  
  # Limiting the dataset by individuals having both a juvenile measure of FGM and an adult immune value
  {
    data2 <- subset(data, ageannee >= 2) # 1277 obs / 572 ind
    unique(factor(data2$numind))
    
    data1 <- subset(data, ageannee == 1)
    data1 <- data1[!is.na(data1$FGMngg),]
    unique(factor(data1$numind)) # 450 obs / 449 ind
    data1[duplicated(data1$numind),]
    data1 <- data1[!data1$idkit=="2014-79-031",] # removing an obs for an individual with two captures the same year
    unique(factor(data1$numind)) # 449 obs / 449 ind
    
    data1 <- data1[,c(3,5,7,8,41)]
    colnames(data1)[5] <- "fgm1" # FGM level during the 1st year
    
    data <- merge(data1, data2,  by="numind", all=F)
    unique(factor(data$numind))
    colnames(data)[c(2,3,4,8,10,12)] <- c("annee_juv", "datejulienne_juv", "dureemin_juv", "annee_adult", "datejulienne_adult", "dureemin_adult")
    
    data$age_diff <- data$ageannee-1
    
    data$age_factor_full <- as.factor(data$ageannee)
    
    data$cohorte <- factor(data$cohorte)
    data$annee_juv <- factor(data$annee_juv)
    data$annee_adult <- factor(data$annee_adult)
  }
}

# FGM CORRECTION ----
# We test whether FGM values are impacted by the julian date of capture (linear and quadratic relationships), the delay between capture and sampling
# and the time between sampling and freezing (immediately upon collection or within 24 hours)
# For the manuscript, we will only use early-life FGM, so we will only correct those data
data$annee_juv <- as.numeric(data$annee_juv)+2009
data$freezing <- ifelse(data$pop=="TF" & data$annee_juv<=2016, "<24h", "immediate")
data$freezing <- as.factor(data$freezing)
data$annee_juv <- as.factor(data$annee_juv)
  
data_unique <- data[!duplicated(data$numind),]
unique(factor(data_unique$numind)) # 162 ind/162 obs
  
mod.fgm <- lm(log(fgm1) ~ datejulienne_juv + I(datejulienne_juv^2) + dureemin_juv + freezing, 
              data=data_unique)
summary(mod.fgm)
  
options(na.action="na.fail")
mod.fgm.dredge <- dredge(mod.fgm, evaluate=T, rank="AICc", trace=2); mod.fgm.dredge
# The null model is kept

