### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## Early-life glucocorticoids accelerate the senescence rate of lymphocyte count in roe deer
##
## Test for FGM effect on immunosenescence patterns
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
  
  # Age of last observation/capture for Selective disappearance
  {
    # last year obs / last year capt / age
    # Calculate the age of last capture or observation
    # lastyearcapt - cohort = age at last capture
    # Same with lastobs - cohort = age at last obs
    # Keep the most recent of both, or the one available when there is a NA
    
    # In case there is multiple value for a same individual, we keep only the max
    # /!\ issue when NA for an individual --> max does not account for NAs ...
    data$lastyearcapt[is.na(data$lastyearcapt)] <- 0 # replace NA by O
    data$lastcaptmax <- ave(data$lastyearcapt, data$numind, FUN=function(x) max(x, na.rm=T)) # max last year capt for each indiv
    data$lastcaptmax[data$lastcaptmax==0] <- NA # replace back 0 by NA
    
    data$lastobs[is.na(data$lastobs)] <- 0 # replace NA by O
    data$lastobsmax <- ave(data$lastobs, data$numind, FUN=function(x) max(x))
    data$lastobsmax[data$lastobsmax==0] <- NA # replace back 0 by NA
    
    # In a new column we keep the most recent value (obs or capture), and if we have only one of the two we keep that one
    data$lastcaptobs <- with(data, 
                             ifelse(is.na(lastcaptmax), lastobsmax, 
                                    ifelse(is.na(lastobsmax), lastcaptmax, 
                                           ifelse(lastcaptmax > lastobsmax, lastcaptmax, lastobsmax))))
    
    # Age at last capture/observation
    data$agelast <- data$lastcaptobs - data$cohorte # agelast will be the variable used to test for selective disappearance
    data$cohorte <- as.factor(data$cohorte)
    data <- data[!is.na(data$agelast),]
  }
  
  # Include residuals of the HL ~ serum color linear regression before limiting the dataset by juvenile FGMs
  {
    dataHL <- data[!is.na(data$HL),] 
    dataHL <- dataHL[!is.na(dataHL$couleur_HAHL),] 
    dataHL$resHL <- residuals(lm(dataHL$HL~dataHL$couleur_HAHL))
    dataHL <- dataHL[,c(1,47)]
    data <- merge(data, dataHL, by="idkit", all=T)
  }

  # Limiting the dataset by individuals having both a juvenile measure of FGM and an adult immune value
  {
    data2 <- subset(data, ageannee >= 2) # 1277 obs / 572 ind
    unique(factor(data2$numind))
    
    data1 <- subset(data, ageannee == 1)
    data1 <- data1[!is.na(data1$FGMngg),]
    unique(factor(data1$numind)) # 450 obs / 449 ind
    data1[duplicated(data1$numind),]
    data1 <- data1[!data1$idkit=="IDkit788",] # removing an obs for an individual with two captures the same year
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

  # Separate cohorts between two cohorts class (good quality v. poor quality)
  {
    # Given than cohort quality do not seem normally distributed we will split cohort into high and low quality according to pop-specific median cohort qual
    # We saw in the previous analysis, on the larger dataset that medians were 13.61 kg in CH and 16.16 kg in TF
    # MEDIAN TF: 16.16
    # MEDIAN CH: 13.61
    
    dataCH <- data[data$pop=="CH",]
    dataTF <- data[data$pop=="TF",]
    dataCH$cohort_class <- ifelse(dataCH$qualite_cohorte >= 13.61, "good", "poor")
    table(dataCH$cohort_class)
    # Good quality: 94 obs
    # Poor quality: 73 obs
    dataTF$cohort_class <- ifelse(dataTF$qualite_cohorte >= 16.16, "good", "poor")
    table(dataTF$cohort_class)
    # Good quality: 100 obs
    # Poor quality: 85 obs
    
    data <- rbind(dataCH, dataTF)
    
    data$cohort_class <- factor(data$cohort_class)
  }
}

# I - IMMUNITY -----
# 1A./ NEUTROPHIL ----
dataneutro <- data[!is.na(data$neutro),]
dataneutro <- dataneutro[!is.na(dataneutro$masse),] 
unique(factor(dataneutro$numind)) # 267 obs/142 ind

# Retained model: QUADRATIC
# Age² + pop + delay
dataneutro$age2 <- I(dataneutro$ageannee^2)
mod.neutro <- lmer(neutro ~ age2*log(fgm1) + pop*log(fgm1) + dureemin_adult
                 + (1|numind) + (1|annee_adult:pop),
                 data=dataneutro, control=lmerControl("bobyqa"), REML=F)
summary(mod.neutro)

options(na.action="na.fail")
dneutro <- dredge(mod.neutro, rank="AICc", evaluate=T, trace=2, fixed=~age2 + pop + dureemin_adult); dneutro
# We retain the same model than previously, no effect of FGM

# 1B./ MONOCYTE ----
datamono <- data[!is.na(data$mono),]
datamono <- datamono[!is.na(datamono$masse),] 
unique(factor(datamono$numind)) # 267 obs/142 ind

# Retained model: Threshold (slope-constant) (5yo)
# Age.1
datamono$age.1 <- datamono$ageannee
datamono$age.1 <- ifelse(datamono$age.1<5, datamono$age.1-5, 0)

mod.mono <- lmer(mono ~ age.1*log(fgm1)
                 + (1|numind) + (1|annee_adult:pop), 
                 data=datamono, control=lmerControl("bobyqa"), REML=F)
summary(mod.mono)

options(na.action="na.fail")
dmono <- dredge(mod.mono, rank="AICc", evaluate=T, trace=2, fixed=~age.1); dmono
# We retain the same model than previously, no effect of FGM

# 1C./ BASOPHIL ----
databaso <- data[!is.na(data$baso),]
databaso <- databaso[!is.na(databaso$masse),] 
unique(factor(databaso$numind)) # 267 obs/142 ind

# Retained model: CONSTANT
# Dureemin (AICc = -2081.26, df = 5)
mod.baso <- lmer(baso ~ dureemin_adult + log(fgm1)
                   + (1|numind) + (1|annee_adult:pop), 
                   data=databaso, control=lmerControl("bobyqa"), REML=F)
summary(mod.baso)

dbaso <- dredge(mod.baso, rank="AICc", evaluate=T, trace=2, fixed=~dureemin_adult); dbaso
# We retain the same model than previously, no effect of FGM

# 1D./ EOSINOPHIL ----
dataeosino <- data[!is.na(data$eosino),]
dataeosino <- dataeosino[!is.na(dataeosino$masse),]
unique(factor(dataeosino$numind)) # 267 obs/142 ind

# Retained model: SEX
# Dureemin + Sex
mod.eosino <- lmer(eosino ~ dureemin_adult + sexe*log(fgm1)
                   + (1|numind) + (1|annee_adult:pop), 
                   data=dataeosino, control=lmerControl("bobyqa"), REML=F)
summary(mod.eosino)

deosino <- dredge(mod.eosino, rank="AICc", evaluate=T, trace=2, fixed=~dureemin_adult + sexe); deosino
# We retain the same model than previously, no effect of FGM

# 1E./ HEMAGLUTINATION ----
dataHA <- data[!is.na(data$HA),]
dataHA <- dataHA[!is.na(dataHA$masse),]
unique(factor(dataHA$numind)) # 332 obs/159 ind

# Retained model: NULL MODEL
# Constant
mod.HA <- lmer(HA ~ log(fgm1)
                     + (1|numind) + (1|annee_adult:pop), 
                   data=dataHA, control=lmerControl("bobyqa"), REML=F)
summary(mod.HA)

dHA <- dredge(mod.HA, rank="AICc", evaluate=T, trace=2); dHA
# We retain the same model than previously, no effect of FGM

# 1F./ HEMOLYSIS ----
dataHL <- data[!is.na(data$HL),]
dataHL <- dataHL[!is.na(dataHL$masse),]
unique(factor(dataHL$numind))  # 328 obs/159 ind

# Retained model: NULL
# colorHAHL + dureemin
mod.HL <- lmer(resHL ~ dureemin_adult + log(fgm1)
                    + (1|numind) + (1|annee_adult:pop), 
                    data=dataHL, control=lmerControl("bobyqa"), REML=F) # AICc = 3358.34, df = 5
summary(mod.HL)

dHL <- dredge(mod.HL, rank="AICc", evaluate=T, trace=2, fixed=~dureemin_adult); dHL
# We retain the same model than previously, no effect of FGM

# 2A./ ALPHA 1-GLOBULIN ----
dataalpha1 <- data[!is.na(data$alpha1),]
dataalpha1 <- dataalpha1[!is.na(dataalpha1$masse),] 
unique(factor(dataalpha1$numind)) # 249 obs/136 ind

# Retained model: QUADRATIC
# Age^2*sex + agelast + mass + datejulienne
dataalpha1$age2 <- I(dataalpha1$ageannee^2)
mod.alpha1 <- lmer(alpha1 ~ age2*sexe + masse + datejulienne_adult + agelast + log(fgm1)*age2 + log(fgm1)*sexe
                 + (1|numind) + (1|annee_adult:pop),
                 data=dataalpha1, control=lmerControl("bobyqa"), REML=F)
summary(mod.alpha1)

options(na.action="na.fail")
dredge(mod.alpha1, rank="AICc", evaluate=T, trace=2, fixed=~age2*sexe + masse + datejulienne_adult + agelast)
# We retain the same model than previously, no effect of FGM

# ID creates singularity
# Create a dataset with only one obs per ID
dataalpha1unique <- dataalpha1 %>% group_by(numind) %>% filter(n()== 1) %>% ungroup() # group individual with only 1 obs
dataalpha1repeti <- dataalpha1 %>% group_by(numind) %>% filter(n()>= 2) %>% ungroup() # group individual with >= 2 obs
unique(factor(dataalpha1repeti$numind)) # 59 ind

mean(dataalpha1unique$ageannee) # 2.84 yo on average for individuals with 1 obs
dataalpha1repeti$diff <- abs(dataalpha1repeti$ageannee - 2.84) # calculate difference between age at observation and the mean age of individuals with 1 obs
dataalpha1repet <-  dataalpha1repeti %>% 
  group_by(numind) %>% 
  slice(which.min(diff)) # group obs for which the age difference between obs and the mean age of unique obs is the lowest

dataalpha1repet <- dataalpha1repet[,-53] # remove "diff" column so dataalpha1unique and dataalpha1repet have the same number of column
dataalpha1 <- rbind(dataalpha1unique, dataalpha1repet)
unique(factor(dataalpha1$numind)) # 136 individuals

dataalpha1$age2 <- I(dataalpha1$ageannee^2)
mod.alpha1 <- lmer(alpha1 ~ age2*sexe + masse + datejulienne_adult + agelast + log(fgm1)*age2 + log(fgm1)*sexe
                   + (1|annee_adult:pop),
                   data=dataalpha1, control=lmerControl("bobyqa"), REML=F)
summary(mod.alpha1)

options(na.action="na.fail")
dalpha1 <- dredge(mod.alpha1, rank="AICc", evaluate=T, trace=2, fixed=~age2*sexe + masse + datejulienne_adult + agelast); dalpha1
# We retain the same model than previously, no effect of FGM

# 2B./ ALPHA 2-GLOBULIN ----
dataalpha2 <- data[!is.na(data$alpha2),]
dataalpha2 <- dataalpha2[!is.na(dataalpha2$masse),]
unique(factor(dataalpha2$numind)) # 249 obs/136 ind

# Retained model: CONSTANT
# Sex
mod.alpha2 <- lmer(alpha2 ~ sexe*log(fgm1)
                  + (1|numind) + (1|annee_adult:pop), 
                  data=dataalpha2, control=lmerControl("bobyqa"), REML=F)
summary(mod.alpha2)

options(na.action="na.fail")
dredge(mod.alpha2, rank="AICc", evaluate=T, trace=2, fixed=~sexe)
# We retain the same model than previously, no effect of FGM

# ID creates singularity
# Create a dataset with only one obs per ID
dataalpha2unique <- dataalpha2 %>% group_by(numind) %>% filter(n()== 1) %>% ungroup() # group individual with only 1 obs
dataalpha2repeti <- dataalpha2 %>% group_by(numind) %>% filter(n()>= 2) %>% ungroup() # group individual with >= 2 obs
unique(factor(dataalpha2repeti$numind)) # 59 ind

mean(dataalpha2unique$ageannee) # 2.84 yo on average for individuals with 1 obs
dataalpha2repeti$diff <- abs(dataalpha2repeti$ageannee - 2.84) # calculate difference between age at observation and the mean age of individuals with 1 obs
dataalpha2repet <-  dataalpha2repeti %>% 
  group_by(numind) %>% 
  slice(which.min(diff)) # group obs for which the age difference between obs and the mean age of unique obs is the lowest
dataalpha2repet <- dataalpha2repet[,-53] # remove "diff" column so dataalpha2unique and dataalpha2repet have the same number of column
dataalpha2 <- rbind(dataalpha2unique, dataalpha2repet)
unique(factor(dataalpha2$numind)) # 136 individuals

mod.alpha2 <- lmer(alpha2 ~ sexe*log(fgm1)
                   + (1|annee_adult:pop), 
                   data=dataalpha2, control=lmerControl("bobyqa"), REML=F)
summary(mod.alpha2)

options(na.action="na.fail")
dalpha2 <- dredge(mod.alpha2, rank="AICc", evaluate=T, trace=2, fixed=~sexe); dalpha2
# We retain the same model than previously

# 2C./ BETAGLOBULIN ----
databeta <- data[!is.na(data$beta),]
databeta <- databeta[!is.na(databeta$masse),]
unique(factor(databeta$numind)) # 249 obs/136 ind

# Retained model: LINEAR 
# Age + sex + mass
mod.beta <- lmer(beta ~ ageannee*log(fgm1) + sexe*log(fgm1) + masse
                   + (1|numind) + (1|annee_adult:pop), 
                   data=databeta, control=lmerControl("bobyqa"), REML=F)
summary(mod.beta)

options(na.action="na.fail")
dbeta <- dredge(mod.beta, rank="AICc", evaluate=T, trace=2, fixed=~ageannee + sexe + masse); dbeta
# We retain the same model than previously, no effect of FGM

# 2D./ HAPTOGLOBIN ----
dataHAP <- data[!is.na(data$HAP),]
dataHAP <- dataHAP[!is.na(dataHAP$masse),]
unique(factor(dataHAP$numind)) # 254 obs/138 ind

# Retained model: QUADRATIC
# Age² + sex
dataHAP$age2 <- I(dataHAP$ageannee^2)
mod.HAP <- lmer(HAP ~ age2*log(fgm1) + sexe*log(fgm1)
                 + (1|numind) + (1|annee_adult:pop), 
                 data=dataHAP, control=lmerControl("bobyqa"), REML=F)
summary(mod.HAP)

options(na.action="na.fail")
dHAP <- dredge(mod.HAP, rank="AICc", evaluate=T, trace=2, fixed=~age2 + sexe); dHAP
# We retain the same model than previously, no effect of FGM

# 3A./ GAMMAGLOBULIN ----
datagamma <- data[!is.na(data$gamma),]
datagamma <- datagamma[!is.na(datagamma$masse),]
unique(factor(datagamma$numind)) # 249 obs/136 ind

# Retained model: LINEAR
# Age + mass + pop
mod.gamma <- lmer(gamma ~ ageannee*log(fgm1) + pop*log(fgm1) + masse
                  + (1|numind) + (1|annee_adult:pop),
                  data=datagamma, control=lmerControl("bobyqa"), REML=F)
summary(mod.gamma)

options(na.action="na.fail")
dgamma <- dredge(mod.gamma, rank="AICc", evaluate=T, trace=2, fixed=~ageannee + pop + masse); dgamma
# We retain the same model than previously, no effect of FGM

# 3B./ LYMPHOCYTE ----
datalympho <- data[!is.na(data$lympho),]
datalympho <- datalympho[!is.na(datalympho$masse),]
unique(factor(datalympho$numind)) # 267 obs/142 ind

# Retained model : 2-SLOPES THRESHOLD (4yo)
# Age.1 + Age.2*pop + sex + Delay
{datalympho$age.1 <- datalympho$ageannee
  datalympho$age.2 <- datalympho$ageannee
  
  datalympho$age.1 <- datalympho$age.2 <- datalympho$ageannee
  datalympho$age.1 [datalympho$age.1 > 4] <- 4
  datalympho$age.2 [datalympho$age.2 <= 4] <- 0
  datalympho$age.2 [datalympho$age.2 > 4] <- datalympho$age.2 [datalympho$age.2 > 4] - 4
  mod.lympho <- lmer(lympho ~ age.1*log(fgm1) + age.2*pop + age.2*log(fgm1) + pop*log(fgm1) + sexe*log(fgm1) + dureemin_adult
                           + (1|numind) + (1|annee_adult:pop),
                           data=datalympho, control=lmerControl("bobyqa"), REML=F)}
summary(mod.lympho)

options(na.action="na.fail")
dlympho <- dredge(mod.lympho, rank="AICc", evaluate=T, trace=2, fixed=~age.1 + age.2*pop + sexe + dureemin_adult); dlympho
# Age.2*FGM retained

{datalympho$age.1 <- datalympho$ageannee
  datalympho$age.2 <- datalympho$ageannee
  datalympho$logfgm1 <- log(datalympho$fgm1)
  
  datalympho$age.1 <- datalympho$age.2 <- datalympho$ageannee
  datalympho$age.1 [datalympho$age.1 > 4] <- 4
  datalympho$age.2 [datalympho$age.2 <= 4] <- 0
  datalympho$age.2 [datalympho$age.2 > 4] <- datalympho$age.2 [datalympho$age.2 > 4] - 4
  mod.lymphoFGM <- lmer(lympho ~ age.1 + age.2*pop + age.2*logfgm1 + sexe + dureemin_adult
                     + (1|numind) + (1|annee_adult:pop),
                     data=datalympho, control=lmerControl("bobyqa"), REML=F)}
summary(mod.lymphoFGM)

# plot
dataFGM <- interplot(mod.lymphoFGM, var1="age.2", var2="logfgm1", plot=F)

plotLYMPHO <- ggplot(dataFGM, aes(x=logfgm1, y=coef)) +
  geom_line(aes(x=logfgm1, y=coef-0.45, colour="darkgreen"), size=1) +
  geom_ribbon(aes(ymin=lb-0.45, ymax=ub-0.45, fill="darkgreen"), alpha=0.2, linetype=0) +
  geom_line(aes(x=logfgm1, y=coef, colour="darkolivegreen3"), size=1) +
  geom_ribbon(aes(ymin=lb, ymax=ub, fill="darkolivegreen3"), alpha=0.2, linetype=0) +
  geom_hline(yintercept=0, linetype="dashed") +
  ylab(expression(Estimated~coefficient~(slope)~of~age~after~the~threshold)) +
  labs(title = "Lymphocyte", tag="A") +
  scale_color_manual(name="Population", values=c("darkolivegreen3", "darkgreen"), labels=c("CH", "TF")) +
  scale_fill_manual(name="Population", values=c("darkolivegreen3", "darkgreen"), labels=c("CH", "TF")) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        plot.tag = element_text(face="bold", size=18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position=c(0.15, 0.23),
        legend.background=element_blank(),
        strip.text.x = element_text(size=13, face="bold"),
        strip.text.y = element_text(size=13, face="bold"),
        strip.background = element_blank()) #change legend text font size
plotLYMPHO

ggsave(filename = "results/02_fgm_immunosenescence_FIGURES/3A_lymphoFGM.tiff", plotLYMPHO, width = 6, height = 5, dpi = 300, units = "in", device='tiff')

# II - PARASITISM ----
# 1./ GASTRO-INTESTINAL STRONGYLES ----
datastrongyles <- data[!is.na(data$SD),]
datastrongyles <- datastrongyles[!is.na(datastrongyles$masse),]
unique(factor(datastrongyles$numind)) # 274 obs/140 ind


# Retained model: 2-SLOPES THRESHOLD (10 years old)
# age.1 + age.2*sexe + agelast + masse
{datastrongyles$age.1 <- datastrongyles$ageannee
  datastrongyles$age.2 <- datastrongyles$ageannee
  
  datastrongyles$age.1 <- datastrongyles$age.2 <- datastrongyles$ageannee
  datastrongyles$age.1 [datastrongyles$age.1 > 10] <- 10
  datastrongyles$age.2 [datastrongyles$age.2 <= 10] <- 0
  datastrongyles$age.2 [datastrongyles$age.2 > 10] <- datastrongyles$age.2 [datastrongyles$age.2 > 10] - 10
  mod.strong <- lmer(log(SD+1) ~ age.1*log(fgm1) + age.2*sexe + age.2*log(fgm1) + sexe*log(fgm1) + agelast + masse
                           + (1|numind) + (1|annee_adult:pop),
                           data=datastrongyles, control=lmerControl("bobyqa"), REML=F)}
summary(mod.strong)

options(na.action="na.fail")
dstrong <- dredge(mod.strong, rank="AICc", evaluate=T, trace=2, fixed=~age.1 + age2*sexe + agelast + masse); dstrong
# No FGM effects

# 2./ TRICHURIS SP. ----
datatrich <- data[!is.na(data$Trich),]
datatrich <- datatrich[!is.na(datatrich$masse),]
unique(factor(datatrich$numind)) # 252 obs/129 ind

# Retained model: 2-SLOPES THRESHOLD (10 years old)
# Age.1*pop + age.1*sex + age.2*pop + age.2*sexe + mass
{datatrich$age.1 <- datatrich$ageannee
  datatrich$age.2 <- datatrich$ageannee
  
  datatrich$age.1 <- datatrich$age.2 <- datatrich$ageannee
  datatrich$age.1 [datatrich$age.1 > 10] <- 10
  datatrich$age.2 [datatrich$age.2 <= 10] <- 0
  datatrich$age.2 [datatrich$age.2 > 10] <- datatrich$age.2 [datatrich$age.2 > 10] - 10
  mod.trich <- lmer(log(Trich+1) ~ age.1*sexe + age.2*sexe + age.1*pop + age.2*pop + masse
                           + age.1*log(fgm1) + age.2*log(fgm1) + sexe*log(fgm1) + pop*log(fgm1)
                           + (1|numind) + (1|annee_adult:pop),
                           data=datatrich, control=lmerControl("bobyqa"), REML=F)}
summary(mod.trich)

options(na.action="na.fail")
dtrich <- dredge(mod.trich, rank="AICc", evaluate=T, trace=2, fixed=~age.1*sexe + age.2*sexe + age.1*pop + age.2*pop + masse); dtrich
# No FGM effects

# 3./ PROTOSTRONGYLIDS ----
dataproto <- data[!is.na(data$Proto),]
dataproto <- dataproto[!is.na(dataproto$masse),]
unique(factor(dataproto$numind)) # 258 obs/133 ind

# Retained model: QUADRATIC
# Age*cohort + Age²*cohort + + Age²*sex + agelast + masse + pop
mod.proto <- lmer(log(Proto+1) ~ ageannee*cohort_class + I(ageannee^2)*sexe + I(ageannee^2)*cohort_class + pop + masse + agelast
                 + ageannee*log(fgm1) + I(ageannee^2)*log(fgm1) + cohort_class*log(fgm1) + sexe*log(fgm1) + pop*log(fgm1)
                 + (1|numind) + (1|annee_adult:pop), 
                 data=dataproto, control=lmerControl("bobyqa"), REML=F)
summary(mod.proto)

options(na.action="na.fail")
dproto <- dredge(mod.proto, rank="AICc", evaluate=T, trace=2, fixed=~ageannee*cohort_class + I(ageannee^2)*sexe + I(ageannee^2)*cohort_class + pop + masse + agelast);dproto
# FGM*cohort_class retained

mod.protoFGM <- lmer(log(Proto+1) ~ ageannee*cohort_class + I(ageannee^2)*sexe + I(ageannee^2)*cohort_class + pop + masse + agelast
                  + log(fgm1)
                  + (1|numind) + (1|annee_adult:pop), 
                  data=dataproto, control=lmerControl("bobyqa"), REML=F)
summary(mod.protoFGM)

# Plot
dataproto$logProto <- log(dataproto$Proto+1)
dataproto$logfgm1 <- log(dataproto$fgm1)

dataproto$res <- residuals(lmer(logProto ~ ageannee*cohort_class + I(ageannee^2)*sexe + I(ageannee^2)*cohort_class + pop + masse + agelast
                      + (1|numind) + (1|annee_adult:pop), 
                      data=dataproto, control=lmerControl("bobyqa"), REML=F))

mod.protoFGM.res <- lmer(res ~ cohort_class*logfgm1 
                         + (1|numind) + (1|annee_adult:pop), 
                         data=dataproto, control=lmerControl("bobyqa"), REML=F)
summary(mod.protoFGM.res)

pred <- ggpredict(mod.protoFGM.res, terms = c("logfgm1 [all]", "cohort_class")) # to have intercept and slope for each sex
colnames(pred)[c(1,6)] <- c("logfgm1", "cohort_class")
pred$cohort_class <- ifelse(pred$cohort_class=="good", "Good", "Poor")
dataproto$cohort_class <- ifelse(dataproto$cohort_class=="good", "Good", "Poor")

plotPROTO <- ggplot(dataproto, aes(x=logfgm1, y=res, color=cohort_class, fill=cohort_class)) +
  geom_point(aes(x=logfgm1, y=res, color=cohort_class), size=2, alpha=0.5) +
  geom_line(data=pred, aes(x=logfgm1, y=predicted), linewidth=1) +
  geom_ribbon(data=pred, aes(x=logfgm1, y=predicted, ymin=conf.low, ymax=conf.high), alpha=0.2, linetype=0) +
  ylab("Residuals of protostrongylid abundance") +
  labs(title="Protostrongylids", tag="B") +
  scale_x_continuous(breaks = seq(5, 9, by = 1)) +
  scale_color_manual(name="Cohort quality", values=c("peachpuff3", "peachpuff4")) +
  scale_fill_manual(name="Cohort quality", values=c("peachpuff3", "peachpuff4")) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        plot.tag = element_text(face="bold", size=18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position=c(0.15, 0.83)) #change legend text font size
plotPROTO

ggsave(filename = "results/02_fgm_immunosenescence_FIGURES/3B_protoFGM.tiff", plotPROTO, width = 6, height = 5, dpi = 300, units = "in", device='tiff')

# Combine figures ----
figure3 <- arrangeGrob(plotLYMPHO, plotPROTO, ncol=2, nrow=1)
grid.newpage()
grid.draw(figure3)

figure3 <- annotate_figure(figure3, bottom = text_grob("Juvenile faecal glucocorticoid metabolites (FGMs, log-transformed)", size=16))
figure3
ggsave(filename = "results/02_fgm_immunosenescence_FIGURES/FIGURE3.tiff", figure3, width = 12, height = 6, dpi = 300, units = "in", device='tiff')

# 4./ COCCIDIA ----
datacoccidia <- data[!is.na(data$Coc),]
datacoccidia <- datacoccidia[!is.na(datacoccidia$masse),]
unique(factor(datacoccidia$numind)) # 272 obs/140 ind

# Retained model: CONSTANT
# NULL
datacoccidia$logCoc <- log(datacoccidia$Coc+1)
mod.coc <- lmer(logCoc ~ log(fgm1)
                     + (1|numind) + (1|annee_adult:pop), 
                   data=datacoccidia, control=lmerControl("bobyqa"), REML=F)
summary(mod.coc)

options(na.action="na.fail")
dcoc <- dredge(mod.coc, rank="AICc", evaluate=T, trace=2); dcoc
# No FGM effects