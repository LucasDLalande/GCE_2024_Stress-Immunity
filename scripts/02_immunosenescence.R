### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## Early-life glucocorticoids accelerate the senescence rate of lymphocyte counts in wild roe deer
##
## Immunosenescence patterns determination
## 
## Lucas Lalande et al. 2024
## May 2024
##
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
rm(list=ls())


source("scripts/immunosenescence_function.R")
# To use the functions cohort_trajectory() and cohort_trajectory_parasitism()
# that automatically run constant, linear, quadratic and the 3 threshold trajectories

# DATASET PREPARATION -----
{
  # Data cleaning and keeping only individuals aged 2 or more years of age (see Cheynel et al., 2017 Scientific Reports)
  {
    data <- read.csv("dataset/dataset_anonym/00_dataset-immunosenescence-FGM_anonym.csv", h=T, sep=";", dec=".")
    data <- data[!is.na(data$ageannee),]
    data <- data[!is.na(data$masse),]
    data$pop <- as.factor(data$pop)
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
    
    # Keeping only individuals aged 2 years of age or more
    data <- data[data$ageannee>=2,]
    data$age_factor_full <- as.factor(data$ageannee)
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
  
  
  # Separate cohorts between two cohorts class (good quality v. poor quality)
  {
    # Given than cohort quality do not seem normally distributed we will split cohort into high and low quality according to pop-specific median cohort qual
    cohort <- data[!duplicated(data[c(2,11,12)]),]
    cohort <- cohort[,c(2,11,12)]
    
    ggplot(cohort, aes(x=qualite_cohorte, fill=pop)) + 
      geom_density(size=1, alpha=0.4) +
      scale_fill_manual(values=wes_palette(n=2, name="Moonrise2"))
    
    hist(cohort$qualite_cohorte[cohort$pop=="CH"])
    hist(cohort$qualite_cohorte[cohort$pop=="TF"])
    
    tapply(cohort$qualite_cohorte, cohort$pop, median)
    # MEDIAN TF: 16.16
    # MEDIAN CH: 13.61
    
    dataCH <- data[data$pop=="CH",]
    dataTF <- data[data$pop=="TF",]
    dataCH$cohort_class <- ifelse(dataCH$qualite_cohorte >= 13.61, "good", "poor")
    table(dataCH$cohort_class)
    # Good quality: 287 obs
    # Poor quality: 362 obs
    dataTF$cohort_class <- ifelse(dataTF$qualite_cohorte >= 16.16, "good", "poor")
    table(dataTF$cohort_class)
    # Good quality: 328 obs
    # Poor quality: 300 obs
    
    data <- rbind(dataCH, dataTF)
    
  
    data$cohort_class <- factor(data$cohort_class)
  }
  
  table(data$ageannee, data$pop, data$sexe, data$cohort_class)
  # 10 years of age is the last age for which we have at least 1 obs after that age for each pop-sex-cohort group
}

# I - IMMUNITY ----
# 1A./ NEUTROPHIL ----
dataneutro <- data[!is.na(data$neutro),]
dataneutro <- dataneutro[!is.na(dataneutro$masse),] 
unique(factor(dataneutro$numind)) # 1021 obs/484 ind

immuno_trajectory(dataneutro, dataneutro$neutro)

# Retained model: QUADRATIC
# Age² + pop + delay (AICc = 4293.620, df = 13)
mod.quad <- lmer(neutro ~ I(ageannee^2) + pop + dureemin
                         + (1|numind) + (1|annee:pop),
                         data=dataneutro, control=lmerControl("bobyqa"), REML=F)
AICc(mod.quad)
summary(mod.quad)
r.squaredGLMM(mod.quad)

# Plot
dataneutro$age2 <- I(dataneutro$ageannee^2)
mod.quad <- lmer(neutro ~ age2 + pop
                 + (1|numind) + (1|annee:pop),
                 data=dataneutro, control=lmerControl("bobyqa"), REML=F)
# same model without confounding variables for visual representation

dataneutro$pop <- as.factor(dataneutro$pop)
dataneutro$age_factor_full <- as.factor(dataneutro$age_factor_full)
dataneutro$mean <- ave(dataneutro$neutro, dataneutro$age_factor_full:dataneutro$pop)
dataneutro$se <- ave(dataneutro$neutro, dataneutro$age_factor_full:dataneutro$pop, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(dataneutro[,c("pop", "ageannee", "age_factor_full", "mean", "se")])

pred <- ggpredict(mod.quad, terms = c("age2", "pop")) # to have intercept and slope for each pop
pred$x <- sqrt(pred$x)
colnames(pred)[c(1,6)] <- c("ageannee","pop")
pred

dataerror <- merge(dataerror, pred, by=c("ageannee", "pop"))

plotNEUTRO <- ggplot(dataerror, aes(x=ageannee, y=mean, color=pop, fill=pop)) +
  geom_point(aes(x=ageannee, y=mean, colour=pop), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=pop), width=.3, size=1.5) +
  geom_line(aes(x=ageannee, y=predicted), size=1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2, linetype=0) +
  ylab(expression(Neutrophil~count~(10^3~cells/mL))) +
  xlab(" ") +
  labs(title = "Neutrophil", tag="A") +
  scale_color_manual(name="Population", values=c("darkolivegreen3", "darkgreen")) +
  scale_fill_manual(name="Population", values=c("darkolivegreen3", "darkgreen")) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        axis.title.x = element_text(size=14),
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
plotNEUTRO

ggsave(filename = "results/01_immunosenescence_FIGURES/1A_NEUTRO.tiff", plotNEUTRO, width = 5, height = 5, dpi = 300, units = "in", device='tiff')

# 1B./ MONOCYTE ----
datamono <- data[!is.na(data$mono),]
datamono <- datamono[!is.na(datamono$masse),] 
unique(factor(datamono$numind)) # 1021 obs/484 ind

immuno_trajectory(datamono, datamono$mono)

# Retained model: Threshold (slope-constant) (5yo)
# Age.1 (AICc = 151.66, df = 6)
datamono$age.1 <- datamono$ageannee
datamono$age.1 <- ifelse(datamono$age.1<5, datamono$age.1-5, 0)

mod.quad <- lmer(mono ~ age.1
                 + (1|numind) + (1|annee:pop), 
                 data=datamono, control=lmerControl("bobyqa"), REML=F)

AICc(mod.quad)
summary(mod.quad)
r.squaredGLMM(mod.quad)

# Plot 
datamono$age_factor_full <- as.factor(datamono$age_factor_full)
datamono$mean <- ave(datamono$mono, datamono$age_factor_full)
datamono$se <- ave(datamono$mono, datamono$age_factor_full, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(datamono[,c("ageannee", "age_factor_full", "mean", "se")])

pred <- ggpredict(mod.quad, terms = c("age.1 [all]")) # to have intercept and slope
colnames(pred)[1] <- "ageannee"
pred$ageannee <- pred$ageannee+5
predint <- data.frame("ageannee"=seq(5,16,1), "predicted"=rep(0.2527781,12), "std.error"=rep(0.04152051,12),
                      "conf.low"=rep(0.1713024,12), "conf.high"=rep(0.3342539,12), "group"=rep(1,12))

pred <- rbind(pred, predint)

plotMONO <- ggplot(data=dataerror, aes(x=ageannee, y=mean)) +
  geom_line(data=pred, aes(x=ageannee, y=predicted), color="black", linewidth=1) +
  geom_ribbon(data=pred, aes(x=ageannee, y=predicted, ymax=conf.high, ymin=conf.low), alpha=0.2) +
  geom_point(data=dataerror, aes(x=ageannee, y=mean), size=2.5) +
  geom_errorbar(data=dataerror, aes(ymin=mean-se, ymax=mean+se), width=.3, linewidth=1.5) +
  geom_vline(xintercept=5, linetype='dashed', linewidth=0.75) +
  ylab(expression(Monocyte~count~(10^3~cells/mL))) +
  xlab(" ") +
  labs(title = "Monocyte", tag="B") +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        plot.tag = element_text(face = "bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position=c(0.15, 0.83)) #change legend text font size
plotMONO

ggsave(filename = "results/01_immunosenescence_FIGURES/1B_MONO.tiff", plotMONO, width = 5, height = 5, dpi = 300, units = "in", device='tiff')


# 1C./ BASOPHIL ----
databaso <- data[!is.na(data$baso),]
databaso <- databaso[!is.na(databaso$masse),] 
unique(factor(databaso$numind)) # 1021 obs/484 ind

immuno_trajectory(databaso, databaso$baso)

# Retained model: CONSTANT
# Dureemin (AICc = -2081.26, df = 5)
mod.linear <- lmer(baso ~ dureemin
                   + (1|numind) + (1|annee:pop), 
                   data=databaso, control=lmerControl("bobyqa"), REML=F)
summary(mod.linear)
AICc(mod.linear)
r.squaredGLMM(mod.linear)

# Plot
mod.linear <- lmer(baso ~
                     + (1|numind) + (1|annee:pop), 
                   data=databaso, control=lmerControl("bobyqa"), REML=F)
# same model without confounding variables for visual representation

databaso$age_factor_full <- as.factor(databaso$age_factor_full)
databaso$mean <- ave(databaso$baso, databaso$age_factor_full)
databaso$se <- ave(databaso$baso, databaso$age_factor_full, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(databaso[,c("ageannee", "age_factor_full", "mean", "se")])

plotBASO <- ggplot(dataerror, aes(x=ageannee, y=mean)) +
  geom_smooth(aes(x=ageannee, y=0.059598), color="black") +
  geom_ribbon(aes(x=ageannee, y=mean, ymax=0.059598+1.96*0.008912, ymin=0.059598-1.96*0.008912), alpha=0.2, linetype=0) +
  geom_point(aes(x=ageannee, y=mean), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, size=1.5) +
  ylab(expression(Basophil~count~(10^3~cells/mL))) +
  xlab(" ") +
  labs(title="Basophil", tag="A") +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        plot.tag = element_text(face = "bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plotBASO

ggsave(filename = "results/01_immunosenescence_FIGURES/S1A_BASO.tiff", plotBASO, width = 5, height = 5, dpi = 300, units = "in", device='tiff')

# 1D./ EOSINOPHIL ----
dataeosino <- data[!is.na(data$eosino),]
dataeosino <- dataeosino[!is.na(dataeosino$masse),] 
unique(factor(dataeosino$numind)) # 1021 obs/484 ind


immuno_trajectory(dataeosino, dataeosino$eosino)

# Retained model: SEX
# Dureemin + Sex (AICc = -909.68, df = 6)
mod.linear <- lmer(eosino ~ dureemin + sexe 
                   + (1|numind) + (1|annee:pop), 
                   data=dataeosino, control=lmerControl("bobyqa"), REML=F)
summary(mod.linear)
AICc(mod.linear, mod.quad)
r.squaredGLMM(mod.linear)

# Plot
mod.linear <- lmer(eosino ~ sexe 
                   + (1|numind) + (1|annee:pop), 
                   data=dataeosino, control=lmerControl("bobyqa"), REML=F)
# same model without confounding variables for visual representation

dataeosino$age_factor_full <- as.factor(dataeosino$age_factor_full)
dataeosino$sexe <- as.factor(dataeosino$sexe)
dataeosino$mean <- ave(dataeosino$eosino, dataeosino$age_factor_full:dataeosino$sexe)
dataeosino$se <- ave(dataeosino$eosino, dataeosino$age_factor_full:dataeosino$sexe, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(dataeosino[,c("ageannee", "age_factor_full", "sexe", "mean", "se")])
dataerror$predicted <- ifelse(dataerror$sexe=="M", 0.10285, 0.14455)
dataerror$conf.high <- ifelse(dataerror$sexe=="M", 0.10285+1.96*0.01077, 0.14455+1.96*0.01077)
dataerror$conf.low <- ifelse(dataerror$sexe=="M", 0.10285-1.96*0.01077, 0.14455-1.96*0.01077)

plotEOSINO <- ggplot(dataerror, aes(x=ageannee, y=mean, color=sexe, fill=sexe)) +
  geom_line(aes(x=ageannee, y=predicted), size=1) +
  geom_ribbon(aes(x=ageannee, y=mean, ymax=conf.high, ymin=conf.low), alpha=0.2, linetype=0) +
  geom_point(aes(x=ageannee, y=mean, color=sexe), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, color=sexe), width=.3, size=1.5) +
  scale_color_manual("Sex", values=c("darkorchid4", "darkorchid1")) +
  scale_fill_manual("Sex", values=c("darkorchid4", "darkorchid1")) +
  ylab(expression(Eosinophil~count~(10^3~cells/mL))) +
  xlab(" ") +
  labs(title="Eosinophil", tag="B") +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        plot.tag = element_text(face = "bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position=c(0.15, 0.87))
plotEOSINO

ggsave(filename = "results/01_immunosenescence_FIGURES/S1B_EOSINO.tiff", plotEOSINO, width = 5, height = 5, dpi = 300, units = "in", device='tiff')

# 1E./ HAEMAGLUTINATION ----
dataHA <- data[!is.na(data$HA),]
dataHA <- dataHA[!is.na(dataHA$masse),] 
dataHA <- dataHA[!is.na(dataHA$dureemin),]
unique(factor(dataHA$numind)) # 1201 obs/556 ind

immuno_trajectory(dataHA, dataHA$HA)

# Retained model: NULL MODEL
# Constant (AICc = 3883.18, df = 4)
mod.linear <- lmer(HA ~ 
                     + (1|numind) + (1|annee:pop), 
                   data=dataHA, control=lmerControl("bobyqa"), REML=F)
AICc(mod.linear)
summary(mod.linear)
r.squaredGLMM(mod.linear)

# Plot
dataHA$age_factor_full <- as.factor(dataHA$age_factor_full)
dataHA$mean <- ave(dataHA$HA, dataHA$age_factor_full)
dataHA$se <- ave(dataHA$HA, dataHA$age_factor_full, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(dataHA[,c("ageannee", "age_factor_full", "mean", "se")])

plotHA <- ggplot(dataerror, aes(x=ageannee, y=mean)) +
  geom_smooth(aes(x=ageannee, y=3.8953), color="black") +
  geom_ribbon(aes(x=ageannee, y=mean, ymax=3.8953+1.96*0.1453, ymin=3.8953-1.96*0.1453), alpha=0.2, linetype=0) +
  geom_point(aes(x=ageannee, y=mean), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, size=1.5) +
  ylab(expression(Hemagglutination~score~(titer))) +
  xlab("Age (years)") +
  labs(title="Haemagglutination", tag="C") +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        plot.tag = element_text(face="bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plotHA

ggsave(filename = "results/01_immunosenescence_FIGURES/S1C_HA.tiff", plotHA, width = 5, height = 5, dpi = 300, units = "in", device='tiff')

# 1F./ HEMOLYSIS ----
dataHL <- data[!is.na(data$HL),]
dataHL <- dataHL[!is.na(dataHL$masse),] 
dataHL <- dataHL[!is.na(dataHL$dureemin),] 
unique(factor(dataHL$numind)) # 1191 obs/555 ind

summary(lm(dataHL$HL~dataHL$couleur_HAHL))
plot(dataHL$HL~dataHL$couleur_HAHL)
abline(lm(dataHL$HL~dataHL$couleur_HAHL), col="red", lwd=2)

dataHL$resHL <- residuals(lm(dataHL$HL~dataHL$couleur_HAHL)) # to obtain residuals of the linear relationship between haemolysis score and serum coloration

# We'll prefer not to transform the response variable and thus use this model including serum coloration as fixed effect
immuno_trajectory_HL(dataHL, dataHL$HL)

# Retained model: NULL
# colorHAHL + dureemin (AICc = 3353.09, df = 6)
mod.linearA <- lmer(HL ~ couleur_HAHL + dureemin
                    + (1|numind) + (1|annee:pop), 
                    data=dataHL, control=lmerControl("bobyqa"), REML=F)
summary(mod.linearA)

# However we will plot the residuals of the HL ~ color regression on the y axis to account for both effects
mod.linearB <- lmer(resHL ~ dureemin
                    + (1|numind) + (1|annee:pop), 
                    data=dataHL, control=lmerControl("bobyqa"), REML=F) # AICc = 3358.34, df = 5
summary(mod.linearB)
r.squaredGLMM(mod.linearB)

# Plot
mod.linearCnocorr <- lmer(resHL ~
                            + (1|numind) + (1|annee:pop), 
                          data=dataHL, control=lmerControl("bobyqa"), REML=F)
# same model without confounding variables for visual representation

dataHL$age_factor_full <- as.factor(dataHL$age_factor_full)
dataHL$mean <- ave(dataHL$resHL, dataHL$age_factor_full)
dataHL$se <- ave(dataHL$resHL, dataHL$age_factor_full, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(dataHL[,c("ageannee", "age_factor_full", "mean", "se")])

text <- grobTree(textGrob("Residuals of haemolysis score according to serum color", x=0.5, y=0.985, hjust=0.5,
                          gp=gpar(fontsize=10.85, fontface="bold")))
plotHL <- ggplot(dataerror, aes(x=ageannee, y=mean)) +
  geom_smooth(aes(x=ageannee, y=0.0463), color="black") +
  geom_ribbon(aes(x=ageannee, y=mean, ymax=0.0463+1.96*0.1617, ymin=0.0463-1.96*0.1617), alpha=0.2, linetype=0) +
  geom_point(aes(x=ageannee, y=mean), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, size=1.5) +
  ylab(expression(Residuals~of~Haemolysis~score)) +
  xlab("Age (years)") +
  annotation_custom(text) +
  labs(title="Haemolysis", tag="D") +
  theme(plot.title = element_text(face="bold", hjust=0.5, size=18),
        plot.tag = element_text(face="bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plotHL

ggsave(filename = "results/01_immunosenescence_FIGURES/S1D_HL.tiff", plotHL, width = 5, height = 5, dpi = 300, units = "in", device='tiff')

# 2A./ ALPHA 1-GLOBULIN ----
dataalpha1 <- data[!is.na(data$alpha1),]
dataalpha1 <- dataalpha1[!is.na(dataalpha1$masse),] 
unique(factor(dataalpha1$numind)) # 1029 obs/485 ind

immuno_trajectory(dataalpha1, dataalpha1$alpha1)

# Retained model: QUADRATIC
# Age^2*sex + agelast + mass + datejulienne (AICc = 1672.62, df = 10)
dataalpha1$age2 <- I(dataalpha1$ageannee^2)
mod.quad <- lmer(alpha1 ~ age2*sexe + masse + datejulienne + agelast
                 + (1|numind) + (1|annee:pop),
                 data=dataalpha1, control=lmerControl("bobyqa"), REML=F)
summary(mod.quad)
r.squaredGLMM(mod.quad)

# Plot
dataalpha1$age2 <- I(dataalpha1$ageannee^2)
mod.quad <- lmer(alpha1 ~ age2*sexe
                 + (1|numind) + (1|annee:pop),
                 data=dataalpha1, control=lmerControl("bobyqa"), REML=F)
# same model without confounding variables for visual representation

dataalpha1$sexe <- as.factor(dataalpha1$sexe)
dataalpha1$age_factor_full <- as.factor(dataalpha1$age_factor_full)
dataalpha1$mean <- ave(dataalpha1$alpha1, dataalpha1$age_factor_full:dataalpha1$sexe)
dataalpha1$se <- ave(dataalpha1$alpha1, dataalpha1$age_factor_full:dataalpha1$sexe, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(dataalpha1[,c("sexe", "ageannee", "age_factor_full", "mean", "se")])

pred <- ggpredict(mod.quad, terms = c("age2", "sexe")) # to have intercept and slope for each pop
colnames(pred)[c(1,6)] <- c("ageannee", "sexe")
pred$ageannee <- sqrt(pred$ageannee)

dataerror <- merge(dataerror, pred, by=c("ageannee", "sexe"))

plotALPHA1 <- ggplot(dataerror, aes(x=ageannee, y=mean, color=sexe, fill=sexe)) +
  geom_point(aes(x=ageannee, y=mean, colour=sexe), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=sexe), width=.3, size=1.5) +
  geom_line(aes(x=ageannee, y=predicted), size=1) +
  geom_ribbon(aes(x=ageannee, y=predicted, ymin=conf.low, ymax=conf.high), alpha=0.2, linetype=0) +
  ylab(expression(Alpha1-globulin~(mg/mL))) +
  xlab(" ") +
  labs(title="Alpha1-globulin", tag="C") +
  scale_color_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  scale_fill_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        plot.tag=element_text(face="bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position=c(0.15, 0.83)) #change legend text font size
plotALPHA1

ggsave(filename = "results/01_immunosenescence_FIGURES/1C_ALPHA1.tiff", plotALPHA1, width = 5, height = 5, dpi = 300, units = "in", device='tiff')

# 2B./ ALPHA 2-GLOBULIN ----
dataalpha2 <- data[!is.na(data$alpha2),]
dataalpha2 <- dataalpha2[!is.na(dataalpha2$masse),] 
unique(factor(dataalpha2$numind)) # 1029 obs/485 ind

immuno_trajectory(dataalpha2, dataalpha2$alpha2)

# Retained model: CONSTANT
# Sex (AIC = 3702.95, df = 5)
mod.const <- lmer(alpha2 ~ sexe 
                 + (1|numind) + (1|annee:pop), 
                 data=dataalpha2, control=lmerControl("bobyqa"), REML=F)
summary(mod.const)
r.squaredGLMM(mod.const)

# Plot
dataalpha2$sexe <- as.factor(dataalpha2$sexe)
dataalpha2$age_factor_full <- as.factor(dataalpha2$age_factor_full)
dataalpha2$mean <- ave(dataalpha2$alpha2, dataalpha2$age_factor_full:dataalpha2$sexe)
dataalpha2$se <- ave(dataalpha2$alpha2, dataalpha2$age_factor_full:dataalpha2$sexe, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(dataalpha2[,c("sexe", "ageannee", "age_factor_full", "mean", "se")])

pred <- ggpredict(mod.const, terms = c("sexe"))
colnames(pred)[1] <- c("sexe")

dataerror <- merge(dataerror, pred, by=c("sexe"))

plotALPHA2 <- ggplot(dataerror, aes(x=ageannee, y=mean, color=sexe, fill=sexe)) +
  geom_point(aes(x=ageannee, y=mean, colour=sexe), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=sexe), width=.3, size=1.5) +
  geom_line(aes(x=ageannee, y=predicted), size=1) +
  geom_ribbon(aes(x=ageannee, y=predicted, ymin=conf.low, ymax=conf.high), alpha=0.2, linetype=0) +
  ylab(expression(Alpha2-globulin~(mg/mL))) +
  xlab("Age (years)") +
  labs(title="Alpha2-globulin", tag="E") +
  scale_color_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  scale_fill_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        plot.tag=element_text(face="bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position=c(0.15, 0.83)) #change legend text font size
plotALPHA2

ggsave(filename = "results/01_immunosenescence_FIGURES/S1E_ALPHA2_final.tiff", plotALPHA2, width = 5, height = 5, dpi = 300, units = "in", device='tiff')


# 2C./ BETAGLOBULIN ----
databeta <- data[!is.na(data$beta),]
databeta <- databeta[!is.na(databeta$masse),] 
unique(factor(databeta$numind)) # 1029 obs/485 ind

immuno_trajectory(databeta, databeta$beta)

# Retained model: LINEAR 
# Age + sex + mass (AICc = 4036.15, df = 7)
mod.linear <- lmer(beta ~ ageannee + sexe + masse
                   + (1|numind) + (1|annee:pop), 
                   data=databeta, control=lmerControl("bobyqa"), REML=F)
AICc(mod.linear)
summary(mod.linear)
r.squaredGLMM(mod.linear)

# Plot
mod.linear <- lmer(beta ~ ageannee + sexe
                   + (1|numind) + (1|annee:pop), 
                   data=databeta, control=lmerControl("bobyqa"), REML=F)
# same model without confounding variables for visual representation

databeta$sexe <- as.factor(databeta$sexe)
databeta$age_factor_full <- as.factor(databeta$age_factor_full)
databeta$mean <- ave(databeta$beta, databeta$age_factor_full:databeta$sexe)
databeta$se <- ave(databeta$beta, databeta$age_factor_full:databeta$sexe, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(databeta[,c("sexe", "ageannee", "age_factor_full", "mean", "se")])

pred <- ggpredict(mod.linear, terms = c("ageannee", "sexe")) # to have intercept and slope for each sex
colnames(pred)[c(1,6)] <- c("ageannee", "sexe")

dataerror <- merge(dataerror, pred, by=c("ageannee", "sexe"))

plotBETA <- ggplot(dataerror, aes(x=ageannee, y=mean, color=sexe, fill=sexe)) +
  geom_point(aes(x=ageannee, y=mean, colour=sexe), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=sexe), width=.3, size=1.5) +
  geom_line(aes(x=ageannee, y=predicted), linewidth=1) +
  geom_ribbon(aes(x=ageannee, y=predicted, ymin=conf.low, ymax=conf.high), alpha=0.2, linetype=0) +
  ylab(expression(Betaglobulin~(mg/mL))) +
  xlab(" ") +
  labs(title="Betaglobulin", tag="D") +
  scale_color_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  scale_fill_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        plot.tag = element_text(face="bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position=c(0.15, 0.83)) #change legend text font size
plotBETA

ggsave(filename = "results/01_immunosenescence_FIGURES/1D_BETA.tiff", plotBETA, width = 5, height = 5, dpi = 300, units = "in", device='tiff')


# 2D./ HAPTOGLOBIN ----
dataHAP <- data[!is.na(data$HAP),]
dataHAP <- dataHAP[!is.na(dataHAP$masse),] 
unique(factor(dataHAP$numind)) # 1034 obs/483 ind

immuno_trajectory(dataHAP, dataHAP$HAP)

# Retained model: QUADRATIC
# Age² + sex (AICc = 3173.39, df = 6)
dataHAP$age2 <- I(dataHAP$ageannee^2)
mod.quad <- lmer(HAP ~ age2 + sexe
                   + (1|numind) + (1|annee:pop), 
                   data=dataHAP, control=lmerControl("bobyqa"), REML=F)
AICc(mod.quad)
summary(mod.quad)
r.squaredGLMM(mod.quad)

# Plot
dataHAP$sexe <- as.factor(dataHAP$sexe)
dataHAP$age_factor_full <- as.factor(dataHAP$age_factor_full)
dataHAP$mean <- ave(dataHAP$HAP, dataHAP$age_factor_full:dataHAP$sexe)
dataHAP$se <- ave(dataHAP$HAP, dataHAP$age_factor_full:dataHAP$sexe, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(dataHAP[,c("sexe", "ageannee", "age_factor_full", "mean", "se")])

pred <- ggpredict(mod.quad, terms = c("age2", "sexe")); pred 
colnames(pred)[c(1,6)] <- c("ageannee", "sexe")
pred$ageannee <- sqrt(pred$ageannee)

dataerror <- merge(dataerror, pred, by=c("ageannee", "sexe"))

plotHAPTO <- ggplot(dataerror, aes(x=ageannee, y=mean, color=sexe, fill=sexe)) +
  geom_point(aes(x=ageannee, y=mean, colour=sexe), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=sexe), width=.3, linewidth=1.5) +
  geom_line(aes(x=ageannee, y=predicted), linewidth=1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2, linetype=0) +
  ylab(expression(Haptoglobin~(mg/mL))) +
  xlab(" ") +
  labs(title="Haptoglobin", tag="E") +
  scale_color_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  scale_fill_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        plot.tag = element_text(face="bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position=c(0.15, 0.83),) #change legend text font size
plotHAPTO

ggsave(filename = "results/01_immunosenescence_FIGURES/1E_HAPTO.tiff", plotHAPTO, width = 5, height = 5, dpi = 300, units = "in", device='tiff')

# 3A./ GAMMAGLOBULIN ----
datagamma <- data[!is.na(data$gamma),]
datagamma <- datagamma[!is.na(datagamma$masse),] 
unique(factor(datagamma$numind)) # 1029 obs/485 ind

immuno_trajectory(datagamma, datagamma$gamma)

# Retained model: LINEAR
# Age + mass + pop (AICc = 5665.51, df = 7)
mod.linear <- lmer(gamma ~ ageannee + pop + masse
                 + (1|numind) + (1|annee:pop),
                 data=datagamma, control=lmerControl("bobyqa"), REML=F)
summary(mod.linear)
r.squaredGLMM(mod.linear)

# Plot
mod.linear <- lmer(gamma ~ ageannee + pop
                 + (1|numind) + (1|annee:pop),
                 data=datagamma, control=lmerControl("bobyqa"), REML=F)
# same model without confounding variables for visual representation

datagamma$pop <- as.factor(datagamma$pop)
datagamma$age_factor_full <- as.factor(datagamma$age_factor_full)
datagamma$mean <- ave(datagamma$gamma, datagamma$age_factor_full:datagamma$pop)
datagamma$se <- ave(datagamma$gamma, datagamma$age_factor_full:datagamma$pop, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(datagamma[,c("pop", "ageannee", "age_factor_full", "mean", "se")])

pred <- ggpredict(mod.linear, terms = c("ageannee", "pop")) # to have intercept and slope for each pop
colnames(pred)[c(1,6)] <- c("ageannee", "pop")

dataerror <- merge(dataerror, pred, by=c("ageannee", "pop"))

plotGAMMA <- ggplot(dataerror, aes(x=ageannee, y=mean, color=pop, fill=pop)) +
  geom_point(aes(x=ageannee, y=mean, colour=pop), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=pop), width=.3, size=1.5) +
  geom_line(aes(x=ageannee, y=predicted), size=1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2, linetype=0) +
  ylab(expression(Gammaglobulin~(mg/mL))) +
  xlab("Age (years)") +
  labs(title="Gammaglobulin", tag="F") +
  scale_color_manual(name="Population", values=c("darkolivegreen3", "darkgreen")) +
  scale_fill_manual(name="Population", values=c("darkolivegreen3", "darkgreen")) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        plot.tag = element_text(face="bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position=c(0.15, 0.83)) #change legend text font size
plotGAMMA

ggsave(filename = "results/01_immunosenescence_FIGURES/1F_GAMMA.tiff", plotGAMMA, width = 5, height = 5, dpi = 300, units = "in", device='tiff')

# 3A./ LYMPHOCYTE ----
datalympho <- data[!is.na(data$lympho),]
datalympho <- datalympho[!is.na(datalympho$masse),] 
unique(factor(datalympho$numind)) # 1021 obs/484 ind

immuno_trajectory(datalympho, datalympho$lympho)

# Retained model : 2-SLOPES THRESHOLD (4yo)
# Age.1 + Age.2*pop + sex + Delay (AICc = 2659.99, df = 11)
{datalympho$age.1 <- datalympho$ageannee
  datalympho$age.2 <- datalympho$ageannee
  
  datalympho$age.1 <- datalympho$age.2 <- datalympho$ageannee
  datalympho$age.1 [datalympho$age.1 > 4] <- 4
  datalympho$age.2 [datalympho$age.2 <= 4] <- 0
  datalympho$age.2 [datalympho$age.2 > 4] <- datalympho$age.2 [datalympho$age.2 > 4] - 4
  mod.threshold_SS <- lmer(lympho ~ age.1 + age.2*pop + sexe + dureemin
                           + (1|numind) + (1|annee:pop),
                           data=datalympho, control=lmerControl("bobyqa"), REML=F)}
AICc(mod.threshold_SS)
summary(mod.threshold_SS)
r.squaredGLMM(mod.threshold_SS)

# Plot
{datalympho$age.1 <- datalympho$ageannee
  datalympho$age.2 <- datalympho$ageannee
  
  datalympho$age.1 <- datalympho$age.2 <- datalympho$ageannee
  datalympho$age.1 [datalympho$age.1 > 4] <- 4
  datalympho$age.2 [datalympho$age.2 <= 4] <- 0
  datalympho$age.2 [datalympho$age.2 > 4] <- datalympho$age.2 [datalympho$age.2 > 4] - 4
  mod.threshold_SS <- lmer(lympho ~ age.1 + age.2*pop + sexe
                           + (1|numind) + (1|annee:pop),
                           data=datalympho, control=lmerControl("bobyqa"), REML=F)}
# same model without confounding variables for visual representation

datalympho$sexe <- as.factor(datalympho$sexe)
datalympho$pop <- as.factor(datalympho$pop)
datalympho$age_factor_full <- as.factor(datalympho$age_factor_full)
datalympho$mean <- ave(datalympho$lympho, datalympho$age_factor_full:datalympho$pop:datalympho$sexe)
datalympho$se <- ave(datalympho$lympho, datalympho$age_factor_full:datalympho$pop:datalympho$sexe, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(datalympho[,c("ageannee", "age_factor_full", "mean", "se", "sexe", "pop")])

pred1 <- ggpredict(mod.threshold_SS, terms = c("age.1", "sexe", "pop"), condition = c("age.2"=0)); pred1 # to have intercept and slope for each pop
colnames(pred1)[c(1,6,7)] <- c("ageannee", "sexe", "pop")
pred2 <- ggpredict(mod.threshold_SS, terms = c("age.2", "sexe", "pop"), condition = c("age.1"=4)); pred1 # to have intercept and slope for each pop
colnames(pred2)[c(1,6,7)] <- c("ageannee", "sexe", "pop")
pred2$ageannee <- pred2$ageannee+4

pred <- rbind(pred1, pred2)

dataerror <- merge(dataerror, pred, by=c("ageannee", "sexe", "pop"))

dataerror$sexe <- ifelse(dataerror$sexe=="F", "Females", "Males")

plotLYMPHO <- ggplot(dataerror, aes(x=ageannee, y=mean, color=pop, fill=pop)) +
  geom_point(aes(x=ageannee, y=mean, colour=pop), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=pop), width=.3, linewidth=1.5) +
  geom_line(aes(x=ageannee, y=predicted), size=1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2, linetype=0) +
  geom_vline(xintercept=4, linetype='dashed', linewidth=0.75) +
  ylab(expression(Lymphocyte~count~(10^3~cells/mL))) +
  xlab("Age (years)") +
  scale_x_continuous(breaks = seq(4, 16, by = 4)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1)) +
  labs(title = "Lymphocyte", tag="G") +
  scale_color_manual(name="Population", values=c("darkolivegreen3", "darkgreen")) +
  scale_fill_manual(name="Population", values=c("darkolivegreen3", "darkgreen")) +
  facet_rep_grid(.~sexe, scales="fixed",
                 repeat.tick.labels = T) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        axis.title.x = element_text(size=14),
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
        legend.position=c(0.15, 0.83),
        legend.background=element_blank(),
        strip.text.x = element_text(size=13, face="bold"),
        strip.text.y = element_text(size=13, face="bold"),
        strip.background = element_blank()) #change legend text font size
plotLYMPHO

ggsave(filename = "results/01_immunosenescence_FIGURES/1G_LYMPHO.tiff", plotLYMPHO, width = 10, height = 5, dpi = 300, units = "in", device='tiff')

# COMBINE ALL FIGURES ----
# Figure 1
neutro <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/1A_NEUTRO.tiff"))
mono <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/1B_MONO.tiff"))
alpha1 <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/1C_ALPHA1.tiff"))
beta <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/1D_BETA.tiff"))
hapto <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/1E_HAPTO.tiff"))
gamma <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/1F_GAMMA.tiff"))
lympho <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/1G_LYMPHO.tiff"))

figure1 <- grid.arrange(neutro, mono, alpha1, 
                        beta, hapto, gamma,
                        lympho,
                        layout_matrix = matrix(c(1, 2, 3, 4, 5, 6, 7, 7, NA), 
                                               byrow = TRUE, ncol = 3))

grid.newpage()
grid.draw(figure1)
ggsave(filename = "results/01_immunosenescence_FIGURES/FIGURE1.tiff", figure1, width = 13, height = 13, dpi = 300, units = "cm", device='tiff')


# Figure S1
baso <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/S1A_BASO.tiff"))
eosino <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/S1B_EOSINO.tiff"))
HA <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/S1C_HA.tiff"))
HL <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/S1D_HL.tiff"))
alpha2 <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/S1E_ALPHA2.tiff"))

figureS1 <- grid.arrange(baso, eosino, HA,
                         HL, alpha2,
                         layout_matrix = matrix(c(1,2,3,4,5,NA), 
                                                byrow = TRUE, ncol = 3))

grid.newpage()
grid.draw(figureS1)
ggsave(filename = "results/01_immunosenescence_FIGURES/FIGURES1.tiff", figureS1, width = 13, height = 8, dpi = 300, units = "cm", device='tiff')


# II - PARASITISM ----
# 1./ GASTRO-INTESTINAL STRONGYLES ----
datastrongyles <- data[!is.na(data$SD),]
datastrongyles <- datastrongyles[!is.na(datastrongyles$masse),] 
unique(factor(datastrongyles$numind)) # 976 obs/469 ind

immuno_trajectory_parasitism(datastrongyles, datastrongyles$SD)

# Retained model: 2-SLOPES THRESHOLD (10 years old)
# age.1 + age.2*sexe + agelast + masse (AICc = 3619.93, df = 11)
{datastrongyles$age.1 <- datastrongyles$ageannee
  datastrongyles$age.2 <- datastrongyles$ageannee
  
  datastrongyles$age.1 <- datastrongyles$age.2 <- datastrongyles$ageannee
  datastrongyles$age.1 [datastrongyles$age.1 > 10] <- 10
  datastrongyles$age.2 [datastrongyles$age.2 <= 10] <- 0
  datastrongyles$age.2 [datastrongyles$age.2 > 10] <- datastrongyles$age.2 [datastrongyles$age.2 > 10] - 10
  mod.threshold_SS <- lmer(log(SD+1) ~ age.1 + age.2*sexe + agelast + masse
                           + (1|numind) + (1|annee:pop),
                           data=datastrongyles, control=lmerControl("bobyqa"), REML=F)}
AICc(mod.threshold_SS)
summary(mod.threshold_SS)
r.squaredGLMM(mod.threshold_SS)


# Plot
datastrongyles$logSD <- log(datastrongyles$SD+1)
{datastrongyles$age.1 <- datastrongyles$ageannee
  datastrongyles$age.2 <- datastrongyles$ageannee
  
  datastrongyles$age.1 <- datastrongyles$age.2 <- datastrongyles$ageannee
  datastrongyles$age.1 [datastrongyles$age.1 > 10] <- 10
  datastrongyles$age.2 [datastrongyles$age.2 <= 10] <- 0
  datastrongyles$age.2 [datastrongyles$age.2 > 10] <- datastrongyles$age.2 [datastrongyles$age.2 > 10] - 10
  mod.threshold_SS <- lmer(logSD ~ age.1 + age.2*sexe
                           + (1|numind) + (1|annee:pop),
                           data=datastrongyles, control=lmerControl("bobyqa"), REML=F)}
# same model without confounding variables for visual representation

datastrongyles$sexe <- as.factor(datastrongyles$sexe)
datastrongyles$age_factor_full <- as.factor(datastrongyles$age_factor_full)
datastrongyles$mean <- ave(datastrongyles$logSD, datastrongyles$age_factor_full:datastrongyles$sexe)
datastrongyles$se <- ave(datastrongyles$logSD, datastrongyles$age_factor_full:datastrongyles$sexe, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(datastrongyles[,c("sexe", "ageannee", "age_factor_full", "mean", "se")])

pred1 <- ggpredict(mod.threshold_SS, terms = c("age.1", "sexe"), condition = c(age.2=0)) # to have intercept and slope for each sex
colnames(pred1)[c(1,6)] <- c("ageannee", "sexe")
pred2 <- ggpredict(mod.threshold_SS, terms = c("age.2", "sexe"), condition = c(age.1=10)) # to have intercept and slope for each sex
colnames(pred2)[c(1,6)] <- c("ageannee", "sexe")
pred2$ageannee <- pred2$ageannee+10

pred <- rbind(pred1, pred2)

dataerror <- merge(dataerror, pred, by=c("ageannee", "sexe"))

plotSTRONGYLES <- ggplot(dataerror, aes(x=ageannee, y=mean, color=sexe, fill=sexe)) +
  geom_point(aes(x=ageannee, y=mean, colour=sexe), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=sexe), width=.3, size=1.5) +
  geom_line(aes(x=ageannee, y=predicted), linewidth=1) +
  geom_ribbon(aes(x=ageannee, y=predicted, ymin=conf.low, ymax=conf.high), alpha=0.2, linetype=0) +
  geom_vline(xintercept=10, linetype="dashed", linewidth=0.75) +
  ylab(expression(paste("Gastro-intestinal Strongyles [EPG, log(n+1)]"))) +
  xlab(" ") +
  labs(title="Gastro-intestinal strongyles", tag="A") + 
  scale_color_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  scale_fill_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        plot.tag = element_text(face="bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position=c(0.15, 0.83)) #change legend text font size
plotSTRONGYLES

ggsave(filename = "results/01_immunosenescence_FIGURES/2A_STRONGYLES.tiff", plotSTRONGYLES, width = 5, height = 5, dpi = 300, units = "in", device='tiff')


# 2./ TRICHURIS SP. ----
datatrich <- data[!is.na(data$Trich),]
datatrich <- datatrich[!is.na(datatrich$masse),] 
unique(factor(datatrich$numind)) # 938 obs/450 ind

immuno_trajectory_parasitism(datatrich, datatrich$Trich)

# Retained model: 2-SLOPES THRESHOLD (10 years old)
# Age.1*pop + age.1*sex + age.2*pop + age.2*sexe + masse  (AICc = 3478.44, df = 14)
{datatrich$age.1 <- datatrich$ageannee
  datatrich$age.2 <- datatrich$ageannee
  
  datatrich$age.1 <- datatrich$age.2 <- datatrich$ageannee
  datatrich$age.1 [datatrich$age.1 > 10] <- 10
  datatrich$age.2 [datatrich$age.2 <= 10] <- 0
  datatrich$age.2 [datatrich$age.2 > 10] <- datatrich$age.2 [datatrich$age.2 > 10] - 10
  mod.threshold_SS <- lmer(log(Trich+1) ~ age.1*sexe + age.2*sexe + age.1*pop + age.2*pop + masse
                           + (1|numind) + (1|annee:pop),
                           data=datatrich, control=lmerControl("bobyqa"), REML=F)}
AICc(mod.threshold_SS)
summary(mod.threshold_SS)
r.squaredGLMM(mod.threshold_SS)

# Plot
datatrich$logTrich <- log(datatrich$Trich+1)

{datatrich$age.1 <- datatrich$ageannee
  datatrich$age.2 <- datatrich$ageannee
  
  datatrich$age.1 <- datatrich$age.2 <- datatrich$ageannee
  datatrich$age.1 [datatrich$age.1 > 10] <- 10
  datatrich$age.2 [datatrich$age.2 <= 10] <- 0
  datatrich$age.2 [datatrich$age.2 > 10] <- datatrich$age.2 [datatrich$age.2 > 10] - 10
  mod.threshold_SS <- lmer(logTrich ~ age.1*sexe + age.2*sexe + age.1*pop + age.2*pop
                           + (1|numind) + (1|annee:pop),
                           data=datatrich, control=lmerControl("bobyqa"), REML=F)}
# same model without confounding variables for visual representation

datatrich$sexe <- as.factor(datatrich$sexe)
datatrich$pop <- as.factor(datatrich$pop)
datatrich$age_factor_full <- as.factor(datatrich$age_factor_full)
datatrich$mean <- ave(datatrich$logTrich, datatrich$age_factor_full:datatrich$sexe:datatrich$pop)
datatrich$se <- ave(datatrich$logTrich, datatrich$age_factor_full:datatrich$sexe:datatrich$pop, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(datatrich[,c("sexe", "pop", "ageannee", "age_factor_full", "mean", "se")])

pred1 <- ggpredict(mod.threshold_SS, terms = c("age.1", "sexe", "pop"), condition = c(age.2=0)) # to have intercept and slope for each sex
colnames(pred1)[c(1,6,7)] <- c("ageannee", "sexe", "pop")
pred2 <- ggpredict(mod.threshold_SS, terms = c("age.2", "sexe", "pop"), condition = c(age.1=10)) # to have intercept and slope for each sex
colnames(pred2)[c(1,6,7)] <- c("ageannee", "sexe", "pop")
pred2$ageannee <- pred2$ageannee+10

pred <- rbind(pred1, pred2)

dataerror <- merge(dataerror, pred, by=c("ageannee", "sexe", "pop"))

dataerror$pop <- ifelse(dataerror$pop=="TF", "Trois-Fontaines", "Chizé")

plotTRICH <- ggplot(dataerror, aes(x=ageannee, y=mean, color=sexe, fill=sexe)) +
  geom_point(aes(x=ageannee, y=mean, colour=sexe), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=sexe), width=.3, linewidth=1.5) +
  geom_line(aes(x=ageannee, y=predicted), linewidth=1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2, linetype=0) +
  geom_vline(xintercept=10, linetype='dashed', linewidth=0.75) +
  ylab(expression(italic(Trichuris)~"sp. [EPG, log(n+1)]")) +
  xlab(" ") +
  scale_x_continuous(breaks = seq(4, 16, by = 4)) +
  labs(title = expression(bolditalic(Trichuris)~bold(sp.)), tag="B") +
  scale_color_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  scale_fill_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  facet_rep_grid(.~pop, scales="fixed",
                 repeat.tick.labels = T) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        axis.title.x = element_text(size=14),
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
        legend.position=c(0.075, 0.91),
        legend.background=element_blank(),
        strip.text.x = element_text(size=13, face="bold"),
        strip.background = element_blank()) #change legend text font size
plotTRICH

ggsave(filename = "results/01_immunosenescence_FIGURES/2B_TRICHURIS.tiff", plotTRICH, width = 10, height = 5, dpi = 300, units = "in", device='tiff')


# 3./ PROTOSTRONGYLIDS ----
dataproto <- data[!is.na(data$Proto),]
dataproto <- dataproto[!is.na(dataproto$masse),] 
unique(factor(dataproto$numind)) # 915 obs/449 ind


immuno_trajectory_parasitism(dataproto, dataproto$Proto)

# Retained model: QUADRATIC
# Age + Age²*cohort + + Age²*sex + agelast + masse + pop (AICc= 1924.77, df = 14)
mod.quad <- lmer(log(Proto+1) ~ ageannee*cohort_class + I(ageannee^2)*sexe + I(ageannee^2)*cohort_class + pop + masse + agelast
                       + (1|numind) + (1|annee:pop), 
                       data=dataproto, control=lmerControl("bobyqa"), REML=F)
AICc(mod.quad)
summary(mod.quad)
r.squaredGLMM(mod.quad)

# Plot
dataproto$logProto <- log(dataproto$Proto+1)
mod.quad <- lmer(logProto ~ ageannee*cohort_class + I(ageannee^2)*sexe + I(ageannee^2)*cohort_class + pop
                 + (1|numind) + (1|annee:pop), 
                 data=dataproto, control=lmerControl("bobyqa"), REML=F)
# same model without confounding variables for visual representation

dataproto$sexe <- as.factor(dataproto$sexe)
dataproto$pop <- as.factor(dataproto$pop)
dataproto$cohort_class <- as.factor(dataproto$cohort_class)
dataproto$age_factor_full <- as.factor(dataproto$age_factor_full)
dataproto$mean <- ave(dataproto$logProto, dataproto$age_factor_full:dataproto$sexe:dataproto$pop:dataproto$cohort_class)
dataproto$se <- ave(dataproto$logProto, dataproto$age_factor_full:dataproto$sexe:dataproto$pop:dataproto$cohort_class, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(dataproto[,c("ageannee", "sexe" , "pop", "cohort_class", "age_factor_full", "mean", "se")])

pred1 <- ggpredict(mod.quad, terms = c("ageannee [all]", "sexe", "pop", "cohort_class")) # to have intercept and slope for each pop
colnames(pred1)[c(1,6,7,8)] <- c("ageannee", "sexe", "pop", "cohort_class")

dataerror <- merge(dataerror, pred1, by=c("ageannee", "sexe", "pop", "cohort_class"))

dataerror$pop <- ifelse(dataerror$pop=="CH", "Chizé", "Trois-Fontaines")
dataerror$cohort_class <- ifelse(dataerror$cohort_class=="good", "Good cohort quality", "Poor cohort quality")

plotPROTO <- ggplot(dataerror, aes(x=ageannee, y=mean, color=sexe, fill=sexe)) +
  geom_point(aes(x=ageannee, y=mean, colour=sexe), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=sexe), width=.3, linewidth=1.5) +
  geom_line(aes(x=ageannee, y=predicted), size=1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2, linetype=0) +
  ylab(expression("Protostrongylids [LPG, log(n+1)]")) +
  xlab("Age (years)") +
  scale_x_continuous(breaks = seq(4, 16, by = 4)) +
  labs(title = "Protostrongylids", tag="C") +
  scale_y_continuous(breaks=seq(0,4,1)) +
  scale_color_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  scale_fill_manual(name="Sex", values=c("darkorchid4", "darkorchid1")) +
  facet_rep_grid(cohort_class~pop, scales="fixed",
                 repeat.tick.labels = T) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        axis.title.x = element_text(size=14),
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
        legend.position=c(0.075, 0.91),
        legend.background=element_blank(),
        strip.text.x = element_text(size=13, face="bold"),
        strip.text.y = element_text(size=13, face="bold"),
        strip.background = element_blank()) #change legend text font size
plotPROTO

ggsave(filename = "results/01_immunosenescence_FIGURES/2C_PROTO.tiff", plotPROTO, width = 10, height = 10, dpi = 300, units = "in", device='tiff')

# 4./ COCCIDIA ----
datacoccidia <- data[!is.na(data$Coc),]
datacoccidia <- datacoccidia[!is.na(datacoccidia$masse),] 
unique(factor(datacoccidia$numind)) # 974 obs/469 ind

immuno_trajectory_parasitism(datacoccidia, datacoccidia$Coc)

# Retained model: CONSTANT
# NULL (AICc = 4060.780, df = 4)

# Plot
datacoccidia$logCoc <- log(datacoccidia$Coc+1)
mod.linear <- lmer(logCoc ~
                     + (1|numind) + (1|annee:pop), 
                   data=datacoccidia, control=lmerControl("bobyqa"), REML=F)

summary(mod.linear)
AICc(mod.linear)
r.squaredGLMM(mod.linear)

datacoccidia$age_factor_full <- as.factor(datacoccidia$age_factor_full)
datacoccidia$mean <- ave(datacoccidia$logCoc, datacoccidia$age_factor_full)
datacoccidia$se <- ave(datacoccidia$logCoc, datacoccidia$age_factor_full, FUN = function(x) sqrt(var(x) / length(x)))
dataerror <- unique(datacoccidia[,c("ageannee", "age_factor_full", "mean", "se")])

plotCOCCIDIA <- ggplot(dataerror, aes(x=ageannee, y=mean)) +
  geom_smooth(aes(x=ageannee, y=1.1692), color="black") +
  geom_ribbon(aes(x=ageannee, y=mean, ymax=1.1692+1.96*0.0717, ymin=1.1692-1.96*0.0717), alpha=0.2, linetype=0) +
  geom_point(aes(x=ageannee, y=mean), size=2.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, size=1.5) +
  ylab(expression("Coccidia [OPG, log(n+1)]")) +
  xlab("Age (years)") +
  labs(title = expression(bold(Coccidia)~bolditalic(Eimeria)~bold(sp.)), tag="D") +
  scale_x_continuous(breaks = seq(4, 16, by = 4)) +
  theme(plot.title = element_text(face = "bold", hjust=0.5, size=18),
        plot.tag = element_text(face="bold", size=18),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=11.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plotCOCCIDIA

ggsave(filename = "results/01_immunosenescence_FIGURES/2D_COCCIDIA.tiff", plotCOCCIDIA, width = 5, height = 5, dpi = 300, units = "in", device='tiff')

# COMBINE ALL FIGURES ----
# Figure 2
SD <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/2A_STRONGYLES.tiff"))
Trich <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/2B_TRICHURIS.tiff"))
Proto <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/2C_PROTO.tiff"))
Coc <- rasterGrob(readTIFF("results/01_immunosenescence_FIGURES/2D_COCCIDIA.tiff"))

figure2 <- grid.arrange(SD, Trich,
                        Proto, Coc,
                        layout_matrix = matrix(c(1, 2, 2, 3, 3, 4, 3, 3, NA), 
                                              byrow = 2, ncol = 3))

grid.newpage()
grid.draw(figure2)
ggsave(filename = "results/01_immunosenescence_FIGURES/FIGURE2.tiff", figure2, width = 13, height = 13, dpi = 300, units = "cm", device='tiff')

