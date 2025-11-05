### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## Early-life glucocorticoids accelerate the senescence rate of lymphocyte counts in wild roe deer
##
## Immunosenescence patterns determination function
## 
## Lucas Lalande et al. 2024
## May 2024
##
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# IMMUNOSENESCENCE: IMMUNE TRAITS ----
# Function to determine the ageing trajectory of immune traits
immuno_trajectory <- function(subset, trait) {
  
  options(na.action="na.fail")
  
  
  # constant ----
  mod.constant <- lmer(trait ~ pop + sexe + cohort_class + 
                       + agelast + masse + datejulienne + dureemin
                     + (1|numind) + (1|annee:pop),
                     data=subset, control=lmerControl("bobyqa"), REML=F)
  
  cat("\n\nConstant model is running:\n")
  d.constant <- dredge(mod.constant, evaluate=T, rank="AICc", trace=2)
  d1.constant <- subset(d.constant, delta <=2, recalc.weights = FALSE)
  d2.constant <- subset(d1.constant, which(df==min(df)), recalc.weights = FALSE)
  d3.constant <- subset(d2.constant, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # linear ----
  mod.linear <- lmer(trait ~ ageannee*pop + ageannee*sexe + ageannee*cohort_class + 
                   + agelast + masse + datejulienne + dureemin
                   + (1|numind) + (1|annee:pop),
                   data=subset, control=lmerControl("bobyqa"), REML=F)
  
  cat("\n\nLinear model is running:\n")
  d.linear <- dredge(mod.linear, evaluate=T, rank="AICc", trace=2, fixed=~ageannee)
  d1.linear <- subset(d.linear, delta <=2, recalc.weights = FALSE)
  d2.linear <- subset(d1.linear, which(df==min(df)), recalc.weights = FALSE)
  d3.linear <- subset(d2.linear, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # quadratic ----
  mod.quad <- lmer(trait ~ ageannee*pop + ageannee*sexe + ageannee*cohort_class + I(ageannee^2)*pop + I(ageannee^2)*sexe + I(ageannee^2)*cohort_class 
                   + agelast + masse + datejulienne + dureemin
                   + (1|numind) + (1|annee:pop),
                   data=subset, control=lmerControl("bobyqa"), REML=F)
  
  cat("\n\nQuadratic model is running:\n")
  d.quad <- dredge(mod.quad, evaluate=T, rank="AICc", trace=2, fixed=~I(ageannee^2))
  d1.quad <- subset(d.quad, delta <=2, recalc.weights = FALSE)
  d2.quad <- subset(d1.quad, which(df==min(df)), recalc.weights = FALSE)
  d3.quad <- subset(d2.quad, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # threshold-CS: constant - slope ----
  {subset$age.2 <- subset$ageannee
    
    threshCS <- data.frame()
    sCS <- seq(3,10,1)
    
    cat("\n\nTest breakpoint ages\n\n")
    
    for (i in sCS) {
      subset$age.2 <- subset$ageannee
      print(paste('Fitting threshold model with breakpoint age', i))
      subset$age.2 <- ifelse(subset$age.2 < i, 0, subset$age.2 - i)
      
      mod.threshCS <- lmer(trait ~ age.2*pop + age.2*sexe + age.2*cohort_class 
                          + agelast + masse + datejulienne + dureemin
                          + (1|numind) + (1|annee:pop), 
                          data=subset, control=lmerControl("bobyqa"), REML=F)
      
      k <- AICc(mod.threshCS)
      
      threshCS <- rbind(threshCS, k)
    }
    threshCS <- cbind(sCS, threshCS)
    colnames(threshCS) <- c("age","AICc")
    par(mfrow=c(1,1))
    plot(AICc ~ age, threshCS,
         xlab="age", ylab="AICc")
  }
  
  print(paste('BEST BREAKPOINT AGE IS', threshCS$age[which.min(threshCS$AICc)]))
  
  {subset$age.2 <- subset$ageannee
    subset$age.2 <- ifelse(subset$age.2 <  threshCS$age[which.min(threshCS$AICc)], 0, subset$age.2 -  threshCS$age[which.min(threshCS$AICc)])
    
    mod.thresholdCS <- lmer(trait ~ age.2*pop + age.2*sexe + age.2*cohort_class 
                           + agelast + masse + datejulienne + dureemin
                           + (1|numind) + (1|annee:pop), 
                           data=subset, control=lmerControl("bobyqa"), REML=F)}
  
  cat("\n\n1-slope Threshold model (constant/slope) is running:\n")
  d.thresholdCS <- dredge(mod.thresholdCS, evaluate=T, rank="AICc", trace=2, fixed=~age.2)
  d1.thresholdCS <- subset(d.thresholdCS, delta <=2, recalc.weights = FALSE)
  d2.thresholdCS <- subset(d1.thresholdCS, which(df==min(df)), recalc.weights = FALSE)
  d3.thresholdCS <- subset(d2.thresholdCS, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # threshold-SC: slope - constant ----
  {subset$age.1 <- subset$ageannee
  
  threshSC <- data.frame()
  sSC <- seq(3,10,1)
  
  cat("\n\nTest breakpoint ages\n\n")
  
  for (i in sSC) {
    subset$age.1 <- subset$ageannee
    print(paste('Fitting threshold model with breakpoint age', i))
    subset$age.1 <- ifelse(subset$age.1 < i, subset$age.1 - i, 0)
    
    mod.threshSC <- lmer(trait ~ age.1*pop + age.1*sexe + age.1*cohort_class 
                        + agelast + masse + datejulienne + dureemin
                        + (1|numind) + (1|annee:pop), 
                        data=subset, control=lmerControl("bobyqa"), REML=F)
    
    k <- AICc(mod.threshSC)
    
    threshSC <- rbind(threshSC, k)
  }
  threshSC <- cbind(sSC, threshSC)
  colnames(threshSC) <- c("age","AICc")
  par(mfrow=c(1,1))
  plot(AICc ~ age, threshSC,
       xlab="age", ylab="AICc")
  }
  
  print(paste('BEST BREAKPOINT AGE IS', threshSC$age[which.min(threshSC$AICc)]))
  
  {subset$age.1 <- subset$ageannee
    subset$age.1 <- ifelse(subset$age.1 <  threshSC$age[which.min(threshSC$AICc)], subset$age.1 -  threshSC$age[which.min(threshSC$AICc)], 0)
    
    mod.thresholdSC <- lmer(trait ~ age.1*pop + age.1*sexe + age.1*cohort_class 
                           + agelast + masse + datejulienne + dureemin
                           + (1|numind) + (1|annee:pop), 
                           data=subset, control=lmerControl("bobyqa"), REML=F)}
  
  cat("\n\n1-slope Threshold model (slope/constant) is running:\n")
  d.thresholdSC <- dredge(mod.thresholdSC, evaluate=T, rank="AICc", trace=2, fixed=~age.1)
  d1.thresholdSC <- subset(d.thresholdSC, delta <=2, recalc.weights = FALSE)
  d2.thresholdSC <- subset(d1.thresholdSC, which(df==min(df)), recalc.weights = FALSE)
  d3.thresholdSC <- subset(d2.thresholdSC, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # threshold-SS: slope - slope ----
  {subset$age.1 <- subset$ageannee
  subset$age.2 <- subset$ageannee
  
  threshSS <- data.frame()
  sSS <- seq(3,10,1)
  
  cat("\n\nTest breakpoint ages\n\n")
  
  for (i in sSS) {
    subset$age.1 <- subset$age.2 <- subset$ageannee
    print(paste('Fitting threshold model with breakpoint age', i))
    subset$age.1 [subset$age.1 > i] <- i
    subset$age.2 [subset$age.2 <= i] <- 0
    subset$age.2 [subset$age.2 > i] <- subset$age.2 [subset$age.2 > i] - i
    
    mod.threshSS <- lmer(trait ~ age.1*pop + age.1*sexe + age.1*cohort_class + age.2*pop + age.2*sexe + age.2*cohort_class 
                        + agelast + masse + datejulienne + dureemin
                        + (1|numind) + (1|annee:pop), 
                        data=subset, control=lmerControl("bobyqa"), REML=F)
    
    k <- AICc(mod.threshSS)
    
    threshSS <- rbind(threshSS, k)
  }
  threshSS <- cbind(sSS, threshSS)
  colnames(threshSS) <- c("age","AICc")
  par(mfrow=c(1,1))
  plot(AICc ~ age, threshSS,
       xlab="age", ylab="AICc")
  }
  
  print(paste('BEST BREAKPOINT AGE IS', threshSS$age[which.min(threshSS$AICc)]))
  
  {subset$age.1 <- subset$ageannee
    subset$age.2 <- subset$ageannee
    
    subset$age.1 <- subset$age.2 <- subset$ageannee
    subset$age.1 [subset$age.1 > threshSS$age[which.min(threshSS$AICc)]] <- threshSS$age[which.min(threshSS$AICc)]
    subset$age.2 [subset$age.2 <= threshSS$age[which.min(threshSS$AICc)]] <- 0
    subset$age.2 [subset$age.2 > threshSS$age[which.min(threshSS$AICc)]] <- subset$age.2 [subset$age.2 > threshSS$age[which.min(threshSS$AICc)]] - threshSS$age[which.min(threshSS$AICc)]
    
    mod.thresholdSS <- lmer(trait ~ age.1*pop + age.1*sexe + age.1*cohort_class + age.2*pop + age.2*sexe + age.2*cohort_class 
                           + agelast + masse + datejulienne + dureemin
                           + (1|numind) + (1|annee:pop), 
                           data=subset, control=lmerControl("bobyqa"), REML=F)}
  
  cat("\n\n2-slopes Threshold model (slope/slope) is running:\n")
  d.thresholdSS <- dredge(mod.thresholdSS, evaluate=T, rank="AICc", trace=2, fixed=~age.1 + age.2)
  d1.thresholdSS <- subset(d.thresholdSS, delta <=2, recalc.weights = FALSE)
  d2.thresholdSS <- subset(d1.thresholdSS, which(df==min(df)), recalc.weights = FALSE)
  d3.thresholdSS <- subset(d2.thresholdSS, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  cat("CONSTANT MODEL\n\n")
  print(data.frame(d1.constant))
  cat("\n") 
  print(data.frame(d3.constant))
  cat("\n\n\n")
  cat("LINEAR MODEL\n\n")
  print(data.frame(d1.linear))
  cat("\n")
  print(data.frame(d3.linear))
  cat("\n\n\n")
  cat("QUADRATIC MODEL\n\n")
  print(data.frame(d1.quad))
  cat("\n")
  print(data.frame(d3.quad))
  cat("\n\n\n")
  cat("1-SLOPE THRESHOLD MODEL; Constant/slope", "\n", paste("The retained threshold is", threshCS$age[which.min(threshCS$AICc)], "years old\n\n"))
  print(data.frame(d1.thresholdCS))
  cat("\n")
  print(data.frame(d3.thresholdCS))
  cat("\n\n\n")
  cat("1-SLOPE THRESHOLD MODEL; Slope/constant", "\n", paste("The retained threshold is", threshSC$age[which.min(threshSC$AICc)], "years old\n\n"))
  print(data.frame(d1.thresholdSC))
  cat("\n")
  print(data.frame(d3.thresholdSC))
  cat("\n\n\n")
  cat("2-SLOPES THRESHOLD MODEL", "\n", paste("The retained threshold is", threshSS$age[which.min(threshSS$AICc)], "years old\n\n"))
  print(data.frame(d1.thresholdSS))
  cat("\n")
  print(data.frame(d3.thresholdSS))
  cat("\n\n\n")
  
  recap <- data.frame("Model type"=c("Constant",    
                                     "Linear",
                                     "Quadratic",
                                     "Threshold-CS",
                                     "Threshold-SC",
                                     "Threshold-SS"),
                      "Threshold"=c("NA",
                                    "NA",
                                    "NA",
                                    threshCS$age[which.min(threshCS$AICc)],
                                    threshSC$age[which.min(threshSC$AICc)],
                                    threshSS$age[which.min(threshSS$AICc)]),
                      "AICc"=c(d3.constant$AICc,
                               d3.linear$AICc,
                               d3.quad$AICc, 
                               d3.thresholdCS$AICc+2,
                               d3.thresholdSC$AICc+2,
                               d3.thresholdSS$AICc+2),
                      "DF"=c(d3.constant$df,
                             d3.linear$df,
                             d3.quad$df, 
                             d3.thresholdCS$df+1,
                             d3.thresholdSC$df+1,
                             d3.thresholdSS$df+1))
  
  
  recap$delta <- recap$AICc-min(recap$AICc)
  recap$weights <- exp(-0.5*recap$delta)/sum(exp(-0.5*recap$delta))
  cat("\n Summary of the model AICcs and DFs according to model type\n\n")
  print(recap)
  cat("\n\n")
}

# IMMUNOSENESCENCE: HAEMOLYSIS SCORE ----
# Function to determine the ageing trajectory of haemolysis score
immuno_trajectory_HL <- function(subset, trait) {
  
  options(na.action="na.fail")
  
  # constant ----
  mod.constant <- lmer(trait ~ pop + sexe + cohort_class + 
                         + agelast + masse + datejulienne + dureemin + couleur_HAHL
                       + (1|numind) + (1|annee:pop),
                       data=subset, control=lmerControl("bobyqa"), REML=F)
  
  cat("\n\nConstant model is running:\n")
  d.constant <- dredge(mod.constant, evaluate=T, rank="AICc", trace=2)
  d1.constant <- subset(d.constant, delta <=2, recalc.weights = FALSE)
  d2.constant <- subset(d1.constant, which(df==min(df)), recalc.weights = FALSE)
  d3.constant <- subset(d2.constant, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # linear ----
  mod.linear <- lmer(trait ~ ageannee*pop + ageannee*sexe + ageannee*cohort_class + 
                       + agelast + masse + datejulienne + dureemin + couleur_HAHL
                     + (1|numind) + (1|annee:pop),
                     data=subset, control=lmerControl("bobyqa"), REML=F)
  
  cat("\n\nLinear model is running:\n")
  d.linear <- dredge(mod.linear, evaluate=T, rank="AICc", trace=2, fixed=~ageannee)
  d1.linear <- subset(d.linear, delta <=2, recalc.weights = FALSE)
  d2.linear <- subset(d1.linear, which(df==min(df)), recalc.weights = FALSE)
  d3.linear <- subset(d2.linear, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # quadratic ----
  mod.quad <- lmer(trait ~ ageannee*pop + ageannee*sexe + ageannee*cohort_class + I(ageannee^2)*pop + I(ageannee^2)*sexe + I(ageannee^2)*cohort_class 
                   + agelast + masse + datejulienne + dureemin + couleur_HAHL
                   + (1|numind) + (1|annee:pop),
                   data=subset, control=lmerControl("bobyqa"), REML=F)
  
  cat("\n\nQuadratic model is running:\n")
  d.quad <- dredge(mod.quad, evaluate=T, rank="AICc", trace=2, fixed=~I(ageannee^2))
  d1.quad <- subset(d.quad, delta <=2, recalc.weights = FALSE)
  d2.quad <- subset(d1.quad, which(df==min(df)), recalc.weights = FALSE)
  d3.quad <- subset(d2.quad, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # threshold-CS: constant - slope ----
  {subset$age.2 <- subset$ageannee
  
  threshCS <- data.frame()
  sCS <- seq(3,10,1)
  
  cat("\n\nTest breakpoint ages\n\n")
  
  for (i in sCS) {
    subset$age.2 <- subset$ageannee
    print(paste('Fitting threshold model with breakpoint age', i))
    subset$age.2 <- ifelse(subset$age.2 < i, 0, subset$age.2 - i)
    
    mod.threshCS <- lmer(trait ~ age.2*pop + age.2*sexe + age.2*cohort_class 
                         + agelast + masse + datejulienne + dureemin + couleur_HAHL
                         + (1|numind) + (1|annee:pop), 
                         data=subset, control=lmerControl("bobyqa"), REML=F)
    
    k <- AICc(mod.threshCS)
    
    threshCS <- rbind(threshCS, k)
  }
  threshCS <- cbind(sCS, threshCS)
  colnames(threshCS) <- c("age","AICc")
  par(mfrow=c(1,1))
  plot(AICc ~ age, threshCS,
       xlab="age", ylab="AICc")
  }
  
  print(paste('BEST BREAKPOINT AGE IS', threshCS$age[which.min(threshCS$AICc)]))
  
  {subset$age.2 <- subset$ageannee
    subset$age.2 <- ifelse(subset$age.2 <  threshCS$age[which.min(threshCS$AICc)], 0, subset$age.2 -  threshCS$age[which.min(threshCS$AICc)])
    
    mod.thresholdCS <- lmer(trait ~ age.2*pop + age.2*sexe + age.2*cohort_class 
                            + agelast + masse + datejulienne + dureemin + couleur_HAHL
                            + (1|numind) + (1|annee:pop), 
                            data=subset, control=lmerControl("bobyqa"), REML=F)}
  
  cat("\n\n1-slope Threshold model (constant/slope) is running:\n")
  d.thresholdCS <- dredge(mod.thresholdCS, evaluate=T, rank="AICc", trace=2, fixed=~age.2)
  d1.thresholdCS <- subset(d.thresholdCS, delta <=2, recalc.weights = FALSE)
  d2.thresholdCS <- subset(d1.thresholdCS, which(df==min(df)), recalc.weights = FALSE)
  d3.thresholdCS <- subset(d2.thresholdCS, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # threshold-SC: slope - constant ----
  {subset$age.1 <- subset$ageannee
  
  threshSC <- data.frame()
  sSC <- seq(3,10,1)
  
  cat("\n\nTest breakpoint ages\n\n")
  
  for (i in sSC) {
    subset$age.1 <- subset$ageannee
    print(paste('Fitting threshold model with breakpoint age', i))
    subset$age.1 <- ifelse(subset$age.1 < i, subset$age.1 - i, 0)
    
    mod.threshSC <- lmer(trait ~ age.1*pop + age.1*sexe + age.1*cohort_class 
                         + agelast + masse + datejulienne + dureemin + couleur_HAHL
                         + (1|numind) + (1|annee:pop), 
                         data=subset, control=lmerControl("bobyqa"), REML=F)
    
    k <- AICc(mod.threshSC)
    
    threshSC <- rbind(threshSC, k)
  }
  threshSC <- cbind(sSC, threshSC)
  colnames(threshSC) <- c("age","AICc")
  par(mfrow=c(1,1))
  plot(AICc ~ age, threshSC,
       xlab="age", ylab="AICc")
  }
  
  print(paste('BEST BREAKPOINT AGE IS', threshSC$age[which.min(threshSC$AICc)]))
  
  {subset$age.1 <- subset$ageannee
    subset$age.1 <- ifelse(subset$age.1 <  threshSC$age[which.min(threshSC$AICc)], subset$age.1 -  threshSC$age[which.min(threshSC$AICc)], 0)
    
    mod.thresholdSC <- lmer(trait ~ age.1*pop + age.1*sexe + age.1*cohort_class 
                            + agelast + masse + datejulienne + dureemin + couleur_HAHL
                            + (1|numind) + (1|annee:pop), 
                            data=subset, control=lmerControl("bobyqa"), REML=F)}
  
  cat("\n\n1-slope Threshold model (slope/constant) is running:\n")
  d.thresholdSC <- dredge(mod.thresholdSC, evaluate=T, rank="AICc", trace=2, fixed=~age.1)
  d1.thresholdSC <- subset(d.thresholdSC, delta <=2, recalc.weights = FALSE)
  d2.thresholdSC <- subset(d1.thresholdSC, which(df==min(df)), recalc.weights = FALSE)
  d3.thresholdSC <- subset(d2.thresholdSC, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # threshold-SS: slope - slope ----
  {subset$age.1 <- subset$ageannee
  subset$age.2 <- subset$ageannee
  
  threshSS <- data.frame()
  sSS <- seq(3,10,1)
  
  cat("\n\nTest breakpoint ages\n\n")
  
  for (i in sSS) {
    subset$age.1 <- subset$age.2 <- subset$ageannee
    print(paste('Fitting threshold model with breakpoint age', i))
    subset$age.1 [subset$age.1 > i] <- i
    subset$age.2 [subset$age.2 <= i] <- 0
    subset$age.2 [subset$age.2 > i] <- subset$age.2 [subset$age.2 > i] - i
    
    mod.threshSS <- lmer(trait ~ age.1*pop + age.1*sexe + age.1*cohort_class + age.2*pop + age.2*sexe + age.2*cohort_class 
                         + agelast + masse + datejulienne + dureemin + couleur_HAHL
                         + (1|numind) + (1|annee:pop), 
                         data=subset, control=lmerControl("bobyqa"), REML=F)
    
    k <- AICc(mod.threshSS)
    
    threshSS <- rbind(threshSS, k)
  }
  threshSS <- cbind(sSS, threshSS)
  colnames(threshSS) <- c("age","AICc")
  par(mfrow=c(1,1))
  plot(AICc ~ age, threshSS,
       xlab="age", ylab="AICc")
  }
  
  print(paste('BEST BREAKPOINT AGE IS', threshSS$age[which.min(threshSS$AICc)]))
  
  {subset$age.1 <- subset$ageannee
    subset$age.2 <- subset$ageannee
    
    subset$age.1 <- subset$age.2 <- subset$ageannee
    subset$age.1 [subset$age.1 > threshSS$age[which.min(threshSS$AICc)]] <- threshSS$age[which.min(threshSS$AICc)]
    subset$age.2 [subset$age.2 <= threshSS$age[which.min(threshSS$AICc)]] <- 0
    subset$age.2 [subset$age.2 > threshSS$age[which.min(threshSS$AICc)]] <- subset$age.2 [subset$age.2 > threshSS$age[which.min(threshSS$AICc)]] - threshSS$age[which.min(threshSS$AICc)]
    
    mod.thresholdSS <- lmer(trait ~ age.1*pop + age.1*sexe + age.1*cohort_class + age.2*pop + age.2*sexe + age.2*cohort_class 
                            + agelast + masse + datejulienne + dureemin + couleur_HAHL
                            + (1|numind) + (1|annee:pop), 
                            data=subset, control=lmerControl("bobyqa"), REML=F)}
  
  cat("\n\n2-slopes Threshold model (slope/slope) is running:\n")
  d.thresholdSS <- dredge(mod.thresholdSS, evaluate=T, rank="AICc", trace=2, fixed=~age.1 + age.2)
  d1.thresholdSS <- subset(d.thresholdSS, delta <=2, recalc.weights = FALSE)
  d2.thresholdSS <- subset(d1.thresholdSS, which(df==min(df)), recalc.weights = FALSE)
  d3.thresholdSS <- subset(d2.thresholdSS, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  
  cat("CONSTANT MODEL\n\n")
  print(data.frame(d1.constant))
  cat("\n") 
  print(data.frame(d3.constant))
  cat("\n\n\n")
  cat("LINEAR MODEL\n\n")
  print(data.frame(d1.linear))
  cat("\n")
  print(data.frame(d3.linear))
  cat("\n\n\n")
  cat("QUADRATIC MODEL\n\n")
  print(data.frame(d1.quad))
  cat("\n")
  print(data.frame(d3.quad))
  cat("\n\n\n")
  cat("1-SLOPE THRESHOLD MODEL; Constant/slope", "\n", paste("The retained threshold is", threshCS$age[which.min(threshCS$AICc)], "years old\n\n"))
  print(data.frame(d1.thresholdCS))
  cat("\n")
  print(data.frame(d3.thresholdCS))
  cat("\n\n\n")
  cat("1-SLOPE THRESHOLD MODEL; Slope/constant", "\n", paste("The retained threshold is", threshSC$age[which.min(threshSC$AICc)], "years old\n\n"))
  print(data.frame(d1.thresholdSC))
  cat("\n")
  print(data.frame(d3.thresholdSC))
  cat("\n\n\n")
  cat("2-SLOPES THRESHOLD MODEL", "\n", paste("The retained threshold is", threshSS$age[which.min(threshSS$AICc)], "years old\n\n"))
  print(data.frame(d1.thresholdSS))
  cat("\n")
  print(data.frame(d3.thresholdSS))
  cat("\n\n\n")
  
  recap <- data.frame("Model type"=c("Constant",    
                                     "Linear",
                                     "Quadratic",
                                     "Threshold-CS",
                                     "Threshold-SC",
                                     "Threshold-SS"),
                      "Threshold"=c("NA",
                                    "NA",
                                    "NA",
                                    threshCS$age[which.min(threshCS$AICc)],
                                    threshSC$age[which.min(threshSC$AICc)],
                                    threshSS$age[which.min(threshSS$AICc)]),
                      "AICc"=c(d3.constant$AICc,
                               d3.linear$AICc,
                               d3.quad$AICc, 
                               d3.thresholdCS$AICc+2,
                               d3.thresholdSC$AICc+2,
                               d3.thresholdSS$AICc+2),
                      "DF"=c(d3.constant$df,
                             d3.linear$df,
                             d3.quad$df, 
                             d3.thresholdCS$df+1,
                             d3.thresholdSC$df+1,
                             d3.thresholdSS$df+1))
  
  
  recap$delta <- recap$AICc-min(recap$AICc)
  recap$weights <- exp(-0.5*recap$delta)/sum(exp(-0.5*recap$delta))
  cat("\n Summary of the model AICcs and DFs according to model type\n\n")
  print(recap)
  cat("\n\n")
}

# IMMUNOSENESCENCE: PARASITISM ----
# Function to determine the ageing trajectory of parasitic traits
immuno_trajectory_parasitism <- function(subset, trait) {

  options(na.action="na.fail")
  
  # constant ----
  mod.constant <- lmer(log(trait+1) ~ pop + sexe + cohort_class + 
                       + agelast + masse + datejulienne
                       + (1|numind) + (1|annee:pop),
                       data=subset, control=lmerControl("bobyqa"), REML=F)
  
  cat("\n\nConstant model is running:\n")
  d.constant <- dredge(mod.constant, evaluate=T, rank="AICc", trace=2)
  d1.constant <- subset(d.constant, delta <=2, recalc.weights = FALSE)
  d2.constant <- subset(d1.constant, which(df==min(df)), recalc.weights = FALSE)
  d3.constant <- subset(d2.constant, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # linear ----
  mod.linear <- lmer(log(trait+1) ~ ageannee*pop + ageannee*sexe + ageannee*cohort_class + 
                       + agelast + masse + datejulienne
                     + (1|numind) + (1|annee:pop),
                     data=subset, control=lmerControl("bobyqa"), REML=F)
  
  cat("\n\nLinear model is running:\n")
  d.linear <- dredge(mod.linear, evaluate=T, rank="AICc", trace=2, fixed=~ageannee)
  d1.linear <- subset(d.linear, delta <=2, recalc.weights = FALSE)
  d2.linear <- subset(d1.linear, which(df==min(df)), recalc.weights = FALSE)
  d3.linear <- subset(d2.linear, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # quadratic ----
  mod.quad <- lmer(log(trait+1) ~ ageannee*pop + ageannee*sexe + ageannee*cohort_class 
                   + I(ageannee^2)*pop + I(ageannee^2)*sexe + I(ageannee^2)*cohort_class 
                   + agelast + masse + datejulienne
                   + (1|numind) + (1|annee:pop),
                   data=subset, control=lmerControl("bobyqa"), REML=F)
  
  cat("\n\nQuadratic model is running:\n")
  d.quad <- dredge(mod.quad, evaluate=T, rank="AICc", trace=2, fixed=~I(ageannee^2))
  d1.quad <- subset(d.quad, delta <=2, recalc.weights = FALSE)
  d2.quad <- subset(d1.quad, which(df==min(df)), recalc.weights = FALSE)
  d3.quad <- subset(d2.quad, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # threshold-CS: constant - slope ----
  {subset$age.2 <- subset$ageannee
  
  threshCS <- data.frame()
  sCS <- seq(3,10,1)
  
  cat("\n\nTest breakpoint ages\n\n")
  
  for (i in sCS) {
    subset$age.2 <- subset$ageannee
    print(paste('Fitting threshold model with breakpoint age', i))
    subset$age.2 <- ifelse(subset$age.2 < i, 0, subset$age.2 - i)
    
    mod.threshCS <- lmer(log(trait+1) ~ age.2*pop + age.2*sexe + age.2*cohort_class 
                         + agelast + masse + datejulienne
                         + (1|numind) + (1|annee:pop), 
                         data=subset, control=lmerControl("bobyqa"), REML=F)
    
    k <- AICc(mod.threshCS)
    
    threshCS <- rbind(threshCS, k)
  }
  threshCS <- cbind(sCS, threshCS)
  colnames(threshCS) <- c("age","AICc")
  par(mfrow=c(1,1))
  plot(AICc ~ age, threshCS,
       xlab="age", ylab="AICc")
  }
  
  print(paste('BEST BREAKPOINT AGE IS', threshCS$age[which.min(threshCS$AICc)]))
  
  {subset$age.2 <- subset$ageannee
    subset$age.2 <- ifelse(subset$age.2 <  threshCS$age[which.min(threshCS$AICc)], 0, subset$age.2 -  threshCS$age[which.min(threshCS$AICc)])
    
    mod.thresholdCS <- lmer(log(trait+1) ~ age.2*pop + age.2*sexe + age.2*cohort_class 
                            + agelast + masse + datejulienne
                            + (1|numind) + (1|annee:pop), 
                            data=subset, control=lmerControl("bobyqa"), REML=F)}
  
  cat("\n\n1-slope Threshold model (constant/slope) is running:\n")
  d.thresholdCS <- dredge(mod.thresholdCS, evaluate=T, rank="AICc", trace=2, fixed=~age.2)
  d1.thresholdCS <- subset(d.thresholdCS, delta <=2, recalc.weights = FALSE)
  d2.thresholdCS <- subset(d1.thresholdCS, which(df==min(df)), recalc.weights = FALSE)
  d3.thresholdCS <- subset(d2.thresholdCS, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # threshold-SC: slope - constant ----
  {subset$age.1 <- subset$ageannee
  
  threshSC <- data.frame()
  sSC <- seq(3,10,1)
  
  cat("\n\nTest breakpoint ages\n\n")
  
  for (i in sSC) {
    subset$age.1 <- subset$ageannee
    print(paste('Fitting threshold model with breakpoint age', i))
    subset$age.1 <- ifelse(subset$age.1 < i, subset$age.1 - i, 0)
    
    mod.threshSC <- lmer(log(trait+1) ~ age.1*pop + age.1*sexe + age.1*cohort_class 
                         + agelast + masse + datejulienne
                         + (1|numind) + (1|annee:pop), 
                         data=subset, control=lmerControl("bobyqa"), REML=F)
    
    k <- AICc(mod.threshSC)
    
    threshSC <- rbind(threshSC, k)
  }
  threshSC <- cbind(sSC, threshSC)
  colnames(threshSC) <- c("age","AICc")
  par(mfrow=c(1,1))
  plot(AICc ~ age, threshSC,
       xlab="age", ylab="AICc")
  }
  
  print(paste('BEST BREAKPOINT AGE IS', threshSC$age[which.min(threshSC$AICc)]))
  
  {subset$age.1 <- subset$ageannee
    subset$age.1 <- ifelse(subset$age.1 <  threshSC$age[which.min(threshSC$AICc)], subset$age.1 -  threshSC$age[which.min(threshSC$AICc)], 0)
    
    mod.thresholdSC <- lmer(log(trait+1) ~ age.1*pop + age.1*sexe + age.1*cohort_class 
                            + agelast + masse + datejulienne
                            + (1|numind) + (1|annee:pop), 
                            data=subset, control=lmerControl("bobyqa"), REML=F)}
  
  cat("\n\n1-slope Threshold model (slope/constant) is running:\n")
  d.thresholdSC <- dredge(mod.thresholdSC, evaluate=T, rank="AICc", trace=2, fixed=~age.1)
  d1.thresholdSC <- subset(d.thresholdSC, delta <=2, recalc.weights = FALSE)
  d2.thresholdSC <- subset(d1.thresholdSC, which(df==min(df)), recalc.weights = FALSE)
  d3.thresholdSC <- subset(d2.thresholdSC, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  # threshold-SS: slope - slope ----
  {subset$age.1 <- subset$ageannee
  subset$age.2 <- subset$ageannee
  
  threshSS <- data.frame()
  sSS <- seq(3,10,1)
  
  cat("\n\nTest breakpoint ages\n\n")
  
  for (i in sSS) {
    subset$age.1 <- subset$age.2 <- subset$ageannee
    print(paste('Fitting threshold model with breakpoint age', i))
    subset$age.1 [subset$age.1 > i] <- i
    subset$age.2 [subset$age.2 <= i] <- 0
    subset$age.2 [subset$age.2 > i] <- subset$age.2 [subset$age.2 > i] - i
    
    mod.threshSS <- lmer(log(trait+1) ~ age.1*pop + age.1*sexe + age.1*cohort_class + age.2*pop + age.2*sexe + age.2*cohort_class 
                         + agelast + masse + datejulienne
                         + (1|numind) + (1|annee:pop), 
                         data=subset, control=lmerControl("bobyqa"), REML=F)
    
    k <- AICc(mod.threshSS)
    
    threshSS <- rbind(threshSS, k)
  }
  threshSS <- cbind(sSS, threshSS)
  colnames(threshSS) <- c("age","AICc")
  par(mfrow=c(1,1))
  plot(AICc ~ age, threshSS,
       xlab="age", ylab="AICc")
  }
  
  print(paste('BEST BREAKPOINT AGE IS', threshSS$age[which.min(threshSS$AICc)]))
  
  {subset$age.1 <- subset$ageannee
    subset$age.2 <- subset$ageannee
    
    subset$age.1 <- subset$age.2 <- subset$ageannee
    subset$age.1 [subset$age.1 > threshSS$age[which.min(threshSS$AICc)]] <- threshSS$age[which.min(threshSS$AICc)]
    subset$age.2 [subset$age.2 <= threshSS$age[which.min(threshSS$AICc)]] <- 0
    subset$age.2 [subset$age.2 > threshSS$age[which.min(threshSS$AICc)]] <- subset$age.2 [subset$age.2 > threshSS$age[which.min(threshSS$AICc)]] - threshSS$age[which.min(threshSS$AICc)]
    
    mod.thresholdSS <- lmer(log(trait+1) ~ age.1*pop + age.1*sexe + age.1*cohort_class + age.2*pop + age.2*sexe + age.2*cohort_class 
                            + agelast + masse + datejulienne
                            + (1|numind) + (1|annee:pop), 
                            data=subset, control=lmerControl("bobyqa"), REML=F)}
  
  cat("\n\n2-slopes Threshold model (slope/slope) is running:\n")
  d.thresholdSS <- dredge(mod.thresholdSS, evaluate=T, rank="AICc", trace=2, fixed=~age.1 + age.2)
  d1.thresholdSS <- subset(d.thresholdSS, delta <=2, recalc.weights = FALSE)
  d2.thresholdSS <- subset(d1.thresholdSS, which(df==min(df)), recalc.weights = FALSE)
  d3.thresholdSS <- subset(d2.thresholdSS, which(AICc==min(AICc)), recalc.weights = FALSE)
  
  
  cat("CONSTANT MODEL\n\n")
  print(data.frame(d1.constant))
  cat("\n")
  print(data.frame(d3.constant))
  cat("\n\n\n")
  cat("LINEAR MODEL\n\n")
  print(data.frame(d1.linear))
  cat("\n")
  print(data.frame(d3.linear))
  cat("\n\n\n")
  cat("QUADRATIC MODEL\n\n")
  print(data.frame(d1.quad))
  cat("\n")
  print(data.frame(d3.quad))
  cat("\n\n\n")
  cat("1-SLOPE THRESHOLD MODEL; Constant/slope", "\n", paste("The retained threshold is", threshCS$age[which.min(threshCS$AICc)], "years old\n\n"))
  print(data.frame(d1.thresholdCS))
  cat("\n")
  print(data.frame(d3.thresholdCS))
  cat("\n\n\n")
  cat("1-SLOPE THRESHOLD MODEL; Slope/constant", "\n", paste("The retained threshold is", threshSC$age[which.min(threshSC$AICc)], "years old\n\n"))
  print(data.frame(d1.thresholdSC))
  cat("\n")
  print(data.frame(d3.thresholdSC))
  cat("\n\n\n")
  cat("2-SLOPES THRESHOLD MODEL", "\n", paste("The retained threshold is", threshSS$age[which.min(threshSS$AICc)], "years old\n\n"))
  print(data.frame(d1.thresholdSS))
  cat("\n")
  print(data.frame(d3.thresholdSS))
  cat("\n\n\n")
  
  recap <- data.frame("Model type"=c("Constant",    
                                     "Linear",
                                     "Quadratic",
                                     "Threshold-CS",
                                     "Threshold-SC",
                                     "Threshold-SS"),
                      "Threshold"=c("NA",
                                    "NA",
                                    "NA",
                                    threshCS$age[which.min(threshCS$AICc)],
                                    threshSC$age[which.min(threshSC$AICc)],
                                    threshSS$age[which.min(threshSS$AICc)]),
                      "AICc"=c(d3.constant$AICc,
                               d3.linear$AICc,
                               d3.quad$AICc, 
                               d3.thresholdCS$AICc+2,
                               d3.thresholdSC$AICc+2,
                               d3.thresholdSS$AICc+2),
                      "DF"=c(d3.constant$df,
                             d3.linear$df,
                             d3.quad$df, 
                             d3.thresholdCS$df+1,
                             d3.thresholdSC$df+1,
                             d3.thresholdSS$df+1))
  
  
  recap$delta <- recap$AICc-min(recap$AICc)
  recap$weights <- exp(-0.5*recap$delta)/sum(exp(-0.5*recap$delta))
  cat("\n Summary of the model AICcs and DFs according to model type\n\n")
  print(recap)
  cat("\n\n")
}

