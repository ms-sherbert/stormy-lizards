#=== Analysis of trends in body size and condition before and after inundation ===#
# Written by S Herbert
# For R version 4.3.1
# Last tested: 09 June 2026

#=== Preamble (Dependencies & set working directory) ===#

rm(list=ls())
setwd("C:/Repositories/stormy-lizards") 
#NB repositories in D: on uni desktop, C: on personal laptop

#Use these bits of code to install latest versions of Matrix and lme4 if lmer models are throwing errors
#options(repos = c(CRAN = "https://cloud.r-project.org"))
#utils::instALL.packages("Matrix")
#utils::instALL.packages("lme4")

library(tidyverse)
library(MASS)
library(lme4)
library(lmerTest) #displays p values in lmer summary tables
library(DHARMa) #assess fit of many model types, including lmer models
library(cAIC4)
library(rcompanion)
library(gridExtra)
library(blme)
source("Code/Model-selection-functions.R")
#library(lmtest) #Previously used for Breusch-Pagan tests; using DHARMa instead for model fit assessment

#=== Read in required files ===#

captures<-read.csv("Data/CC-grid-captures_S1-S8.csv")

#=== Compute additional variables and tidy data frame ===#

#---Replace season number with the midpoint dates for accurate time effect---#
season.midpoint<- c("5/12/2017","23/01/2018","27/03/2018","7/12/2018",
                    "25/10/2019","2/03/2020","1/12/2020","25/03/2021")
SDates <- strptime(season.midpoint,format="%d/%m/%Y")
SDates <- decimal_date(SDates)

captures$Season[captures$Season == 1] <- SDates[1]
captures$Season[captures$Season == 2] <- SDates[2]
captures$Season[captures$Season == 3] <- SDates[3]
captures$Season[captures$Season == 4] <- SDates[4]
captures$Season[captures$Season == 5] <- SDates[5]
captures$Season[captures$Season == 6] <- SDates[6]
captures$Season[captures$Season == 7] <- SDates[7]
captures$Season[captures$Season == 8] <- SDates[8]

#---Compute other variables---#

captures$Site <- substr(captures$Grid,1,2)
captures$inundated <- rep(0, times = length(captures$uID))
captures$inundated[captures$Grid != "MP1" | captures$Grid != "WP3"] <- 0
captures$inundated[captures$Grid == "MP1" | captures$Grid == "WP3"] <- 1
captures$BA <- rep(0, times = length(captures$uID))
captures$BA[captures$Season < SDates[7]] <- 0
captures$BA[captures$Season >= SDates[7]] <- 1
captures$Season.centred <- captures$Season - SDates[7]
#captures$Wt_lizard_g[captures$Wt_lizard_g == 0] <- NA #find and change any data entry mistakes to NA
captures$SVL_mm[captures$SVL_mm == 0] <- NA #find and change any data entry mistakes to NA
captures <- subset(captures,!is.na(SVL_mm)) #remove missing values of SVL

#=== Create new data frames for OP and WM ===#

OPdf <- subset(captures, Species == "op")
WMdf <- subset(captures, Species == "wm")

#Export for graphing later

OPWMdf <- subset(captures, Species == "op" | Species == "wm")
write.csv(OPWMdf,"Outputs/SVLs.csv")


#=== Examine underlying distribution of SVL ===#

#O. polychroma
hist(OPdf$SVL_mm) #left skewed (i.e. longer tail to left)
qqnorm(OPdf$SVL_mm)
qqline(OPdf$SVL_mm) #yep, definitely skewed although fit through centre ok
ks.test(OPdf$SVL_mm,"pnorm", alternative = "two.sided",exact=TRUE)
# D = 1, p-value = 5.218e-15 - very not normal

SVL_transOP <- transformTukey(OPdf$SVL_mm,plot=FALSE) #try Tukey's Ladder of Powers transformation
hist(SVL_transOP) # this looks much more normal

#W. maculata
hist(WMdf$SVL_mm) #very left skewed (i.e. longer tail to left)
qqnorm(WMdf$SVL_mm)
qqline(WMdf$SVL_mm) #yep, definitely skewed
ks.test(WMdf$SVL_mm,"pnorm", alternative = "two.sided",exact=TRUE)
# D = 1, p-value = 5.107e-15 - very not normal

SVL_transWM <- transformTukey(WMdf$SVL_mm,plot=FALSE) #try Tukey's Ladder of Powers transformation
hist(SVL_transWM) #better, but still not great

#check for covariance among numeric fixed covariates
cor(OPdf[,c(12:13)],method="spearman",use="complete.obs")  # rho = -0.017
cor(WMdf[,c(12:13)],method="spearman",use="complete.obs")  # rho = 0.006

#=== SVL models ===#

#--- OP SVL ---#
M1a.svlOP <- lmer(SVL_transOP ~ BA*inundated + (Season.centred|Grid), data = OPdf) #fit is singular, grid intercept explains zero variance
M1.svlOP <- lmer(SVL_transOP ~ BA*inundated + (1|Grid), data = OPdf)
M2.svlOP <- lm(SVL_transOP ~ BA*inundated, data = OPdf)

anova(M1a.svlOP,M1.svlOP) #Inclusion of random term doesn't improve fit at all
anova(M1a.svlOP,M2.svlOP) #Inclusion of random term doesn't improve fit at all

summary(M1a.svlOP)

par(mfrow=c(2,2))
plot(M2.svlOP) #assess model fit
par(mfrow=c(1,1))

write.csv(coef(summary(M1a.svlOP)), "Outputs/FixedEff_M1a_svlOP.csv")


#--- WM SVL ---#
M1.svlWM <- lmer(SVL_transWM ~ BA*inundated + (Season.centred|Grid), data = WMdf) #potential convergence issues, but allFit checks seem to be ok for all optimisers
M2.svlWM <- lmer(SVL_transWM ~ BA*inundated + (1|Grid), data = WMdf) #potential convergence issues, but allFit checks seem to be ok for all optimisers

#Try decreasing stopping tolerances
strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
if (all(M1.svlWM@optinfo$optimizer=="nloptwrap")) {
  fm1.tol <- update(M1.svlWM, control=strict_tol)
}

strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
if (all(M2.svlWM@optinfo$optimizer=="nloptwrap")) {
  fm1.tol <- update(M2.svlWM, control=strict_tol)
}

## try recomputing gradient and Hessian with Richardson extrapolation

# Model 1
devfun <- update(M1.svlWM, devFunOnly=TRUE)
if (isLMM(M1.svlWM)) {
  pars <- getME(M1.svlWM,"theta")
} else {
  ## GLMM: requires both random and fixed parameters
  pars <- getME(M1.svlWM, c("theta","fixef"))
}
if (require("numDeriv")) {
  cat("hess:\n"); print(hess <- hessian(devfun, unlist(pars)))
  cat("grad:\n"); print(grad <- grad(devfun, unlist(pars)))
  cat("scaled gradient:\n")
  print(scgrad <- solve(chol(hess), grad))
}
## compare with internal calculations:
M1.svlWM@optinfo$derivs

## compute reciprocal condition number of Hessian
H <- M1.svlWM@optinfo$derivs$Hessian
Matrix::rcond(H) #0.4776378

# Model 2

devfun <- update(M2.svlWM, devFunOnly=TRUE)
if (isLMM(M2.svlWM)) {
  pars <- getME(M2.svlWM,"theta")
} else {
  ## GLMM: requires both random and fixed parameters
  pars <- getME(M2.svlWM, c("theta","fixef"))
}
if (require("numDeriv")) {
  cat("hess:\n"); print(hess <- hessian(devfun, unlist(pars)))
  cat("grad:\n"); print(grad <- grad(devfun, unlist(pars)))
  cat("scaled gradient:\n")
  print(scgrad <- solve(chol(hess), grad))
}
## compare with internal calculations:
M2.svlWM@optinfo$derivs

## compute reciprocal condition number of Hessian
H <- M2.svlWM@optinfo$derivs$Hessian
Matrix::rcond(H) #0.4636483

## 4. restart the fit from the original value (or
## a slightly perturbed value):
# Lines 183-193 are throwing errors as at last test run; have commented out
# but left in script to document previous process

#M1.svlWM.restart <- update(M1.svlWM, start=pars) 
#set.seed(101)
#pars_x <- runif(length(pars),pars/1.01,pars*1.01)
#M1.svlWM.restart2 <- update(M1.svlWM, start=pars_x, control=strict_tol)

#M2.svlWM.restart <- update(M2.svlWM, start=pars)
#set.seed(101)
#pars_x <- runif(length(pars),pars/1.01,pars*1.01)
#M2.svlWM.restart2 <- update(M2.svlWM, start=pars_x, control=strict_tol)

#---Try all available optimizers

M1.svlWM.fitcheck <- allFit(M1.svlWM)
M2.svlWM.fitcheck <- allFit(M2.svlWM)

summary(M1.svlWM.fitcheck)
summary(M2.svlWM.fitcheck) #Actually, these summaries seem to indicate that convergence is ok with all optimizers
# Will go ahead and assume that convergence warnings don't mean the fit is incorrect in this case
# see: https://rdrr.io/cran/lme4/man/convergence.html

summary(M1.svlWM)

# This is weird; the DHARMa diagnostic plots from the June 2025 test of this 
# script suggest a lack of fit
# I don't recall seeing this during the original analysis
simulationOutput.M1.svlWM <- simulateResiduals(fittedModel = M1.svlWM, plot = F) #assess fit of model with lowest AIC
plotQQunif(simulationOutput.M1.svlWM) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M1.svlWM) # right plot in plot.DHARMa()

write.csv(coef(summary(M1.svlWM)), "Outputs/FixedEff_M1_svlWM.csv")