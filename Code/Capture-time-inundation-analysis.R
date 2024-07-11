#=== Analysis of trends in CPUE and Mt+1 before and after inundation ===#
# Written by S Herbert
# For R version 4.3.1
# Last tested: 11/07/2024

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
library(gridExtra)
source("Code/Model-selection-functions.R")
#library(lmtest) #Previously used for Breusch-Pagan tests; using DHARMa instead for model fit assessment

#=== Read in required files ===#

CPUE<-read.csv("Outputs/CPUElizards.csv")
Mt1<-read.csv("Outputs/Mtlizards.csv")
disturbance<-read.csv("Data/Site_inundation_and_trap_replacement.csv")
temperatures<-read.csv("Data/CMR_check_schedule.csv")

#=== Compute additional variables and tidy data frame ===#

names(CPUE)[names(CPUE)=="n.Nsessions"] <- "CPUE"
names(Mt1)[names(Mt1)=="n.Nsessions"] <- "Mt1"

Site <- substr(CPUE$Grid,1,2)
DF <- cbind(CPUE[,2:4],CPUE[,11],Mt1[,11],Site)

names(DF)[names(DF)=="CPUE[, 11]"] <- "CPUE"
names(DF)[names(DF)=="Mt1[, 11]"] <- "Mt1"

#=== Prepare trap replacement by time period and site vector ===#

Replacements <- disturbance %>%
                group_by(Site.number) %>%
                summarise("ACO.replaced" = sum(ACO.replaced),
                          "Pitfall.replaced" = sum(Pitfall.replaced)
                          )

#write.csv(Replacements, "Outputs/Replaced-stations.csv")

#Create new empty covariate columns and append to CPUE and Mt+1 dataframes

ACO.repl <- rep(0, times = 112)
Pit.repl <- rep(0, times = 112)
BA <- rep(0, times = 112)
Inundated <- rep(0, times = 112)

DF <- cbind(DF,ACO.repl,Pit.repl,BA,Inundated)

#Update relevant cells of dataframe with values from Replacements
DF$ACO.repl[DF$Grid == "MP1" & DF$Season == 7] <- 6
DF$Pit.repl[DF$Grid == "MP1" & DF$Season == 7] <- 17
DF$ACO.repl[DF$Grid == "WP2" & DF$Season == 7] <- 1
DF$Pit.repl[DF$Grid == "WP2" & DF$Season == 7] <- 1
DF$ACO.repl[DF$Grid == "WP3" & DF$Season == 7] <- 9
DF$Pit.repl[DF$Grid == "WP3" & DF$Season == 7] <- 10

DF$Inundated[DF$Grid == "MP1" | DF$Grid == "WP3"] <- 1
DF$BA[DF$Season >= 7] <- 1

#--- Replace season number with the midpoint dates for accurate time effect ---#
season.midpoint<- c("5/12/2017","23/01/2018","27/03/2018","7/12/2018",
                    "25/10/2019","2/03/2020","1/12/2020","25/03/2021")
SDates <- strptime(season.midpoint,format="%d/%m/%Y")
SDates <- decimal_date(SDates)

DF$Season[DF$Season == 1] <- SDates[1]
DF$Season[DF$Season == 2] <- SDates[2]
DF$Season[DF$Season == 3] <- SDates[3]
DF$Season[DF$Season == 4] <- SDates[4]
DF$Season[DF$Season == 5] <- SDates[5]
DF$Season[DF$Season == 6] <- SDates[6]
DF$Season[DF$Season == 7] <- SDates[7]
DF$Season[DF$Season == 8] <- SDates[8]
 
#=== Calculate average daily maximum and minimum temperatures for each season ===#

TMAX <- temperatures %>%
         group_by(Season) %>%
         summarise("MP1" = sum(Tmax*MP1,na.rm=TRUE)/sum(MP1,na.rm=TRUE),
                   "MP2" = sum(Tmax*MP2,na.rm=TRUE)/sum(MP2,na.rm=TRUE),
                   "RE1" = sum(Tmax*RE1,na.rm=TRUE)/sum(RE1,na.rm=TRUE),
                   "RE2" = sum(Tmax*RE2,na.rm=TRUE)/sum(RE2,na.rm=TRUE),
                   "WP2" = sum(Tmax*WP2,na.rm=TRUE)/sum(WP2,na.rm=TRUE),
                   "WP3" = sum(Tmax*WP3,na.rm=TRUE)/sum(WP3,na.rm=TRUE))

TMIN <- temperatures %>%
  group_by(Season) %>%
  summarise("MP1" = sum(Tmin*MP1,na.rm=TRUE)/sum(MP1,na.rm=TRUE),
            "MP2" = sum(Tmin*MP2,na.rm=TRUE)/sum(MP2,na.rm=TRUE),
            "RE1" = sum(Tmin*RE1,na.rm=TRUE)/sum(RE1,na.rm=TRUE),
            "RE2" = sum(Tmin*RE2,na.rm=TRUE)/sum(RE2,na.rm=TRUE),
            "WP2" = sum(Tmin*WP2,na.rm=TRUE)/sum(WP2,na.rm=TRUE),
            "WP3" = sum(Tmin*WP3,na.rm=TRUE)/sum(WP3,na.rm=TRUE))       

T.max <- c(TMAX$RE1,TMAX$RE2,
          TMAX$MP1,TMAX$MP2,TMAX$RE1,TMAX$RE2,TMAX$WP2,TMAX$WP3,
          TMAX$MP1,TMAX$MP2,TMAX$RE1,TMAX$RE2,TMAX$WP2,TMAX$WP3)

T.min <- c(TMIN$RE1,TMIN$RE2,
          TMIN$MP1,TMIN$MP2,TMIN$RE1,TMIN$RE2,TMIN$WP2,TMIN$WP3,
          TMIN$MP1,TMIN$MP2,TMIN$RE1,TMIN$RE2,TMIN$WP2,TMIN$WP3)

DF <- cbind(DF,T.max,T.min) 

DF$Season.centred <- DF$Season - 2020.915

#Create subsets of the data by species
OP <- subset(DF, Species == "op")
WM <- subset(DF, Species == "wm")

#=== Examine underlying distributions of CPUE and Mt1 ===#

par(mfrow=c(2,2))
hist(DF$CPUE) #looks negative binomial
hist(sqrt(DF$CPUE))
hist(log(DF$CPUE))
qqnorm(sqrt(DF$CPUE))
qqline(sqrt(DF$CPUE)) #Approx normal
par(mfrow=c(1,1))

ks.test(sqrt(DF$CPUE),"pnorm", alternative = "two.sided", exact = TRUE)
# D = 0.5348, p-value = 2.776e-15 
# KS test says very non-normal, but we know this test can be over-sensitive

par(mfrow=c(2,2))
hist(DF$Mt1) #looks negative binomial
hist(sqrt(DF$Mt1))
hist(log(DF$Mt1))
qqnorm(sqrt(DF$Mt1))
qqline(sqrt(DF$Mt1)) #Approx normal
par(mfrow=c(1,1))

ks.test(sqrt(DF$Mt1),"pnorm", alternative = "two.sided", exact = TRUE)
#D = 0.5, p-value = 2.776e-15

cor(OP[,c(3,7:12)],method="spearman") #check for covariance among potential fixed factors
# Correlations between T.max and T.min, Season and BA, Pit.rep and ACO.rep >|0.7|, 
# So don't use both in same model
# Consistently more pitfalls needed to be replaced than ACOs among grids, so have chosen this variable for disturbance

#=== Create models for data and assess fits ===#
#https://stats.oarc.ucla.edu/r/dae/negative-binomial-regression/
#https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

#--- Oligosoma polychroma CPUE ---#

#No warning messages from any of the following eight models
#fit is singular if spatial random component is Site/Grid

M1.cOP <- lmer(log(CPUE+0.01) ~ Season.centred*BA*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = OP)
M2.cOP <- lmer(log(CPUE+0.01) ~ Season.centred*BA*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = OP)
M3.cOP <- lmer(log(CPUE+0.01) ~ Season.centred*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = OP)
M4.cOP <- lmer(log(CPUE+0.01) ~ Season.centred*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = OP)

models.cOP <- list(M1.cOP,M2.cOP,M3.cOP,M4.cOP)
Avg.models.cOP <- modelAvg(models.cOP, opt = TRUE)

cAIC.M1.cOP <- cAIC(M1.cOP)
cAIC.M2.cOP <- cAIC(M2.cOP)
cAIC.M3.cOP <- cAIC(M3.cOP)
cAIC.M4.cOP <- cAIC(M4.cOP)

AIC.cOP <- unlist(c(cAIC.M1.cOP[5],cAIC.M2.cOP[5],cAIC.M3.cOP[5],cAIC.M4.cOP[5]))
df.cOP <- unlist(c(cAIC.M1.cOP[2],cAIC.M2.cOP[2],cAIC.M3.cOP[2],cAIC.M4.cOP[2]))
DAIC.cOP <- Delta.AIC(AIC.cOP)

Cand.model.cOP <- as.data.frame(cbind(AIC.cOP, DAIC.cOP, df.cOP))
Model.selection.cOP <- Cand.model.cOP %>%
  reframe("Model" = c(1:4),
          "cAIC" = round(AIC.cOP,2), 
          "DcAIC" = round(DAIC.cOP,2),
          "Estimated DF" = round(df.cOP,2)
  ) %>%
  arrange(cAIC)

Model.selection.cOP

summary(M1.cOP) #view summary of model with lowest AIC

simulationOutput.M1.cOP <- simulateResiduals(fittedModel = M1.cOP, plot = F) #assess fit of model with lowest AIC
plotQQunif(simulationOutput.M1.cOP) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M1.cOP) # right plot in plot.DHARMa()

write.csv(coef(summary(M1.cOP)), "Outputs/FixedEff_M1_cOP.csv")
write.csv(Model.selection.cOP,"Outputs/ModSel_cOP.csv")

#--- Woodworthia maculata CPUE ---#

#No warning messages produced by any of these models
M1.cWM <- lmer(log(CPUE+0.01) ~ Season.centred*BA*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = WM)
M2.cWM <- lmer(log(CPUE+0.01) ~ Season.centred*BA*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = WM)
M3.cWM <- lmer(log(CPUE+0.01) ~ Season.centred*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = WM)
M4.cWM <- lmer(log(CPUE+0.01) ~ Season.centred*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = WM)

models.cWM <- list(M1.cWM,M2.cWM,M3.cWM,M4.cWM)
Avg.models.cWM <- modelAvg(models.cWM, opt = TRUE)

cAIC.M1.cWM <- cAIC(M1.cWM)
cAIC.M2.cWM <- cAIC(M2.cWM)
cAIC.M3.cWM <- cAIC(M3.cWM)
cAIC.M4.cWM <- cAIC(M4.cWM)

AIC.cWM <- unlist(c(cAIC.M1.cWM[5],cAIC.M2.cWM[5],cAIC.M3.cWM[5],cAIC.M4.cWM[5]))
df.cWM <- unlist(c(cAIC.M1.cWM[2],cAIC.M2.cWM[2],cAIC.M3.cWM[2],cAIC.M4.cWM[2]))
DAIC.cWM <- Delta.AIC(AIC.cOP)

Cand.model.cWM <- as.data.frame(cbind(AIC.cWM, DAIC.cWM, df.cWM))
Model.selection.cWM <- Cand.model.cWM %>%
  reframe("Model" = c(1:4),
          "cAIC" = round(AIC.cWM,2), 
          "DcAIC" = round(DAIC.cWM,2),
          "Estimated DF" = round(df.cWM,2)
  ) %>%
  arrange(cAIC)

Model.selection.cWM

summary(M3.cWM) #view summary of model with lowest AIC

simulationOutput.M3.cWM <- simulateResiduals(fittedModel = M3.cWM, plot = F) #assess fit of model with lowest AIC
plotQQunif(simulationOutput.M3.cWM) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M3.cWM) # right plot in plot.DHARMa()

write.csv(coef(summary(M3.cWM)), "Outputs/FixedEff_M3_cWM.csv")
write.csv(Model.selection.cWM,"Outputs/ModSel_cWM.csv")

#--- All lizards CPUE ---#

#No warning messages produced by these models
M1.cALL <- lmer(sqrt(CPUE+0.01) ~ Season.centred*BA*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = DF) #fit is singular if spatial random component is Site/Grid
M2.cALL <- lmer(sqrt(CPUE+0.01) ~ Season.centred*BA*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = DF)
M3.cALL <- lmer(sqrt(CPUE+0.01) ~ Season.centred*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = DF) 
M4.cALL <- lmer(sqrt(CPUE+0.01) ~ Season.centred*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = DF)

models.cALL <- list(M1.cALL,M2.cALL,M3.cALL,M4.cALL)
Avg.models.cALL <- modelAvg(models.cALL, opt = TRUE)

cAIC.M1.cALL <- cAIC(M1.cALL)
cAIC.M2.cALL <- cAIC(M2.cALL)
cAIC.M3.cALL <- cAIC(M3.cALL)
cAIC.M4.cALL <- cAIC(M4.cALL)

AIC.cALL <- unlist(c(cAIC.M1.cALL[5],cAIC.M2.cALL[5],cAIC.M3.cALL[5],cAIC.M4.cALL[5]))
df.cALL <- unlist(c(cAIC.M1.cALL[2],cAIC.M2.cALL[2],cAIC.M3.cALL[2],cAIC.M4.cALL[2]))
DAIC.cALL <- Delta.AIC(AIC.cALL)

Cand.model.cALL <- as.data.frame(cbind(AIC.cALL, DAIC.cALL, df.cALL))
Model.selection.cALL <- Cand.model.cALL %>%
  reframe("Model" = c(1:4),
          "cAIC" = round(AIC.cALL,2), 
          "DcAIC" = round(DAIC.cALL,2),
          "Estimated DF" = round(df.cALL,2)
  ) %>%
  arrange(cAIC)

Model.selection.cALL

summary(M3.cALL) #view summary of model with lowest AIC

simulationOutput.M3.cALL <- simulateResiduals(fittedModel = M3.cALL, plot = F) #assess fit of model with lowest AIC
plotQQunif(simulationOutput.M3.cALL) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M3.cALL) # right plot in plot.DHARMa()

write.csv(coef(summary(M3.cALL)), "Outputs/FixedEff_M3_cALL.csv")
write.csv(Model.selection.cALL,"Outputs/ModSel_cALL.csv")

#--- Oligosoma polychroma Mt+1 ---#

M1.mOP <- lmer(log(Mt1+0.01) ~ Season.centred*BA*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = OP) #singular
M2.mOP <- lmer(log(Mt1+0.01) ~ Season.centred*BA*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = OP)
M3.mOP <- lmer(log(Mt1+0.01) ~ Season.centred*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = OP) #singular
M4.mOP <- lmer(log(Mt1+0.01) ~ Season.centred*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = OP) #singular
M5.mOP <- lmer(log(Mt1+0.01) ~ Season.centred*BA*Inundated + T.min + Pit.repl + (1|Grid), data = OP)
M6.mOP <- lmer(log(Mt1+0.01) ~ Season.centred*BA*Inundated + T.max + Pit.repl + (1|Grid), data = OP)
M7.mOP <- lmer(log(Mt1+0.01) ~ Season.centred*Inundated + T.min + Pit.repl + (1|Grid), data = OP)
M8.mOP <- lmer(log(Mt1+0.01) ~ Season.centred*Inundated + T.max + Pit.repl + (1|Grid), data = OP)

#Test whether simplifying random term creates a significant difference 
#in any pair of models with same fixed factor structure
anova(M1.mOP,M5.mOP) #p = 0.9854
anova(M2.mOP,M6.mOP) #p = 0.9943
anova(M3.mOP,M7.mOP) #p = 0.9989
anova(M4.mOP,M8.mOP) #p = 0.9998
#No, so decided to use models with random intercept term (1|Grid) in candidate set
models.mOP <- list(M5.mOP,M6.mOP,M7.mOP,M8.mOP)
Avg.models.mOP <- modelAvg(models.mOP, opt = TRUE)

cAIC.M5.mOP <- cAIC(M5.mOP)
cAIC.M6.mOP <- cAIC(M6.mOP)
cAIC.M7.mOP <- cAIC(M7.mOP)
cAIC.M8.mOP <- cAIC(M8.mOP)

AIC.mOP <- unlist(c(cAIC.M5.mOP[5],cAIC.M6.mOP[5],cAIC.M7.mOP[5],cAIC.M8.mOP[5]))
df.mOP <- unlist(c(cAIC.M5.mOP[2],cAIC.M6.mOP[2],cAIC.M7.mOP[2],cAIC.M8.mOP[2]))
DAIC.mOP <- Delta.AIC(AIC.mOP)

Cand.model.mOP <- as.data.frame(cbind(AIC.mOP, DAIC.mOP, df.mOP))
Model.selection.mOP <- Cand.model.mOP %>%
  reframe("Model" = c(1:4),
          "cAIC" = round(AIC.mOP,2), 
          "DcAIC" = round(DAIC.mOP,2),
          "Estimated DF" = round(df.mOP,2)
  ) %>%
  arrange(cAIC)

Model.selection.mOP

summary(M7.mOP) #view summary of model with lowest AIC (in this set model 3 is model M7)

simulationOutput.M7.mOP <- simulateResiduals(fittedModel = M7.mOP, plot = F) #assess fit of model with lowest AIC
plotQQunif(simulationOutput.M7.mOP) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M7.mOP) # right plot in plot.DHARMa()

write.csv(coef(summary(M7.mOP)), "Outputs/FixedEff_M7_mOP.csv")
write.csv(Model.selection.mOP,"Outputs/ModSel_mOP.csv")

#--- Woodworthia maculata Mt+1 ---#

M1.mWM <- lmer(log(Mt1+0.01) ~ Season.centred*BA*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = WM)
M2.mWM <- lmer(log(Mt1+0.01) ~ Season.centred*BA*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = WM)
M3.mWM <- lmer(log(Mt1+0.01) ~ Season.centred*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = WM)
M4.mWM <- lmer(log(Mt1+0.01) ~ Season.centred*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = WM)

models.mWM <- list(M1.mWM,M2.mWM,M3.mWM,M4.mWM)
Avg.models.mWM <- modelAvg(models.mWM, opt = TRUE)

cAIC.M1.mWM <- cAIC(M1.mWM)
cAIC.M2.mWM <- cAIC(M2.mWM)
cAIC.M3.mWM <- cAIC(M3.mWM)
cAIC.M4.mWM <- cAIC(M4.mWM)

AIC.mWM <- unlist(c(cAIC.M1.mWM[5],cAIC.M2.mWM[5],cAIC.M3.mWM[5],cAIC.M4.mWM[5]))
df.mWM <- unlist(c(cAIC.M1.mWM[2],cAIC.M2.mWM[2],cAIC.M3.mWM[2],cAIC.M4.mWM[2]))
DAIC.mWM <- Delta.AIC(AIC.mWM)

Cand.model.mWM <- as.data.frame(cbind(AIC.mWM, DAIC.mWM, df.mWM))
Model.selection.mWM <- Cand.model.mWM %>%
  reframe("Model" = c(1:4),
          "cAIC" = round(AIC.mWM,2), 
          "DcAIC" = round(DAIC.mWM,2),
          "Estimated DF" = round(df.mWM,2)
  ) %>%
  arrange(cAIC)

Model.selection.mWM

summary(M1.mWM) #view summary of model with lowest AIC

simulationOutput.M1.mWM <- simulateResiduals(fittedModel = M1.mWM, plot = F) #assess fit of model with lowest AIC
plotQQunif(simulationOutput.M1.mWM) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M1.mWM) # right plot in plot.DHARMa()

write.csv(coef(summary(M1.mWM)), "Outputs/FixedEff_M1_mWM.csv")
write.csv(Model.selection.mWM,"Outputs/ModSel_mWM.csv")

#--- All lizards Mt+1 ---#

M1.mALL <- lmer(sqrt(Mt1+0.01) ~ Season.centred*BA*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = DF)
M2.mALL <- lmer(sqrt(Mt1+0.01) ~ Season.centred*BA*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = DF)
M3.mALL <- lmer(sqrt(Mt1+0.01) ~ Season.centred*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = DF)
M4.mALL <- lmer(sqrt(Mt1+0.01) ~ Season.centred*Inundated + T.max + Pit.repl + (Season.centred|Grid), data = DF)

models.mALL <- list(M1.mALL,M2.mALL,M3.mALL,M4.mALL)
Avg.models.mALL <- modelAvg(models.mALL, opt = TRUE)

cAIC.M1.mALL <- cAIC(M1.mALL)
cAIC.M2.mALL <- cAIC(M2.mALL)
cAIC.M3.mALL <- cAIC(M3.mALL)
cAIC.M4.mALL <- cAIC(M4.mALL)

AIC.mALL <- unlist(c(cAIC.M1.mALL[5],cAIC.M2.mALL[5],cAIC.M3.mALL[5],cAIC.M4.mALL[5]))
df.mALL <- unlist(c(cAIC.M1.mALL[2],cAIC.M2.mALL[2],cAIC.M3.mALL[2],cAIC.M4.mALL[2]))
DAIC.mALL <- Delta.AIC(AIC.mALL)

Cand.model.mALL <- as.data.frame(cbind(AIC.mALL, DAIC.mALL, df.mALL))
Model.selection.mALL <- Cand.model.mALL %>%
  reframe("Model" = c(1:4),
          "cAIC" = round(AIC.mALL,2), 
          "DcAIC" = round(DAIC.mALL,2),
          "Estimated DF" = round(df.mALL,2)
  ) %>%
  arrange(cAIC)

Model.selection.mALL

summary(M4.mALL) #view summary of model with lowest AIC

simulationOutput.M4.mALL <- simulateResiduals(fittedModel = M4.mALL, plot = F) #assess fit of model with lowest AIC
plotQQunif(simulationOutput.M4.mALL) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M4.mALL) # right plot in plot.DHARMa()

write.csv(coef(summary(M4.mALL)), "Outputs/FixedEff_M4_mALL.csv")
write.csv(Model.selection.mALL,"Outputs/ModSel_mALL.csv")

#=== Plot best model outcomes ===#

#Referenced: https://www.azandisresearch.com/2022/12/31/visualize-mixed-effect-regressions-in-r-with-ggplot2/

#--- CPUE ---#

#OP

OP$labelsOP <- rep("Inundated",times=length(OP$CPUE))
OP$labelsOP[OP$Inundated == 0] <- "Not affected"
OP$Labels <- paste(OP$Grid,sep=": ",OP$labelsOP)

OP <- OP %>% 
  mutate(fit.m = predict(M1.cOP, re.form = NA), # marginal fits (i.e. fixed effects only)
         fit.c = predict(M1.cOP, re.form = NULL)) #conditional fits (i.e. fixed + random effects)

OPplot <- 
  ggplot(OP, aes(y = log(CPUE+0.01),x=Season,color=labelsOP)) +
  geom_smooth(data = subset(OP,Season < 2020.915), aes(y = fit.m),method="lm",se=FALSE) +
  geom_smooth(data = subset(OP,Season > 2020.915), aes(y = fit.m),method="lm",se=FALSE) +
  geom_point() +
  #geom_point(data = subset(OP,Season < 2020.915), aes(y = fit.m),pch=3) +
  #geom_point(data = subset(OP,Season > 2020.915), aes(y = fit.m),pch=3) +
  scale_colour_manual(values=c('Not affected'="#ffac4c",'Inundated'= "#CC5500")) +
  geom_vline(xintercept = 2020.33,color="#89c8f4",linewidth=1) +
  #facet_wrap(~factor(Labels,
  #                   levels=c('MP1: Inundated','MP2: Not affected','RE1: Not affected',
  #                            'WP3: Inundated','WP2: Not affected','RE2: Not affected'))) +
  #facet_wrap(~factor(labelsOP)) +
  scale_y_continuous(name="Log-transformed CPUE") +
  scale_x_continuous(name="Date") +
  ggtitle(label = "Oligosoma polychroma") + 
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),legend.title=element_blank())

#WM

WM$labelsWM <- rep("Inundated",times=length(WM$CPUE))
WM$labelsWM[WM$Inundated == 0] <- "Not affected"
WM$Labels <- paste(WM$Grid,sep=": ",WM$labelsWM)

WM <- WM %>% 
  mutate(fit.m = predict(M3.cWM, re.form = NA), # marginal fits (i.e. fixed effects only)
         fit.c = predict(M3.cWM, re.form = NULL)) # conditional fits (i.e. fixed + random effects)

WMplot <- 
  ggplot(WM, aes(y = log(Mt1+0.01),x=Season,color=labelsWM)) +
  geom_smooth(data = WM,aes(y = fit.m),method="lm",se = FALSE) + #No BA term in best model
  geom_point() + 
  #geom_point(data = subset(WM,Season < 2020.915), aes(y = fit.m),pch=3) +
  #geom_point(data = subset(WM,Season > 2020.915), aes(y = fit.m),pch=3) +
  scale_colour_manual(values=c('Not affected'="#92d050",'Inundated'= "#228B22")) +
  geom_vline(xintercept = 2020.33,color="#89c8f4",linewidth=1) +
  #facet_wrap(~factor(Labels,
  #                   levels=c('MP1: Inundated','MP2: Not affected','RE1: Not affected',
  #                            'WP3: Inundated','WP2: Not affected','RE2: Not affected'))) +
  #facet_wrap(~factor(labelsWM)) +
  scale_y_continuous(name="Log-transformed CPUE") +
  scale_x_continuous(name="Date") +
  ggtitle(label = "Woodworthia maculata") + 
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),legend.title=element_blank())

#All

DF$labelsDF <- rep("Inundated",times=length(DF$CPUE))
DF$labelsDF[DF$Inundated == 0] <- "Not affected"
DF$Labels <- paste(DF$Grid,sep=": ",DF$labelsDF)

DF <- DF %>% 
  mutate(fit.m = predict(M3.cALL, re.form = NA), # marginal fits (i.e. fixed effects only)
         fit.c = predict(M3.cALL, re.form = NULL)) #conditional fits (i.e. fixed + random effects)

DFplot <- 
  ggplot(DF, aes(y = sqrt(Mt1+0.01),x=Season,color=labelsDF)) +
  geom_smooth(data = DF,aes(y = fit.m),method="lm",se = FALSE) + #no BA term
  geom_point() + 
  #geom_point(data = subset(DF,Season < 2020.915), aes(y = fit.m),pch=3) +
  #geom_point(data = subset(DF,Season > 2020.915), aes(y = fit.m),pch=3) +
  scale_colour_manual(values=c('Not affected'="grey70",'Inundated'= "grey40")) +
  geom_vline(xintercept = 2020.33,color="#89c8f4",linewidth=1) +
  #facet_wrap(~factor(Labels,
  #                   levels=c('MP1: Inundated','MP2: Not affected','RE1: Not affected',
  #                            'WP3: Inundated','WP2: Not affected','RE2: Not affected'))) +
  #facet_wrap(~factor(labelsDF)) +
  scale_y_continuous(name="Square root-transformed CPUE") +
  scale_x_continuous(name="Date") +
  ggtitle(label = "All lizards") + 
  theme_bw() +
  theme(legend.title=element_blank())

#--- Mt+1 ---#

OP <- OP %>% 
  mutate(fit.m = predict(M7.mOP, re.form = NA), #marginal fits (i.e. fixed effects only)
         fit.c = predict(M7.mOP, re.form = NULL)) #conditional fits (i.e. fixed + random effects)

OPplot <- 
  ggplot(OP, aes(y = log(Mt1+0.01),x=Season,color=labelsOP)) +
  #geom_smooth(data = OP, aes(y = fit.m),method="lm",se=FALSE) + #No BA term, no significant effects of time
  geom_point() +
  #geom_point(data = subset(OP,Season < 2020.915), aes(y = fit.m),pch=3) +
  #geom_point(data = subset(OP,Season > 2020.915), aes(y = fit.m),pch=3) +
  scale_colour_manual(values=c('Not affected'="#ffac4c",'Inundated'= "#CC5500")) +
  geom_vline(xintercept = 2020.33,color="#89c8f4",linewidth=1) +
  #facet_wrap(~factor(Labels,
  #                   levels=c('MP1: Inundated','MP2: Not affected','RE1: Not affected',
  #                            'WP3: Inundated','WP2: Not affected','RE2: Not affected'))) +
  #facet_wrap(~factor(labelsOP)) +
  scale_y_continuous(name="Log-transformed Mt+1 / day") +
  scale_x_continuous(name="Date") +
  ggtitle(label = "Oligosoma polychroma") + 
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),legend.title=element_blank())

#WM

WM <- WM %>% 
  mutate(fit.m = predict(M1.cWM, re.form = NA), # marginal fits (i.e. fixed effects only)
         fit.c = predict(M1.cWM, re.form = NULL)) #conditional fits (i.e. fixed + random effects)

WMplot <- 
  ggplot(WM, aes(y = log(Mt1+0.01),x=Season,color=labelsWM)) +
  geom_smooth(data = subset(WM,Season < 2020.915),aes(y = fit.m),method="lm",se = FALSE) + 
  geom_smooth(data = subset(WM,Season > 2020.915),aes(y = fit.m),method="lm", se = FALSE) +
  geom_point() + 
  #geom_point(data = subset(WM,Season < 2020.915), aes(y = fit.m),pch=3) +
  #geom_point(data = subset(WM,Season > 2020.915), aes(y = fit.m),pch=3) +
  scale_colour_manual(values=c('Not affected'="#92d050",'Inundated'= "#228B22")) +
  geom_vline(xintercept = 2020.33,color="#89c8f4",linewidth=1) +
  #facet_wrap(~factor(Labels,
  #                   levels=c('MP1: Inundated','MP2: Not affected','RE1: Not affected',
  #                            'WP3: Inundated','WP2: Not affected','RE2: Not affected'))) +
  #facet_wrap(~factor(labelsWM)) +
  scale_y_continuous(name="Log-transformed Mt+1 / day") +
  scale_x_continuous(name="Date") +
  ggtitle(label = "Woodworthia maculata") + 
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),legend.title=element_blank())

#All

DF <- DF %>% 
  mutate(fit.m = predict(M4.mALL, re.form = NA), # marginal fits (i.e. fixed effects only)
         fit.c = predict(M4.mALL, re.form = NULL)) #conditional fits (i.e. fixed + random effects)

DFplot <- 
  ggplot(DF, aes(y = sqrt(Mt1+0.01),x=Season,color=labelsDF)) +
  #geom_smooth(data = subset(DF,Season < 2020.915),aes(y = fit.m),method="lm",se = FALSE) + 
  #geom_smooth(data = subset(DF,Season > 2020.915),aes(y = fit.m),method="lm", se = FALSE) +
  geom_point() + 
  #geom_point(data = subset(DF,Season < 2020.915), aes(y = fit.m),pch=3) +
  #geom_point(data = subset(DF,Season > 2020.915), aes(y = fit.m),pch=3) +
  scale_colour_manual(values=c('Not affected'="grey70",'Inundated'= "grey40")) +
  geom_vline(xintercept = 2020.33,color="#89c8f4",linewidth=1) +
  #facet_wrap(~factor(Labels,
  #                   levels=c('MP1: Inundated','MP2: Not affected','RE1: Not affected',
  #                            'WP3: Inundated','WP2: Not affected','RE2: Not affected'))) +
  #facet_wrap(~factor(labelsDF)) +
  scale_y_continuous(name="Square root-transformed Mt+1 / day") +
  scale_x_continuous(name="Date") +
  ggtitle(label = "All lizards") + 
  theme_bw() +
  theme(legend.title=element_blank())

#--- Export multiplots ---#

#This is generic code for plotting either the CPUE or Mt+1 panel plots; change as needed
LotsaPlots <- grid.arrange(OPplot, WMplot, DFplot, ncol=1)
ggsave("Outputs/Mt1plots.png",plot=LotsaPlots,
       units="cm",width=16.5,height=25.5)
