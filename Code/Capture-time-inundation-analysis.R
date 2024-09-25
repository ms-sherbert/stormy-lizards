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


M.cOP <- lmer(log(CPUE+0.01) ~ BA*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = OP)
summary(M.cOP) #view summary of model

write.csv(coef(summary(M.cOP)), "Outputs/FixedEff_M_cOP.csv")

simulationOutput.M.cOP <- simulateResiduals(fittedModel = M.cOP, plot = F) #assess fit of model
plotQQunif(simulationOutput.M.cOP) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M.cOP) # right plot in plot.DHARMa()

#--- Woodworthia maculata CPUE ---#

#No warning messages produced by any of these models
M.cWM <- lmer(log(CPUE+0.01) ~ BA*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = WM)
summary(M.cWM) #view summary of model

simulationOutput.M.cWM <- simulateResiduals(fittedModel = M.cWM, plot = F) #assess fit of model with lowest AIC
plotQQunif(simulationOutput.M.cWM) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M.cWM) # right plot in plot.DHARMa()

write.csv(coef(summary(M.cWM)), "Outputs/FixedEff_M_cWM.csv")


#--- Oligosoma polychroma Mt+1 ---#

M1.mOP <- lmer(log(Mt1+0.01) ~ BA*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = OP)
M2.mOP <- lmer(log(Mt1+0.01) ~ BA*Inundated + T.min + Pit.repl + (1|Grid), data = OP)


#Test whether simplifying random term creates a significant difference 
#in any pair of models with same fixed factor structure
anova(M1.mOP,M2.mOP) #p = 1

#No, so decided to use models with random intercept term (1|Grid)\

summary(M2.mOP) #view summary of model

simulationOutput.M2.mOP <- simulateResiduals(fittedModel = M2.mOP, plot = F) #assess fit of model with lowest AIC
plotQQunif(simulationOutput.M2.mOP) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M2.mOP) # right plot in plot.DHARMa()

write.csv(coef(summary(M2.mOP)), "Outputs/FixedEff_M2_mOP.csv")


#--- Woodworthia maculata Mt+1 ---#

M.mWM <- lmer(log(Mt1+0.01) ~ BA*Inundated + T.min + Pit.repl + (Season.centred|Grid), data = WM)
summary(M.mWM) #view summary of model

simulationOutput.M.mWM <- simulateResiduals(fittedModel = M.mWM, plot = F) #assess fit of model with lowest AIC
plotQQunif(simulationOutput.M.mWM) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M.mWM) # right plot in plot.DHARMa()

write.csv(coef(summary(M.mWM)), "Outputs/FixedEff_M_mWM.csv")

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


#--- Export multiplots ---#

#This is generic code for plotting either the CPUE or Mt+1 panel plots; change as needed
LotsaPlots <- grid.arrange(OPplot, WMplot, DFplot, ncol=1)
ggsave("Outputs/Mt1plots.png",plot=LotsaPlots,
       units="cm",width=16.5,height=25.5)
