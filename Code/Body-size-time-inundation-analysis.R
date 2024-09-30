#=== Analysis of trends in body size and condition before and after inundation ===#
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
cor(OPdf[,c(2,11)],method="spearman",use="complete.obs")  # rho = 0.107
cor(WMdf[,c(2,11)],method="spearman",use="complete.obs")  # rho = 0.218

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
M1.svlWM.restart <- update(M1.svlWM, start=pars)
set.seed(101)
pars_x <- runif(length(pars),pars/1.01,pars*1.01)
M1.svlWM.restart2 <- update(M1.svlWM, start=pars_x,
                       control=strict_tol)

M2.svlWM.restart <- update(M2.svlWM, start=pars)
set.seed(101)
pars_x <- runif(length(pars),pars/1.01,pars*1.01)
M2.svlWM.restart2 <- update(M2.svlWM, start=pars_x,
                            control=strict_tol)

#---Try all available optimizers

M1.svlWM.fitcheck <- allFit(M1.svlWM)
M2.svlWM.fitcheck <- allFit(M2.svlWM)

summary(M1.svlWM.fitcheck)
summary(M2.svlWM.fitcheck) #Actually, these summaries seem to indicate that convergence is ok with all optimizers
# Will go ahead and assume that convergence warnings don't mean the fit is incorrect in this case
# see: https://rdrr.io/cran/lme4/man/convergence.html

summary(M1.svlWM)

simulationOutput.M1.svlWM <- simulateResiduals(fittedModel = M1.svlWM, plot = F) #assess fit of model with lowest AIC
plotQQunif(simulationOutput.M1.svlWM) # left plot in plot.DHARMa()
plotResiduals(simulationOutput.M1.svlWM) # right plot in plot.DHARMa()

write.csv(coef(summary(M1.svlWM)), "Outputs/FixedEff_M1_svlWM.csv")


#=== Plot best model outcomes ===#

#Referenced: https://www.azandisresearch.com/2022/12/31/visualize-mixed-effect-regressions-in-r-with-ggplot2/

#--- OP ---#

OPdf$labelsOP <- rep("Inundated",times=length(OPdf$uID))
OPdf$labelsOP[OPdf$inundated == 0] <- "Not affected"

OPdf <- OPdf %>% 
  mutate(fit.m = predict(M1.svlOP, re.form = NA), #marginal fits (i.e. fixed effects only)
         fit.c = predict(M1.svlOP, re.form = NULL)) #conditional fits (i.e. fixed + random effects)

OPplot <- 
  ggplot(OPdf, aes(y = SVL_transOP,x=Season,color=labelsOP)) +
  geom_jitter(width=0.02) +
  geom_smooth(data = OPdf, aes(y = fit.m),method="lm",se=FALSE) + 
  #geom_point(data = subset(OP,Season < 2020.915), aes(y = fit.m),pch=3) +
  #geom_point(data = subset(OP,Season > 2020.915), aes(y = fit.m),pch=3) +
  scale_colour_manual(values=c('Not affected'="#ffac4c",'Inundated'= "#CC5500")) +
  geom_vline(xintercept = 2020.33,color="#89c8f4",linewidth=1) +
  #facet_wrap(~factor(Labels,
  #                   levels=c('MP1: Inundated','MP2: Not affected','RE1: Not affected',
  #                            'WP3: Inundated','WP2: Not affected','RE2: Not affected'))) +
  #facet_wrap(~factor(labelsOP)) +
  scale_y_continuous(name="SVL (mm) ^ 2.475") +
  scale_x_continuous(name="Date") +
  ggtitle(label = "Oligosoma polychroma") + 
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),legend.title=element_blank())

#--- WM ---#

WMdf$labelsWM <- rep("Inundated",times=length(WMdf$uID))
WMdf$labelsWM[WMdf$inundated == 0] <- "Not affected"

WMdf <- WMdf %>% 
  mutate(fit.m = predict(M1.svlWM, re.form = NA), #marginal fits
         fit.c = predict(M1.svlWM, re.form = NULL)) #conditional fits

WMplot <- 
  ggplot(WMdf, aes(y = SVL,x=Season,color=labelsWM)) +
  geom_jitter(width=0.02) + 
  geom_smooth(data = subset(WMdf,Season < 2020.915),aes(y = fit.m),method="lm",se = FALSE) + 
  geom_smooth(data = subset(WMdf,Season > 2020.915),aes(y = fit.m),method="lm",se = FALSE) + 
  #geom_point(data = subset(WMdf,Season < 2020.915), aes(y = fit.m),pch=3) +
  #geom_point(data = subset(WMdf,Season > 2020.915), aes(y = fit.m),pch=3) +
  scale_colour_manual(values=c('Not affected'="#92d050",'Inundated'= "#228B22")) +
  geom_vline(xintercept = 2020.33,color="#89c8f4",linewidth=1) +
  #facet_wrap(~factor(Labels,
  #                   levels=c('MP1: Inundated','MP2: Not affected','RE1: Not affected',
  #                            'WP3: Inundated','WP2: Not affected','RE2: Not affected'))) +
  #facet_wrap(~factor(labelsWM)) +
  scale_y_continuous(name="SVL (mm) ^ 4.5") +
  scale_x_continuous(name="Date") +
  ggtitle(label = "Woodworthia maculata") + 
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),legend.title=element_blank())

multiplot <- grid.arrange(OPplot, WMplot, ncol=1)
ggsave("Outputs/SVLplots.png",plot=multiplot,
       units="cm",width=16.5,height=16.5)

CPUEplot <- ggplot(data=CPUElizards, aes(x=Date,y=n/Nsessions,color=Species)) +
  geom_point(stat="identity") +
  #geom_point(data=TotalCPUE,aes(x=Date,y=nNsessions),stat="identity",
  #           color="black", pch=1) +
  facet_wrap(~factor(Labels, 
                     levels=c('MP1: Inundated','WP3: Inundated',
                              'MP2: Not affected','WP2: Not affected',
                              'RE1: Not affected','RE2: Not affected')),
             nrow=3,ncol=2)+
  scale_color_manual(values=c('wm'="#92d050",'op'= "#ffac4c",'oa'="#815f46"),
                     labels=c('wm' = "Woodworthia maculata", 'op' = "Oligosoma polychroma", 'oa' = "Oligosoma aeneum")) +
  geom_vline(xintercept = inundation,color="#89c8f4",linewidth=1) +
  scale_y_continuous(name = "CPUE") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"))

ggsave("Outputs/Fig3.png",plot=CPUEplot,device="png",
       height=20, width = 16, units = "cm",dpi="print")

#=== Simple panel plots of SVL vs time ===#

inundation<-2020.333

OPsummary <- OPdf %>%
      group_by(Season,Grid,inundated) %>%
      summarize("Mean" = mean(SVL_mm),
                "Min" = min(SVL_mm),
                "Max"=max(SVL_mm))

OPsummary$Inundated <- "Not affected"
OPsummary$Inundated[OPsummary$inundated == 1] <- "Inundated"

OPsummary$Labels <- paste(OPsummary$Grid,sep=": ",OPsummary$Inundated)

WMsummary <- WMdf %>%
  group_by(Season,Grid,inundated) %>%
  summarize("Mean" = mean(SVL_mm),
            "Min" = min(SVL_mm),
            "Max"=max(SVL_mm))

WMsummary$Inundated <- "Not affected"
WMsummary$Inundated[WMsummary$inundated == 1] <- "Inundated"

WMsummary$Labels <- paste(WMsummary$Grid,sep=": ",WMsummary$Inundated)

SVLplot <- ggplot() +
  geom_pointrange(stat="identity",data=WMsummary, aes(x=Season+0.01,y=Mean,ymax=Max,ymin=Min),color="#92d050") +
  geom_pointrange(stat="identity",data=OPsummary, aes(x=Season-0.01,y=Mean,ymax=Max,ymin=Min),color="#ffac4c") +
  facet_wrap(~factor(Labels, 
                     levels=c('MP1: Inundated','WP3: Inundated',
                              'MP2: Not affected','WP2: Not affected',
                              'RE1: Not affected','RE2: Not affected')),
             nrow=3,ncol=2) +
  geom_vline(xintercept = inundation,color="#89c8f4",linewidth=1) +
  scale_y_continuous(name = "SVL (mm)") +
  scale_x_continuous(name = "Season") +
  theme_bw() 
  #theme(legend.text = element_text(face = "italic"))

ggsave("Outputs/Fig4.png",plot=SVLplot,device="png",
       height=20, width = 16, units = "cm",dpi="print")
