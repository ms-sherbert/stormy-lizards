# Functions to quickly calculate Delta AIC and AIC weight values from candidate model sets
# Written by S Herbert
# For R version 4.3.1
# Last tested: 09 June 2025

Delta.AIC <- function(AIC.values){
  DAIC <- AIC.values - min(AIC.values)
  DAIC              
}

AIC.weight <- function (DIC.values){ 
  Llik <- exp(-0.5*DIC.values)
  Wts <- Llik/sum(Llik)
  Wts
}