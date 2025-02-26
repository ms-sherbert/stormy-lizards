#--- Simplified graphics ---#
# Written by Sarah Herbert
# for R version 4.3.1
# Last test date: 26/02/2025

#--- Preamble ---#

rm(list=ls())
setwd("C:/Repositories/stormy-lizards") 
#NB repositories in D: on uni desktop, C: on personal laptop


library(ggplot2)
library(lubridate)
library(tidyverse)


CPUElizards <- read.csv("Outputs/CPUElizards.csv")
inundation<-strptime("15/04/2020",format="%d/%m/%Y")

# Simplified CPUE graph

CPUEsummary <- CPUElizards %>%
  group_by(Species,Inundated,Date) %>%
  summarise('Treatment_mean' = mean(n/Nsessions),
            'SE' = sd(n/Nsessions,na.rm=TRUE)/sum(n/Nsessions,na.rm=TRUE))

CPUEsummary$Date<-strptime(CPUEsummary$Date,format="%Y-%m-%d")

CPUEsummaryplot<-ggplot(data=subset(CPUEsummary,Species != "oa"),
                        aes(x=Date,y=Treatment_mean,ymin=Treatment_mean-SE,ymax=Treatment_mean+SE,color=Inundated)) +
  geom_pointrange(stat="identity") +
  facet_wrap(~factor(Species)) +
  #scale_color_manual(values=c('wm'="#440154FF",'op'= "#20A387FF"),
  #                  labels=c('wm' = "Woodworthia maculata", 'op' = "Oligosoma polychroma")) +
  geom_vline(xintercept = inundation,color="#89c8f4",linewidth=1) +
  scale_y_continuous(name = "CPUE") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"))