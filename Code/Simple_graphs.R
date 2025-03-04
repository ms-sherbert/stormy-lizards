#--- Simplified graphics ---#
# Written by Sarah Herbert
# for R version 4.3.1
# Last test date: 26/02/2025

#--- Preamble ---#

rm(list=ls())
setwd("D:/Repositories/stormy-lizards") 
#NB repositories in D: on uni desktop, C: on personal laptop


library(ggplot2)
library(lubridate)
library(tidyverse)
library(gridExtra)

CPUElizards <- read.csv("Outputs/CPUElizards.csv")
SVLs <- read.csv("outputs/SVLs.csv")
Mtlizards <- read.csv("Outputs/Mtlizards.csv")

inundation<-strptime("15/04/2020",format="%d/%m/%Y")

#Make nicer species display names
CPUElizards["Species"][CPUElizards["Species"] == "wm"] <- "Woodworthia maculata"
CPUElizards["Species"][CPUElizards["Species"] == "op"] <- "Oligosoma polychroma"

Mtlizards["Species"][Mtlizards["Species"] == "wm"] <- "Woodworthia maculata"
Mtlizards["Species"][Mtlizards["Species"] == "op"] <- "Oligosoma polychroma"

SVLs["Species"][SVLs["Species"] == "wm"] <- "Woodworthia maculata"
SVLs["Species"][SVLs["Species"] == "op"] <- "Oligosoma polychroma"

#Add n plots to labels
CPUElizards["Inundated"][CPUElizards["Inundated"] == "Inundated"] <- "Inundated (N = 2)"
CPUElizards["Inundated"][CPUElizards["Inundated"] == "Not affected"] <- "Not affected (N = 4)"

Mtlizards["Inundated"][Mtlizards["Inundated"] == "Inundated"] <- "Inundated (N = 2)"
Mtlizards["Inundated"][Mtlizards["Inundated"] == "Not affected"] <- "Not affected (N = 4)"

SVLs["inundated"][SVLs["inundated"] == 1] <- "Inundated (N = 2)"
SVLs["inundated"][SVLs["inundated"] == 0] <- "Not affected (N = 4)"

# Simplified CPUE graph

CPUEsummary <- CPUElizards %>%
  group_by(Species,Inundated,Date) %>%
  summarise('Treatment_mean' = mean(n/Nsessions),
            'SE' = sd(n/Nsessions,na.rm=TRUE)/sqrt(sum(n/Nsessions,na.rm=TRUE)))

CPUEsummary$Date<-as.POSIXct(strptime(CPUEsummary$Date,format="%Y-%m-%d"))

CPUEsummaryplot<-ggplot(data=subset(CPUEsummary,Species != "oa"),
                        aes(x=Date,y=Treatment_mean,ymin=Treatment_mean-SE,ymax=Treatment_mean+SE,color=Inundated)) +
  geom_pointrange(stat="identity", position="jitter") +
  facet_wrap(~factor(Species),ncol=2) +
  scale_color_manual(values=c('Inundated (N = 2)'="#023E8A",'Not affected (N = 4)'= "#FFC20A")) +
  labs(color = "Storm surge impact", title="A") +
  geom_vline(xintercept = inundation,color="#89c8f4",linewidth=1) +
  scale_y_continuous(name = "CPUE") +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"))
  
#Simplified SVL graph

SVLsummary <- SVLs %>%
  group_by(Species,inundated,Season) %>%
  summarise('Treatment_mean' = mean(SVL_mm),
            'SE' = sd(SVL_mm,na.rm=TRUE)/sqrt(sum(SVL_mm,na.rm=TRUE)),
            'SD' = sd(SVL_mm,na.rm=TRUE),
            'Min' = min(SVL_mm),
            'Max' = max(SVL_mm))

SVLsummaryplot<-ggplot(data=SVLsummary,
                        aes(x=Season,y=Treatment_mean,ymin=Min,ymax=Max,color=inundated)) +
  geom_pointrange(stat="identity", position="jitter") +
  facet_wrap(~factor(Species),ncol=2) +
  scale_color_manual(values=c('Inundated (N = 2)'="#023E8A",'Not affected (N = 4)'= "#FFC20A")) +
  labs(color = "Storm surge impact",title="B") +
  geom_vline(xintercept = decimal_date(inundation),color="#89c8f4",linewidth=1) +
  scale_y_continuous(name = "SVL (mm)") +
  scale_x_continuous(name = "Date") +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"))

# Simplified CPUE graph

Mtsummary <- Mtlizards %>%
  group_by(Species,Inundated,Date) %>%
  summarise('Treatment_mean' = mean(n.Nsessions),
            'SE' = sd(n.Nsessions,na.rm=TRUE)/sqrt(sum(n.Nsessions,na.rm=TRUE)))

Mtsummary$Date<-as.POSIXct(strptime(Mtsummary$Date,format="%Y-%m-%d"))

Mtsummaryplot<-ggplot(data=subset(Mtsummary,Species != "oa"),
                        aes(x=Date,y=Treatment_mean,ymin=Treatment_mean-SE,ymax=Treatment_mean+SE,color=Inundated)) +
  geom_pointrange(stat="identity", position="jitter") +
  facet_wrap(~factor(Species),ncol=2) +
  scale_color_manual(values=c('Inundated (N = 2)'="#023E8A",'Not affected (N = 4)'= "#FFC20A")) +
  labs(color = "Storm surge impact", y = expression(paste('M'[t+1]*'day'^-1))) +
  geom_vline(xintercept = inundation,color="#89c8f4",linewidth=1) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"))


#Export graphs

Fig2 <- grid.arrange(grobs=list(CPUEsummaryplot,SVLsummaryplot),nrow=2)
ggsave("Outputs/Fig2ab.tif",plot=Fig2,device="tiff",
       height=16, width = 16, units = "cm",dpi="print")

ggsave("Outputs/Mtsummaryplot.tif",plot=Mtsummaryplot,device="tiff",
       height=10, width = 16, units = "cm",dpi="print")
