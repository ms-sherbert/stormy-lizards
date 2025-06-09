##=== CPUE and Mt+1 (weighted by the number of check days) calculations ===##
# Written by Sarah Herbert
# for R version 4.3.1
# Last test date: 09 June 2025

#--- Preamble ---#

rm(list=ls())
setwd("C:/Repositories/stormy-lizards") 
#NB repositories in D: on uni desktop, C: on personal laptop

library(tidyverse)

caps<-read.csv("Data/CC-grid-captures_S1-S8.csv")
sessions<-read.csv("Data/N-sessions.csv")

caps["Capture_type"][caps["Capture_type"] == "Marked "] <- "Marked"
levels(as.factor(caps$Capture_type))

Date<-strptime(sessions$date,format="%d/%m/%Y")
sessions<-cbind(sessions,Date)

#--- Mt+1 ---#

Nlizards1<-subset(caps, Capture_type == "Marked" & Species == "op" 
                  | Capture_type == "Marked" & Species == "wm"
                  | Capture_type == "Marked" & Species == "oa") %>%
  count(Grid,Species,Season,ID,sort=FALSE,.drop = FALSE)


Nlizards2<-Nlizards1 %>%
  count(Grid,Species,Season,sort=FALSE,.drop = FALSE)

missing<-cbind(c("MP1","MP1","MP2","WP2","RE1","RE1","RE1","RE1","RE1","RE2","RE1","RE2","WP2"),
               c("wm","wm","wm","op","wm","wm","wm","wm","wm","wm","oa","wm","wm"),
               c(1,2,1,6,1,2,3,4,5,1,3,7,7),
               rep(0,13)
)
colnames(missing) <- c("Grid","Species","Season","n")

Nlizards2<-rbind(Nlizards2,missing)
Nlizards2$Season<-as.integer(Nlizards2$Season)
Nlizards2$n<-as.numeric(Nlizards2$n)

Nlizards3<-full_join(Nlizards2,sessions,join_by(Season,Grid))

#Add column with Mt+1 (n individuals) weighted by the number of sessions
Mtlizards<-Nlizards3 %>%
  group_by(Grid,Season) %>%
  mutate(n/Nsessions)

Mtlizards <- Mtlizards %>% arrange(Species, Grid, Season)
Site <- substr(Mtlizards$Grid,1,2)
Inundated <- rep("Not affected",times=length(Mtlizards$Grid))
Mtlizards <- cbind(Mtlizards,Site,Inundated)
names(Mtlizards)[names(Mtlizards) == "...11"] <- "Site" 
names(Mtlizards)[names(Mtlizards) == "...12"] <- "Inundated"  

Mtlizards["Inundated"][Mtlizards["Grid"] == "WP3" | Mtlizards["Grid"] == "MP1"] <- "Inundated"
Mtlizards$Labels <- paste(Mtlizards$Grid,sep=": ",Mtlizards$Inundated)

write.csv(Mtlizards,file="Outputs/Mtlizards.csv")

# Graph Mt+1 values by grid (Figure 2)

TotalMt<-Mtlizards %>%
  group_by(Labels,Date) %>%
  summarise('nNsessions' = sum(n/Nsessions))

TotalMt$Species <- rep("All",times=length(TotalMt$Date)) 

inundation<-strptime("15/04/2020",format="%d/%m/%Y")

Mtplot <- ggplot(data=Mtlizards, aes(x=Date,y=n/Nsessions,color=Species)) +
  geom_point(stat="identity") +
  #geom_point(data=TotalMt,aes(x=Date,y=nNsessions),stat="identity",
  #           color="black", pch=1) +
  facet_wrap(~factor(Labels, 
             levels=c('MP1: Inundated','WP3: Inundated',
                      'MP2: Not affected','WP2: Not affected',
                      'RE1: Not affected','RE2: Not affected')),
             nrow=3,ncol=2) +
  scale_color_manual(values=c('wm'="#440154FF",'op'= "#20A387FF",'oa'="#B8DE29FF"),
                     labels=c('wm' = "Woodworthia maculata", 
                              'op' = "Oligosoma polychroma", 
                              'oa' = "Oligosoma aeneum")) +
  geom_vline(xintercept = inundation,color="#89c8f4",linewidth=1) +
  #scale_y_continuous(name = "Mt+1 / day") +
  labs(y = expression(paste('M'[t+1]*'day'^-1))) +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"))

ggsave("Outputs/Fig2.png",plot=Mtplot,device="png",
       height=20, width = 16, units = "cm",dpi="print")

#--- CPUE ---#

CPUElizards1<-subset(caps, Capture_type != "Dead" & Species == "op" 
                  | Capture_type != "Dead" & Species == "wm"
                  | Capture_type != "Dead" & Species == "oa") %>%
  count(Grid,Species,Season,ID,sort=FALSE,.drop = FALSE)


CPUElizards2<-CPUElizards1 %>%
  count(Grid,Species,Season,sort=FALSE,.drop = FALSE)

CPUEmissing<-cbind(c("MP1","MP1","MP2","RE1","RE1","RE1","RE1","RE1","RE2","RE2"),
               c("wm","wm","wm","wm","wm","wm","wm","wm","wm","wm"),
               c(1,2,1,1,2,3,4,5,1,7),
               rep(0,10)
)
colnames(CPUEmissing) <- c("Grid","Species","Season","n")

CPUElizards2<-rbind(CPUElizards2,CPUEmissing)
CPUElizards2$Season<-as.integer(CPUElizards2$Season)
CPUElizards2$n<-as.numeric(CPUElizards2$n)

CPUElizards3<-full_join(CPUElizards2,sessions,join_by(Season,Grid))

#Add column with Mt+1 (n individuals) weighted by the number of sessions
CPUElizards<-CPUElizards3 %>%
  group_by(Grid,Season) %>%
  mutate(n/Nsessions)

CPUElizards <- CPUElizards %>% arrange(Species, Grid, Season)

Site <- substr(CPUElizards$Grid,1,2)
Inundated <- rep("Not affected",times=length(CPUElizards$Grid))
CPUElizards <- cbind(CPUElizards,Site,Inundated)
names(CPUElizards)[names(CPUElizards) == "...11"] <- "Site" 
names(CPUElizards)[names(CPUElizards) == "...12"] <- "Inundated"  

CPUElizards["Inundated"][CPUElizards["Grid"] == "WP3" | CPUElizards["Grid"] == "MP1"] <- "Inundated"
CPUElizards$Labels <- paste(CPUElizards$Grid,sep=": ",CPUElizards$Inundated)

write.csv(CPUElizards,file="Outputs/CPUElizards.csv")

TotalCPUE<-CPUElizards %>%
  group_by(Labels,Date) %>%
  summarise('nNsessions' = sum(n/Nsessions))

TotalCPUE$Species <- rep("All",times=length(TotalCPUE$Date)) 

# Graph CPUE values by grid (Figure 3)


CPUEplot <- ggplot(data=CPUElizards, aes(x=Date,y=n/Nsessions,color=Species)) +
  geom_point(stat="identity") +
  #geom_point(data=TotalCPUE,aes(x=Date,y=nNsessions),stat="identity",
  #           color="black", pch=1) +
  facet_wrap(~factor(Labels, 
                     levels=c('MP1: Inundated','WP3: Inundated',
                              'MP2: Not affected','WP2: Not affected',
                              'RE1: Not affected','RE2: Not affected')),
                     nrow=3,ncol=2)+
  scale_color_manual(values=c('wm'="#440154FF",'op'= "#20A387FF",'oa'="#B8DE29FF"),
                     labels=c('wm' = "Woodworthia maculata", 'op' = "Oligosoma polychroma", 'oa' = "Oligosoma aeneum")) +
  geom_vline(xintercept = inundation,color="#89c8f4",linewidth=1) +
  scale_y_continuous(name = "CPUE") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"))

ggsave("Outputs/Fig3.png",plot=CPUEplot,device="png",
       height=20, width = 16, units = "cm",dpi="print")

