##########################################################################################################
## Run trait model
## Author Daniel Anstett
## 
##
## Last Modified June 19, 2023
##########################################################################################################
#Import libraries
library(tidyverse)
library(RColorBrewer)
library(car)


slope.trait <-read_csv("data/trait_data/slope.trait.csv")

offset_pop <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/snp_change/data/binomial_data/time_cumul_beagle.csv")

#Join
evol <- left_join(offset_pop,slope.trait,by="Paper_ID")

##########################################################################################################
#Run models

lm.sla.wet.cumul <- lm(SLA_Wet~cumul_all,data=evol)
Anova(lm.sla.wet.cumul,type="III")
summary(lm.sla.wet.cumul)



# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(evol$Paper_ID))
color.list <- lat_cols(n.sites)

##########################################################################################################
#SLA Wet
ggplot(evol, aes(x=cumul_pos, y=SLA_Wet)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA_Wet")+
  scale_x_continuous(name="Cumulative Positive Evolution")+
                     #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("Graphs/traits/1_cumul_pos_SLA_wet.pdf",width=8, height = 6, units = "in")

#2012-2015 offset plotted against lambda
ggplot(evol, aes(x=cumul_all, y=SLA_Wet)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA_Wet")+
  scale_x_continuous(name="Cumulative Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("Graphs/traits/2_cumul_all_SLA_wet.pdf",width=8, height = 6, units = "in")



##########################################################################################################
#SLA_Dry
ggplot(evol, aes(x=cumul_pos, y=SLA_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA_Dry")+
  scale_x_continuous(name="Cumulative Positive Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("Graphs/traits/3_cumul_pos_SLA_dry.pdf",width=8, height = 6, units = "in")

#2012-2015 offset plotted against lambda
ggplot(evol, aes(x=cumul_all, y=SLA_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA_Dry")+
  scale_x_continuous(name="Cumulative Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("Graphs/traits/4_cumul_all_SLA_dry.pdf",width=8, height = 6, units = "in")



##########################################################################################################
#2012-2015 offset plotted against lambda
ggplot(evol, aes(x=cumul_pos, y=fl_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="fl_wet")+
  scale_x_continuous(name="Cumulative Positive Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("Graphs/traits/5_cumul_pos_fl_Dry.pdf",width=8, height = 6, units = "in")

#2012-2015 offset plotted against lambda
ggplot(evol, aes(x=cumul_all, y=fl_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="fl_wet")+
  scale_x_continuous(name="Cumulative Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("Graphs/traits/6_cumul_all_fl_Dry.pdf",width=8, height = 6, units = "in")



##########################################################################################################
#fl_Dry
ggplot(evol, aes(x=cumul_pos, y=fl_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="fl_Dry")+
  scale_x_continuous(name="Cumulative Positive Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("Graphs/traits/7_cumul_pos_fl_dry.pdf",width=8, height = 6, units = "in")

#2012-2015 offset plotted against lambda
ggplot(evol, aes(x=cumul_all, y=fl_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="fl_Dry")+
  scale_x_continuous(name="Cumulative Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("Graphs/traits/8_cumul_all_fl_dry.pdf",width=8, height = 6, units = "in")




#lm.1215 <- lm(lambda.slope.trunc~offset_1215,data=offset_pop)
