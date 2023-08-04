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
library(cowplot)


slope.trait <-read_csv("data/trait_data/slope.trait.csv")
offset_pop <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/snp_change/data/binomial_data_half/time_cumul_beagle.csv")
pi.df <- read_csv("data/genomic_data/raw_pi.csv")

#Join
evol <- left_join(offset_pop,slope.trait,by="Paper_ID") %>% left_join(pi.df,by=c("Paper_ID"="Site"))


##########################################################################################################
#Run models
#SLA Wet
lm1 <- lm(SLA_Wet~cumul_pos,data=evol)
Anova(lm1,type="III")
summary(lm1)

lm2 <- lm(SLA_Wet~cumul_all,data=evol)
Anova(lm2,type="III")
summary(lm2)

#SLA Dry
lm3 <- lm(SLA_Dry~cumul_pos,data=evol)
Anova(lm3,type="III")
summary(lm3)

lm4 <- lm(SLA_Dry~cumul_all,data=evol)
Anova(lm4,type="III")
summary(lm4)

#Fl Wet
lm5 <- lm(fl_Wet~cumul_pos,data=evol)
Anova(lm5,type="III")
summary(lm5)

lm6 <- lm(fl_Wet~cumul_all,data=evol)
Anova(lm6,type="III")
summary(lm6)

#Fl Dry
lm7 <- lm(fl_Dry~cumul_pos,data=evol)
Anova(lm7,type="III")
summary(lm7)

lm8 <- lm(fl_Dry~cumul_all,data=evol)
Anova(lm8,type="III")
summary(lm8)



# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(evol$Paper_ID))
color.list <- lat_cols(n.sites)

##########################################################################################################
#SLA Wet
plot_1 <- ggplot(evol, aes(x=cumul_pos, y=SLA_Wet)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA Slope Wet")+
  scale_x_continuous(name="Cumulative Positive Evolution")+
                     #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )


#Positive and negative rapid evolution
plot_2 <- ggplot(evol, aes(x=cumul_all, y=SLA_Wet)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA Slope Wet")+
  scale_x_continuous(name="Cumulative Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )




##########################################################################################################
#SLA_Dry
plot_3 <- ggplot(evol, aes(x=cumul_pos, y=SLA_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA Slope Dry")+
  scale_x_continuous(name="Cumulative Positive Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )


#Positive and negative rapid evolution
plot_4 <- ggplot(evol, aes(x=cumul_all, y=SLA_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA Slope Dry")+
  scale_x_continuous(name="Cumulative Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )




##########################################################################################################
#fl_Wet
plot_5 <- ggplot(evol, aes(x=cumul_pos, y=fl_Wet)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Flowering Slope Wet")+
  scale_x_continuous(name="Cumulative Positive Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )


#Positive and negative rapid evolution
plot_6 <- ggplot(evol, aes(x=cumul_all, y=fl_Wet)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Flowering Slope Wet")+
  scale_x_continuous(name="Cumulative Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )




##########################################################################################################
#fl_Dry
plot_7 <- ggplot(evol, aes(x=cumul_pos, y=fl_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Flowering Slope Dry")+
  scale_x_continuous(name="Cumulative Positive Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )


#Positive and negative rapid evolution
plot_8 <- ggplot(evol, aes(x=cumul_all, y=fl_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Flowering Slope Dry")+
  scale_x_continuous(name="Cumulative Evolution")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )



############################################################################################################
#Cowplot
plot_grid(plot_1,plot_2,plot_3,plot_4,plot_5,plot_6,plot_7,plot_8, labels = "AUTO",ncol = 2,label_x = 0.23) #export at 7 X 14 


