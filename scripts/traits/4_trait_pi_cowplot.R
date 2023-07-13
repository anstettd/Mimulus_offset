##########################################################################################################
## Corelate trait evolution to pi
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
offset_pop <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/snp_change/data/binomial_data/time_cumul_beagle.csv")
pi.df <- read_csv("data/genomic_data/raw_pi.csv")

#Join
evol <- left_join(offset_pop,slope.trait,by="Paper_ID") %>% left_join(pi.df,by=c("Paper_ID"="Site"))


##########################################################################################################
#Run models
#SLA Wet
lm1 <- lm(SLA_Wet~pi_snp_set,data=evol)
Anova(lm1,type="III")
summary(lm1)

lm2 <- lm(SLA_Wet~pi_all_snps,data=evol)
Anova(lm2,type="III")
summary(lm2)

#SLA Dry
lm3 <- lm(SLA_Dry~pi_snp_set,data=evol)
Anova(lm3,type="III")
summary(lm3)

lm4 <- lm(SLA_Dry~pi_all_snps,data=evol)
Anova(lm4,type="III")
summary(lm4)

#Fl Wet
lm5 <- lm(fl_Wet~pi_snp_set,data=evol)
Anova(lm5,type="III")
summary(lm5)

lm6 <- lm(fl_Wet~pi_all_snps,data=evol)
Anova(lm6,type="III")
summary(lm6)

#Fl Dry
lm7 <- lm(fl_Dry~pi_snp_set,data=evol)
Anova(lm7,type="III")
summary(lm7)

lm8 <- lm(fl_Dry~pi_all_snps,data=evol)
Anova(lm8,type="III")
summary(lm8)



# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(evol$Paper_ID))
color.list <- lat_cols(n.sites)

##########################################################################################################
#SLA Wet
plot_1 <- ggplot(evol, aes(x=pi_snp_set, y=SLA_Wet)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA Wet")+
  scale_x_continuous(name="pi SNP Set")+
                     #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
legend.position = "none",
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
#ggsave("Graphs/traits/pi/1_pi_snp_set_SLA_wet.pdf",width=8, height = 6, units = "in")

#Positive and negative rapid evolution
plot_2 <- ggplot(evol, aes(x=pi_all_snps, y=SLA_Wet)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA Wet")+
  scale_x_continuous(name="pi Genome-Wide")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
legend.position = "none",
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
#ggsave("Graphs/traits/pi/2_pi_all_snps_SLA_wet.pdf",width=8, height = 6, units = "in")



##########################################################################################################
#SLA_Dry
plot_3 <- ggplot(evol, aes(x=pi_snp_set, y=SLA_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA Dry")+
  scale_x_continuous(name="pi SNP Set")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
legend.position = "none",
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
#ggsave("Graphs/traits/pi/3_pi_snp_set_SLA_dry.pdf",width=8, height = 6, units = "in")

#Positive and negative rapid evolution
plot_4 <- ggplot(evol, aes(x=pi_all_snps, y=SLA_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="SLA Dry")+
  scale_x_continuous(name="pi Genome-Wide")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
legend.position = "none",
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
#ggsave("Graphs/traits/pi/4_pi_all_snps_SLA_dry.pdf",width=8, height = 6, units = "in")



##########################################################################################################
#fl_Wet
plot_5 <- ggplot(evol, aes(x=pi_snp_set, y=fl_Wet)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Date of FLowering Wet")+
  scale_x_continuous(name="pi SNP Set")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
legend.position = "none",
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
#ggsave("Graphs/traits/pi/5_pi_snp_set_fl_Dry.pdf",width=8, height = 6, units = "in")

#Positive and negative rapid evolution
plot_6 <- ggplot(evol, aes(x=pi_all_snps, y=fl_Wet)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Date of Flowering Wet")+
  scale_x_continuous(name="pi Genome-Wide")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
legend.position = "none",
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
#ggsave("Graphs/traits/pi/6_pi_all_snps_fl_Dry.pdf",width=8, height = 6, units = "in")



##########################################################################################################
#fl_Dry
plot_7 <- ggplot(evol, aes(x=pi_snp_set, y=fl_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Date of Flowering Dry")+
  scale_x_continuous(name="pi SNP set")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
legend.position = "none",
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
#ggsave("Graphs/traits/pi/7_pi_snp_set_fl_dry.pdf",width=8, height = 6, units = "in")

#Positive and negative rapid evolution
plot_8 <- ggplot(evol, aes(x=pi_all_snps, y=fl_Dry)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Date of Flowering Dry")+
  scale_x_continuous(name="pi Genome-Wide")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
legend.position = "none",
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
#ggsave("Graphs/traits/pi/8_pi_all_snps_fl_dry.pdf",width=8, height = 6, units = "in")

############################################################################################################
#Cowplot
plot_grid(plot_1,plot_2,plot_3,plot_4,plot_5,plot_6,plot_7,plot_8, labels = "AUTO",ncol = 2,label_x = 0.23) #export at 7 X 14 




