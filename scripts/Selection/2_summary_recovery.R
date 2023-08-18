##################################################################################
## Plot lambda mean against slope summaries
## Author Daniel Anstett
## 
## 
## Last Modified Aug 14, 2023
###################################################################################
#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)
library(RColorBrewer)

#Import data & Prepare data frame
offset_pop <- read_csv("data/genomic_data/offset_pop_beagle.csv") %>% dplyr::select(Site, Paper_ID) #just to get translation of pop names <--> numbers
offset_pop[20,1] <- "Mill Creek"
offset_pop[20,2] <- 12
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")
demog_recovery <- left_join(demog_recovery,offset_pop,by=c("Site"="Site")) %>% rename(Site_Name=Site)

slope.summary <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/snp_change/data/binomial_data_half/slope_summary.csv")

demo_pop <- left_join(slope.summary,demog_recovery,by=c("Site"="Paper_ID")) %>% filter(Site!=10) %>% filter(Site!=12)


#stats
lm1 <- lm(lambda.mean.recovery~median,data=demo_pop)
summary(lm1)
Anova(lm1,type="III")

lm2 <- lm(lambda.mean.recovery~pos_slope,data=demo_pop)
summary(lm2)
Anova(lm2,type="III")

lm3 <- lm(lambda.mean.recovery~abs_slope,data=demo_pop)
summary(lm3)
Anova(lm3,type="III")

lm4 <- lm(lambda.slope.recovery~median,data=demo_pop)
summary(lm4)
Anova(lm4,type="III")

lm5 <- lm(lambda.slope.recovery~pos_slope,data=demo_pop)
summary(lm5)
Anova(lm5,type="III")

lm6 <- lm(lambda.slope.recovery~abs_slope,data=demo_pop)
summary(lm6)
Anova(lm6,type="III")

###########################################################################################################
# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(10,"Spectral"))
n.sites <- length(unique(demo_pop$Site))
color.list <- lat_cols(n.sites)


#Median slope vs. lambda.mean.recovery
ggplot(demo_pop, aes(x=median, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Population Growth Rate")+
  scale_x_continuous(name="Median Slope")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Selection_demo/Slope_summary/1_median_mean_recovery_lambda.pdf",width=8, height = 6, units = "in")


#Positive slope vs. lambda.mean.recovery
ggplot(demo_pop, aes(x=pos_slope, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Population Growth Rate")+
  scale_x_continuous(name="Positive Selection")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Selection_demo/Slope_summary/2_pos_mean_recovery_lambda.pdf",width=8, height = 6, units = "in")

#Total slope vs. lambda.mean.recovery
ggplot(demo_pop, aes(x=abs_slope, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Population Growth Rate")+
  scale_x_continuous(name="Total Selection")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Selection_demo/Slope_summary/3_abs_mean_recovery_lambda.pdf",width=8, height = 6, units = "in")

###########################################################################################################

#Median slope vs. lambda.slope.recovery
ggplot(demo_pop, aes(x=median, y=lambda.slope.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Rate of Increase in Lambda")+
  scale_x_continuous(name="Median Slope")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Selection_demo/Slope_summary/4_median_slope_recovery_lambda.pdf",width=8, height = 6, units = "in")


#Positive slope vs. lambda.slope.recovery
ggplot(demo_pop, aes(x=pos_slope, y=lambda.slope.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Rate of Increase in Lambda")+
  scale_x_continuous(name="Positive Selection")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Selection_demo/Slope_summary/5_pos_slope_recovery_lambda.pdf",width=8, height = 6, units = "in")

#Total slope vs. lambda.slope.recovery
ggplot(demo_pop, aes(x=abs_slope, y=lambda.slope.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Rate of Increase in Lambda")+
  scale_x_continuous(name="Total Selection")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Selection_demo/Slope_summary/6_abs_slope_recovery_lambda.pdf",width=8, height = 6, units = "in")




