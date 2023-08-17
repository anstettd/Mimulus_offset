##################################################################################
## Plot lambda.slope & lambda mean against cumulative selection
## Author Daniel Anstett
## 
## 
## Last Modified Aug 3, 2023
###################################################################################

#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)
library(RColorBrewer)
library(GGally)

#Import data & Prepare data frame
offset_pop <- read_csv("data/genomic_data/offset_pop_beagle.csv") %>% dplyr::select(Site, Paper_ID) #just to get translation of pop names <--> numbers
offset_pop[20,1] <- "Mill Creek"
offset_pop[20,2] <- 12

demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

#cumul <- read_csv("~/Documents/Git repos/snp_change/data/binomial_data_half/time_cumul_beagle.csv")
cumul <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/snp_change/data/binomial_data_half/time_cumul_beagle.csv")

demog_recovery <- left_join(demog_recovery,offset_pop,by=c("Site"="Site")) %>% rename(Site_Name=Site)

demo_pop <- left_join(cumul,demog_recovery,by="Paper_ID") %>% filter(Paper_ID!=10) %>% filter(Paper_ID!=12)

anoms <- read_csv("data/climate_data/climate_anomaly.csv")

demo_pop <- left_join(demo_pop, anoms)


#stats
lm1 <- lm(lambda.mean.recovery~cumul_pos_env2,data=demo_pop)
lm2 <- lm(lambda.mean.recovery~cumul_all_env2,data=demo_pop)
lm3 <- lm(lambda.mean.recovery~cumul_pos_env9,data=demo_pop)
lm4 <- lm(lambda.mean.recovery~cumul_all_env9,data=demo_pop)

# lam slope ~ pos sel
summary(lm1)
Anova(lm1,type="III")

# lam slope ~ all sel
summary(lm2)
Anova(lm2,type="III")

# lam mean ~ pos sel
summary(lm3)
Anova(lm3,type="III")

# lam mean ~ all sel
summary(lm4)
Anova(lm4,type="III")


###########################################################################################################

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(10,"Spectral"))
n.sites <- length(unique(demo_pop$Paper_ID))
color.list <- lat_cols(n.sites)


#Mean Lambda

#Directional Positive Selection
ggplot(demo_pop, aes(x=cumul_pos, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Directional Selection")+
  # ,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )
#ggsave("Graphs/Selection_demo/Cumul/1_cumul_env2_pos_recovery_mean.pdf",width=8, height = 6, units = "in")


#Total Selection
ggplot(demo_pop, aes(x=cumul_all, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Total Selection")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )
#ggsave("Graphs/Selection_demo/Cumul/1_cumul_env2_all_recovery_mean.pdf",width=8, height = 6, units = "in")


###########################################################################################################
#Mean Lambda

#Directional Positive Selection
ggplot(demo_pop, aes(x=cumul_pos, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Directional Selection")+
                    # ,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )
#ggsave("Graphs/Selection_demo/Cumul/3_cumul_pos_recovery_mean.pdf",width=8, height = 6, units = "in")


#Total Selection
ggplot(demo_pop, aes(x=cumul_all, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Total Selection")+
                     #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )
#ggsave("Graphs/Selection_demo/Cumul/4_cumul_all_recovery_mean.pdf",width=8, height = 6, units = "in")

