##################################################################################
## Plot offset against lambda.mean
## Author Daniel Anstett
## 
## 
## Last Modified April 19, 2023
###################################################################################

#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)
library(RColorBrewer)


#Import data
offset <- read_csv("data/genomic_data/offset_pop_beagle.csv") %>% select(Site,offset_1215,offset_climate,offset_SSP245,offset_SSP585)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

offset_pop <- left_join(offset,demog_recovery,by="Site")


#stats
lm.1215 <- lm(lambda.mean.drought~offset_1215,data=offset_pop)
lm.SSP245 <- lm(lambda.mean.drought~offset_SSP245,data=offset_pop)
lm.SSP585 <- lm(lambda.mean.drought~offset_SSP585,data=offset_pop)
lm.dist <- lm(lambda.mean.drought~offset_climate,data=offset_pop)

summary(lm.1215)
Anova(lm.1215,type="III")

summary(lm.SSP245)
Anova(lm.SSP245,type="III")

summary(lm.SSP585)
Anova(lm.SSP585,type="III")

summary(lm.dist)
Anova(lm.dist,type="III")



###########################################################################################################


# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(offset_pop$Site))
color.list <- lat_cols(n.sites)


#2012-2015 offset plotted against lambda
ggplot(offset_pop, aes(x=offset_1215, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black",lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda")+
  scale_x_continuous(name="2012-2015 Genetic Offset",breaks=c(0.025,0.03,0.035,0.04,0.045))+
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
ggsave("Graphs/lambda_gradient/5_offset_mean_1215.pdf",width=8, height = 6, units = "in")


#SSP 245 offset plotted against lambda
ggplot(offset_pop, aes(x=offset_SSP245, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda")+
  scale_x_continuous(name="2040-2070 SSP245 Genetic Offset",breaks=c(0.04,0.045,0.05,0.055,0.06))+
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
#ggsave("Graphs/lambda_gradient/6_offset_mean_ssp245.pdf",width=8, height = 6, units = "in")

#SSP 585 offset plotted against lambda
ggplot(offset_pop, aes(x=offset_SSP585, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda")+
  scale_x_continuous(name="2040-2070 SSP585 Genetic Offset",breaks=c(0.055,0.06,0.065,0.07,0.075))+
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
#ggsave("Graphs/lambda_gradient/7_offset_mean_ssp585.pdf",width=8, height = 6, units = "in")




###########################

#climate distance plotted against lambda
ggplot(offset_pop, aes(x=offset_climate, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black",lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
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
ggsave("Graphs/lambda_gradient/8_distance_mean.pdf",width=8, height = 6, units = "in")




