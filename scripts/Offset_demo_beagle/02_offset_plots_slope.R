##################################################################################
## Plot offset against lambda.slope
## Author Daniel Anstett
## 
## 
## Last Modified April 19, 2023
###################################################################################

#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)

#Import data
offset_pop <- read_csv("data/genomic_data/offset_pop_beagle.csv")
#offset_pop <- offset_pop %>% filter(Demo_ID!=24)
offset_pop_10 <- offset_pop %>% filter(Paper_ID<13)

#stats
lm.1215 <- lm(lambda.slope.trunc~offset_1215,data=offset_pop)
lm.SSP245 <- lm(lambda.slope.trunc~offset_SSP245,data=offset_pop)
lm.SSP585 <- lm(lambda.slope.trunc~offset_SSP585,data=offset_pop)
lm.dist <- lm(lambda.slope.trunc~offset_climate,data=offset_pop)

lm.1215_10 <- lm(lambda.slope.trunc~offset_1215,data=offset_pop_10)
lm.SSP245_10 <- lm(lambda.slope.trunc~offset_SSP245,data=offset_pop_10)
lm.SSP585_10 <- lm(lambda.slope.trunc~offset_SSP585,data=offset_pop_10)
lm.dist_10 <- lm(lambda.slope.trunc~offset_climate,data=offset_pop_10)

summary(lm.1215)
summary(lm.SSP245)
summary(lm.SSP585)


Anova(lm.1215,type="III")
Anova(lm.SSP245,type="III")
Anova(lm.SSP585,type="III")
Anova(lm.dist,type="III")


Anova(lm.1215_10,type="III")
Anova(lm.SSP245_10,type="III")
Anova(lm.SSP585_10,type="III")
Anova(lm.dist_10,type="III")






###########################################################################################################



#2012-2015 offset plotted against lambda
ggplot(offset_pop, aes(x=offset_1215, y=lambda.slope.trunc, label=Demo_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
#  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
#ggsave("Graphs/lambda_slope/1_offset_lambda_1215.pdf",width=8, height = 6, units = "in")


#SSP 245 offset plotted against lambda
ggplot(offset_pop, aes(x=offset_SSP245, y=lambda.slope.trunc, label=Demo_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2040-2070 SSP245 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/lambda_slope/2_offset_lambda_ssp245.pdf",width=7, height = 5, units = "in")

#SSP 585 offset plotted against lambda
ggplot(offset_pop, aes(x=offset_SSP585, y=lambda.slope.trunc, label=Demo_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2040-2070 SSP585 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/lambda_slope/3_offset_lambda_ssp585.pdf",width=7, height = 5, units = "in")




###########################

#climate distance plotted against lambda
ggplot(offset_pop, aes(x=offset_climate, y=lambda.slope.trunc, label=Demo_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
#ggsave("Graphs/lambda_slope/4_distance_lambda.pdf",width=8, height = 6, units = "in")





##################################################################################################################
#Offset pop with 10 initial populations

#2012-2015 offset plotted against lambda only 10
ggplot(offset_pop_10, aes(x=offset_1215, y=lambda.slope.trunc, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
#ggsave("Graphs/1_offset_lambda_1215_only10.pdf",width=8, height = 6, units = "in")



#SSP 245 offset plotted against lambda
ggplot(offset_pop_10, aes(x=offset_SSP245, y=lambda.slope.trunc, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2040-2070 SSP245 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
#ggsave("Graphs/2_offset_lambda_45_only10.pdf",width=7, height = 5, units = "in")

#SSP 585 offset plotted against lambda
ggplot(offset_pop_10, aes(x=offset_SSP585, y=lambda.slope.trunc, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2040-2070 SSP585 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
#ggsave("Graphs/3_offset_lambda_85_only10.pdf",width=7, height = 5, units = "in")


#climate distance plotted against lambda
ggplot(offset_pop_10, aes(x=offset_climate, y=lambda.slope.trunc, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
#ggsave("Graphs/4_distance_lambda_only10.pdf",width=8, height = 6, units = "in")


##################################################################################################################









#By region
ggplot(offset_pop,aes(x=offset_1215, y=lambda.slope.trunc, color=Region,label=Demo_ID))+
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method="lm")+
  geom_text(hjust=-.15, vjust=-.2,color="black")+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))


