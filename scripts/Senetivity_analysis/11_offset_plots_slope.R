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
offset_pop <- read_csv("data/genomic_data/offset_pop_sense.csv")
offset_pop_10 <- offset_pop %>% filter(Paper_ID<13)

#bf30 
#env9 
#oldWZA 
#old_snp_set 
#new_snp_set 


#stats
lm.bf30 <- lm(lambda.slope~offset_bf30,data=offset_pop)
lm.env9  <- lm(lambda.slope~offset_env9 ,data=offset_pop)
lm.oldWZA <- lm(lambda.slope~offset_oldWZA,data=offset_pop)
lm.old_snp_set  <- lm(lambda.slope~offset_old_snp_set,data=offset_pop)
lm.new_snp_set  <- lm(lambda.slope~offset_new_snp_set ,data=offset_pop)

lm.bf30_10 <- lm(lambda.slope~offset_bf30,data=offset_pop_10)
lm.env9_10 <- lm(lambda.slope~offset_env9,data=offset_pop_10)
lm.oldWZA_10 <- lm(lambda.slope~offset_oldWZA,data=offset_pop_10)
lm.old_snp_set_10 <- lm(lambda.slope~offset_old_snp_set,data=offset_pop_10)
lm.new_snp_set_10 <- lm(lambda.slope~offset_new_snp_set,data=offset_pop_10)


Anova(lm.bf30,type="III")
Anova(lm.env9,type="III")
Anova(lm.oldWZA ,type="III")
Anova(lm.old_snp_set ,type="III")
Anova(lm.new_snp_set,type="III")


Anova(lm.bf30_10,type="III")
Anova(lm.env9_10,type="III")
Anova(lm.oldWZA_10,type="III")
Anova(lm.old_snp_set_10,type="III")
Anova(lm.new_snp_set_10,type="III")






###########################################################################################################
#2012-2015 offset plotted against lambda

#BF30
ggplot(offset_pop, aes(x=offset_bf30, y=lambda.slope, label=Demo_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
#  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="BF>30 Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/Sensetivity/1_offset_lambda_bf30.pdf",width=8, height = 6, units = "in")


#env9 - new WZA 
ggplot(offset_pop, aes(x=offset_env9, y=lambda.slope, label=Demo_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="New WZA Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/Sensetivity/2_offset_lambda_env9.pdf",width=7, height = 5, units = "in")

#Old WZA
ggplot(offset_pop, aes(x=offset_oldWZA , y=lambda.slope, label=Demo_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="Old WZA Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/Sensetivity/3_offset_lambda_oldWZA.pdf",width=7, height = 5, units = "in")

#Old snp set
ggplot(offset_pop, aes(x=offset_old_snp_set  , y=lambda.slope, label=Demo_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="Old WZA + BF Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/Sensetivity/4_offset_lambda_old_snp_set.pdf",width=7, height = 5, units = "in")


#New WZA
ggplot(offset_pop, aes(x=offset_new_snp_set  , y=lambda.slope, label=Demo_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="New WZA + BF Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/Sensetivity/5_offset_lambda_new_snp_set.pdf",width=7, height = 5, units = "in")



##################################################################################################################
#Offset pop with 10 initial populations

#BF30
ggplot(offset_pop_10, aes(x=offset_bf30, y=lambda.slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="BF>30 Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/Sensetivity/1_offset_lambda_bf30_only10.pdf",width=8, height = 6, units = "in")


#env9 - new WZA 
ggplot(offset_pop_10, aes(x=offset_env9, y=lambda.slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="New WZA Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/Sensetivity/2_offset_lambda_env9_only10.pdf",width=7, height = 5, units = "in")

#Old WZA
ggplot(offset_pop_10, aes(x=offset_oldWZA , y=lambda.slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="Old WZA Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/Sensetivity/3_offset_lambda_oldWZA_only10.pdf",width=7, height = 5, units = "in")

#Old snp set
ggplot(offset_pop_10, aes(x=offset_old_snp_set  , y=lambda.slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="Old WZA + BF Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/Sensetivity/4_offset_lambda_old_snp_set_only10.pdf",width=7, height = 5, units = "in")


#New WZA
ggplot(offset_pop_10, aes(x=offset_new_snp_set  , y=lambda.slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="New WZA + BF Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs/Sensetivity/5_offset_lambda_new_snp_set_only10.pdf",width=7, height = 5, units = "in")








#2012-2015 offset plotted against lambda only 10
ggplot(offset_pop_10, aes(x=offset_bf30, y=lambda.slope, label=Paper_ID)) + 
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
#ggsave("Graphs/Sensetivity/1_offset_lambda_bf30_only10.pdf",width=8, height = 6, units = "in")



#SSP 245 offset plotted against lambda
ggplot(offset_pop_10, aes(x=offset_env9, y=lambda.slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2040-2070 env9 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
#ggsave("Graphs/Sensetivity/2_offset_lambda_45_only10.pdf",width=7, height = 5, units = "in")

#SSP 585 offset plotted against lambda
ggplot(offset_pop_10, aes(x=offset_oldWZA, y=lambda.slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2040-2070 oldWZA Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
#ggsave("Graphs/Sensetivity/3_offset_lambda_85_only10.pdf",width=7, height = 5, units = "in")

