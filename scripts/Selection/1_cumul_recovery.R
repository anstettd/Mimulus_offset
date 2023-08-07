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

cumul <- read_csv("~/Documents/Git repos/snp_change/data/binomial_data_half/time_cumul_beagle.csv")

demog_recovery <- left_join(demog_recovery,offset_pop,by=c("Site"="Site")) %>% rename(Site_Name=Site)

demo_pop <- left_join(cumul,demog_recovery,by="Paper_ID") %>% filter(Paper_ID!=10) %>% filter(Paper_ID!=12)

anoms <- read_csv("data/climate_data/climate_anomaly.csv")

demo_pop <- left_join(demo_pop, anoms)


# inspect effect of climate anomalies on demographic trends
drought.period <- demo_pop %>% 
  dplyr::select(cumul_pos, cumul_all, lambda.slope.decline, lambda.mean.drought, MAT_1215, MAP_1215, PAS_1215, CMD_1215, Tave_wt_1215, Tave_sm_1215, PPT_wt_1215, PPT_sm_1215)

ggpairs(drought.period)
# populations with lower mean lambda during drought have higher cumulative positive selection (r=-0.692*); when lambda was lower selection was stronger 
# populations with lower mean lambda during drought had stronger summer precipitation anomalies (r=+-0.667*); when drought was stronger lambda was lower

recovery.period <- demo_pop %>% 
  dplyr::select(cumul_pos, cumul_all, lambda.slope.recovery, lambda.mean.recovery, MAT_1619, MAP_1619, PAS_1619, CMD_1619, Tave_wt_1619, Tave_sm_1619, PPT_wt_1619, PPT_sm_1619)

ggpairs(recovery.period)
# cumulative positive selection positively associated with MAP anomaly: drier sites associated with greater response to selection
# lambda mean during recovery period not associated with any climate variables or cumulative selection

#stats
lm1 <- lm(lambda.slope.recovery~cumul_pos,data=demo_pop)
lm2 <- lm(lambda.slope.recovery~cumul_all,data=demo_pop)
lm3 <- lm(lambda.mean.recovery~cumul_pos,data=demo_pop)
lm3b <- lm(lambda.mean.recovery~cumul_pos + CMD_1619, data=demo_pop)
lm3c <- lm(lambda.mean.recovery~cumul_pos + PPT_sm_1619, data=demo_pop)
lm4 <- lm(lambda.mean.recovery~cumul_all,data=demo_pop)
lm5 <- lm(lambda.mean.recovery~MAT_1619, data=recovery.period)
lm6 <- lm(lambda.mean.recovery~MAP_1619, data=recovery.period)
lm7 <- lm(lambda.mean.recovery~PAS_1619, data=recovery.period)
lm8 <- lm(lambda.mean.recovery~CMD_1619, data=recovery.period)
lm9 <- lm(lambda.mean.recovery~Tave_wt_1619, data=recovery.period)
lm10 <- lm(lambda.mean.recovery~Tave_sm_1619, data=recovery.period)
lm11 <- lm(lambda.mean.recovery~PPT_wt_1619, data=recovery.period)
lm12 <- lm(lambda.mean.recovery~PPT_sm_1619, data=recovery.period)
lm13 <- lm(lambda.mean.drought~PPT_sm_1215, data=demo_pop)

# lam slope ~ pos sel
summary(lm1)
Anova(lm1,type="III")

# lam slope ~ all sel
summary(lm2)
Anova(lm2,type="III")

# lam mean ~ pos sel
summary(lm3)
Anova(lm3,type="III")
summary(lm3b) #w/ CMD covariate
Anova(lm3b,type="III")
summary(lm3c) #w/ PPT_sm covariate
Anova(lm3c,type="III")

# lam mean ~ all sel
summary(lm4)
Anova(lm4,type="III")

summary(lm5)
Anova(lm5,type="III")

summary(lm6)
Anova(lm6,type="III")

summary(lm7)
Anova(lm7,type="III")

summary(lm8)
Anova(lm8,type="III")

summary(lm9)
Anova(lm9,type="III")

summary(lm10)
Anova(lm10,type="III")

summary(lm11)
Anova(lm11,type="III")

summary(lm12)
Anova(lm12,type="III")

summary(lm13)
Anova(lm13,type="III")

###########################################################################################################

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(10,"Spectral"))
n.sites <- length(unique(demo_pop$Paper_ID))
color.list <- lat_cols(n.sites)


#Lambda Slope Recovery

#Directional Positive Selection
ggplot(demo_pop, aes(x=cumul_pos, y=lambda.slope.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Rate of Increase in Lambda")+
  scale_x_continuous(name="Directional Selection")+
                     #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Selection_demo/1_cumul_pos_recovery_lambda.pdf",width=8, height = 6, units = "in")


#Total Selection
ggplot(demo_pop, aes(x=cumul_all, y=lambda.slope.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Rate of Increase in Lambda")+
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
ggsave("Graphs/Selection_demo/2_cumul_all_recovery_lambda.pdf",width=8, height = 6, units = "in")


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
ggsave("Graphs/Selection_demo/3_cumul_pos_recovery_mean.pdf",width=8, height = 6, units = "in")


#Total Selection
ggplot(demo_pop, aes(x=cumul_all, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
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
ggsave("Graphs/Selection_demo/4_cumul_all_recovery_mean.pdf",width=8, height = 6, units = "in")


#Population decline and drought anomaly
ggplot(demo_pop, aes(x=PPT_sm_1215, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda during Drought")+
  scale_x_continuous(name="Summer precipitation anomaly")+
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
ggsave("Graphs/Selection_demo/ppt_anom_decline_mean.pdf",width=8, height = 6, units = "in")


#Population recovery and climate anomalies
a <- ggplot(demo_pop, aes(x=PPT_sm_1619, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm, color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Summer precipitation anomaly")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5),
    axis.title.y = element_text(color="black", size=16,vjust = 2, hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )

b <- ggplot(demo_pop, aes(x=PPT_wt_1619, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm, color="black", lty="dashed", se=FALSE)+
  #scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Winter precipitation anomaly")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5),
    #axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )

c <- ggplot(demo_pop, aes(x=MAP_1619, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm, color="black", lty="dashed", se=FALSE)+
  #scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Annual precipitation anomaly")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5),
    #axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )

d <- ggplot(demo_pop, aes(x=Tave_sm_1619, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm, color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Summer temperature anomaly")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5),
    axis.title.y = element_text(color="black", size=16,vjust = 2, hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )

e <- ggplot(demo_pop, aes(x=Tave_wt_1619, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm, color="black", lty="dashed", se=FALSE)+
  #scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Winter temperature anomaly")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5),
    #axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )

f <- ggplot(demo_pop, aes(x=MAT_1619, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm, color="black", lty="dashed", se=FALSE)+
  #scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Annual temperature anomaly")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5),
    #axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )

g <- ggplot(demo_pop, aes(x=CMD_1619, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm, color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Climate moisture deficit anomaly")+
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
library(patchwork)
(a+b+c)/(d+e+f) + plot_layout(guides='collect')

c + f + plot_layout(guides="collect")
