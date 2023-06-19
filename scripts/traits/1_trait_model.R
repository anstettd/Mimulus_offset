##########################################################################################################
## Run trait model
## Author Daniel Anstett
## 
##
## Last Modified June 19, 2023
##########################################################################################################
#Import libraries
library(tidyverse)
library(lme4)
library(lmtest)
library(visreg)


##########################################################################################################
#Import files
y3 <- read.csv("data/trait_data/y3.csv", header=T) #Import trait data

#Set factors
y3$Block <- as.factor(y3$Block) ; y3$Family <- as.factor(y3$Family) # prep factors

# Set up vectors with treatment and regional information
treatment.v<-c("W", "D")
pop.v<-c("32.9_S02", "34.3_S07", "36.7_S08", "36.2_S10", "34.1_S11", "43.4_S15",
            "41.8_S16", "41.7_S17", "39.7_S18", "39.4_S29","37.5_S32", "42.3_S36")
site.v<-c("1","12","2","3","4","5","6","7","8","9","10","11")
order.row<-1

#Make input data frame
trait_dat <- data.frame()
trait_dat[1:24,1] <- as.data.frame(pop.v)
trait_dat[1:24,2] <- as.data.frame(site.v)
trait_dat[1:12,3] <- "Wet"
trait_dat[13:24,3] <- "Dry"

trait_dat_SLA <- trait_dat
trait_dat_SLA$trait <- "SLA"

trait_dat_fl <- trait_dat
trait_dat_fl$trait <- "fl"


##########################################################################################################
##########################################################################################################
## Run Site*Year*Drought Mixed Models
trait.SLA <- lmer(SLA ~ Site.Lat*Year*Drought + Block + (1|Family),data=y3)
trait.fl <- lmer(Experiment_Date ~ Site.Lat*Year*Drought + Block + (1|Family),data=y3)


for (i in 1:2){
  vis_SLA<-visreg(trait.SLA, xvar="Year", by="Site.Lat", cond=list(Drought=treatment.v[i]))
  Res_SLA<-vis_SLA$res
  for (j in 1:12){
    Ref_SLA_filter<-Res_SLA %>% filter(Site.Lat==pop.v[j])
    Ref_SLA_filter<- Ref_SLA_filter %>% mutate(Res.scale=scale(visregRes))
    lm_SLA<-lm(Res.scale~Year, data=Ref_SLA_filter)
    summary_SLA<-summary(lm_SLA)
    trait_dat_SLA[order.row,5]<-summary_SLA$coefficients[2,1]
    trait_dat_SLA[order.row,6]<-summary_SLA$coefficients[2,2]
    order.row<-order.row+1
  }
}

order.row<-1
for (i in 1:2){
  vis_fl<-visreg(trait.fl, xvar="Year", by="Site.Lat", cond=list(Drought=treatment.v[i]))
  Res_fl<-vis_fl$res
  for (j in 1:12){
    Ref_fl_filter<-Res_fl %>% filter(Site.Lat==pop.v[j])
    Ref_fl_filter<- Ref_fl_filter %>% mutate(Res.scale=scale(visregRes))
    lm_fl<-lm(Res.scale~Year, data=Ref_fl_filter)
    summary_fl<-summary(lm_fl)
    trait_dat_fl[order.row,5]<-summary_fl$coefficients[2,1]
    trait_dat_fl[order.row,6]<-summary_fl$coefficients[2,2]
    order.row<-order.row+1
  }
}

#Join slope tables
slope.trait <- rbind(trait_dat_SLA,trait_dat_fl)
colnames(slope.trait) <- c("Site_ID","Pop","Treatment","Trait","Slope","STDER")


write.csv(slope.trait,"data/trait_data/slope.trait.csv")








