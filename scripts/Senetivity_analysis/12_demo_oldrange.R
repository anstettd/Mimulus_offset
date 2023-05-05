##################################################################################
## Get positive slope info per population
## Author Daniel Anstett
## 
## 
## Last Modified Jan 25, 2023
###################################################################################

#Library install and import
library(tidyverse) 

#Import offset
offset_pop <- read_csv("data/genomic_data/offset_pop_sense.csv")


#Import lambda from demography data for all sites
demo_all <- read_csv("data/demography data/siteYear.lambda_2010-2016.csv")

demo_select <- demo_all %>% filter(Site !="Mill Creek" & Site != "Deer Creek" ) 
demo_select <- demo_select %>% filter(!Year%in%c(2010,2011))
demo_select <- na.omit(demo_select)


demo_order <- demo_select[order(demo_select$Latitude),]
demo_ID <- data.frame(unique(demo_order$Site))
demo_ID$ID<-1:nrow(demo_ID)
names(demo_ID) <- c("Site","ID")
demo_select_ID <- left_join(demo_select,demo_ID,by="Site")
demo_select_raw <- demo_select %>% select(-Year,-lambda,-SiteYear)
demo_ID_join <- left_join(demo_ID,unique(demo_select_raw),by="Site")

###################################################################################
#Calc lambda slopes over time

lambda_pop <- data.frame()
for (i in 1:nrow(demo_ID)){
  demo_pop <- demo_select_ID %>% filter(ID == i)
  lambda_pop[i,1] <- summary(lm(lambda~Year,data=demo_pop))$coefficients[2]
}

lambda_pop_id <- cbind(demo_ID_join,lambda_pop)
colnames(lambda_pop_id)[8] <- "lambda_slope_1215"

lambda_pop_sort <- lambda_pop_id[order(lambda_pop_id$Latitude, decreasing = TRUE),] 

offset_pop_lambda <- cbind(offset_pop,lambda_pop_sort)

write_csv(offset_pop_lambda,"data/genomic_data/offset_pop_sense_oldrange.csv")






