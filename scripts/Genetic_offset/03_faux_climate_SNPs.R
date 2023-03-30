##################################################################################
## Gradient forest data prep for env only analysis
## Author Daniel Anstett
## 
## 
## 
## Last Modified March 30, 2023
###################################################################################
##Libraries
library(tidyverse)

#Function 
#scale values between 0 & 1
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

##Import Data
climate <- read_csv("data/genomic_data/climate_pop.csv")

faux_1 <- as.data.frame(scale_values(climate$MAT))
faux_2 <- as.data.frame(scale_values(climate$MAP))
faux_3 <- as.data.frame(scale_values(climate$PAS))
faux_4 <- as.data.frame(scale_values(climate$EXT))
faux_5 <- as.data.frame(scale_values(climate$CMD))
faux_6 <- as.data.frame(scale_values(climate$Tave_wt))
faux_7 <- as.data.frame(scale_values(climate$Tave_sm))
faux_8 <- as.data.frame(scale_values(climate$PPT_wt))
faux_9 <- as.data.frame(scale_values(climate$PPT_sm))

faux_prep <- cbind(faux_1,faux_2,faux_3,faux_4,faux_5,faux_6,faux_7,faux_8,faux_9)
colnames(faux_prep) <- c("faux1","faux2","faux3","faux4","faux5","faux6","faux7","faux8","faux9")
faux <- cbind(climate,faux_prep)


write_csv(faux,"data/genomic_data/faux_snp.csv")




