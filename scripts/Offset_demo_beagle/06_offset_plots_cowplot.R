##################################################################################
## Plot offset against lambda.mean
## Author Daniel Anstett
## 
## 
## Last Modified April 19, 2023
###################################################################################

#Library install and import
library(tidyverse) 
library(RColorBrewer)
library(cowplot)

#Import data
offset_pop <- read_csv("data/genomic_data/offset_pop_beagle.csv")
#Import data
offset_cumul <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/snp_change/data/binomial_data/time_cumul_beagle.csv")


###########################################################################################################

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(offset_pop$Site))
color.list <- lat_cols(n.sites)


#2012-2015 offset plotted against lambda
plot_1 <- ggplot(offset_pop, aes(x=offset_1215, y=lambda.slope.trunc)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2012-2015 Genetic Offset",breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(1, "lines")  # Increase the size of the legend dots
  )

#SSP 245 offset plotted against lambda
plot_2 <- ggplot(offset_pop, aes(x=offset_SSP245, y=lambda.slope.trunc)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2040-2070 SSP245 Genetic Offset",breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(1, "lines")  # Increase the size of the legend dots
  )

#SSP 585 offset plotted against lambda
plot_3 <- ggplot(offset_pop, aes(x=offset_SSP585, y=lambda.slope.trunc)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2040-2070 SSP585 Genetic Offset",breaks=c(0.055,0.06,0.065,0.07,0.075))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(1, "lines")  # Increase the size of the legend dots
  )

#climate distance plotted against lambda
plot_4 <- ggplot(offset_pop, aes(x=offset_climate, y=lambda.slope.trunc)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(1, "lines")  # Increase the size of the legend dots
  )


##################################################################################################################

#2012-2015 offset plotted against lambda
plot_5 <- ggplot(offset_pop, aes(x=offset_1215, y=lambda.mean)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
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
    legend.key.size = unit(1, "lines")  # Increase the size of the legend dots
  )

#SSP 245 offset plotted against lambda
plot_6 <- ggplot(offset_pop, aes(x=offset_SSP245, y=lambda.mean)) + 
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
    legend.key.size = unit(1, "lines")  # Increase the size of the legend dots
  )

#SSP 585 offset plotted against lambda
plot_7 <- ggplot(offset_pop, aes(x=offset_SSP585, y=lambda.mean)) + 
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
    legend.key.size = unit(1, "lines")  # Increase the size of the legend dots
  )

#climate distance plotted against lambda
plot_8 <- ggplot(offset_pop, aes(x=offset_climate, y=lambda.mean)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
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
    legend.key.size = unit(1, "lines")  # Increase the size of the legend dots
  )
###########################################################################################################
#Cumul against climate distance


#cumul slope plotted against 2012-2015 Climate Distance
plot_cumul1 <- ggplot(offset_cumul, aes(x=offset_climate, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank())



#cumul slope plotted against 2012-2015 Climate Distance
plot_cumul2 <-ggplot(offset_cumul, aes(x=offset_climate, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank())


###########################################################################################################
#Make Cowplots

plot_grid(plot_2,plot_3,plot_5,plot_7, labels = "AUTO",ncol = 2,label_x = 0.23) #export at 6 X 8 

plot_grid(plot_cumul1,plot_cumul2,plot_4,plot_8, labels = "AUTO",ncol = 2,label_x = 0.23,label_size = 0) #export at 6 X 8 







