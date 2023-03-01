#############################################################################################################
## Plot Manhattan for WZA windows
## Author Daniel Anstett
## 
## 
## Modified from Tom Booker WZA Vignette
## Last Modified August 3, 2022
#############################################################################################################
#Import libraries

library(tidyverse)
library(Kendall)

#Import files

wza_win_mat <- read_csv("data/genomic_data/WZA_win_mat.csv")
wza_win_map <- read_csv("data/genomic_data/WZA_win_map.csv")
wza_win_cmd <- read_csv("data/genomic_data/WZA_win_cmd.csv")


###########################################################################################################


#Plot Distribution
ggplot(data = wza_win_mat, aes( x = Z_pVal))+
  geom_histogram( bins = 50)+
  scale_y_continuous("Count")+
  theme_classic() #Approx normal

ggplot(data = wza_win_map, aes( x = Z_pVal))+
  geom_histogram( bins = 50)+
  scale_y_continuous("Count")+
  theme_classic() #Approx normal

ggplot(data = wza_win_cmd, aes( x = Z_pVal))+
  geom_histogram( bins = 50)+
  scale_y_continuous("Count")+
  theme_classic() #Approx normal



#######################################################################################################


# WZA Empirical p-value

#MAT
wza_empri_mat <- ggplot(data = wza_win_mat, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_mat)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "darkgoldenrod"), 22 )) +
  theme_classic()+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_mat
ggsave("Graphs/wza_mat.png", wza_empri_mat, width=10, height = 5, units = "in")


#MAP
wza_empri_map <- ggplot(data = wza_win_map, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_map)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  theme_classic()+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_map
ggsave("Graphs/wza_map.png", wza_empri_map, width=10, height = 5, units = "in")


#CMD
wza_empri_cmd <- ggplot(data = wza_win_cmd, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_cmd)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "magenta3"), 22 )) +
  theme_classic()+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_cmd
ggsave("Graphs/wza_cmd.png", wza_empri_cmd, width=10, height = 5, units = "in")



