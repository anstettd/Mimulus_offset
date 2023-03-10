#############################################################################################################
## Get peak windows from manhattan plots with empirical p-values
## Author Daniel Anstett
## 
## 
## 
## Last Modified Sept 2, 2022
#############################################################################################################
#Import libraries

library(tidyverse)
library(Kendall)
library(plotly)

#Import files
wza_win_mat <- read_csv("data/genomic_data/WZA_win_mat.csv")
wza_win_map <- read_csv("data/genomic_data/WZA_win_map.csv")
wza_win_cmd <- read_csv("data/genomic_data/WZA_win_cmd.csv")

snps_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_mat.csv")
snps_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_map.csv")
snps_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_cmd.csv")

#Filter by Bonferonnii correction alpha critical = 1.29423e-06, aka 5.887988 sigma
mat_bon <- wza_win_mat %>% filter(Z_pVal<1.29423e-06)
map_bon <- wza_win_map %>% filter(Z_pVal<1.29423e-06)
cmd_bon <- wza_win_cmd %>% filter(Z_pVal<1.29423e-06)

#Get peak windows only
#Filter out windows that are not peak windows
mat_bon_peak <- mat_bon %>% filter(!win %in% c(9934,9937,9939,3340)) # 11 windows
map_bon_peak <- map_bon # 2 windows
#%>% filter(!win %in% c())
cmd_bon_peak <- cmd_bon #9 windows
#%>% filter(!win %in% c())

#Filter SNPs for peak
snps_mat_peak <- snps_mat %>% filter(win %in% mat_bon_peak$win)
snps_map_peak <- snps_mat %>% filter(win %in% map_bon_peak$win)
snps_cmd_peak <- snps_mat %>% filter(win %in% cmd_bon_peak$win)

#Add chr_snp column
snps_mat_peak <- snps_mat_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_map_peak <- snps_map_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_cmd_peak <- snps_cmd_peak %>% unite(chr_snp,"chr","snp",sep="_")




#View empirical p-value plots using plotly

#CMD
wza_empri_cmd <- ggplot(data = wza_win_cmd, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(0,410))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
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
#wza_empri_cmd
ggplotly(wza_empri_cmd)

#MAP
wza_empri_map <- ggplot(data = wza_win_map, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(0,410))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
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
#wza_empri_map
ggplotly(wza_empri_map)

#MAT
wza_empri_mat <- ggplot(data = wza_win_mat, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(341,342))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
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
#wza_empri_mat
ggplotly(wza_empri_mat)


#Export peak snps
write_csv(snps_mat_peak,"data/genomic_data/snps_peak_mat.csv")
write_csv(snps_map_peak,"data/genomic_data/snps_peak_map.csv")
write_csv(snps_cmd_peak,"data/genomic_data/snps_peak_cmd.csv")

#Export peak windows
write_csv(mat_bon_peak,"data/genomic_data/peak_window_mat.csv")
write_csv(map_bon_peak,"data/genomic_data/peak_window_map.csv")
write_csv(cmd_bon_peak,"data/genomic_data/peak_window_cmd.csv")






