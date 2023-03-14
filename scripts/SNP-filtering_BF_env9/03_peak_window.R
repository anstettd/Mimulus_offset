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

#Windows
#wza_win_mat <- read_csv("data/genomic_data/WZA_win_mat_bf.csv")

wza_win_env1 <- read_csv("data/genomic_data/WZA_win_env1_bf.csv")
wza_win_env2 <- read_csv("data/genomic_data/WZA_win_env2_bf.csv")
wza_win_env3 <- read_csv("data/genomic_data/WZA_win_env3_bf.csv")
wza_win_env4 <- read_csv("data/genomic_data/WZA_win_env4_bf.csv")
wza_win_env5 <- read_csv("data/genomic_data/WZA_win_env5_bf.csv")

wza_win_env6 <- read_csv("data/genomic_data/WZA_win_env6_bf.csv")
wza_win_env7 <- read_csv("data/genomic_data/WZA_win_env7_bf.csv")
wza_win_env8 <- read_csv("data/genomic_data/WZA_win_env8_bf.csv")
wza_win_env9 <- read_csv("data/genomic_data/WZA_win_env9_bf.csv")

#SNPs
#snps_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_mat_bf.csv")

snps_env1 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env1_bf.csv")
snps_env2 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env2_bf.csv")
snps_env3 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env3_bf.csv")
snps_env4 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env4_bf.csv")
snps_env5 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env5_bf.csv")

snps_env6 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env6_bf.csv")
snps_env7 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env7_bf.csv")
snps_env8 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env8_bf.csv")
snps_env9 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env9_bf.csv")


#Filter by Bonferonnii correction alpha critical = 1.29423e-06, aka 5.887988 sigma
#mat_bon <- wza_win_mat %>% filter(Z_pVal<1.294264e-06)
env1_bon <- wza_win_env1 %>% filter(Z_pVal<1.294264e-06) #21 windows
env2_bon <- wza_win_env2 %>% filter(Z_pVal<1.294264e-06) #16 windows
env3_bon <- wza_win_env3 %>% filter(Z_pVal<1.294264e-06) #32 windows
env4_bon <- wza_win_env4 %>% filter(Z_pVal<1.294264e-06) #28 windows
env5_bon <- wza_win_env5 %>% filter(Z_pVal<1.294264e-06) #12 windows

env6_bon <- wza_win_env6 %>% filter(Z_pVal<1.294264e-06) #12 windows
env7_bon <- wza_win_env7 %>% filter(Z_pVal<1.294264e-06) #27 windows
env8_bon <- wza_win_env8 %>% filter(Z_pVal<1.294264e-06) #18 windows
env9_bon <- wza_win_env9 %>% filter(Z_pVal<1.294264e-06) #19 windows

#Get peak windows only
#Filter out windows that are not peak windows
#mat_bon_peak <- mat_bon %>% filter(!win %in% c(6458,6459)) # 19 windows

env1_bon_peak <- env1_bon %>% filter(!win %in% c(6458,6459)) # 19 windows
env2_bon_peak <- evn2_bon %>% filter(!win %in% c(3881,3882)) # 14 windows
env_bon_peak <- evn_bon %>% filter(!win %in% c()) # windows
env_bon_peak <- evn_bon %>% filter(!win %in% c()) # windows
env5_bon_peak <- env5_bon # 12 windows

env_bon_peak <- evn_bon %>% filter(!win %in% c()) # windows
env_bon_peak <- evn_bon %>% filter(!win %in% c()) # windows
env_bon_peak <- evn_bon %>% filter(!win %in% c()) # windows
env_bon_peak <- evn_bon %>% filter(!win %in% c()) # windows




#Filter SNPs for peak
snps_mat_peak <- snps_mat %>% filter(win %in% mat_bon_peak$win)
snps_map_peak <- snps_mat %>% filter(win %in% map_bon_peak$win)
snps_cmd_peak <- snps_mat %>% filter(win %in% cmd_bon_peak$win)

#Add chr_snp column
snps_mat_peak <- snps_mat_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_map_peak <- snps_map_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_cmd_peak <- snps_cmd_peak %>% unite(chr_snp,"chr","snp",sep="_")



#Export peak snps
#write_csv(snps_mat_peak,"data/genomic_data/snps_peak_mat_bf.csv")
#write_csv(snps_map_peak,"data/genomic_data/snps_peak_map_bf.csv")
#write_csv(snps_cmd_peak,"data/genomic_data/snps_peak_cmd_bf.csv")

#Export peak windows
#write_csv(mat_bon_peak,"data/genomic_data/peak_window_mat_bf.csv")
#write_csv(map_bon_peak,"data/genomic_data/peak_window_map_bf.csv")
#write_csv(cmd_bon_peak,"data/genomic_data/peak_window_cmd_bf.csv")


#View empirical p-value plots using plotly

#Annual

#ENV 1 - MAT
wza_empri_mat <- ggplot(data = wza_win_env1, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(64.5,64.7))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env1)[1])), col = "red", lty = 2, lwd = 1)+
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

#ENV2 - MAP
wza_empri_map <- ggplot(data = wza_win_env2, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(353,354))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env2)[1])), col = "red", lty = 2, lwd = 1)+
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

#ENV3
wza_empri_env3 <- ggplot(data = wza_win_env3, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(0,50))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env3)[1])), col = "red", lty = 2, lwd = 1)+
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
ggplotly(wza_empri_env3)

#ENV4
wza_empri_env4 <- ggplot(data = wza_win_env4, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(0,410))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env4)[1])), col = "red", lty = 2, lwd = 1)+
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
ggplotly(wza_empri_env4)


#ENV 5 - CMD
wza_empri_cmd <- ggplot(data = wza_win_env5, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(0,410))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env5)[1])), col = "red", lty = 2, lwd = 1)+
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


#Seasonal
#ENV6
wza_empri_env6 <- ggplot(data = wza_win_env6, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(0,410))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env6)[1])), col = "red", lty = 2, lwd = 1)+
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
ggplotly(wza_empri_env6)

#ENV7
wza_empri_env7 <- ggplot(data = wza_win_env7, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(0,410))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env7)[1])), col = "red", lty = 2, lwd = 1)+
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
ggplotly(wza_empri_env7)

#ENV8
wza_empri_env8 <- ggplot(data = wza_win_env8, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(0,410))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env8)[1])), col = "red", lty = 2, lwd = 1)+
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
ggplotly(wza_empri_env8)

#ENV9
wza_empri_env9 <- ggplot(data = wza_win_env9, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(0,410))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env9)[1])), col = "red", lty = 2, lwd = 1)+
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
ggplotly(wza_empri_env9)
















