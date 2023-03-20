#############################################################################################################
## Calc WZA for windows
## Get empirical p-value
## Author Daniel Anstett
## 
## 
## Modified from Tom Booker WZA Vignette
## Last Modified August 3, 2022
#############################################################################################################
#Import libraries

library(tidyverse)
#library(Kendall)

#Import files

snps_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_mat.csv")
snps_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_map.csv")
snps_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_cmd.csv")

###########################################################################################################

#Filter needed variables for WZA
#snps_mat_input <- snps_mat %>% select(empirical_p,win,q_bar)
#snps_map_input <- snps_map %>% select(empirical_p,win,q_bar)
#snps_cmd_input <- snps_cmd %>% select(empirical_p,win,q_bar)

#Too large to store on github. Store locally
#write_csv(snps_mat_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_mat_input.csv")
#write_csv(snps_map_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_map_input.csv")     
#write_csv(snps_cmd_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_cmd_input.csv")   

#Run python script run_WZA.txt
#Import WZA scores
WZA_df_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/mat_WZA.csv")
WZA_df_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/map_WZA.csv")
WZA_df_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/cmd_WZA.csv")

###########################################################################################################
#Re-name windows
names(WZA_df_mat)[names(WZA_df_mat) == 'gene'] <- 'win'
names(WZA_df_map)[names(WZA_df_map) == 'gene'] <- 'win'
names(WZA_df_cmd)[names(WZA_df_cmd) == 'gene'] <- 'win'

#Make DF that has chr ID for each window
chr_mat <- snps_mat %>% select(chr,win) %>% distinct()
chr_map <- snps_map %>% select(chr,win) %>% distinct()
chr_cmd <- snps_cmd %>% select(chr,win) %>% distinct()

#Get chromosome information into WZA score
WZA_df_mat_chr <- left_join(WZA_df_mat,chr_mat,by="win")
WZA_df_map_chr <- left_join(WZA_df_map,chr_map,by="win")
WZA_df_cmd_chr <- left_join(WZA_df_cmd,chr_cmd,by="win")

#Input position for each window
WZA_df_mat_chr$pos <- (WZA_df_mat_chr$win+0.5)*10000
WZA_df_map_chr$pos <- (WZA_df_map_chr$win+0.5)*10000
WZA_df_cmd_chr$pos <- (WZA_df_cmd_chr$win+0.5)*10000


write_csv(WZA_df_mat_chr, "data/genomic_data/WZA_win_mat.csv")
write_csv(WZA_df_map_chr, "data/genomic_data/WZA_win_map.csv")
write_csv(WZA_df_cmd_chr, "data/genomic_data/WZA_win_cmd.csv")





