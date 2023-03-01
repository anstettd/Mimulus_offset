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

#snps_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_mat.csv")
#snps_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_map.csv")
#snps_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_cmd.csv")

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
#Calc empirical p-value based on WZA scores

#MAT
WZA_mean_mat <- mean(WZA_df_mat$Z_pVal,na.rm = T) ## Calculate the mean and sd of the WZA distribution
WZA_sd_mat <- sd(WZA_df_mat$Z_pVal,na.rm = T)
WZA_df_mat$approx_p <- 2*pnorm(-abs(WZA_df_mat$Z_pVal), mean = WZA_mean_mat, sd= WZA_sd_mat) ## Calculate an approximate p-value based on the assumption of normality

#MAP
WZA_mean_map <- mean(WZA_df_map$Z_pVal,na.rm = T) ## Calculate the mean and sd of the WZA distribution
WZA_sd_map <- sd(WZA_df_map$Z_pVal,na.rm = T)
WZA_df_map$approx_p <- 2*pnorm(-abs(WZA_df_map$Z_pVal), mean = WZA_mean_map, sd= WZA_sd_map) ## Calculate an approximate p-value based on the assumption of normality

#CMD
WZA_mean_cmd <- mean(WZA_df_cmd$Z_pVal,na.rm = T) ## Calculate the mean and sd of the WZA distribution
WZA_sd_cmd <- sd(WZA_df_cmd$Z_pVal,na.rm = T)
WZA_df_cmd$approx_p <- 2*pnorm(-abs(WZA_df_cmd$Z_pVal), mean = WZA_mean_cmd, sd= WZA_sd_cmd) ## Calculate an approximate p-value based on the assumption of normality

write_csv(WZA_df_mat, "data/genomic_data/WZA_win_mat.csv")
write_csv(WZA_df_map, "data/genomic_data/WZA_win_map.csv")
write_csv(WZA_df_cmd, "data/genomic_data/WZA_win_cmd.csv")





