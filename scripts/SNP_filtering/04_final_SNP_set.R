#############################################################################################################
## Generate Final SNP Set
## BF > 30
## Peak Boferroni window, BF >5
## 
## 
## Last Modified Sept 28, 2022
#############################################################################################################
#Import libraries
library(tidyverse)

#Import Large Loci Win
loci_win <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/loci_win.csv")                                                                                 
loci_win <- loci_win %>% unite(col="chr_snp", c("chr","snp"), sep="_")

#Import files
snps_peak_mat <- read_csv("data/genomic_data/snps_peak_mat.csv")
snps_peak_map <- read_csv("data/genomic_data/snps_peak_map.csv")
snps_peak_cmd <- read_csv("data/genomic_data/snps_peak_cmd.csv")

#Import snp env associations (Baseline)
env1 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
env2 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_2_trim.tsv",header=F, sep=" ")
env5 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_5_trim.tsv",header=F, sep=" ")

#Name Columns
colnames(env1) <- c("Chromosome","SNP","Env","BF")
colnames(env2) <- c("Chromosome","SNP","Env","BF")
colnames(env5) <- c("Chromosome","SNP","Env","BF")

env1_united <- env1 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env2_united <- env2 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env5_united <- env5 %>% unite(chr_snp,"Chromosome","SNP",sep="_")

#############################################################################################################
#Investigate how many SNPs are left
#Filter Bayes factor by peak windows
snp_mat_peakbf <- env1_united %>% #filter(BF<=30) %>% 
  filter(chr_snp %in% as.character(snps_peak_mat$chr_snp))
snp_map_peakbf <- env2_united %>% #filter(BF<=30) %>% 
  filter(chr_snp %in% as.character(snps_peak_map$chr_snp))
snp_cmd_peakbf <- env5_united %>% #filter(BF<=30) %>% 
  filter(chr_snp %in% as.character(snps_peak_cmd$chr_snp))

#Get window IDs for all selected SNPs
snp_mat_peakbf_win <- left_join(snp_mat_peakbf,loci_win, by="chr_snp")
snp_map_peakbf_win <- left_join(snp_map_peakbf,loci_win, by="chr_snp")
snp_cmd_peakbf_win <- left_join(snp_cmd_peakbf,loci_win, by="chr_snp")

colnames(snp_mat_peakbf_win)[5] <- "win"
colnames(snp_map_peakbf_win)[5] <- "win"
colnames(snp_cmd_peakbf_win)[5] <- "win"

#How many windowed SNPs clear BF thresholds
snp_mat_peakbf30_win <- snp_mat_peakbf_win %>% filter(BF>30) 
snp_map_peakbf30_win <- snp_map_peakbf_win %>% filter(BF>30)
snp_cmd_peakbf30_win <- snp_cmd_peakbf_win %>% filter(BF>30)

snp_mat_peakbf20_win <- snp_mat_peakbf_win %>% filter(BF>20) 
snp_map_peakbf20_win <- snp_map_peakbf_win %>% filter(BF>20)
snp_cmd_peakbf20_win <- snp_cmd_peakbf_win %>% filter(BF>20)

snp_mat_peakbf10_win <- snp_mat_peakbf_win %>% filter(BF>10) 
snp_map_peakbf10_win <- snp_map_peakbf_win %>% filter(BF>10)
snp_cmd_peakbf10_win <- snp_cmd_peakbf_win %>% filter(BF>10)

snp_mat_peakbf5_win <- snp_mat_peakbf_win %>% filter(BF>5) 
snp_map_peakbf5_win <- snp_map_peakbf_win %>% filter(BF>5)
snp_cmd_peakbf5_win <- snp_cmd_peakbf_win %>% filter(BF>5)

snp_mat_peakbf2_win <- snp_mat_peakbf_win %>% filter(BF>2) 
snp_map_peakbf2_win <- snp_map_peakbf_win %>% filter(BF>2)
snp_cmd_peakbf2_win <- snp_cmd_peakbf_win %>% filter(BF>2)

#BF>30 for all windows
env1_united_bf30 <- env1_united %>% filter(BF>30)
env2_united_bf30 <- env2_united %>% filter(BF>30)
env5_united_bf30 <- env5_united %>% filter(BF>30)


#Make SNPs in Window Summary Data Frame
snps_in_win <- data.frame()

snps_in_win[1,1] <- ">30 all win"
snps_in_win[1,2] <- dim(env1_united_bf30 )[1]
snps_in_win[1,3] <- dim(env2_united_bf30)[1]
snps_in_win[1,4] <- dim(env5_united_bf30)[1]

snps_in_win[2,1] <- ">30"
snps_in_win[2,2] <- dim(snp_mat_peakbf30_win)[1]
snps_in_win[2,3] <- dim(snp_map_peakbf30_win)[1]
snps_in_win[2,4] <- dim(snp_cmd_peakbf30_win)[1]

snps_in_win[3,1] <- ">20"
snps_in_win[3,2] <- dim(snp_mat_peakbf20_win)[1]
snps_in_win[3,3] <- dim(snp_map_peakbf20_win)[1]
snps_in_win[3,4] <- dim(snp_cmd_peakbf20_win)[1]

snps_in_win[4,1] <- ">10"
snps_in_win[4,2] <- dim(snp_mat_peakbf10_win)[1]
snps_in_win[4,3] <- dim(snp_map_peakbf10_win)[1]
snps_in_win[4,4] <- dim(snp_cmd_peakbf10_win)[1]

snps_in_win[5,1] <- ">5"
snps_in_win[5,2] <- dim(snp_mat_peakbf5_win)[1]
snps_in_win[5,3] <- dim(snp_map_peakbf5_win)[1]
snps_in_win[5,4] <- dim(snp_cmd_peakbf5_win)[1]

snps_in_win[6,1] <- ">2"
snps_in_win[6,2] <- dim(snp_mat_peakbf2_win)[1]
snps_in_win[6,3] <- dim(snp_map_peakbf2_win)[1]
snps_in_win[6,4] <- dim(snp_cmd_peakbf2_win)[1]

colnames(snps_in_win) <- c("log10BF","MAT SNPs", "MAP SNPs", "CMD SNPs")

snps_in_win


#Summarize how many windows are retained
windows_in <- data.frame()

windows_in[1,1] <- "any"
windows_in[1,2] <- 11
windows_in[1,3] <- 2
windows_in[1,4] <- 9

windows_in[2,1] <- ">5"
windows_in[2,2] <- length(unique(snp_mat_peakbf5_win$win))
windows_in[2,3] <- length(unique(snp_map_peakbf5_win$win))
windows_in[2,4] <- length(unique(snp_cmd_peakbf5_win$win))
  
windows_in[3,1] <- ">2"
windows_in[3,2] <- length(unique(snp_mat_peakbf2_win$win))
windows_in[3,3] <- length(unique(snp_map_peakbf2_win$win))
windows_in[3,4] <- length(unique(snp_cmd_peakbf2_win$win))

colnames(windows_in) <- c("log10BF","MAT Windows","MAP Windows","CMD Windows")
  
windows_in

write_csv(snps_in_win,"data/genomic_data/snps_in_win.csv")
write_csv(windows_in,"data/genomic_data/windows_in.csv")




dim(env1_united_bf30)[1]



#############################################################################################################











#Filter baseline by BF=>30

##Filter Bayes factor by peak windows
env1_united_bf30 <- env1_united %>% filter(BF>30)
env2_united_bf30 <- env2_united %>% filter(BF>30)
env5_united_bf30 <- env5_united %>% filter(BF>30)

#Filter Bayes factor by peak windows
snp_mat_peakbf <- env1_united %>% #filter(BF<=30) %>% 
  filter(chr_snp %in% as.character(snps_peak_mat$chr_snp))
snp_map_peakbf <- env2_united %>% #filter(BF<=30) %>% 
  filter(chr_snp %in% as.character(snps_peak_map$chr_snp))
snp_cmd_peakbf <- env5_united %>% #filter(BF<=30) %>% 
  filter(chr_snp %in% as.character(snps_peak_cmd$chr_snp))

#Filter peak windows by BF >10
snp_mat_peakbf10 <- snp_mat_peakbf %>% filter(BF>5)
snp_map_peakbf10 <- snp_map_peakbf %>% filter(BF>5)
snp_cmd_peakbf10 <- snp_cmd_peakbf %>% filter(BF>5)

#Merge peakbf10 with bf30
snp_mat_peakbf30 <- rbind(snp_mat_peakbf10,env1_united_bf30)
snp_map_peakbf30 <- rbind(snp_map_peakbf10,env2_united_bf30)
snp_cmd_peakbf30 <- rbind(snp_cmd_peakbf10,env5_united_bf30)

#Get window IDs for all selected SNPs
snp_mat_peakbf30_win <- left_join(snp_mat_peakbf30,loci_win, by="chr_snp")
snp_map_peakbf30_win <- left_join(snp_map_peakbf30,loci_win, by="chr_snp")
snp_cmd_peakbf30_win <- left_join(snp_cmd_peakbf30,loci_win, by="chr_snp")

colnames(snp_mat_peakbf30_win)[5] <- "win"
colnames(snp_map_peakbf30_win)[5] <- "win"
colnames(snp_cmd_peakbf30_win)[5] <- "win"


#write_csv(snp_mat_peakbf30_win,"data/genomic_data/win_bf_mat30_5.csv")
#write_csv(snp_map_peakbf30_win,"data/genomic_data/win_bf_map30_5.csv")
#write_csv(snp_cmd_peakbf30_win,"data/genomic_data/win_bf_cmd30_5.csv")


