#############################################################################################################
## Format snps for offset sensetivity analysis
##
## 
## 
## Last Modified Mat 2, 2023
#############################################################################################################
#Import libraries
library(tidyverse)

#Import files
snps_peak_mat <- read_csv("data/genomic_data/snps_peak_mat_oldWZA.csv")
snps_peak_map <- read_csv("data/genomic_data/snps_peak_map_oldWZA.csv")
snps_peak_cmd <- read_csv("data/genomic_data/snps_peak_cmd_oldWZA.csv")

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

#Filter Bayes factor by peak windows
snp_mat_peakbf <- env1_united %>% filter(chr_snp %in% as.character(snps_peak_mat$chr_snp))
snp_map_peakbf <- env2_united %>% filter(chr_snp %in% as.character(snps_peak_map$chr_snp))
snp_cmd_peakbf <- env5_united %>% filter(chr_snp %in% as.character(snps_peak_cmd$chr_snp))

#Filter peak windows by BF >5
snp_env1_peakbf5 <- snp_mat_peakbf %>% filter(BF>5)
snp_env2_peakbf5 <- snp_map_peakbf %>% filter(BF>5)
snp_env5_peakbf5 <- snp_cmd_peakbf %>% filter(BF>5)


#Get unique snp set
snps_old_all <- rbind(snp_env1_peakbf5,
                      snp_env2_peakbf5,
                      snp_env5_peakbf5)
snps_old_all_u <- as.data.frame(unique(snps_old_all$chr_snp))
colnames(snps_old_all_u) <- "chr_snp"

#Import full snp table for baseline
climate <- read_csv("data/genomic_data/climate_pop.csv")
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Select BayPass SNPs from merged dataframe
loci_env <- loci_snp %>% filter (chr_snp %in% as.character(snps_old_all_u$chr_snp))

#Generate snpA table for BayPass SNPs
snp_bay<-data.frame()
counter<-1
for (i in seq (2,dim(loci_env)[2]-1,2)){
  for(j in 1:dim(loci_env)[1]){
    tmp_total<-as.numeric(loci_env[j,i]) + as.numeric(loci_env[j,i+1])
    snp_bay[j,counter]<-as.numeric(loci_env[j,i])/tmp_total
  }
  counter<-counter+1
}
rownames(snp_bay)<- loci_env$chr_snp
snp_bay_T <- as.data.frame(t(snp_bay))
snp_clim_bayNA <- cbind(climate,snp_bay_T)

snp_clim_bay_noNA <- snp_clim_bayNA %>% select_if(~ !any(is.na(.)))
write_csv(snp_clim_bay_noNA, "data/genomic_data/old_WZA_snp_set.csv")

