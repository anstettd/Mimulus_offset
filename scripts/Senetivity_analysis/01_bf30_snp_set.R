#############################################################################################################
## Generate BF>30 only snpset and frequency table for gradient forest
##
## 
## 
## Last Modified Mat 2, 2023
#############################################################################################################
#Import libraries
library(tidyverse)

#Import all BF
#Import snp env associations (Baseline)
env1 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
env2 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_2_trim.tsv",header=F, sep=" ")
env3 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_3_trim.tsv",header=F, sep=" ")
env4 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_4_trim.tsv",header=F, sep=" ")
env5 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_5_trim.tsv",header=F, sep=" ")

env6 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_6_trim.tsv",header=F, sep=" ")
env7 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_7_trim.tsv",header=F, sep=" ")
env8 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_8_trim.tsv",header=F, sep=" ")
env9 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_9_trim.tsv",header=F, sep=" ")

#Name Columns
colnames(env1) <- c("Chromosome","SNP","Env","BF")
colnames(env2) <- c("Chromosome","SNP","Env","BF")
colnames(env3) <- c("Chromosome","SNP","Env","BF")
colnames(env4) <- c("Chromosome","SNP","Env","BF")
colnames(env5) <- c("Chromosome","SNP","Env","BF")
colnames(env6) <- c("Chromosome","SNP","Env","BF")
colnames(env7) <- c("Chromosome","SNP","Env","BF")
colnames(env8) <- c("Chromosome","SNP","Env","BF")
colnames(env9) <- c("Chromosome","SNP","Env","BF")

env1_allBF <- env1 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env2_allBF <- env2 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env3_allBF <- env3 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env4_allBF <- env4 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env5_allBF <- env5 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env6_allBF <- env6 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env7_allBF <- env7 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env8_allBF <- env8 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env9_allBF <- env9 %>% unite(chr_snp,"Chromosome","SNP",sep="_")

#######################################################################################################

#BF>30 for all windows
env1_united_bf30 <- env1_allBF %>% filter(BF>30)
env2_united_bf30 <- env2_allBF %>% filter(BF>30)
env3_united_bf30 <- env3_allBF %>% filter(BF>30)
env4_united_bf30 <- env4_allBF %>% filter(BF>30)
env5_united_bf30 <- env5_allBF %>% filter(BF>30)
env6_united_bf30 <- env6_allBF %>% filter(BF>30)
env7_united_bf30 <- env7_allBF %>% filter(BF>30)
env8_united_bf30 <- env8_allBF %>% filter(BF>30)
env9_united_bf30 <- env9_allBF %>% filter(BF>30)

#Remove large files
rm(env1)
rm(env2)
rm(env3)
rm(env4)
rm(env5)
rm(env6)
rm(env7)
rm(env8)
rm(env9)

rm(env1_allBF)
rm(env2_allBF)
rm(env3_allBF)
rm(env4_allBF)
rm(env5_allBF)
rm(env6_allBF)
rm(env7_allBF)
rm(env8_allBF)
rm(env9_allBF)

#Get unique BF>30 snps
snps_bf30_all <- rbind(env1_united_bf30,
                       env2_united_bf30,
                       env3_united_bf30,
                       env4_united_bf30,
                       env5_united_bf30,
                       env6_united_bf30,
                       env7_united_bf30,
                       env8_united_bf30,
                       env9_united_bf30)
snps_bf30_all_u <- as.data.frame(unique(snps_bf30_all$chr_snp))
colnames(snps_bf30_all_u) <- "chr_snp"

###################################################################################################
##Set up frequency table

#Import full snp table for baseline
climate <- read_csv("data/genomic_data/climate_pop.csv")
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Select BayPass SNPs from merged dataframe
loci_env <- loci_snp %>% filter (chr_snp %in% as.character(snps_bf30_all_u$chr_snp))

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
write_csv(snp_clim_bay_noNA, "data/genomic_data/bf30only_snp_set.csv")


