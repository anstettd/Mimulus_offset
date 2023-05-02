#############################################################################################################
## Format snps for offset sensetivity analysis
##
## 
## 
## Last Modified Mat 2, 2023
#############################################################################################################
#Import libraries
library(tidyverse)


#Import old snp set frequnecies
env1_united <- read_csv("data/genomic_data/win_bf_mat30_5.csv") 
env2_united <- read_csv("data/genomic_data/win_bf_map30_5.csv") 
env5_united <- read_csv("data/genomic_data/win_bf_cmd30_5.csv")


#Get unique snp set
snps_old_all <- rbind(env1_united,
                       env2_united,
                       env5_united)
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
write_csv(snp_clim_bay_noNA, "data/genomic_data/oldonly_snp_set.csv")

