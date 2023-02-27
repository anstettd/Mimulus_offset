##################################################################################
## Gradient forest data prep
## Author Daniel Anstett
## Filter climate SNPs for baseline preBaypass SNP table 
## Also merge data from 3 climate SNPs into one table
## Generate snpA table for both Bay and full data sets and merge with climate data
## 
## Last Modified September 19, 2022
###################################################################################
##Libraries
library(tidyverse)

##Import Data
climate <- read_csv("Genomics_scripts/Data/climate_pop.csv")

#Import peak SNP data for each environmental variable
env1_peakbf2 <- read_csv("Genomics_scripts/Data/win_bf_mat30_5.csv")
env2_peakbf2 <- read_csv("Genomics_scripts/Data/win_bf_map30_5.csv")
env5_peakbf2 <- read_csv("Genomics_scripts/Data/win_bf_cmd30_5.csv")

#Merge into 1 dataframe that includes every unique allele
env_merge <- rbind(env1_peakbf2,env2_peakbf2,env5_peakbf2)
env_all <- env_merge[,1]
env_snp <- unique(env_all)

#Import full snp table for baseline
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Select BayPass SNPs from merged dataframe
loci_env <- loci_snp %>% filter (chr_snp %in% as.character(env_snp$chr_snp))

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
#write_csv(snp_clim_bayNA, "Genomics_scripts/Data/snp_clim_peakbf5NA.csv")

#Remove NA's from all 55 populations
snp_clim_bay_noNA <- snp_clim_bayNA %>% select_if(~ !any(is.na(.)))
#write_csv(snp_clim_bay_noNA, "Genomics_scripts/Data/snp_clim_peakbf5_noNA.csv")

294-(dim(snp_clim_bayNA %>% select_if(~ !any(is.na(.))))[2]-14) #number of SNPs lost

#Final dataset has XX populations removed to increase number of SNPs going into gradient forest
#See below for reasoning


#Remove NA's from all 55 populations

snp_clim_bayless <- snp_clim_bayNA %>% filter(Paper_ID!= & Paper_ID!= & Paper_ID!=)
294-(dim(snp_clim_bayless %>% select_if(~ !any(is.na(.))))[2]-14) #number of SNPs lost

snp_clim_bayless_noNA <- snp_clim_bayless %>% select_if(~ !any(is.na(.)))
#write_csv(snp_clim_bayless_noNA, "Genomics_scripts/Data/snp_clim_peakbf5_noNA_less.csv")




###############################################################################################
#Count NA's in each row and column

#Tally NAs in each row
row_nas <- apply(snp_clim_bayNA, 1, function(row) sum(is.na(row)))
row_nas <- row_nas[order(-row_nas)]
print(row_nas)

# tally NAs for each column
col_nas <- apply(snp_clim_bayNA, 2, function(col) sum(is.na(col)))
col_nas <- col_nas[order(-col_nas)]
print(col_nas)

# create histogram of number of NAs per row
ggplot(data.frame(nas = row_nas), aes(x = nas)) + 
  geom_histogram() +
  labs(title = "Number of NAs per Pop", x = "Number of NAs", y = "Frequency")

# create histogram of number of NAs per column
ggplot(data.frame(nas = col_nas), aes(x = nas)) + 
  geom_histogram() +
  labs(title = "Number of NAs per SNP", x = "Number of NAs", y = "Frequency")

###############################################################################################
##Possible populations to remove
#confirm populations have only one sampled individual
print(row_nas)
snp_clim_bayNA[51,15:40] #has 57 NA's, 1 sequenced individual
snp_clim_bayNA[54,15:40] #has 31 NA's, 1 sequenced individual
snp_clim_bayNA[31,15:40] #has 29 NA's, 1 sequenced individual
snp_clim_bayNA[20,15:40] #has 28 NA's, 2 sequenced individuals
snp_clim_bayNA[28,15:40] #has 28 NA's, 1 sequenced individual

snp_clim_bayNA[18,15:40] #has 19 NA's, 3 sequenced individuals
snp_clim_bayNA[19,15:40] #has 19 NA's, 4 sequenced individuals
snp_clim_bayNA[52,15:40] #has 19 NA's, 2 sequenced individuals

#Filter out pops with most NAs, then re-check histogram
snp_clim_bayNA_51 <- snp_clim_bayNA %>% filter(Paper_ID!=51)
294-(dim(snp_clim_bayNA_51 %>% select_if(~ !any(is.na(.))))[2]-14) #number of SNPs lost

col_nas_rm <- apply(snp_clim_bayNA_51, 2, function(col) sum(is.na(col)))
col_nas_rm <- col_nas_rm[order(-col_nas_rm)]
ggplot(data.frame(nas = col_nas_rm), aes(x = nas)) + 
  geom_histogram() +
  labs(title = "Number of NAs per Column", x = "Number of NAs", y = "Frequency")

snp_clim_bayNA_54 <- snp_clim_bayNA_51 %>% filter(Paper_ID!=54)
294-(dim(snp_clim_bayNA_54 %>% select_if(~ !any(is.na(.))))[2]-14) #number of SNPs lost

snp_clim_bayNA_31 <- snp_clim_bayNA_54 %>% filter(Paper_ID!=31)
294-(dim(snp_clim_bayNA_31 %>% select_if(~ !any(is.na(.))))[2]-14) #number of SNPs lost

snp_clim_bayNA_20 <- snp_clim_bayNA_31 %>% filter(Paper_ID!=20)
294-(dim(snp_clim_bayNA_20 %>% select_if(~ !any(is.na(.))))[2]-14) #number of SNPs lost

snp_clim_bayNA_28 <- snp_clim_bayNA_20 %>% filter(Paper_ID!=28)
294-(dim(snp_clim_bayNA_28 %>% select_if(~ !any(is.na(.))))[2]-14) #number of SNPs lost

col_nas_rm <- apply(snp_clim_bayNA_28, 2, function(col) sum(is.na(col)))
col_nas_rm <- col_nas_rm[order(-col_nas_rm)]
ggplot(data.frame(nas = col_nas_rm), aes(x = nas)) + 
  geom_histogram() +
  labs(title = "Number of NAs per Column", x = "Number of NAs", y = "Frequency")

snp_clim_bayNA_3more <- snp_clim_bayNA_28 %>% filter(Paper_ID!=18 & Paper_ID!=19 & Paper_ID!=52)
294-(dim(snp_clim_bayNA_3more %>% select_if(~ !any(is.na(.))))[2]-14) #number of SNPs lost




