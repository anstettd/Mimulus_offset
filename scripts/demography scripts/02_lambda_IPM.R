#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Create data frame of vital rate parameters and build integral projection models to obtain estimates of annual lambdas for each population
#### AUTHOR: Seema Sheth and Amy Angert
#### DATE LAST MODIFIED: 20230207


#*******************************************************************************
#### 0. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Make vector of packages needed
packages_needed <- c("lme4", "glmmTMB", "tidyverse")
require(lme4)
require(glmmTMB)
require(plyr)
require(dplyr)
require(ggplot2)
require(tidyr)
require(gridExtra)

# Install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# Load packages needed
for (i in 1:length(packages_needed)){
  library(packages_needed[i], character.only = TRUE)
}


#*******************************************************************************
#### 2. Read in vital rate data frames ###
#*******************************************************************************

data <- read.csv("data/demography data/Mcard_demog_data_2010-2016_cleanindivs.csv")
data$Site = factor(data$Site)
data$Year = factor(data$Year)
data$Year = factor(data$Year)
data$SiteYear = factor(data$SiteYear)

site_fruit_count_data <- read.csv("data/demography data/Mcard_demog_data_2010-2016_seedinput.csv")
site_fruit_count_data$Site = factor(site_fruit_count_data$Site)
site_fruit_count_data$Year = factor(site_fruit_count_data$Year)
site_fruit_count_data$Year = factor(site_fruit_count_data$Year)
site_fruit_count_data$SiteYear = factor(site_fruit_count_data$SiteYear)

#*******************************************************************************
#### 1. Create global survival, growth and fecundity models using data from all sites ###
#*******************************************************************************

# Create a vector of unique Site x Year for subsetting; note this is sorted by decreasing latitude 
SiteYear <- unique(data$SiteYear)

# Set up data frame of model parameters
params=c()
growth=c()
fruit=c()

  #*******************************************************************************
  ### 1A. Survival ###
  #*******************************************************************************

  # Read in top survival model output (Formula: Surv ~ logSize + (logSize | Year/Site))
  surv.reg <- load("data/demography data/surv.reg.rda")

  # Store model coefficients
  params$SiteYear=rownames(coefficients(s3)$'Site:Year')
  params$surv.int=coefficients(s3)$'Site:Year'[,1] 
  params$surv.slope=coefficients(s3)$'Site:Year'[,2] 

  
  #*******************************************************************************
  ### 1B. Growth ###
  #*******************************************************************************
  
  # Read in top growth model output (Formula: logSizeNext ~ logSize + (logSize|Year/Site))
  growth.reg <- load("data/demography data/growth.reg.rda")
  
  # Store model coefficients
  growth$SiteYear=rownames(coefficients(g3)$'Site:Year')
  growth$growth.int=coefficients(g3)$'Site:Year'[,1] 
  growth$growth.slope=coefficients(g3)$'Site:Year'[,2] 
  growth$growth.sd=rep(sigma(g3), times=length(coefficients(g3)$'Site:Year'[,2])) 
  # TO DO: Address Seema's comment about growth.sd: "WARNING! I'M UNCERTAIN THAT THIS IS BEST METHOD!"
  
  # make a data frame
  growth=data.frame(growth)
  
  #*******************************************************************************
  ### 1C. Flowering ###
  #*******************************************************************************
  
  # Read in top flowering model output (Formula: Fec0 ~ logSize + (logSize|Year/Site))
  flowering.reg=load("data/demography data/flowering.reg.rda")

  # Store model coefficients
  params$flowering.int=coefficients(fl3)$'Site:Year'[,1] 
  params$flowering.slope=coefficients(fl3)$'Site:Year'[,2] 
  
  #*******************************************************************************
  ### 1D. Fruit number (untransformed) using negative binomial regression ###
  #*******************************************************************************
  
  # Read in top model output for fruit.reg (Formula: Fec1 ~ logSize + (logSize|Year/Site))   
  fruit.reg=load("data/demography data/fruit.reg.rda")

  # Store model coefficients (fr3 from glmmTMB)
  fruit$SiteYear=rownames(ranef(fr3)$cond$'Site:Year')
  fruit$fruit.int=fixef(fr3)$cond[1]+ranef(fr3)$cond$'Site:Year'[,1] 
  fruit$fruit.slope=fixef(fr3)$cond[2]+ranef(fr3)$cond$'Site:Year'[,2] 
  
  # make data frame
  fruit=data.frame(fruit) 
  
  #*******************************************************************************
  ### 1E. Size distribution of recruits ###
  #*******************************************************************************
  
  # Calculate mean size in year t+1 of plants that were new that year (unmeasured at time t)
  recruit.size.mean <- tapply(data$logSizeNext[is.na(data$logSize)], data$SiteYear[is.na(data$logSize)], FUN="mean") %>% data.frame()
  
  # Calculate st dev of size in year t+1 of plants that were new that year (unmeasured at time t)
  recruit.size.sd <- tapply(data$logSizeNext[is.na(data$logSize)], data$SiteYear[is.na(data$logSize)], FUN="sd") %>% data.frame()
  
  # Assemble mean and sd into same data frame
  recruit_size <- data.frame(rownames(recruit.size.mean), recruit.size.mean, recruit.size.sd)
  recruit_size <- rename(recruit_size, SiteYear="rownames.recruit.size.mean.", recruit.logSize.mean=".",recruit.logSize.sd="..1")
  
  # Assign 0 values to sites with no recruitment (mean & SD = NA) or only 1 recruit (mean = value but SD= A)
  recruit_size$recruit.logSize.mean[is.na(recruit_size$recruit.logSize.mean)]=0 
  recruit_size$recruit.logSize.sd[is.na(recruit_size$recruit.logSize.sd)]=0 
  # TO DO: Address Seema's comment about assigning 0s: "WARNING! MAY NEED TO MODIFY THIS!"

  #*******************************************************************************
  ### 1F. Create data frame of site-specific parameter estimates and join all estimates ###
  #*******************************************************************************
  
  params <- data.frame(params)
  params <- full_join(params,growth) %>% full_join(fruit) %>% full_join(recruit_size)
  
  #*******************************************************************************
  ### 1G. Number of seeds per fruit ###
  #*******************************************************************************
  
  # Note: these are site-specific, but not site-by-year-specific, values
  
  # obtain mean seed counts per fruit per site
  seeds.per.site <- data %>% 
    group_by(Site) %>% 
    summarize(SeedCt = mean(SeedCt, na.rm=T)) 
  
  #seeds.per.site = tapply(data$SeedCt, data$SiteYear, FUN=mean, na.rm=T) 
  
  # make into a data frame
  seeds.per.site <- data.frame(seeds.per.site, rownames(seeds.per.site)) 
  
  # define column names for merging
  colnames(seeds.per.site)=c("seed.ct","SiteYear") 
  
  # replace missing values with NA
  seeds.per.site$seed.ct[seeds.per.site$seed.ct=="Inf"]=NA
  
  # join seeds per fruit to other parameters
  params <- left_join(params, seeds.per.site)
                      
  #*******************************************************************************
  ### 1H. Establishment probability ###
  #*******************************************************************************
  
  # Obtain number of new recruits per site at year = t+1
  recruit.number <- data %>% 
    dplyr::filter(is.na(logSize)) %>% 
    group_by(SiteYear) %>% 
    summarise(recruit.number = n())
    
  # Obtain total fruit count per site at year = t
  fruits.per.site <- site_fruit_count_data %>% 
    filter(!is.na(Fec1)) %>% 
    group_by(SiteYear) %>% 
    summarise(fruits.per.site = sum(Fec1))
  
  # Join recruit # and fruits per site into one frame
  # Obtain total seed count per site (= # fruits per site * # seeds per fruit per site) at year = t
  total.seeds.per.site <- full_join(recruit.number, fruits.per.site) %>% full_join(seeds.per.site)
    mutate(total.seeds.per.site = )
  
  # 
  total.seeds.per.site=fruits.per.site*seeds.per.site$seed.ct	
  
  # Estimate establishment probability as # of new recruits at year = t+1/# of seeds at year = t
  params$establishment.prob=recruit.number$recruit.number/total.seeds.per.site
  
  # Set establishment probability as 0 for Hauser Creek (was calculated as NA because Hauser creek has 0 new recruits)
  params$establishment.prob[is.na(params$establishment.prob)]=0	

  # Add separate columns for year and site
  params=params %>% separate(SiteYear,c("Site","Year"),":",remove=FALSE)
  
  
#*******************************************************************************
### 2. Create site-specific IPMs parameterized by site-specific parameters derived from global vital rates models 
#*******************************************************************************

  #*******************************************************************************
  ### 2A. Subset data for site f
  #*******************************************************************************
  
  # remove rows of parameters with NA
  params=params[complete.cases(params),]
  
  # Store parameters in .csv file for later use
  #write.csv(params,"R_output/vital_rate_coefficients.csv",row.names=FALSE)
  
  # remove site x year combinations without parameter estimates
  siteYear=siteYear[siteYear %in% params$SiteYear]
  
  # create empty vectors for lambda and site to be filled
  lambda=c()
  SiteYear=character()
  
  for (f in 1:length(siteYear)) {
    data1=subset(data,SiteYear==siteYear[f])
    params1=subset(params,SiteYear==siteYear[f])
    params1=subset(params1,select=-SiteYear)
    
    #*******************************************************************************
    ### 2B. Create survival, growth, and fecundity functions and build IPM by running integral_projection_model.R script
    #*******************************************************************************
    
    source("R_scripts/integral_projection_model.R")
    
    #*******************************************************************************
    ### 2C. Obtain lambda estimate for site f
    #*******************************************************************************
    
    lambda[f] <- Re(eigen(K)$values[1])
    SiteYear[f]=as.character(siteYear[f])
    } # end loop to run IPMs and estimate lambdas for each site
    
    # make data frame of site and lambda
    siteYear.lambda=data.frame(SiteYear,lambda)
    
#*******************************************************************************
### 3. Merge site information with lambda estimates and save to .csv file
#*******************************************************************************

# Create data frame of Site, Latitude, Longitude, Region, and Elevation for hypothesis testing
site.info=subset(data,select=c(Site,Year,SiteYear,Latitude,Longitude,Elevation,Region,RegionRank)) %>% unique() %>% arrange(-Latitude)
    
# merge site info with lambda estimates
site.info=join(site.info,siteYear.lambda)

# Manually add lambda = 0 for Hauser Creek in 2012 & 2013
Hauser=filter(site.info,Site=="Hauser Creek"&Year==2011)
Hauser$Year=2012 %>% factor()
Hauser$SiteYear="Hauser Creek:2012" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info,Hauser) 
Hauser$Year=2013 %>% factor()
Hauser$SiteYear="Hauser Creek:2013" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info,Hauser) %>% mutate(Year=factor(Year),SiteYear=factor(SiteYear))

# view final data frame
str(site.info)

# save to .csv file 
write.csv(site.info,"R_output/siteYear.lambda_2010-2016.csv",row.names=FALSE)


