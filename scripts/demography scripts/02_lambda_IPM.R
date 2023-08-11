#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Create data frame of vital rate parameters and build integral projection models to obtain estimates of annual lambdas for each population
#### AUTHOR: Seema Sheth and Amy Angert
#### DATE LAST MODIFIED: 20230809


#*******************************************************************************
#### 0. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Make vector of packages needed
packages_needed <- c("lme4", "glmmTMB", "tidyverse")

# Install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# if there are errors related to glmmTMB, try: install.packages("glmmTMB", type="source")

# Load packages needed
for (i in 1:length(packages_needed)){
  library(packages_needed[i], character.only = TRUE)
}


#*******************************************************************************
#### 2. Read in vital rate data frames ###
#*******************************************************************************

data <- read.csv("data/demography data/Mcard_demog_data_2010-2019_cleanindivs.csv") %>% 
  mutate(Period = ifelse(Year<2015&Year>2011, "drought",
                         ifelse(Year>2014, "recovery", NA))) %>% 
  filter(!is.na(Period))
data$Site = factor(data$Site)
data$Year = factor(data$Year)
data$SiteYear = factor(data$SiteYear)

data.drought <- data %>% filter(Period=="drought")
data.recovery <- data %>% filter(Period=="recovery")

site_fruit_count_data <- read.csv("data/demography data/Mcard_demog_data_2010-2019_seedinput.csv") %>% 
  mutate(Period = ifelse(Year<2015&Year>2011, "drought",
                         ifelse(Year>2014, "recovery", NA))) %>% 
  filter(!is.na(Period))
site_fruit_count_data$Site = factor(site_fruit_count_data$Site)
site_fruit_count_data$Year = factor(site_fruit_count_data$Year)
site_fruit_count_data$SiteYear = factor(site_fruit_count_data$SiteYear)

site_fruit_count_data.drought <- site_fruit_count_data %>% filter(Period=="drought")
site_fruit_count_data.recovery <- site_fruit_count_data %>% filter(Period=="recovery")


#*******************************************************************************
#### 1. Create global survival, growth and fecundity models using data from all sites ###
#*******************************************************************************

# Set up data frames of model parameters
surv.drought=c()
growth.drought=c()
flower.drought=c()
fruit.drought=c()
surv.recovery=c()
growth.recovery=c()
flower.recovery=c()
fruit.recovery=c()


  #*******************************************************************************
  ### 1A. Survival ###
  #*******************************************************************************

  # Read in top survival model output during drought (Formula: Surv ~ logSize + (1|Year) + (logSize|Site))
  surv.drought.reg <- load("data/demography data/surv.drought.reg.rda")

  # Store model coefficients for each Site:Year
  Site=rownames(coefficients(s7.d)$'Site')
  Year=rownames(coefficients(s7.d)$'Year')
  counter=0
  for (i in 1:length(Site)) {
    for (j in 1:length(Year)) {
      counter=counter+1
      surv.drought$SiteYear[counter]=paste(rownames(coefficients(s7.d)$'Site')[i],":",rownames(coefficients(s7.d)$'Year')[j], sep="")
      surv.drought$surv.int[counter]=coefficients(s7.d)$'Site'[i,1] + coefficients(s7.d)$'Year'[j,1]
      surv.drought$surv.slope[counter]=coefficients(s7.d)$'Site'[i,2] + + coefficients(s7.d)$'Year'[j,2]
    }
  }
  
  # make a data frame
  surv.drought <- data.frame(surv.drought)
  
  # Read in top survival model output during recovery (Formula: Surv ~ logSize + (1|Year/Site)
  surv.recovery.reg <- load("data/demography data/surv.recovery.reg.rda")
  
  # Store model coefficients for each Site:Year
  surv.recovery$SiteYear = rownames(coefficients(s4.r)$'Site:Year')
  surv.recovery$surv.int=coefficients(s4.r)$'Site:Year'[,1] 
  surv.recovery$surv.slope=coefficients(s4.r)$'Site:Year'[,2] 
  
  # make a data frame
  surv.recovery <- data.frame(surv.recovery)

    #*******************************************************************************
  ### 1B. Growth ###
  #*******************************************************************************
  # Note: script "dnorm.R" determines that we do not need to switch to the cdf estimation for the growth function
  
  # Read in top growth model output during drought (Formula: logSizeNext ~ logSize + (logSize|Year) + (logSize|Site))
  growth.drought.reg <- load("data/demography data/growth.drought.reg.rda")
  
  # Store model coefficients
  Site=rownames(coefficients(g5.d)$'Site')
  Year=rownames(coefficients(g5.d)$'Year')
  counter=0
  for (i in 1:length(Site)) {
    for (j in 1:length(Year)) {
      counter=counter+1
      growth.drought$SiteYear[counter]=paste(rownames(coefficients(g5.d)$'Site')[i],":",rownames(coefficients(g5.d)$'Year')[j], sep="")
      growth.drought$growth.int[counter]=coefficients(g5.d)$'Site'[i,1] + coefficients(g5.d)$'Year'[j,1]
      growth.drought$growth.slope[counter]=coefficients(g5.d)$'Site'[i,2] + coefficients(g5.d)$'Year'[j,2]
    }
  }

  # Make a data frame
  growth.drought=data.frame(growth.drought)
  
  growth.drought$growth.sd = rep(sigma(g5.d), dim(growth.drought)[1])
  
  # Read in top growth model output during recovery (Formula: logSizeNext ~ logSize + (1|Year/Site))
  growth.recovery.reg <- load("data/demography data/growth.recovery.reg.rda")
  
  # Store model coefficients
  growth.recovery$SiteYear = rownames(coefficients(g4.r)$'Site:Year')
  growth.recovery$growth.int=coefficients(g4.r)$'Site:Year'[,1] 
  growth.recovery$growth.slope=coefficients(g4.r)$'Site:Year'[,2] 
  
  # Make a data frame
  growth.recovery=data.frame(growth.recovery)
  
  growth.recovery$growth.sd = rep(sigma(g4.r), dim(growth.recovery)[1])
  
  #*******************************************************************************
  ### 1C. Flowering ###
  #*******************************************************************************
  
  # Read in top flowering model output during drought (Formula: Fec0 ~ logSize + (1|Year) + (logSize|Site))
  flowering.drought.reg=load("data/demography data/flowering.drought.reg.rda")

  # Store model coefficients
  Site=rownames(coefficients(fl7.d)$'Site')
  Year=rownames(coefficients(fl7.d)$'Year')
  counter=0
  for (i in 1:length(Site)) {
    for (j in 1:length(Year)) {
      counter=counter+1
      flower.drought$SiteYear[counter]=paste(rownames(coefficients(fl7.d)$'Site')[i],":",rownames(coefficients(fl7.d)$'Year')[j], sep="")
      flower.drought$flowering.int[counter]=coefficients(fl7.d)$'Site'[i,1] + coefficients(fl7.d)$'Year'[j,1]
      flower.drought$flowering.slope[counter]=coefficients(fl7.d)$'Site'[i,2] + coefficients(fl7.d)$'Year'[j,2]
    }
  }
  
  # make a data frame
  flower.drought=data.frame(flower.drought)

  # Read in top flowering model output during recovery (Formula: Fec0 ~ logSize + (logSize|Year) + (logSize|Site))
  flowering.recovery.reg=load("data/demography data/flowering.recovery.reg.rda")
  
  # Store model coefficients
  Site=rownames(coefficients(fl5.r)$'Site')
  Year=rownames(coefficients(fl5.r)$'Year')
  counter=0
  for (i in 1:length(Site)) {
    for (j in 1:length(Year)) {
      counter=counter+1
      flower.recovery$SiteYear[counter]=paste(rownames(coefficients(fl5.r)$'Site')[i],":",rownames(coefficients(fl5.r)$'Year')[j], sep="")
      flower.recovery$flowering.int[counter]=coefficients(fl5.r)$'Site'[i,1] + coefficients(fl5.r)$'Year'[j,1]
      flower.recovery$flowering.slope[counter]=coefficients(fl5.r)$'Site'[i,2] + coefficients(fl5.r)$'Year'[j,2]
    }
  }
  
  # make a data frame
  flower.recovery=data.frame(flower.recovery)
  
  #*******************************************************************************
  ### 1D. Fruit number (untransformed) using negative binomial regression ###
  #*******************************************************************************
  # Read in top model output for fruit.reg during drought (Formula: Fec1 ~ logSize + (logSize|Year/Site)  
  fruit.drought.reg=load("data/demography data/fruit.drought.reg.rda")

  # Store model coefficients (fr3.d from glmmTMB)
  SiteYear = rownames(ranef(fr3.d)$cond$'Site:Year') %>% 
    as.data.frame()
  colnames(SiteYear) = "SiteYear"
  SiteYear <- SiteYear %>% 
    separate_wider_delim(SiteYear, delim=":", names=c("Site", "Year"))
  
  temp = c()
  temp$Year = rownames(ranef(fr3.d)$cond$'Year')
  temp$year.int = fixef(fr3.d)$cond[1] + ranef(fr3.d)$cond$'Year'[1]
  temp$year.slope = fixef(fr3.d)$cond[2] + ranef(fr3.d)$cond$'Year'[2]
  temp = data.frame(temp)
  
  SiteYear <-  left_join(SiteYear, temp)
  
  fruit.drought$SiteYear = rownames(ranef(fr3.d)$cond$'Site:Year')
  fruit.drought$fruit.int = fixef(fr3.d)$cond[1] + ranef(fr3.d)$cond$'Site:Year'[1] + SiteYear$X.Intercept.
  fruit.drought$fruit.slope = fixef(fr3.d)$cond[2] + ranef(fr3.d)$cond$'Site:Year'[2] + SiteYear$logSize
  
  # make data frame
  fruit.drought=data.frame(fruit.drought, row.names=NULL) %>% 
    rename(fruit.int = X.Intercept., fruit.slope = logSize)

  # Read in top model output for fruit.reg during recovery (Formula: Fec1 ~ logSize + (1|Year) + (logSize|Site))  
  fruit.recovery.reg=load("data/demography data/fruit.recovery.reg.rda")
  
  # Store model coefficients (fr7.r from glmmTMB)
  Site=rownames(ranef(fr7.r)$cond$'Site')
  Year=rownames(ranef(fr7.r)$cond$'Year')
  counter=0
  for (i in 1:length(Site)) {
    for (j in 1:length(Year)) {
      counter=counter+1
      fruit.recovery$SiteYear[counter]=paste(Site[i],":",Year[j], sep="")
      fruit.recovery$fruit.int[counter]=fixef(fr7.r)$cond[1] + ranef(fr3.d)$cond$'Site'[i,1] + ranef(fr7.r)$cond$'Year'[j,1]
      fruit.recovery$fruit.slope[counter]=fixef(fr7.r)$cond[2] + ranef(fr3.d)$cond$'Site'[i,2] #no year-specific slopes
    }
  }
  
  # make data frame
  fruit.recovery=data.frame(fruit.recovery) 
  
  #*******************************************************************************
  ### 1E. Size distribution of recruits ###
  #*******************************************************************************
  
  # Calculate mean  and standard deviation of size in year t+1 of plants that were new that year (unrecorded at time t)
  recruit_size.drought <- data.drought %>% 
    filter(is.na(logSize)) %>% 
    group_by(SiteYear) %>% 
    summarise(recruit.size.mean = mean(logSizeNext),
              recruit.size.sd = sd(logSizeNext))
  
  recruit_size.recovery <- data.recovery %>% 
    filter(is.na(logSize)) %>% 
    group_by(SiteYear) %>% 
    summarise(recruit.size.mean = mean(logSizeNext),
              recruit.size.sd = sd(logSizeNext))
  
  #*******************************************************************************
  ### 1F. Create data frame of site-specific parameter estimates and join all estimates ###
  #*******************************************************************************
  
  params.drought <- full_join(surv.drought,growth.drought) %>% full_join(flower.drought) %>% full_join(fruit.drought) %>% full_join(recruit_size.drought) %>% 
    replace_na(list(recruit.size.mean=0, recruit.size.sd=0)) # assign 0s to sites with no recruitment (mean & SD = NA) or only 1 recruit (mean = value but SD=NA)

  params.recovery <- full_join(surv.recovery,growth.recovery) %>% full_join(flower.recovery) %>% full_join(fruit.recovery) %>% full_join(recruit_size.recovery) %>% 
    replace_na(list(recruit.size.mean=0, recruit.size.sd=0)) # assign 0s to sites with no recruitment (mean & SD = NA) or only 1 recruit (mean = value but SD=NA)
  
  #*******************************************************************************
  ### 1G. Number of seeds per fruit ###
  #*******************************************************************************
  
  # Note: these are site-specific, but not site-by-year-specific, values, but we want to copy the site-specific values into each site-year level
  
  # obtain mean seed counts per fruit per site
  seeds.per.siteyear.drought <- data.drought %>% 
    group_by(SiteYear, Site) %>% 
    summarize(SeedCtSY = mean(SeedCt, na.rm=T)) # this yields NaN for some site-year combos

  seeds.per.site.drought <- data.drought %>% 
    group_by(Site) %>% 
    summarize(seed.ct = mean(SeedCt, na.rm=T)) 
  
  seeds.per.site.drought <- left_join(seeds.per.siteyear.drought, seeds.per.site.drought) %>% dplyr::select(SiteYear, seed.ct)
  
  # join seeds per fruit to other parameters
  params.drought <- left_join(params.drought, seeds.per.site.drought)
                      
  # obtain mean seed counts per fruit per site
  seeds.per.siteyear.recovery <- data.recovery %>% 
    group_by(SiteYear, Site) %>% 
    summarize(SeedCtSY = mean(SeedCt, na.rm=T)) # this yields NaN for some site-year combos
  
  seeds.per.site.recovery <- data.recovery %>% 
    group_by(Site) %>% 
    summarize(seed.ct = mean(SeedCt, na.rm=T)) 
  
  seeds.per.site.recovery <- left_join(seeds.per.siteyear.recovery, seeds.per.site.recovery) %>% dplyr::select(SiteYear, seed.ct)
  
  # join seeds per fruit to other parameters
  params.recovery <- left_join(params.recovery, seeds.per.site.recovery)
  
  #*******************************************************************************
  ### 1H. Establishment probability ###
  #*******************************************************************************
  
  # Obtain number of new recruits per site at year = t+1
  recruit.number.drought <- data.drought %>% 
    dplyr::filter(is.na(logSize)) %>% 
    group_by(SiteYear) %>% 
    summarise(recruit.number = n()) 
    
  recruit.number.recovery <- data.recovery %>% 
    dplyr::filter(is.na(logSize)) %>% 
    group_by(SiteYear) %>% 
    summarise(recruit.number = n()) 

    # Obtain total fruit count per site at year = t
  fruits.per.site.drought <- site_fruit_count_data.drought %>% 
    filter(!is.na(Fec1)) %>% 
    group_by(SiteYear) %>% 
    summarise(fruits.per.site = sum(Fec1))
  
  fruits.per.site.recovery <- site_fruit_count_data.recovery %>% 
    filter(!is.na(Fec1)) %>% 
    group_by(SiteYear) %>% 
    summarise(fruits.per.site = sum(Fec1))

    # Join recruit # and fruits per site into one frame
  # Obtain total seed count per site (= # fruits per site at year t * # seeds per fruit per site)
  # Obtain establishment probability (= # of new recruits at year t+1/total seed count per site at year t)
  establishment.drought <- full_join(recruit.number.drought, fruits.per.site.drought) %>% full_join(seeds.per.site.drought) %>% 
    replace_na(list(recruit.number=0, fruits.per.site=0)) %>% # assign 0s to sites with no recruitment or no fertility
    # NOTE: Amy added this in 2023. 
    mutate(total.seeds.per.site = fruits.per.site*seed.ct,
           establishment.prob = ifelse(recruit.number==0, 0, 
                                  ifelse(recruit.number>0 & fruits.per.site==0, recruit.number/seed.ct,recruit.number/total.seeds.per.site))) 
    # NOTE: there are a handful of site-years where recruit.number>0 but fruits.per.site=0, resulting in Inf for establishment.prob. Since establishment is non-zero in these site-years, we need a reasonable number for fruits.per.site. Solution used here is to replace fruits.per.site with value of 1 because fecundity was very low that year.
  
  establishment.recovery <- full_join(recruit.number.recovery, fruits.per.site.recovery) %>% full_join(seeds.per.site.recovery) %>% 
    replace_na(list(recruit.number=0, fruits.per.site=0)) %>% # assign 0s to sites with no recruitment or no fertility
    mutate(total.seeds.per.site = fruits.per.site*seed.ct,
           establishment.prob = ifelse(recruit.number==0, 0, 
                                       ifelse(recruit.number>0 & fruits.per.site==0, recruit.number/seed.ct,recruit.number/total.seeds.per.site))) 
  
  # Join with params frame
    params.drought <- left_join(params.drought, establishment.drought)

    params.recovery <- left_join(params.recovery, establishment.recovery)
    
  # Add separate columns for year and site
  params.drought <- params.drought %>% separate(SiteYear, c("Site","Year"), ":", remove=FALSE)
  
  params.recovery <- params.recovery %>% separate(SiteYear, c("Site","Year"), ":", remove=FALSE)
  
#*******************************************************************************
### 2. Create site-specific IPMs parameterized by site-specific parameters derived from global vital rates models 
#*******************************************************************************

  #*******************************************************************************
  ### 2A. Subset data for site f
  #*******************************************************************************
  
  # Which site-years are missing parameter estimates?
  params.missing.drought <- params.drought[!complete.cases(params.drought),]
  params.missing.drought$SiteYear
  # Buck Meadows:2012 --> site inaccessible in 2013 due to fire; lambda is NA
  # Buck Meadows:2013 --> site inaccessible in 2013 due to fire; lambda is NA
  # Carlon:2012 --> site inaccessible in 2013 due to fire; lambda is NA
  # Carlon:2013 --> site inaccessible in 2013 due to fire; lambda is NA
  ###COAST FORK WILLAMETTE 2012 BECAUSE FRUITS = 0
  # Kitchen Creek:2014 --> all plants dead; manually set lambda to 0 below
  # Mill Creek:2013 --> site inaccessible in 2013 due to flood; lambda is NA
  # Mill Creek:2014 --> site inaccessible in 2013 due to flood; lambda is NA
  # Rainbow Pool:2012 --> site inaccessible in 2013 due to fire; lambda is NA
  # Rainbow Pool:2013 --> site inaccessible in 2013 due to fire; lambda is NA
  # Rock Creek:2012 --> 2013 data folder lost; lambda is NA
  # Rock Creek:2013 --> 2013 data folder lost; lambda is NA
  # Whitewater Canyon:2014 --> all plants died; manually set lambda to 0 below

  ###ENTIRELY MISSING FROM PARAMS.DROUGHT FILE
  # Hauser Creek:2012 --> all plants dead; manually set lambda to 0 below
  # Hauser Creek:2013 --> all plants dead; manually set lambda to 0 below
  # Hauser Creek:2014 --> all plants dead; manually set lambda to 0 below
  
  ### NEWLY APPEARING AS COMPLETE CASES
  # Kitchen Creek: 2013 --> only remaining plant died; manually set lambda to 0 below
  # West Fork Mojave River:2013 --> existing plots 1-5 all dead, so some parameters inestimable. But new plot 6 established in 2014, so the entire site was not dead, only the main area where we were observing 2010-2013. Keep as NA because the entire site was not dead (in contrast to Hauser, Kitchen, Whitewater, where we set lambda to 0 when all plants died).
  
  params.missing.recovery <- params.recovery[!complete.cases(params.recovery),]
  params.missing.recovery$SiteYear
  # Carlon:2016 --> growth slopes, intercepts, and sd are inestimable because of sparse data, but site was visited and censused and all other parameters are estimable; use mean across other years at this site for this year's growth parameters
  params.recovery$growth.int[params.recovery$SiteYear=="Carlon:2016"] = mean(params.recovery$growth.int[params.recovery$Site=="Carlon"], na.rm=T)
  params.recovery$growth.slope[params.recovery$SiteYear=="Carlon:2016"] = mean(params.recovery$growth.slope[params.recovery$Site=="Carlon"], na.rm=T)
  params.recovery$growth.sd[params.recovery$SiteYear=="Carlon:2016"] = mean(params.recovery$growth.sd[params.recovery$Site=="Carlon"], na.rm=T)
  # Kitchen Creek:2016 --> all plants dead; manually set lambda to 0 below
  # Whitewater Canyon:2017 --> site visited but no plants; manually set lambda to 0 below
  # Whitewater Canyon:2018 --> site visited but no plants; manually set lambda to 0 below
  # Canton Creek:2016 --> site inaccessible in 2017 due to high water; lambda is NA
  # Canton Creek 2017 --> site inaccessible in 2017 due to high water; lambda is NA
  # Hauser Creek:2015 --> all plants dead; manually set lambda to 0 below
  # Hauser Creek:2016 --> all plants dead; manually set lambda to 0 below
  # Hauser Creek:2017 --> all plants dead; manually set lambda to 0 below
  # Kitchen Creek:2015 --> all plants dead; manually set lambda to 0 below
  # Kitchen Creek:2017 --> all plants dead; manually set lambda to 0 below
  # Kitchen Creek:2018 --> new recruits appeared but all other transitions inestimable because of no survivors; lambda is NA
  # North Fork Middle Fork Tule:2016 --> fire closure in 2017; lambda is NA
  # North Fork Middle Fork Tule:2017 --> fire closure in 2017; lambda is NA
  # Redwood Creek:2015 --> fire closure in 2016; lambda is NA
  # Redwood Creek:2016 --> fire closure in 2016; lambda is NA
  # South Fork Middle Fork Tule:2016 --> fire closure in 2017; lambda is NA
  # South Fork Middle Fork Tule:2017 --> fire closure in 2017; lambda is NA
  # Whitewater Canyon:2015 --> site visited but no plants; manually set lambda to 0 below
  # Whitewater Canyon:2016 --> site visited but no plants; manually set lambda to 0 below
  
  
# Store parameters in .csv file for later use
write.csv(params.drought,"data/demography data/vital_rate_coefficients_droughtperiod.csv", row.names=FALSE)
write.csv(params.recovery,"data/demography data/vital_rate_coefficients_recoveryperiod.csv", row.names=FALSE)

# Stack parameter estimates now that they are uniquely fit to drought vs recovery periods
  params <- bind_rows(params.drought, params.recovery)

# Remove rows of parameters with NA
  params <- params[complete.cases(params), ]

  siteYear <- params$SiteYear

  # Create empty vectors for lambda and site to be filled
  lambda=c()
  SiteYear=character()

  for (f in 1:length(siteYear)) {
    data1 = subset(data, SiteYear==siteYear[f])
    params1 = subset(params, SiteYear==siteYear[f])
    params1 = subset(params1, select=-SiteYear)
    
    #*******************************************************************************
    ### 2B. Create survival, growth, and fecundity functions and build IPM by running integral_projection_model.R script
    #*******************************************************************************
    
    source("scripts/demography scripts/integral_projection_model.R")
    
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

# Create data frame of Site, Latitude, Longitude, Region, and Elevation 
site.info <- subset(data, select=c(Site,Year,SiteYear,Latitude,Longitude,Elevation,Region,RegionRank)) %>% unique() %>% arrange(-Latitude)
    
# Merge site info with lambda estimates
site.info <- full_join(site.info, siteYear.lambda)

# How many site-year combos are possible? 
20*9 #180
focal.sites <- c(rep("Coast Fork of Williamette",9),
                 rep("Canton Creek",9),
                 rep("Rock Creek",9),
                 rep("O'Neil Creek",9),
                 rep("Deep Creek",9),
                 rep("Little Jameson Creek",9),
                 rep("Oregon Creek",9),
                 rep("Rainbow Pool",9),
                 rep("Carlon",9),
                 rep("Buck Meadows",9),
                 rep("Wawona",9),
                 rep("Redwood Creek",9),
                 rep("North Fork Middle Fork Tule",9),
                 rep("South Fork Middle Fork Tule",9),
                 rep("West Fork Mojave River",9),
                 rep("Mill Creek",9),
                 rep("Whitewater Canyon",9),
                 rep("Sweetwater River",9),
                 rep("Kitchen Creek",9),
                 rep("Hauser Creek",9)) 
years <- rep(c("2010","2011","2012","2013","2014","2015","2016","2017","2018"), 20) 
site.years.max <- as.data.frame(cbind(focal.sites, years)) %>% mutate(SiteYear = paste(focal.sites,":", years, sep="")) %>% dplyr::select(SiteYear)

# How many site-years have some observed data but are missing lambda estimates?
length(unique(data$SiteYear)) #153
site.years.obs <- as.data.frame(unique(data$SiteYear))
colnames(site.years.obs) = "SiteYear"
   
# Which site-years are missing, and why?
missing.site.years <- anti_join(site.years.max, site.years.obs)

# Canton Creek:2016 --> site inaccessible in 2017 due to high water; this is a real NA
# Rock Creek:2012 --> 2013 data folder lost; this is a real NA 
# Rock Creek:2013 --> 2013 data folder lost; this is a real NA 
# Rainbow Pool:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Rainbow Pool:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Carlon:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Carlon:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Buck Meadows:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Buck Meadows:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Redwood Creek:2015 --> site inaccessible in 2016 due to fire; this is a real NA
# North Fork Middle Fork Tule: 2016 --> site inaccessible in 2017 due to fire; this is a real NA
# North Fork Middle Fork Tule: 2017 --> site inaccessible in 2017 due to fire; this is a real NA
# South Fork Middle Fork Tule:2010 --> site inaccessible in 2011 due to fire; this is a real NA
# South Fork Middle Fork Tule:2011 --> site inaccessible in 2011 due to fire; this is a real NA
# South Fork Middle Fork Tule:2016 --> site inaccessible in 2017 due to fire; this is a real NA
# South Fork Middle Fork Tule:2017 --> site inaccessible in 2017 due to fire; this is a real NA
# Mill Creek: 2013 --> site inaccessible in 2013 due to flood; this is a real NA
# Whitewater Canyon:2015 --> site visited but no plants so manually set lambda=0
Whitewater=filter(site.info, Site=="Whitewater Canyon" & Year==2014)
Whitewater$Year=2015 %>% factor()
Whitewater$SiteYear="Whitewater Canyon:2015" %>% factor()
Whitewater$lambda=0
site.info=bind_rows(site.info, Whitewater) 

# Kitchen Creek:2014 --> site visited but no plants so manually set lambda=0
Kitchen=filter(site.info, Site=="Kitchen Creek" & Year==2013)
Kitchen$Year=2014 %>% factor()
Kitchen$SiteYear="Kitchen Creek:2014" %>% factor()
Kitchen$lambda=0
site.info=bind_rows(site.info, Kitchen) 

# Kitchen Creek:2015 --> site visited but no plants so manually set lambda=0
Kitchen$Year=2015 %>% factor()
Kitchen$SiteYear="Kitchen Creek:2015" %>% factor()
Kitchen$lambda=0
site.info=bind_rows(site.info, Kitchen) 

# Kitchen Creek:2017 --> site visited but no plants so manually set lambda=0
Kitchen$Year=2017 %>% factor()
Kitchen$SiteYear="Kitchen Creek:2017" %>% factor()
Kitchen$lambda=0
site.info=bind_rows(site.info, Kitchen) 

# Hauser Creek:2012 --> site visited but no plants so manually set lambda=0
Hauser=filter(site.info, Site=="Hauser Creek" & Year==2011)
Hauser$Year=2012 %>% factor()
Hauser$SiteYear="Hauser Creek:2012" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 

# Hauser Creek:2013 --> site visited but no plants so manually set lambda=0
Hauser$Year=2013 %>% factor()
Hauser$SiteYear="Hauser Creek:2013" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 

# Hauser Creek:2014 --> site visited but no plants so manually set lambda=0
Hauser$Year=2014 %>% factor()
Hauser$SiteYear="Hauser Creek:2014" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 

# Hauser Creek:2015 --> site visited but no plants so manually set lambda=0
Hauser$Year=2015 %>% factor()
Hauser$SiteYear="Hauser Creek:2015" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 

# Hauser Creek:2016 --> site visited but no plants so manually set lambda=0
Hauser$Year=2016 %>% factor()
Hauser$SiteYear="Hauser Creek:2016" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 

# Hauser Creek:2017 --> site visited but no plants so manually set lambda=0
Hauser$Year=2017 %>% factor()
Hauser$SiteYear="Hauser Creek:2017" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 


# Which site-years have lambda=NA, and why?
lambda.calc.failed <- site.info %>% dplyr::filter(is.na(lambda)) %>% dplyr::select(SiteYear)

# Canton Creek:2017 --> site skipped due to high water; real NA
# Redwood Creek:2016 --> site skipped due to fire closure; real NA
# West Fork Mojave River:2013 --> existing plots 1-5 all dead, so some parameters inestimable. But new plot 6 established in 2014, so the entire site was not dead, only the main area where we were observing 2010-2013. Keep as NA because the entire site was not dead (in contrast to Hauser, Kitchen, Whitewater, where we set lambda to 0 when all plants died).
# Mill Creek:2010 --> all 2010 plots washed out and new plots established in 2011; lambda is NA
# Mill Creek:2014 --> site inaccessible in 2013 due to flood; lambda is NA
# Whitewater Canyon:2014 --> all plants died, so set lambda to 0 
site.info$lambda[site.info$SiteYear=="Whitewater Canyon:2014"] = 0
# Whitewater Canyon:2016 --> all plants died, so set lambda to 0 
site.info$lambda[site.info$SiteYear=="Whitewater Canyon:2016"] = 0
# Whitewater Canyon:2017 --> all plants died, so set lambda to 0 
site.info$lambda[site.info$SiteYear=="Whitewater Canyon:2017"] = 0
# Whitewater Canyon:2018 --> all plants died, so set lambda to 0 
site.info$lambda[site.info$SiteYear=="Whitewater Canyon:2018"] = 0
# Kitchen Creek: 2013 --> all plants died, so set lambda to 0
site.info$lambda[site.info$SiteYear=="Kitchen Creek:2013"] = 0
# Kitchen Creek: 2016 --> all plants died, so set lambda to 0
site.info$lambda[site.info$SiteYear=="Kitchen Creek:2016"] = 0
# Kitchen Creek: 2018 --> new recruits, but no survivors to estimate surv and growth functions; keep as NA

# View final data frame
str(site.info)

# Save to lambdas to .csv file 
write.csv(site.info,"data/demography data/siteYear.lambda_2010-2019.csv",row.names=FALSE)

# Merge params and lambdas for supplementary table
full.table <- left_join(site.info, params) %>% mutate_at(9:20, round, 2) %>% arrange(Latitude)
write.csv(full.table,"data/demography data/siteYear.paramslambda_2010-2019.csv",row.names=FALSE)

#*******************************************************************************
### 4. Compare to preliminary estimates 
#*******************************************************************************
# 
old.lambda <- read.csv("data/demography data/siteYear.lambda_2010-2016_old.csv") %>% 
  dplyr::select(SiteYear, lambda.old=lambda)
# Note: this file is copied in from demography_analyses repo

site.info <- site.info %>% 
  mutate(SiteYear = gsub(" ", "", SiteYear))

all <- left_join(site.info, old.lambda) %>% 
  mutate(lambda.diff = lambda - lambda.old,
         check = ifelse(abs(lambda.diff)>1, "YES", "no"),
         check.prop = ifelse(abs(lambda.diff/lambda)>0.50, "YES", "no"))

summary(all$lambda.diff)
table(all$check)
table(all$check.prop)

check.list = all$SiteYear[all$check=="YES"]

ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + geom_abline(x=y) + xlim(0,70)
# Note: Buck Meadows:2015 and Sweetwater:2016 are absurd - what happened??
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(0,0.5) + ylim(0,0.5) + geom_abline(x=y)
# looks great - symmetrical and very close to 1:1
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(0.5,1) + ylim(0.5,1) + geom_abline(x=y)
# looks good - symmetrical and fairly close to 1:1
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(1,2) + ylim(1,2) + geom_abline(x=y)
# looks good - symmetrical and mostly close to 1:1
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(2,3) + ylim(2,3) + geom_abline(x=y)
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(3,10) + ylim(3,10) + geom_abline(x=y)
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(10,50) + ylim(10,50) + geom_abline(x=y)


