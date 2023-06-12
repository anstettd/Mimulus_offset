#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Create data frame of vital rate parameters and build integral projection models to obtain estimates of annual lambdas for each population
#### AUTHOR: Seema Sheth and Amy Angert
#### DATE LAST MODIFIED: 20230612


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

# Set up data frames of model parameters
params=c()
flower=c()
fruit=c()


  #*******************************************************************************
  ### 1A. Survival ###
  #*******************************************************************************

  # Read in top survival model output (Formula: Surv ~ logSize + (logSize|Year) + (logSize|Site))
  surv.reg <- load("data/demography data/surv.reg.rda")

  # Store model coefficients for each Site:Year
  Site=rownames(coefficients(s5)$'Site')
  Year=rownames(coefficients(s5)$'Year')
  counter=0
    for (i in 1:length(Site)) {
      for (j in 1:length(Year)) {
        counter=counter+1
        params$surv.int[counter]=coefficients(s5)$'Site'[i,1] + coefficients(s5)$'Year'[j,1]
        params$surv.slope[counter]=coefficients(s5)$'Site'[i,2] + coefficients(s5)$'Year'[j,2]
      params$SiteYear[counter]=paste(rownames(coefficients(s5)$'Site')[i],":",rownames(coefficients(s5)$'Year')[j], sep="")
        }
      }
  
  
  #*******************************************************************************
  ### 1B. Growth ###
  #*******************************************************************************
  
  # Read in top growth model output (Formula: logSizeNext ~ logSize + (logSize|Year) + (logSize|Site))
  growth.reg <- load("data/demography data/growth.reg.rda")
  
  # Store model coefficients
  Site=rownames(coefficients(g5)$'Site')
  Year=rownames(coefficients(g5)$'Year')
  counter=0
  for (i in 1:length(Site)) {
    for (j in 1:length(Year)) {
      counter=counter+1
      params$growth.int[counter]=coefficients(g5)$'Site'[i,1] + coefficients(g5)$'Year'[j,1]
      params$growth.slope[counter]=coefficients(g5)$'Site'[i,2] + coefficients(g5)$'Year'[j,2]
      params$growth.sd[counter]=sigma(g5)
      params$SiteYear[counter]=paste(rownames(coefficients(g5)$'Site')[i],":",rownames(coefficients(g5)$'Year')[j], sep="")
    }
  }
  
  # Note: script "dnorm.R" determines that we do not need to switch to the cdf estimation for the growth function
  
  #*******************************************************************************
  ### 1C. Flowering ###
  #*******************************************************************************
  
  # Read in top flowering model output (Formula: Fec0 ~ logSize + (1|Year/Site))
  flowering.reg=load("data/demography data/flowering.reg.rda")

  # Store model coefficients
  flower$SiteYear = rownames(coefficients(fl4)$'Site:Year')
  flower$flowering.int=coefficients(fl4)$'Site:Year'[,1] 
  flower$flowering.slope=coefficients(fl4)$'Site:Year'[,2] 
  
  # make a data frame
  flower=data.frame(flower)

    #*******************************************************************************
  ### 1D. Fruit number (untransformed) using negative binomial regression ###
  #*******************************************************************************
  # Read in top model output for fruit.reg (Formula: Fec1 ~ logSize + (1|Year/Site))   
  fruit.reg=load("data/demography data/fruit.reg.rda")

  # Store model coefficients (fr4 from glmmTMB)
  fruit$SiteYear=rownames(ranef(fr4)$cond$'Site:Year')
  fruit$fruit.int=fixef(fr4)$cond[1]+ranef(fr4)$cond$'Site:Year'[,1] 
  fruit$fruit.slope=rep(fixef(fr4)$cond[2],times=length(rownames(ranef(fr4)$cond$'Site:Year'))) 
  
  # make data frame
  fruit=data.frame(fruit) 
  
  #*******************************************************************************
  ### 1E. Size distribution of recruits ###
  #*******************************************************************************
  
  # Calculate mean  and standard deviation of size in year t+1 of plants that were new that year (unrecorded at time t)
  recruit_size <- data %>% 
    filter(is.na(logSize)) %>% 
    group_by(SiteYear) %>% 
    summarise(recruit.size.mean = mean(logSizeNext),
              recruit.size.sd = sd(logSizeNext))
  
   #*******************************************************************************
  ### 1F. Create data frame of site-specific parameter estimates and join all estimates ###
  #*******************************************************************************
  
  params <- data.frame(params)
  params <- full_join(params,flower) %>% full_join(fruit) %>% full_join(recruit_size) %>% 
    replace_na(list(recruit.size.mean=0, recruit.size.sd=0)) # assign 0s to sites with no recruitment (mean & SD = NA) or only 1 recruit (mean = value but SD=NA)

  
  #*******************************************************************************
  ### 1G. Number of seeds per fruit ###
  #*******************************************************************************
  
  # Note: these are site-specific, but not site-by-year-specific, values, but we want to copy the site-specific values into each site-year level
  
  # obtain mean seed counts per fruit per site
  seeds.per.siteyear <- data %>% 
    group_by(SiteYear, Site) %>% 
    summarize(SeedCtSY = mean(SeedCt, na.rm=T)) # this yields NaN for some site-year combos

  seeds.per.site <- data %>% 
    group_by(Site) %>% 
    summarize(seed.ct = mean(SeedCt, na.rm=T)) 
  
  seeds.per.site <- left_join(seeds.per.siteyear, seeds.per.site) %>% dplyr::select(SiteYear, seed.ct)
  
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
  # Obtain total seed count per site (= # fruits per site at year t * # seeds per fruit per site)
  # Obtain establishment probability (= # of new recruits at year t+1/total seed count per site at year t)
  establishment <- full_join(recruit.number, fruits.per.site) %>% full_join(seeds.per.site) %>% 
    replace_na(list(recruit.number=0, fruits.per.site=0)) %>% # assign 0s to sites with no recruitment or no fertility
    # NOTE: Amy added this in 2023. 
    mutate(total.seeds.per.site = fruits.per.site*seed.ct,
           establishment.prob = ifelse(recruit.number==0, 0, 
                                  ifelse(recruit.number>0 & fruits.per.site==0, recruit.number/seed.ct,recruit.number/total.seeds.per.site))) 
    # NOTE: there are a handful of site-years where recruit.number>0 but fruits.per.site=0, resulting in Inf for establishment.prob. Since establishment is non-zero in these site-years, we need a reasonable number for fruits.per.site. Solution used here is to replace fruits.per.site with value of 1 because fecundity was very low that year.
  
  # Join with params frame
    params <- left_join(params, establishment)

  # Add separate columns for year and site
  params <- params %>% separate(SiteYear, c("Site","Year"), ":", remove=FALSE)
  
  
#*******************************************************************************
### 2. Create site-specific IPMs parameterized by site-specific parameters derived from global vital rates models 
#*******************************************************************************

  #*******************************************************************************
  ### 2A. Subset data for site f
  #*******************************************************************************
  
  # Which site-years are missing parameter estimates?
  params.missing <- params[!complete.cases(params),]
  params.missing$SiteYear
  # Buck Meadows:2012 --> site inaccessible in 2013 due to fire; lambda is NA
  # Buck Meadows:2013 --> site inaccessible in 2013 due to fire; lambda is NA
  # Carlon:2012 --> site inaccessible in 2013 due to fire; lambda is NA
  # Carlon:2013 --> site inaccessible in 2013 due to fire; lambda is NA
  # Coast Fork of Willamette:2012 --> no fruits in 2012 so fruit slopes & intercepts are inestimable because of sparse data, but site was visited and censused and all other parameters are estimable; use mean across other years at this site for this year's fruit parameters or set lambda to NA?
  # Hauser Creek: 2012 --> sites visited but no plants; manually set lambda to 0 below
  # Hauser Creek: 2013 --> sites visited but no plants; manually set lambda to 0 below
  # Hauser Creek: 2014 --> sites visited but no plants; manually set lambda to 0 below
  # Hauser Creek: 2015 --> sites visited but no plants; manually set lambda to 0 below
  # Kitchen Creek: 2014 --> sites visited but no plants; manually set lambda to 0 below
  # Kitchen Creek: 2015 --> sites visited but no plants; manually set lambda to 0 below
  # Rainbow Pool: 2012 --> site inaccessible in 2013 due to fire; lambda is NA
  # Rainbow Pool: 2013 --> site inaccessible in 2013 due to fire; lambda is NA
  # Redwood Creek: 2015 --> ??
  # Rock Creek:2012 --> 2013 data folder lost; lambda is NA
  # Rock Creek:2013 --> 2013 data folder lost; lambda is NA
  # South Fork Middle Fork Tule:2010 --> site inaccessible in 2011; lambda is NA
  # South Fork Middle Fork Tule:2011 --> site inaccessible in 2011; lambda is NA
  # Whitewater Canyon:2014 --> all plants died; manually lambda to 0 below
  # Whitewater Canyon:2015 --> all plants died; manually lambda to 0 below

  
# Store parameters in .csv file for later use
write.csv(params,"data/demography data/vital_rate_coefficients.csv", row.names=FALSE)

# Remove rows of parameters with NA
  params <- params[complete.cases(params), ]
  
  
  # Remove site x year combinations without parameter estimates
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
19*6 #114
focal.sites <- c(rep("Coast Fork of Williamette",6),
                 rep("Canton Creek",6),
                 rep("Rock Creek",6),
                 rep("O'Neil Creek",6),
                 rep("Deep Creek",6),
                 rep("Little Jameson Creek",6),
                 rep("Oregon Creek",6),
                 rep("Rainbow Pool",6),
                 rep("Carlon",6),
                 rep("Buck Meadows",6),
                 rep("Wawona",6),
                 rep("Redwood Creek",6),
                 rep("North Fork Middle Fork Tule",6),
                 rep("South Fork Middle Fork Tule",6),
                 rep("West Fork Mojave River",6),
                 rep("Whitewater Canyon",6),
                 rep("Sweetwater River",6),
                 rep("Kitchen Creek",6),
                 rep("Hauser Creek",6)) 
years <- rep(c("2010","2011","2012","2013","2014","2015"), 19) 
site.years.max <- as.data.frame(cbind(focal.sites, years)) %>% mutate(SiteYear = paste(focal.sites,":", years, sep="")) %>% dplyr::select(SiteYear)

# How many site-years have some observed data but are missing lambda estimates?
length(unique(data$SiteYear)) #96
site.years.obs <- as.data.frame(unique(data$SiteYear))
colnames(site.years.obs) = "SiteYear"
   
# Which site-years are missing, and why?
missing.site.years <- anti_join(site.years.max, site.years.obs)

# Rock Creek:2012 --> 2013 data folder lost; this is a real NA 
# Rock Creek:2013 --> 2013 data folder lost; this is a real NA 
# Rainbow Pool:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Rainbow Pool:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Carlon:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Carlon:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Buck Meadows:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Buck Meadows:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Redwood Creek:2015 --> site inaccessible in 2016 due to fire; this is a real NA
# South Fork Middle Fork Tule:2010 --> site inaccessible in 2011; this is a real NA
# South Fork Middle Fork Tule:2011 --> site inaccessible in 2011; this is a real NA
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
Kitchen=filter(site.info, Site=="Kitchen Creek" & Year==2014)
Kitchen$Year=2015 %>% factor()
Kitchen$SiteYear="Kitchen Creek:2015" %>% factor()
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

# Which site-years have lambda=NA, and why?
lambda.calc.failed <- site.info %>% dplyr::filter(is.na(lambda)) %>% dplyr::select(SiteYear)

# Coast Fork Willamette:2012 --> no fruits; keep as NA
# Whitewater Canyon:2014 --> all plants died, so set lambda to 0 
site.info$lambda[site.info$SiteYear=="Whitewater Canyon:2014"] = 0

# View final data frame
str(site.info)

# Save to lambdas to .csv file 
write.csv(site.info,"data/demography data/siteYear.lambda_2010-2016.csv",row.names=FALSE)

# Merge params and lambdas for supplementary table
full.table <- left_join(site.info, params) %>% mutate_at(9:20, round, 2) %>% arrange(Latitude)
write.csv(full.table,"data/demography data/siteYear.paramslambda_2010-2016.csv",row.names=FALSE)

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
# Note: Deer Creek:2012 has come down from >30 to <20 but is unreasonable regardless
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(0,0.5) + ylim(0,0.5) + geom_abline(x=y)
# Note: many more below 1:1 than above
# Note: Whitewater:2010 went from ~0.0 to ~0.3
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(0.5,1) + ylim(0.5,1) + geom_abline(x=y)
# Note: symmetrical around 1:1 
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(1,2) + ylim(1,2) + geom_abline(x=y)
# Note: many more above 1:1 line than below. Deep Creek:2014 has come up a lot due to exclusion of plot 4 where all plants were erroneously coded as D instead of E
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(2,3) + ylim(2,3) + geom_abline(x=y)
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(3,10) + ylim(3,10) + geom_abline(x=y)
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(10,50) + ylim(10,50) + geom_abline(x=y)
# Note: Deer Creek:2012 has come down from >30 to ~20 but is unreasonable regardless


