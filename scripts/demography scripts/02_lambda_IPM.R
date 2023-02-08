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
#SiteYear <- unique(data$SiteYear)

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
  params <- full_join(params,growth) %>% full_join(fruit) %>% full_join(recruit_size) %>% 
    replace_na(list(recruit.size.mean=0, recruit.size.sd=0)) # assign 0s to sites with no recruitment (mean & SD = NA) or only 1 recruit (mean = value but SD=NA)
  # TO DO: decide if there any other situations where param values of NA should be re-coded as 0.
    
  
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
  # Obtain total seed count per site (= # fruits per site * # seeds per fruit per site) at year = t
  # Obtain establishment probability as # of new recruits at year = t+1/# of seeds at year = t
  establishment <- full_join(recruit.number, fruits.per.site) %>% full_join(seeds.per.site) %>% 
    replace_na(list(recruit.number=0, fruits.per.site=0)) %>% # assign 0s to sites with no recruitment or no fertility
    # NOTE: Amy added this in 2023. Ok?
    mutate(total.seeds.per.site = fruits.per.site*seed.ct,
           establishment.prob = ifelse(recruit.number==0, 0, 
                                  ifelse(recruit.number>0 & fruits.per.site==0, recruit.number/seed.ct,                                          recruit.number/total.seeds.per.site))) 
    # TO DO: there are a handful of site-years where recruit.number>0 but fruits.per.site=0, resulting in Inf for establishment.prob. Since establishment is non-zero in these site-years, we need a reasonable number for fruits.per.site. Ideas: (a) replace fruits.per.site with mean across other years or (b) replace fruits.per.site with value of 1 because fecundity was very low that year? Currently using (b)
  
  # Join with params frame
    params <- left_join(params, establishment)

  # Add separate columns for year and site
  params <- params %>% separate(SiteYear, c("Site","Year"), ":", remove=FALSE)
  
  # Store parameters in .csv file for later use
  write.csv(params,"data/demography data/vital_rate_coefficients.csv", row.names=FALSE)
  
#*******************************************************************************
### 2. Create site-specific IPMs parameterized by site-specific parameters derived from global vital rates models 
#*******************************************************************************

  #*******************************************************************************
  ### 2A. Subset data for site f
  #*******************************************************************************
  
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
21*6 #126
focal.sites <- c(rep("Coast Fork of Williamette",6),
                 rep("Canton Creek",6),
                 rep("Rock Creek",6),
                 rep("Deer Creek",6),
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
                 rep("Mill Creek",6),
                 rep("Whitewater Canyon",6),
                 rep("Sweetwater River",6),
                 rep("Kitchen Creek",6),
                 rep("Hauser Creek",6)) 
years <- rep(c("2010","2011","2012","2013","2014","2015"), 21) 

site.years.max <- as.data.frame(cbind(focal.sites, years)) %>% mutate(SiteYear = paste(focal.sites,":", years, sep="")) %>% dplyr::select(SiteYear)

# How many site-years have some observed data but are missing lambda estimates?
length(unique(data$SiteYear)) #106
site.years.obs <- as.data.frame(unique(data$SiteYear))
colnames(site.years.obs) = "SiteYear"
   
# Which are missing, and why?
missing.site.years <- anti_join(site.years.max, site.years.obs)

# Deer Creek:2010 --> site established in 2011; this is a real NA
# Mill Creek:2013 --> site inaccessible in 2014 due to flood; this is a real NA
### TO DO: Figure out why Mill Creek:2014 doesn't show up on this scan for inestimable site-years
# Note: Mill Creek:2014 has lambda=NA in site.info file
# South Fork Middle Fork Tule:2010 --> site inaccessible in 2011; this is a real NA
# South Fork Middle Fork Tule:2011 --> site inaccessible in 2011; this is a real NA
# Carlon:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Carlon:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Rainbow Pool:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Rainbow Pool:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Buck Meadows:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Buck Meadows:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Rock Creek:2012 --> 2013 data folder lost; this is a real NA 
# Rock Creek:2013 --> 2013 data folder lost; this is a real NA 

# Redwood Creek:2015

# Sites visited but no plants --> manually set lambda=0
# Hauser Creek:2012
# Hauser Creek:2013
# Hauser Creek:2014
# Hauser Creek:2015
# Kitchen Creek:2014
# Kitchen Creek:2015
# Whitewater Canyon:2015

Hauser=filter(site.info, Site=="Hauser Creek" & Year==2011)
Hauser$Year=2012 %>% factor()
Hauser$SiteYear="Hauser Creek:2012" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 
Hauser$Year=2013 %>% factor()
Hauser$SiteYear="Hauser Creek:2013" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 
Hauser$Year=2014 %>% factor()
Hauser$SiteYear="Hauser Creek:2014" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 
Hauser$Year=2015 %>% factor()
Hauser$SiteYear="Hauser Creek:2015" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) %>% mutate(Year=factor(Year), SiteYear=factor(SiteYear))
Kitchen=filter(site.info, Site=="Kitchen Creek" & Year==2013)
Kitchen$Year=2014 %>% factor()
Kitchen$SiteYear="Kitchen Creek:2014" %>% factor()
Kitchen$lambda=0
site.info=bind_rows(site.info, Kitchen) 
Kitchen$Year=2015 %>% factor()
Kitchen$SiteYear="Kitchen Creek:2015" %>% factor()
Kitchen$lambda=0
site.info=bind_rows(site.info, Kitchen) %>% mutate(Year=factor(Year), SiteYear=factor(SiteYear))
Whitewater=filter(site.info, Site=="Whitewater Canyon" & Year==2014)
Whitewater$Year=2015 %>% factor()
Whitewater$SiteYear="Whitewater Canyon:2015" %>% factor()
Whitewater$lambda=0
site.info=bind_rows(site.info, Whitewater) %>% mutate(Year=factor(Year), SiteYear=factor(SiteYear))

# Which sites have lambda=NA, and why?
lambda.calc.failed <- site.info %>% dplyr::filter(is.na(lambda)) %>% dplyr::select(SiteYear)

# Coast Fork of Willamette:2012
# Canton Creek:2011
# West Fork Mojave River:2013
# Mill Creek:2010 --> all 2010 plots washed out and new plots established in 2011; keep as NA
# Mill Creek:2014 --> site not visited in 2013, so this makes sense; keep as NA
# Whitewater Canyon:2014 --> all plants died, so set lambda to 0
site.info$lambda[site.info$SiteYear=="Whitewater Canyon:2014"] = 0
# Kitchen Creek:2013 --> only remaining plant died, so set lambda to 0
site.info$lambda[site.info$SiteYear=="Kitchen Creek:2013"] = 0
# Hauser Creek:2011 --> all plants dead on plots 1-2, many plants on plot 3 but individuals indistinguishable; keep as NA


# view final data frame
str(site.info)

# save to .csv file 
write.csv(site.info,"R_output/siteYear.lambda_2010-2016.csv",row.names=FALSE)


