#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Tidy raw demographic vital rate data in preparation for IPM analyses 
#### AUTHOR: Seema Sheth and Amy Angert
#### DATE LAST MODIFIED: 20230202

#*******************************************************************************
#### 1. Load required packages
#*******************************************************************************

# Easy code for installing packages in R (if not installed) and calling their libraries
# From: https://gist.github.com/DrK-Lo/a945a29d6606b899022d0f03109b9483

# Make vector of packages needed
packages_needed <- c("tidyverse")

# Install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# Load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

#*******************************************************************************
#### 2. Bring in and combine M. cardinalis vital rate data from 2010-2016
#*******************************************************************************

# Read in seed count per fruit data, select and rename relevant columns, and round to nearest integer
seed.ct <- read.csv("data/demography data/fec2.seed.per.fruit.2010.2011.2012.csv")
seed.ct <- subset(seed.ct,select=c(site,newgrandmean))
colnames(seed.ct) = c("Site","SeedCt")
seed.ct$SeedCt = round(seed.ct$SeedCt,digits=0) # round seed count to nearest integer
seed.ct$Site = factor(seed.ct$Site) # make site column a factor to streamline joining
# TO DO: update with more recent years' collections

# Read in vital rate data for 2010-11, 2011-12, 2012-13, 2013-14 transitions 
# Note: PY=previous year (time t), CY=current year (time t+1); ignore PPY
data_2010.2014=read.csv("data/demography data/Mcard_demog_data_2010-2014.csv") %>% dplyr::select(-Reasoning, -Reasoning.1) #remove unwanted columns
# Note: this file was created for the analyses published in Sheth and Angert 2018 PNAS 
# It results from Amy Angert's work in July 2016 (original file: "Mcard_demog_data_2010-2013_ALA.xlsx") to scan datasheet notes to identify individuals to exclude, based on these columns:
# Column 'NotAnIndividual': 
#   0 = ok, definitely include (includes "scattered" as long as nothing else noted)
#   1 = not ok, definitely exclude from survival, growth, and fecundity but ok for seed input denominator for recruitment (history of lumping/splitting/relumping; redundant IDs)
#   2 = maybe (notes about difficult to distinguish, or merged once)
# Column 'NotARecruit':
#   0 = ok, definitely include
#   1 = wrong, definitely exclude (reasons include new plot, site not visited in prior year, ID within prior years' ranges, coordinates well outside of prior year's search
#   2 = plant noted as possibly missed in prior year (looks old, missed?, J-20xx?, could be [old ID], etc)
#   3 = plant not noted as unusually old-looking, but is within size range of or larger than most plants that are in category 2
#   NA = size measures at year t

# Read in vital rate data for 2014-15 transition 
# Note: PY=previous year (time t), CY=current year (time t+1); ignore PPY
data_2014.2015=read.csv("data/demography data/SS_Horizontal_2015_Notes.csv") %>% 
  rename(Site = SiteID, #rename to match 2010-2014 data
         Class = Class_PY,
         ClassNext = Class_CY,
         Fec1 = TotFr_PY,
         Size = TotStLn_PY,
         SizeNext = TotStLn_CY,
         Surv = SurvPYCY) %>% 
  mutate(Year = 2014) #add year at time t (=2014) column
# Note: this file was queried from the database on 2023-02-01
# We will tidy it below based on the decision rules as above for 2010-14 data, but scripted here instead of done manually in excel. 

# Read in vital rate data for 2015-16 transition
# Note: PY=previous year (time t), CY=current year (time t+1); ignore PPY
data_2015.2016 = read.csv("data/demography data/SS_Horizontal_2016_Notes.csv") %>% 
  rename(Site = SiteID, #rename to match 2010-2014 data
         Class = Class_PY,
         ClassNext = Class_CY,
         Fec1 = TotFr_PY,
         Size = TotStLn_PY,
         SizeNext = TotStLn_CY,
         Surv = SurvPYCY) %>% 
  mutate(Year = 2015) #add year at time t (=2015) column)
# Note: this file was queried from the database on 2023-02-01
# We script the creation of NotARecruit and NotAnIndividual columns below (instead of  manually in excel)

# Combine 2014-16 data into one data frame 
data_2014.2016 = rbind(data_2014.2015,data_2015.2016)

# Add 'NotAnIndividual' column
data_2014.2016 <- data_2014.2016 %>% 
  mutate(NotAnIndividual = ifelse(str_detect(OtherNotesCY, "lump"), 1, 
                           ifelse(str_detect(OtherNotesCY, "Lump"), 1,
                           ifelse(str_detect(OtherNotesCY, "split"), 1,
                           ifelse(str_detect(OtherNotesCY, "Split"), 1,
                           ifelse(str_detect(OtherNotesCY, "same as"), 1,
                           ifelse(str_detect(OtherNotesCY, "merge"), 1,
                           ifelse(str_detect(OtherNotesCY, "Merge"), 1,
                           ifelse(str_detect(OtherNotesPY, "part of"), 1, 
                           ifelse(str_detect(OtherNotesPY, "lump"), 1, 
                           ifelse(str_detect(OtherNotesPY, "Lump"), 1,
                           ifelse(str_detect(OtherNotesPY, "split"), 1,
                           ifelse(str_detect(OtherNotesPY, "Split"), 1,
                           ifelse(str_detect(OtherNotesPY, "same as"), 1,
                           ifelse(str_detect(OtherNotesPY, "merge"), 1,
                           ifelse(str_detect(OtherNotesPY, "Merge"), 1,
                           ifelse(str_detect(OtherNotesPY, "part of"), 1,       
                                0))))))))))))))))) 
# Note: this is lacking level 2 (=maybe), so is possibly more restrictive than 2010-14 filter
# TO DO: Inspect resulting data frame to make sure this is working as expected
# TO DO: repeat this automatic coding for 2010-2014 as a sensitivity analysis

# Add 'NotARecruit' column
data_2014.2016 <- data_2014.2016 %>% 
  mutate(NotARecruit = ifelse(str_detect(OtherNotesCY, "old"), 2, 
                       ifelse(str_detect(OtherNotesCY, "Old"), 2,
                       ifelse(str_detect(OtherNotesCY, "missed"), 2,
                       ifelse(str_detect(OtherNotesCY, "14?"), 2, 
                       ifelse(str_detect(OtherNotesCY, "15?"), 2, 
                       ifelse(!is.na(Size), NA, 
                       ifelse(NewPlot_CY==TRUE, 1, 0))))))))
# TO DO: Consult other queries (e.g., skipped in, exclusion areas) to identify rows that should be coded as level 1
# TO DO: Inspect resulting data frame to make sure this is working as expected
# Note: this is lacking level 3 (=size range of other recruits), which is not reliable

# Create columns of log-transformed sizes
data_2014.2016$logSize = log(data_2014.2016$Size)
data_2014.2016$logSizeNext = log(data_2014.2016$SizeNext)

# Add a column ranking regions from south to north
data_2014.2016$RegionRank[data_2014.2016$Region=="S1"]=1
data_2014.2016$RegionRank[data_2014.2016$Region=="S2"]=2
data_2014.2016$RegionRank[data_2014.2016$Region=="C1"]=3
data_2014.2016$RegionRank[data_2014.2016$Region=="C2"]=4
data_2014.2016$RegionRank[data_2014.2016$Region=="C3"]=5
data_2014.2016$RegionRank[data_2014.2016$Region=="N1"]=6
data_2014.2016$RegionRank[data_2014.2016$Region=="N2"]=7

# Create column for probability of flowering
data_2014.2016$Fec0 = (ifelse(data_2014.2016$Class=="A", 1, 
                              ifelse(data_2014.2016$Class=="J",0 , NA))) 

# Merge transition data with seed count data
data_2014.2016 <- merge(data_2014.2016,seed.ct,by="Site",all.x=TRUE,all.y=FALSE)

data_2014.2016 <- data_2014.2016 %>% select(colnames(data_2010.2014))

# Combine all years
data <- rbind(data_2010.2014, data_2014.2016)

# Make site x year variable
data$SiteYear = paste(data$Site, data$Year, sep=":") %>% factor()

#*******************************************************************************
#### 3. Remove unwanted data
#*******************************************************************************

# Remove plants that do not have site info
data <- filter(data, Site!="") 

# Remove plants that were dead at time t (Class=="D") 
data <- subset(data, Class!="D" | is.na(Class))

# Remove plants with Class=? or Class=excluded
data <- subset(data, Class != "E" & Class != "?" | is.na(Class))

# Remove rows for which size at time t AND size at t+1 is NA
data <- data[!(is.na(data$logSize) & is.na(data$logSizeNext)),]

# Remove plants that were recorded as having a class at time t but have no size measurements in that year
data <- subset(data, !(!is.na(Class) & is.na(logSize)))

# Remove plants that were recorded as having a class at time t but have no survival recorded from t to t+1 (were either excluded at time t+1, recorded as "?" in Class field, or recorded as "NA" in Class field) 
data <- subset(data, !(!is.na(Class) & is.na(Surv)))
# TO DO: check annotation; should it say ?/NA in Class or ClassNext? We have already removed ? and E at time t (lines 159-160). 

# Make Fec1 numeric and round fruit # to nearest integer
data$Fec1 <- round(as.numeric(data$Fec1, digits=0)) 

# Remove data where fruit # has errors resulting from 0s in denominator
data <- subset(data, Fec1 != "#Num!" | Fec1 != "#Div/0!" | is.na(Fec1)) 
# TO DO: Remove? I think this is unnecessary because they are forced to NA by conversion to numeric.

# Only include seed counts for plants that produced at least one fruit
data$SeedCt[data$Fec1 < 1 | is.na(data$Fec1)]=NA

# For most vital rate transitions, remove monster plants where individuals were not distinguishable 
# Note: these should still be included when calculating seed input denominator for recruitment probability
data.indivs = subset(data, NotAnIndividual != 1 | is.na(NotAnIndividual))

# For estimating size distribution of new recruits, also remove individuals that should not be scored as new recruits
data.indivs = subset(data.indivs, NotARecruit != 1 | is.na(NotARecruit))
# Note: because these rows are NA at time t, their exclusion has no bearing on survival, growth, and fecundity transitions. Confirmed:
summary(data.indivs$logSize[data.indivs$NotARecruit==1])


#*******************************************************************************
#### 4. Prepare data for IPMs and write to new .csv files
#*******************************************************************************

# Obtain total fruit and seed counts for each individual at each site in each year, including monster plants
site_fruit_count_data = subset(data, select=c(Site,Year,SiteYear,Region,Fec1,SeedCt)) 

# Examine column names and classes of data
names(data.indivs)
str(data.indivs)

# Check site names to be sure no discrepancies
unique(site_fruit_count_data$Site)
unique(data.indivs$Site)

# Check classes to be sure no discrepancies
unique(data.indivs$Class)
unique(data.indivs$ClassNext)

# Make appropriate columns of data frame a factor
site_fruit_count_data$Site=factor(site_fruit_count_data$Site)
site_fruit_count_data$Year=factor(site_fruit_count_data$Year)
site_fruit_count_data$Region=factor(site_fruit_count_data$Region)
data.indivs$Site=factor(data.indivs$Site)
data.indivs$Year=factor(data.indivs$Year)
data.indivs$Region=factor(data.indivs$Region)
data.indivs$ID=factor(data.indivs$ID)
data.indivs$NotARecruit=factor(data.indivs$NotARecruit)
data.indivs$NotAnIndividual=factor(data.indivs$NotAnIndividual)

# Examine data
names(data.indivs)
head(data.indivs)
tail(data.indivs)

# Sort data by latitude
data.indivs=data.indivs[order(-data.indivs$Latitude,data.indivs$Year),]

# write to .csv
write.csv(data.indivs,"data/demography data/Mcard_demog_data_2010-2016_cleanindivs.csv",row.names=FALSE)
write.csv(site_fruit_count_data,"data/demography data/Mcard_demog_data_2010-2016_seedinput.csv",row.names=FALSE)


#*******************************************************************************
#### 5. Description of columns in _cleanindivs .csv file
#*******************************************************************************

# Site: population
# ID: unique identifier for each individual
# Region: latitudinal region that population is nested within
# Latitude: latitude of population
# Longitude: longitude of population
# Elevation: elevation of population in meters
# Class: stage class (juvenile, adult, or NA) of plant at time = t (PY)
# Fec1: TotFr (Total number of fruits per individual, rounded to nearest integer)   
# logSize: total stem length of the individual, in log-transformed cm
# ClassNext: stage class (juvenile, adult, dead, or NA) of plant at time = t+1 (CY)
# logSizeNext: same as "size" above, for t+1
# Surv: survival (1) or not (0) of individuals between time = t (PY) and time = t+1 (CY)
# Year: annual transition of the long-term data at time = t (2010-2015)
# Fec0: Probability of flowering (1 if Class=="A" for adult, 0 if Class=="J" for juvenile)
# RegionRank: ordinal rank of regions from south to north
# SeedCt: mean seed count for each site
# NotAnIndividual: see above
# NotARecruit: see above
