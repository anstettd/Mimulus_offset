#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Prepare raw demographic vital rate data for IPM analyses 
#### AUTHOR: Seema Sheth and Amy Angert
#### DATE LAST MODIFIED: 20230201

# Easy code for installing packages in R (if not installed) and calling their libraries
# From: https://gist.github.com/DrK-Lo/a945a29d6606b899022d0f03109b9483

# make vector of packages needed
packages_needed <- c("tidyverse")

# install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

#*******************************************************************************
#### 1. Bring in and combine M. cardinalis demography data from 2010-2014; in all files, PY=previous year, CY=current year; ignore PPY
#*******************************************************************************

# read in transition data for 2010-11, 2011-12, 2012-13, 2013-14  
data_2010.2014=read.csv("data/demography data/Mcard_demog_data_2010-2014.csv")
# Note: this file was created for the analyses published in Sheth and Angert 2018 PNAS paper
# It results from Amy Angert's work in July 2016 (original file: "Mcard_demog_data_2010-2013_ALA.xlsx") to scan datasheet notes to identify individuals to exclude, based on these columns:
# Column 'NotAnIndividual': 
# 0 = ok, definitely include (includes "scattered" as long as nothing else noted)
# 1 = not ok, definitely exclude from survival, growth, and fecundity but ok for seed input denominator for recruitment (history of lumping/splitting/relumping; redundant IDs)
# 2 = maybe (notes about difficult to distinguish, or merged once)
# Column 'NotARecruit':
# 0 = ok, definitely include
# 1 = wrong, definitely exclude (reasons include new plot, site not visited in prior year, ID within prior years' ranges, coordinates well outside of prior year's search
# 2 = plant noted as possibly missed in prior year (looks old, missed?, J-20xx?, could be [old ID], etc)
# 3 = plant not noted as unusually old-looking, but is within size range of or larger than most plants that are in category 2
# NA = size measures at year t


# read in transition data for 2014-15 and add year at time = t (2014) column
data.2014.2015=read.csv("data/demography data/SS_Horizontal_2015_Notes.csv")
data.2014.2015$Year=rep("2014",times=length(data.2014.2015$ID))
# note: this file was queried from the database on 2023-02-01
# we will tidy it based on the decision rules as above for 2010-14 data, but scripted here instead of done manually in excel. 

# add 'NotAnIndividual' column
data.2014.2015 <- data.2014.2015 %>% 
  mutate(NotAnIndividual=ifelse(str_detect(OtherNotesCY, "lump"), 1, 
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
# Note: this is lacking level=2 (=maybe), so is possibly more restrictive than 2010-14 filter

# add 'NotARecruit' column
data.2014.2015 <- data.2014.2015 %>% 
  mutate(NotARecruit=ifelse(str_detect(OtherNotesCY, "old"), 1, 
                     ifelse(str_detect(OtherNotesCY, "Old"), 1,
                     ifelse(str_detect(OtherNotesCY, "missed"), 1,
                     ifelse(str_detect(OtherNotesCY, "14?"), 1, 
                     ifelse(TotStLn_PY, 1)
                            0)))))

# read in 2015-2016 data and year at time = t (2015) column
data.2015.2016=read.csv("data/demography data/SS_Horizontal_2016_Notes.csv")
data.2015.2016$Year=rep("2015",times=length(data.2015.2016$ID))
# note: this file was queried from the database on 2023-02-01
# we will tidy it based on the same decision rules as above for 2010-14 data, but scripted here instead of done manually in excel

# add 'NotAnIndividual' column
data.2015.2016 <- data.2015.2016 %>% 
  mutate(NotAnIndividual=ifelse(str_detect(OtherNotesCY, "lump"), 1, 
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
# Note: this is lacking the level=2 (=maybe), so is possibly more restrictive than 2010-14 filter


# combine data from all years into one data frame 
data=rbind(data.2014.2015,data.2015.2016,data.2016.2017,data.2017.2018, data.2018.2019)

# Read in seed count per fruit data, select and rename relevant columns, and round to nearest integer
seed.ct=read.csv("data/vital_rates/fec2.seed.per.fruit.2010.2011.2012.csv")
seed.ct=subset(seed.ct,select=c(site,newgrandmean))
colnames(seed.ct)=c("Site","SeedCt")
seed.ct$SeedCt=round(seed.ct$SeedCt,digits=0) # round seed count to nearest integer
seed.ct$Site=factor(seed.ct$Site) # make site column a factor to streamline joining

#*******************************************************************************
#### 2. Remove unwanted data
#*******************************************************************************

# omit plants that do not have site info
data=filter(data,SiteID!="")

# omit plants that were dead in previous year (Class_PY=="D") 
data=subset(data,Class_PY!="D"|is.na(Class_PY))

# remove data for which size in previous year (PY) AND size in current year (CY) is NA
data=data[!(is.na(data$TotStLn_PY)&is.na(data$TotStLn_CY)),]

# remove data from plots that were new in current year (CY)
data=subset(data,NewPlot_CY==FALSE) 

# get rid of plants with Class_PY=? or Class_PY= excluded
data=subset(data, Class_PY!="E"&Class_PY!="?"|is.na(Class_PY))

# get rid of plants that were recorded as having a class  in previous year (PY) but have no size measurements in previous year (PY)
data=subset(data,!(!is.na(Class_PY)&is.na(TotStLn_PY)))

# get rid of plants that were recorded as having a class in previous year (PY) but have no survival recorded from previous to current year (were either excluded in current year, recorded as "?" in Class_CY field, or recorded as "NA" in Class_CY field)
data=subset(data,!(!is.na(Class_PY)&is.na(SurvPYCY)))

# ID 75412 from Rainbow Pool in 2017 has TotStLn_CY = 0. This must be a typo and we will exclude for now
data=subset(data,TotStLn_CY!=0|is.na(TotStLn_CY))

# Remove data whereTotFr_PY either equal “#Num!” or “#Div/0!”
data=subset(data,TotFr_PY!="#Num!"|TotFr_PY!="#Div/0!"|is.na(TotFr_PY))
## WARNING: THIS NEEDS TO BE FIXED!!! arises when Per_PY = 0
# there are similar issues with TotFr_CY but this column is irrelevant

#*******************************************************************************
#### 3. Prepare data for IPMs and write to new .csv file
#*******************************************************************************

# examine column names and classes of data
names(data)
str(data)

# select certain columns
data=subset(data,select=c("SiteID","ID","Region","Latitude","Longitude","Elevation","Class_PY","TotFr_PY","TotStLn_PY","Class_CY","TotStLn_CY","SurvPYCY","Year"))

# rename column names to IPM friendly names
colnames(data)=c("Site","ID","Region","Latitude","Longitude","Elevation","Class","Fec1","logSize","ClassNext","logSizeNext","Surv","Year")

# create additional  vectors for analysis (see descriptions below)
data$Fec0=(ifelse(data$Class=="A",1,ifelse(data$Class=="J",0,NA))) #create vector for probability of flowering, where if adult, probability of flowering is 1, if juvenile, probability of flowering is 0, otherwise probability of flowering is NA

# make Fec1 numeric
data$Fec1=round(as.numeric(data$Fec1,digits=0)) #round total fruit # to nearest integer
## WARNING: I ORIGINALLY GOT AN ERROR HERE BECAUSE SOME VALUES = "#Num!" AND OTHERS = "#Div/0!" SO COLUMN WAS READ AS CHARACTER

# log transform size and sizeNext
data$logSize=log(data$logSize)
data$logSizeNext=log(data$logSizeNext)

# Add a column ranking regions from south to north
data$RegionRank[data$Region=="S1"]=1
data$RegionRank[data$Region=="S2"]=2
data$RegionRank[data$Region=="C1"]=3
data$RegionRank[data$Region=="C2"]=4
data$RegionRank[data$Region=="C3"]=5
data$RegionRank[data$Region=="N1"]=6
data$RegionRank[data$Region=="N2"]=7

# 2016 and 2017 and 2018 and 2019 data do not have spaces between words in site names; fix this across the dataset
unique(data$Site)
data$Site[data$Year==2016|data$Year==2017|data$Year==2018|data$Year==2019]=gsub("(?!^)(?=[[:upper:]])", " ", data$Site[data$Year==2016|data$Year==2017|data$Year==2018|data$Year==2019], perl=T)
unique(data$Site) # Coast Forkof Williamette and O' Neil Creek need to be corrected
data$Site[data$Site=="Coast Forkof Williamette"]="Coast Fork of Williamette" 
data$Site[data$Site=="O' Neil Creek"]="O'Neil Creek" 

# Make appropriate columns of data frame a factor
data$Site=factor(data$Site)
data$Region=factor(data$Region)

# Merge data with seed count data
data=merge(data,seed.ct,by="Site",all.x=TRUE,all.y=FALSE)

# Only include seed counts for plants that produced at least one fruit
data$SeedCt[data$Fec1<1|is.na(data$Fec1)]=NA

# Add columns so that 2014-2019 data can be combined with 2010-2013 data
data$NotAnIndividual=""
data$Reasoning=""
data$NotARecruit=""
data$Reasoning.1=""

# examine data
names(data)
head(data)
tail(data)

# combine 2010-2013 data with 2014-2019 data
data=rbind(data,data_2010.2013)

# sort data by latitude
data=data[order(-data$Latitude,data$Year),]

# check site names one more time
unique(data$Site)

# write to .csv
write.csv(data,"output/vital_rates/Mcard_demog_data_2010-2018.csv",row.names=FALSE)


#*******************************************************************************
#### 4. Filter to 21 focal sites for LTREB
#******************************************************************************

focal.sites <- c("Coast Fork of Williamette",
                 "Canton Creek",
                 "Rock Creek",
                 "Deer Creek",
                 "O'Neil Creek",
                 "Deep Creek",
                 "Little Jameson Creek",
                 "Oregon Creek",
                 "Rainbow Pool",
                 "Carlon",
                 "Buck Meadows",
                 "Wawona",
                 "Redwood Creek",
                 "North Fork Middle Fork Tule",
                 "South Fork Middle Fork Tule",
                 "West Fork Mojave River",
                 "Mill Creek",
                 "Whitewater Canyon",
                 "Sweetwater River",
                 "Kitchen Creek",
                 "Hauser Creek")

data.focal <- data %>% filter(Site%in%focal.sites) %>% droplevels()
unique(data.focal$Site) #21

write.csv(data.focal,"output/vital_rates/Mcard_demog_data_2010-2018_focal.csv",row.names=FALSE)




#*******************************************************************************
#### 4. Description of columns in new .csv file
#*******************************************************************************

#Site: population
#ID: unique identifier for each individual
#Region: latitudinal region that population is nested within
#Latitude: latitude of population
#Longitude: longitude of population
#Elevation: elevation of population
#Class: stage class (juvenile, adult, or NA) of plant at time = t (PY)
#Fec1: TotFr (Total number of fruits per individual)   
#logSize: total stem length of the individual
#ClassNext: stage class (juvenile, adult, dead, or NA) of plant at time = t+1 (CY)
#logSizeNext: same as "size" above, for t+1
#Surv: survival (1) or not (0) of individuals between time = t (PY) and time = t+1 (CY)
#Year: annual transition of the long-term data at time = t (2010-2013)
#Fec0: Probability of flowering (1 if Class=="A" for adult, 0 if Class=="J" for juvenile)
#RegionRank: ordinal rank of regions from south to north
#SeedCt: mean seed count for each site
#NotAnIndividual: see above
#Reasoning: justification for NotAnIndividual
#NotARecruit: see above
#Reasoning.1: justification for NotARecruit