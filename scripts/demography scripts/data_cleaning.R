#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Prepare tidied Mimulus cardinalis demography data for IPM analyses by filtering out problematic individuals
#### AUTHOR: Seema Sheth and Amy Angert
#### DATE LAST MODIFIED: 20230202

#*******************************************************************************
#### 1) Read in tidied data; convert certain columns to factors and sort in order of descending latitude
#*******************************************************************************

# read in data
data=read.csv("data/demography data/Mcard_demog_data_2010-2016.csv")


#*******************************************************************************
#### 2) Omit data that should not have been recorded as new recruits, based on Amy Angert's cleaning in July 2016 (for 2010-11, 2011-12, 2012-13, 2013-14 transitions) and script 00_data_prep.R (for 2014-15, 2015-16 transitions)
#*******************************************************************************

# Remove plants that should not have been recorded as new recruits
#### NOTE: these are plants that A. Angert noted as "wrong, definitely exclude (reasons include new plot, site not visited in prior year, ID within prior years' ranges, coordinates well outside of prior year's search)"

data = subset(data, NotARecruit != 1 | is.na(NotARecruit))

unique(data$NotARecruit)

length(data$Site) # 34760 rows; NOTE: there are 8 rows in which NotAnIndividual=1 & NotARecruit=1

#*******************************************************************************
#### 3) Save object of all remaining data that should go towards total fruit count and seed set per site and year
#*******************************************************************************

# Obtain total fruit and seed counts for each indivdiual at each site in each year, including monster plants
site_fruit_count_data = subset(data, select=c(Site,Year,SiteYear,Region,Fec1,SeedCt)) 
length(site_fruit_count_data$Site) # 34760 rows

#*******************************************************************************
#### 4) Remove data corresponding to plants that did not represent single unique individuals based on Amy Angert's cleaning in July 2016
#*******************************************************************************

# Remove monster plants where individuals were not distinguished
#### NOTE: these are plants that A. Angert noted as "not ok, definitely exclude from survival, growth, and fecundity but ok for seed input denominator for recruitment (history of lumping/splitting/relumping; redundant IDs)"

data = subset(data, NotAnIndividual != 1 | is.na(NotAnIndividual))

unique(data$NotAnIndividual)

length(data$Site) # 34699 rows

data = subset(data, select = -c(NotAnIndividual,NotARecruit))
