#### PROJECT: Mimulus cardinalis demography 2010-2014
#### PURPOSE: Perform model selection for each vital rate for subsequent use in IPMs
############# Vital rates include survival, growth, flowering, and fruit count
############# Fixed effect: size; Random effects: site, year
#### AUTHOR: Seema Sheth
#### DATE LAST MODIFIED: 20180613

# remove objects and clear workspace
rm(list = ls(all=TRUE))

# Install GLMMADMB package following instructions here: http://glmmadmb.r-forge.r-project.org/
# install.packages("R2admb")
# install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos",getOption("repos")),type="source")

#require packages
require(lme4)
#require(glmmADMB)
require(MuMIn)
require(MASS)
require(pscl)
require(tidyverse)
require(glmmTMB)
require(bbmle)

# set working directory
setwd("/Users/ssheth/Google Drive/2018_demography_climate")

#*******************************************************************************
#### 1. run data_prep.R script to clean up data ###
#*******************************************************************************

source("R_scripts/data_prep.R")

# Variables are: 

# Site: population
# ID: unique identifier for each individual
# Region: latitudinal region that population is nested within
# Latitude: latitude of population
# Longitude: longitude of population
# Elevation: elevation of population
# Class: stage class (juvenile, adult, or NA) of plant at time = t 
# Fec1: Total number of fruits per individual   
# logSize: total stem length of the individual
# ClassNext: stage class (juvenile, adult, dead, or NA) of plant at time = t+1 
# logSizeNext: same as "logSize" above, for t+1
# Surv: survival (1) or not (0) of individuals between time = t and time = t+1 
# Year: annual transition of the long-term data at time = t (2010-2013)
# Fec0: Probability of flowering (1 if Class=="A" for adult, 0 if Class=="J" for juvenile)
# RegionRank: ordinal rank of regions from south to north
# SeedCt: mean seed count, rounded to the nearest integer, for each site

#*******************************************************************************
#### 2. Summarize sample sizes per site per year for each vital rate
#*******************************************************************************

# survival
data_summary_surv=data %>% group_by(SiteYear) %>% filter(!is.na(Surv)) %>% summarise(nSurv=n())
# growth
data_summary_growth=data %>% group_by(SiteYear) %>% filter(!is.na(logSize)&!is.na(logSizeNext))  %>% summarise(nGrowth=n())
# flowering
data_summary_fl=data %>% group_by(SiteYear) %>% filter(!is.na(Fec0)) %>% summarise(nFl=n())
# fruit #
data_summary_fr=data %>% group_by(SiteYear) %>% filter(!is.na(Fec1))  %>% summarise(nFr=n())

# create data frame of site info to join to sample size summary
site.info=data %>% select(Site,Latitude,Year,SiteYear) %>% distinct()

# join summaries together, sort by decreasing latitude, and write to .csv file
data_summary=full_join(data_summary_surv,data_summary_growth,by="SiteYear") %>% full_join(data_summary_fl,by="SiteYear") %>% full_join(data_summary_fr,by="SiteYear") %>% full_join(site.info,by="SiteYear") %>% arrange (-Latitude,Year) %>% select(-Latitude)
# write.csv(data_summary,"R_output/data_summary.csv",row.names = FALSE)

#*******************************************************************************
#### 3. Survival ###
#*******************************************************************************

# fixed effects model w/ and w/out size
s1=glm(Surv~logSize,data=data,family=binomial)
s2=glm(Surv~1,data=data,family=binomial)
model.sel(s1,s2) # model w/ size is preferred

# random intercepts and slopes for Year; random intercepts & slopes for Site nested within Year
s3=glmer(Surv~logSize+(logSize|Year/Site),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# random intercepts and slopes for Year; random intercepts & constant slopes for Site nested within Year
s4=glmer(Surv~logSize+(1|Year/Site),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# Compare models
anova(s3, s4)

model.sel(s3, s4) 

# PREFERRED MODEL IS s3
r.squaredGLMM(s3) 

# Save top survival model to .rda file
save(s3, file='R_output/surv.reg.rda')   

#*******************************************************************************
#### 4. Growth ###
#*******************************************************************************

# fixed effects model w/ and w/out size
g1=glm(logSizeNext~logSize,data=data[!is.na(data$logSize),])
g2=glm(logSizeNext~1,data=data[!is.na(data$logSize),])
model.sel(g1,g2) # model w/ size is preferred

# A. random intercepts and slopes for Year; random intercepts & random slopes for Site nested within Year
g3=lmer(logSizeNext~logSize+(logSize|Year/Site),data=data,control=lmerControl(optimizer = "bobyqa")) 

# B. random intercepts and slopes for Year; random intercepts & constant slopes for Site nested within Year
g4=lmer(logSizeNext~logSize+(1|Year/Site),data=data,control=lmerControl(optimizer = "bobyqa")) 

# Compare models
anova(g3, g4)
model.sel(g3, g4)

# PREFERRED MODEL IS g3
r.squaredGLMM(g3) 

# Save top growth model to .rda file
save(g3, file='R_output/growth.reg.rda')   

#*******************************************************************************
#### 5. Flowering ###
#*******************************************************************************

# fixed effects model w/ and w/out size
fl1=glm(Fec0~logSize,data=data,family=binomial)
fl2=glm(Fec0~1,data=data,family=binomial)
model.sel(fl1,fl2) # model w/ size is preferred

# A. random intercepts and slopes for Year; random intercepts & slopes for Site nested within Year
fl3=glmer(Fec0~logSize+(logSize|Year/Site),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# B. random intercepts and constant slopes for Year; random intercepts & constant slopes for Site nested within Year
fl4=glmer(Fec0~logSize+(1|Year/Site),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# Compare models
anova(fl3, fl4)
model.sel(fl3, fl4) 

# PREFERRED MODEL IS fl3
r.squaredGLMM(fl3) 
r.squaredGLMM(fl4)

# Save top flowering model to .rda file
save(fl3, file='R_output/flowering.reg.rda')   

#*******************************************************************************
#### 6. Fruit number ###
#*******************************************************************************
  
	   #*******************************************************************************
	   ####6A. Fit fixed effects models (GLMs) only for initial exploratory model selection of variance structure (poisson vs. negative binomial vs. 0- inflation)
	   #*******************************************************************************
	  
	   fr1_glm=glm(Fec1~logSize,data=data,na.action=na.omit,family=poisson) # poisson without 0-inflation 
	   fr2_glm=glm.nb(Fec1~logSize,data=data,na.action=na.omit) # negative binomial without 0-inflation
	   fr3_glm=zeroinfl(Fec1~logSize,data=data,na.action=na.omit,dist="poisson") # poisson with 0-inflation
	   fr4_glm=zeroinfl(Fec1~logSize,data=data,na.action=na.omit,dist="negbin") # negative binomial with 0-inflation
	   model.sel(fr1_glm,fr2_glm,fr3_glm,fr4_glm) 
	   
	   # PREFERRED MODEL IS fr2 (negative binomial w/out 0-inflation)
	   # fixed effects model w/ and w/out size
	  	fr5_glm=glm.nb(Fec1~1,data=data,na.action=na.omit)
	  	model.sel(fr2_glm,fr5_glm) # model w/ size is preferred
	  	
		#*******************************************************************************
		####6B. Model selection of random effects structure
	    ####### NOTE: when I ran same models with glmmadmb() function and family="poisson" none of the models except for fr4 converged
	    ####### NOTE: I also ran same models with glmer.nb; computationally faster but less widely used in literature
		#*******************************************************************************
		
	  	# A. random intercepts and slopes for Year; random intercepts & slopes for Site nested within Year
	  	fr3=glmmadmb(Fec1~logSize+(logSize|Year/Site),data=data[!is.na(data$Fec1),],family="nbinom1",link="log") 
	  	
	  	# B. random intercepts and constant slopes for Year; random intercepts & constant slopes for Site nested within Year
	  	fr4=glmmadmb(Fec1~logSize+(1|Year/Site),data=data[!is.na(data$Fec1),],family="nbinom",link="log") 
      
	  	# Compare models
	  	anova(fr3, fr4)
	  	
	  	AICc(fr3, fr4) 
	  	
	  	model.sel(fr3, fr4) 
	  	
		# PREFERRED MODEL IS fr3
		
		# Save top fruit # model to .rda file because it takes a long time to run
		# save(fr3, file='R_output/fruit.reg.rda')   
    
		# Try glmmTMB, which is must faster than glmmADMB
		fr5=glmmTMB(Fec1~logSize+(logSize|Year/Site),data=data[!is.na(data$Fec1),],family=nbinom1()) 
		fr6=glmmTMB(Fec1~logSize+(logSize|Year/Site),data=data[!is.na(data$Fec1),],family=nbinom2()) 
		fr7=glmmTMB(Fec1~logSize+(1|Year/Site),data=data[!is.na(data$Fec1),],family=nbinom1()) # many warnings
		fr8=glmmTMB(Fec1~logSize+(1|Year/Site),data=data[!is.na(data$Fec1),],family=nbinom2()) 
		
		# compare AIC (seems unnecessary given other models did not converge)
		AICtab(fr5,fr6,fr7,fr8) 
		# fr6 is preferred; it's worth trying this since it'll make bootstrapping go much faster
		
		# PREFERRED MODEL IS fr6
		# run script to create function for obtaining R^2 from glmmTMB object
		source("R_scripts/rsqglmmTMB.R")
		my_rsq(fr6) 
		summary(fr6)
		
		# Save top fruit # model to .rda file because it takes a long time to run
		save(fr6, file='R_output/fruit.reg.rda')   
		
		
		
		
		
		
		
		
		
		
		
		
		
				
### JUNK CODE
		# random intercepts and random slopes for Site; random intercepts and slopes for Year; random intercepts & random slopes for Site nested within Year
		s5=glmer(Surv~logSize+(logSize|Year/Site)+(logSize|Site),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 
		
		# random intercepts and constant slopes for Site; random intercepts and slopes for Year; random intercepts & random slopes for Site nested within Year
		s6=glmer(Surv~logSize+(logSize|Year/Site)+(1|Site),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 
		
		# random intercepts and slopes for Site; random intercepts and constant slopes for Year; random intercepts & constant slopes for Site nested within Year
		s7=glmer(Surv~logSize+(1|Year/Site)+(logSize|Site),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) # did not converge
		
		# random intercepts and constant slopes for Site; random intercepts and constant slopes for Year; random intercepts & constant slopes for Site nested within Year
		s8=glmer(Surv~logSize+(1|Year/Site)+(1|Site),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 
		