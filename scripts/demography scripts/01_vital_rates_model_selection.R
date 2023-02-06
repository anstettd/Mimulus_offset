#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Perform model selection for each vital rate for subsequent use in IPMs
#### AUTHOR: Seema Sheth and Amy Angert
#### DATE LAST MODIFIED: 20230206

#*******************************************************************************
#### 0. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Make vector of packages needed
packages_needed <- c("lme4", "MuMIn", "MASS", "pscl", "tidyverse", "glmmTMB", "bbmle", "R2admb", "glmmADMB")

# Install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}
# Note: might need to install from source:
# install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos",getOption("repos")),type="source")

# Load packages needed
for (i in 1:length(packages_needed)){
  library(packages_needed[i], character.only = TRUE)
}



#*******************************************************************************
#### 2. Read in vital rate data ###
#*******************************************************************************
data <- read.csv("data/demography data/Mcard_demog_data_2010-2016_.csv")

#*******************************************************************************
#### 2. Survival models ###
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

# random intercepts and slopes for Year; random intercepts & random slopes for Site nested within Year
g3=lmer(logSizeNext~logSize+(logSize|Year/Site),data=data,control=lmerControl(optimizer = "bobyqa")) 

# random intercepts and slopes for Year; random intercepts & constant slopes for Site nested within Year
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

# random intercepts and slopes for Year; random intercepts & slopes for Site nested within Year
fl3=glmer(Fec0~logSize+(logSize|Year/Site),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# random intercepts and constant slopes for Year; random intercepts & constant slopes for Site nested within Year
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
		
	  	# random intercepts and slopes for Year; random intercepts & slopes for Site nested within Year
	  	fr3=glmmadmb(Fec1~logSize+(logSize|Year/Site),data=data[!is.na(data$Fec1),],family="nbinom1",link="log") 
	  	
	  	# random intercepts and constant slopes for Year; random intercepts & constant slopes for Site nested within Year
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
