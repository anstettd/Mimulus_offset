#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Perform model selection for each vital rate for subsequent use in IPMs
#### AUTHOR: Seema Sheth and Amy Angert
#### DATE LAST MODIFIED: 20230809

#*******************************************************************************
#### 0. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Make vector of packages needed
packages_needed <- c("lme4", "MuMIn", "MASS", "pscl", "tidyverse", "glmmTMB", "R2admb") 

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
#### 2. Read in vital rate data ###
#*******************************************************************************

data <- read.csv("data/demography data/Mcard_demog_data_2010-2019_cleanindivs.csv") %>% 
  mutate(Period = ifelse(Year<2015&Year>2011, "drought",
                         ifelse(Year>2014, "recovery", NA))) %>% 
  filter(!is.na(Period))
# make factor for drought and recovery time periods
data$Site = factor(data$Site)
data$Year = factor(data$Year)

data.drought <- data %>% filter(Period=="drought")
data.recovery <- data %>% filter(Period=="recovery")

#*******************************************************************************
#### 2. Survival models ###
#*******************************************************************************

# Fixed effects model w/ and w/out size
s1 <- glm(Surv ~ logSize, data=data, family=binomial)
s2 <- glm(Surv ~ 1, data=data, family=binomial)
model.sel(s1, s2) # model w/ size is preferred

# Test whether vital rate dependence on size differs by period
sp <- glm(Surv ~ logSize*Period, data=data, family=binomial)
model.sel(s1, sp) #dependence on size differs by period
ggplot(data=data, aes(x=logSize, y=Surv, color=Period)) + 
  geom_point(position="jitter") +
  stat_smooth(method="glm", method.args = list(family=binomial))

## Drought period
# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
s3.d <- glmer(Surv  ~logSize + (logSize|Year/Site), data=data.drought, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: this model has a singularity warning

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
s4.d <- glmer(Surv ~ logSize  +(1|Year/Site), data=data.drought, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: this model has a singularity warning

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
s5.d <- glmer(Surv ~ logSize + (logSize|Year) + (logSize|Site), data=data.drought, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
s6.d <- glmer(Surv ~ logSize + (logSize|Year) + (1|Site), data=data.drought, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
s7.d <- glmer(Surv ~ logSize + (1|Year) + (logSize|Site), data=data.drought, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Compare models
anova(s3.d, s4.d, s5.d, s6.d, s7.d)
model.sel(s3.d, s4.d, s5.d, s6.d, s7.d) 

# PREFERRED MODEL IS s3.d, then s4.d, but due to singularity issues, we are going with the next best model that doesn't have singularity issues, s7.d
r.squaredGLMM(s7.d) 

# Save top survival model to .rda file
save(s7.d, file='data/demography data/surv.drought.reg.rda')   

## Recovery period
# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
s3.r <- glmer(Surv  ~logSize + (logSize|Year/Site), data=data.recovery, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: this model has a singularity warning

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
s4.r <- glmer(Surv ~ logSize  +(1|Year/Site), data=data.recovery, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
s5.r <- glmer(Surv ~ logSize + (logSize|Year) + (logSize|Site), data=data.recovery, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
s6.r <- glmer(Surv ~ logSize + (logSize|Year) + (1|Site), data=data.recovery, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
#NOTE: this model has a singularity warning

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
s7.r <- glmer(Surv ~ logSize + (1|Year) + (logSize|Site), data=data.recovery, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Compare models
anova(s3.r, s4.r, s5.r, s6.r, s7.r)
model.sel(s3.r, s4.r, s5.r, s6.r, s7.r) 

# PREFERRED MODEL IS s3.r, but due to singularity issues, we are going with the next best model that doesn't have singularity issues, s4.r
r.squaredGLMM(s4.r) 

# Save top survival model to .rda file
save(s4.r, file='data/demography data/surv.recovery.reg.rda')   

#*******************************************************************************
#### 4. Growth ###
#*******************************************************************************

# Fixed effects model w/ and w/out size
g1 <- glm(logSizeNext ~ logSize, data=data[!is.na(data$logSize),])
g2 <- glm(logSizeNext ~ 1, data=data[!is.na(data$logSize),])
model.sel(g1,g2) # model w/ size is preferred

# Test whether vital rate dependence on size differs by period
gp <- glm(logSizeNext ~ logSize*Period, data=data[!is.na(data$logSize),])
model.sel(g1, gp) #dependence on size differs by period
ggplot(data=data, aes(x=logSize, y=logSizeNext, color=Period)) + 
  geom_point(position="jitter") +
  stat_smooth(method="lm")

## Drought period
# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
g3.d <- lmer(logSizeNext ~ logSize + (logSize|Year/Site), data=data.drought, control=lmerControl(optimizer = "bobyqa")) 
#NOTE: singularity warning

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
g4.d <- lmer(logSizeNext ~ logSize + (1|Year/Site), data=data.drought, control=lmerControl(optimizer = "bobyqa")) 
#NOTE: singularity warning

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
g5.d <- lmer(logSizeNext ~ logSize + (logSize|Year) + (logSize|Site), data=data.drought, control=lmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
g6.d <- lmer(logSizeNext ~ logSize + (logSize|Year) + (1|Site), data=data.drought, control=lmerControl(optimizer = "bobyqa")) 

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
g7.d <- lmer(logSizeNext ~ logSize + (1|Year) + (logSize|Site), data=data.drought, control=lmerControl(optimizer = "bobyqa")) 

# Compare models
anova(g3.d, g4.d, g5.d, g6.d, g7.d)
model.sel(g3.d, g4.d, g5.d, g6.d, g7.d)

# # PREFERRED MODEL IS g3.d, then g4.d,  but these have  singularity issues so the next best model is g5.d
r.squaredGLMM(g5.d) 

# Save top growth model to .rda file
save(g5.d, file='data/demography data/growth.drought.reg.rda')   

## Recovery period
# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
g3.r <- lmer(logSizeNext ~ logSize + (logSize|Year/Site), data=data.recovery, control=lmerControl(optimizer = "bobyqa")) 
#NOTE: singularity warning

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
g4.r <- lmer(logSizeNext ~ logSize + (1|Year/Site), data=data.recovery, control=lmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
g5.r <- lmer(logSizeNext ~ logSize + (logSize|Year) + (logSize|Site), data=data.recovery, control=lmerControl(optimizer = "bobyqa")) 
#NOTE: singularity warning

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
g6.r <- lmer(logSizeNext ~ logSize + (logSize|Year) + (1|Site), data=data.recovery, control=lmerControl(optimizer = "bobyqa")) 

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
g7.r <- lmer(logSizeNext ~ logSize + (1|Year) + (logSize|Site), data=data.recovery, control=lmerControl(optimizer = "bobyqa")) 

# Compare models
anova(g3.r, g4.r, g5.r, g6.r, g7.r)
model.sel(g3.r, g4.r, g5.r, g6.r, g7.r)

# # PREFERRED MODEL IS g3.r, but this has a singularity issue so the next best model is g4.r
r.squaredGLMM(g4.r) 

# Save top growth model to .rda file
save(g4.r, file='data/demography data/growth.recovery.reg.rda')   

#*******************************************************************************
#### 5. Flowering ###
#*******************************************************************************

# fixed effects model w/ and w/out size
fl1 <- glm(Fec0 ~ logSize, data=data, family=binomial)
fl2 <- glm(Fec0 ~ 1, data=data, family=binomial)
model.sel(fl1, fl2) # model w/ size is preferred

# Test whether vital rate dependence on size differs by period
flp <- glm(Fec0 ~ logSize*Period, data=data, family=binomial)
model.sel(g=fl1, flp) #dependence on size differs by period
ggplot(data=data, aes(x=logSize, y=Fec0, color=Period)) + 
  geom_point(position="jitter") +
  stat_smooth(method="glm", method.args =list(family=binomial))

## Drought period
# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
fl3.d <- glmer(Fec0 ~ logSize + (logSize|Year/Site), data=data.drought, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: singularity warning

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
fl4.d <- glmer(Fec0 ~ logSize + (1|Year/Site), data=data.drought, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: this model has a singularity warning 

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fl5.d <- glmer(Fec0 ~ logSize + (logSize|Year) + (logSize|Site), data=data.drought, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: this model has a singularity warning 

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
fl6.d <- glmer(Fec0 ~ logSize + (logSize|Year) + (1|Site), data=data.drought, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: this model has a singularity warning 

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fl7.d <- glmer(Fec0 ~ logSize + (1|Year) + (logSize|Site), data=data.drought, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Compare models
anova(fl3.d, fl4.d, fl5.d, fl6.d, fl7.d)
model.sel(fl3.d, fl4.d, fl5.d, fl6.d, fl7.d) 

# # PREFERRED MODEL IS fl7.d
r.squaredGLMM(fl7.d) 

# Save top flowering model to .rda file
save(fl7.d, file='data/demography data/flowering.drought.reg.rda')   
## Recovery period
# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
fl3.r <- glmer(Fec0 ~ logSize + (logSize|Year/Site), data=data.recovery, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: singularity warning

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
fl4.r <- glmer(Fec0 ~ logSize + (1|Year/Site), data=data.recovery, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: this model has a singularity warning 

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fl5.r <- glmer(Fec0 ~ logSize + (logSize|Year) + (logSize|Site), data=data.recovery, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
fl6.r <- glmer(Fec0 ~ logSize + (logSize|Year) + (1|Site), data=data.recovery, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: this model has a singularity warning 

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fl7.r <- glmer(Fec0 ~ logSize + (1|Year) + (logSize|Site), data=data.recovery, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Compare models
anova(fl3.r, fl4.r, fl5.r, fl6.r, fl7.r)
model.sel(fl3.r, fl4.r, fl5.r, fl6.r, fl7.r) 

# # PREFERRED MODEL IS fl3.r or fl4.r, but due to singularity we go with fl5.r
r.squaredGLMM(fl5.r) 

# Save top flowering model to .rda file
save(fl5.r, file='data/demography data/flowering.recovery.reg.rda')  

#*******************************************************************************
#### 6. Fruit number ###
#*******************************************************************************
  
##6A. Fit fixed effects models (GLMs) only for initial exploratory model selection of variance structure (poisson vs. negative binomial vs. 0- inflation)
fr1_glm <- glm(Fec1 ~ logSize, data=data, na.action=na.omit, family=poisson) # poisson without 0-inflation 
fr2_glm <- glm.nb(Fec1 ~ logSize, data=data, na.action=na.omit) # negative binomial without 0-inflation
fr3_glm <- zeroinfl(Fec1 ~ logSize, data=data, na.action=na.omit, dist="poisson") # poisson with 0-inflation
fr4_glm <- zeroinfl(Fec1 ~ logSize, data=data, na.action=na.omit, dist="negbin") # negative binomial with 0-inflation

model.sel(fr1_glm, fr2_glm, fr3_glm, fr4_glm) 
	   
# PREFERRED MODEL IS fr2 (negative binomial w/out 0-inflation)

# Fixed effects model w/ and w/out size
fr5_glm <- glm.nb(Fec1 ~ 1, data=data, na.action=na.omit)
model.sel(fr2_glm, fr5_glm) # model w/ size is preferred
	  	
# Test whether vital rate dependence on size differs by period
frp <- glm.nb(Fec1 ~ logSize*Period, data=data, , na.action=na.omit)
model.sel(fr2_glm, frp) #dependence on size differs by time period
ggplot(data=data, aes(x=logSize, y=Fec1, color=Period)) + 
  geom_point(position="jitter") +
  stat_smooth(method="glm.nb")

##6B. Model selection of random effects structure

## Drought period
# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
fr3.d <- glmmTMB(Fec1 ~ logSize + (logSize|Year/Site), data=data.drought[!is.na(data.drought$Fec1),], family=nbinom1()) 

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
fr4.d <- glmmTMB(Fec1 ~ logSize + (1|Year/Site), data=data.drought[!is.na(data.drought$Fec1),], family=nbinom1()) 
      
# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fr5.d <- glmmTMB(Fec1 ~ logSize + (logSize|Year) + (logSize|Site), data=data.drought[!is.na(data.drought$Fec1),], family=nbinom1) 

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
fr6.d <- glmmTMB(Fec1 ~ logSize + (logSize|Year) + (1|Site), data=data.drought[!is.na(data.drought$Fec1),], family=nbinom1()) 

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fr7.d <- glmmTMB(Fec1 ~ logSize + (1|Year) + (logSize|Site), data=data.drought[!is.na(data.drought$Fec1),], family=nbinom1) 

# No random effects
fr8.d <- glmmTMB(Fec1 ~ logSize, data=data.drought[!is.na(data.drought$Fec1),], family=nbinom1) 

# Intercept only
fr9.d <- glmmTMB(Fec1 ~ 1, data=data.drought[!is.na(data.drought$Fec1),], family=nbinom1) 

# Compare models
anova(fr3.d, fr4.d, fr5.d, fr6.d, fr7.d)
model.sel(fr3.d, fr4.d, fr5.d, fr6.d, fr7.d, fr8.d, fr9.d) 

# PREFERRED MODEL IS fr3.d (no singularity issues)
r.squaredGLMM(fr3.d)

# Save top fruit # model to .rda file 
save(fr3.d, file='data/demography data/fruit.drought.reg.rda')   

## Recovery period
# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
fr3.r <- glmmTMB(Fec1 ~ logSize + (logSize|Year/Site), data=data.recovery[!is.na(data.recovery$Fec1),], family=nbinom1()) 
# NOTE: singularity warning

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
fr4.r <- glmmTMB(Fec1 ~ logSize + (1|Year/Site), data=data.recovery[!is.na(data.recovery$Fec1),], family=nbinom1()) 

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fr5.r <- glmmTMB(Fec1 ~ logSize + (logSize|Year) + (logSize|Site), data=data.recovery[!is.na(data.recovery$Fec1),], family=nbinom1) 
# NOTE: singularity warning

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
fr6.r <- glmmTMB(Fec1 ~ logSize + (logSize|Year) + (1|Site), data=data.recovery[!is.na(data.recovery$Fec1),], family=nbinom1()) 

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fr7.r <- glmmTMB(Fec1 ~ logSize + (1|Year) + (logSize|Site), data=data.recovery[!is.na(data.recovery$Fec1),], family=nbinom1) 

# No random effects
fr8.r <- glmmTMB(Fec1 ~ logSize, data=data.recovery[!is.na(data.recovery$Fec1),], family=nbinom1) 

# Intercept only
fr9.r <- glmmTMB(Fec1 ~ 1, data=data.recovery[!is.na(data.recovery$Fec1),], family=nbinom1) 

# Compare models
anova(fr3.r, fr4.r, fr5.r, fr6.r, fr7.r)
model.sel(fr3.r, fr4.r, fr5.r, fr6.r, fr7.r, fr8.r, fr9.r) 

# PREFERRED MODEL IS fr3.r then fr5.r, but due to singularity issues, going with fr7.r
r.squaredGLMM(fr7.r)

# Save top fruit # model to .rda file 
save(fr7.r, file='data/demography data/fruit.recovery.reg.rda')   
