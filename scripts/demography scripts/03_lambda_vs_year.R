#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Calculate slopes of lambda versus year as a metric of the rate of demographic decline during drought
#### AUTHOR: Amy Angert
#### DATE LAST MODIFIED: 20230210


#*******************************************************************************
#### 0. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Make vector of packages needed
packages_needed <- c("tidyverse", "RColorBrewer")

# Install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# Load packages needed
for (i in 1:length(packages_needed)){
  library(packages_needed[i], character.only = TRUE)
}


#*******************************************************************************
### 1. Read in lambda estimates for each site and year
#*******************************************************************************
dat <- read.csv("data/demography data/siteYear.lambda_2010-2016.csv")

focal.sites <- c("CoastForkofWilliamette",
                 "CantonCreek",
                 "RockCreek",
                 "DeerCreek",
                 "O'NeilCreek",
                 "DeepCreek",
                 "LittleJamesonCreek",
                 "OregonCreek",
                 "RainbowPool",
                 "Carlon",
                 "BuckMeadows",
                 "Wawona",
                 "RedwoodCreek",
                 "NorthForkMiddleForkTule",
                 "SouthForkMiddleForkTule",
                 "WestForkMojaveRiver",
                 "MillCreek",
                 "WhitewaterCanyon",
                 "SweetwaterRiver",
                 "KitchenCreek",
                 "HauserCreek")

dat.old <- read.csv("data/demography data/siteYear.lambda_2010-2016_old.csv") %>% 
  filter(Site %in% focal.sites)


#*******************************************************************************
### 2. Visualize estimates over time for each site
#*******************************************************************************

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(dat$Site))
color.list <- lat_cols(n.sites)

ggplot(dat, aes(x=Year, y=lambda, color=as.factor(round(Latitude, 1)))) +
  geom_point() +
  geom_smooth(method="lm") +
  scale_color_manual(values=color.list) +
  ylab("Lambda") +
  #ylim(0,2) +
  geom_hline(yintercept=1, linetype="dotted") +
  facet_wrap(~Latitude, scale="free", nrow=3) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank())

ggplot(dat.old, aes(x=Year, y=lambda, color=as.factor(round(Latitude, 1)))) +
  geom_point() +
  geom_smooth(method="lm") +
  scale_color_manual(values=color.list) +
  ylab("Lambda") +
  #ylim(0,2) +
  geom_hline(yintercept=1, linetype="dotted") +
  facet_wrap(~Latitude, scale="free", nrow=3) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank())

#*******************************************************************************
### 3. Calculate slopes of lambda over time for each site
#*******************************************************************************

# with cleaned lambdas
site.vec <- unique(dat$Site)
slopes.lam <- c()
site.lam <- c()

for (i in 1:length(site.vec)) {
  dat.site <- dat %>% filter(Site==site.vec[i])
  mod <- lm(lambda ~ Year, dat.site)
  slopes.lam[i] <- coefficients(mod)[2]
  site.lam[i] <- site.vec[i]
  }

slopes.lambda <- bind_cols(site.lam, slopes.lam) %>% 
  select(Site=...1, Lambda.Slope.New=...2) %>% 
  mutate(Site = gsub(" ", "", Site))

# compare to preliminary older estimates
site.vec.old <- unique(dat.old$Site)
slopes.lam.old <- c()
site.lam.old <- c()

for (i in 1:length(site.vec.old)) {
  dat.site.old <- dat.old %>% filter(Site==site.vec.old[i])
  mod <- lm(lambda ~ Year, dat.site.old)
  slopes.lam.old[i] <- coefficients(mod)[2]
  site.lam.old[i] <- site.vec.old[i]
  }

slopes.lambda.old <- bind_cols(site.lam.old, slopes.lam.old) %>% 
  select(Site=...1, Lambda.Slope.Old=...2) 
  

slopes.all <- left_join(slopes.lambda, slopes.lambda.old)

ggplot(data=slopes.all, aes(x=Lambda.Slope.New, y=Lambda.Slope.Old, label=Site)) +
  geom_point() + 
  geom_text() + 
  geom_abline() +
  ylim(-3, 3)
  #xlim(-3, 3)
# TO DO ***High Priority*** Figure out what is driving Buck Meadows and Mill Creek estimates to be so different after cleaning
