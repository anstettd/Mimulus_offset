#*******************************************************************************
### 6. Plot latitude by lambda, with separate colors for each year
#*******************************************************************************

# regress lambda on latitude x year
model=lm(lambda~Latitude*Year,data=site.info)
summary(model)

# set graphing theme
theme_set(theme_minimal())

# make plot
ggplot(data = filter(site.info,SiteYear!="Deer Creek:2012"), aes(x = Latitude, y = lambda)) + geom_point(aes(color=Year))

# make list for plots of lambda
plot.lambda=list()

# make vector of unique sites
siteID=unique(site.info$Site) 

# make abbreviated year
site.info$Year_short=ifelse(site.info$Year==2010,10,ifelse(site.info$Year==2011,11,ifelse(site.info$Year==2012,12,ifelse(site.info$Year==2013,13,ifelse(site.info$Year==2014,14,ifelse(site.info$Year==2015,15,16))))))
site.info$Year_short=as.integer(site.info$Year_short)

# make plot
for (i in 1:length(siteID)) {
  site.lam=filter(site.info,Site==siteID[i]&!is.na(lambda))
  plot.lambda[[i]]=ggplot(data=site.lam,aes(x=Year_short,y=lambda)) + 
    geom_point(color="black",fill="grey",shape=21,size=2) + 
    ggtitle(paste(siteID[i])) +
    ylab(expression(lambda)) +
    xlab("Year") +
    theme(plot.title = element_text(size = 7),axis.title = element_text(size=8)) +
    xlim(10, 16)
}

pdf("Figures/lambda_site_year.pdf",width=11,height=8)
do.call(grid.arrange,plot.lambda)
dev.off()
