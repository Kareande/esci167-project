##################****** 1992-2011 ******##################
setwd("/home/kareande/lchl-mhw-events")
library("quantreg") #QR modeling
library("dplyr") #tidymodels
library("tidyr")
library("lubridate") #date and time manipulation
library("lmtest") #DW test for independence
#install.packages()

################## Data Processing ##################
# Get data from senior thesis project (OISST, NCEP, CMEMS)
#cmp_df = read.csv("cmpndData/cmpnd_trn_qrm.csv",header=T) #balanced df
cmp_df = read.csv("cmpndData/cmpnd_trn_unbal_qrm.csv",header=T) #entire df
cmp_df$date <- as.Date(cmp_df$ds, origin = "1992-01-01") #convert days since to date
cmp_df <- cmp_df %>% dplyr::mutate(year = lubridate::year(cmp_df$date), #seperate date into year,
                         month = lubridate::month(cmp_df$date),         #month, and
                         day = lubridate::day(cmp_df$date))             #day columns
cmp_df <- cmp_df[ , ! names(cmp_df) %in% c("yr","lon","lat","ds","cmpCat","date")] #remove columns
cmp_df <- cmp_df %>% relocate(c(day, month, year), .before = doy) #move date vars to front
head(cmp_df)

# Calculate monthly means
cmp_mo = aggregate(cmp_df[,6:17],
                   by=list(cmp_df$month, cmp_df$year, cmp_df$location), #aggregate by month, year, and location
                   FUN=mean, data=cmp_df, na.rm=TRUE)
cmp_mo <- cmp_mo[order(cmp_mo$Group.2),] #reorder dataset by year
cmp_mo <- cmp_mo[order(cmp_mo$Group.1),] #reorder dataset by month
names(cmp_mo)[names(cmp_mo) == "Group.1"] <- "mo" #rename column
names(cmp_mo)[names(cmp_mo) == "Group.2"] <- "yr" #rename column
names(cmp_mo)[names(cmp_mo) == "Group.3"] <- "location" #rename column
cmp_mo$date <- gsub(" ", "", paste(cmp_mo$yr, "-", cmp_mo$mo, "-01")) #merge month and year
cmp_mo$date <- ymd(cmp_mo$date) #get date column
cmp_mo$date <- decimal_date(cmp_mo$date) #convert date to decimal date
cmp_mo <- cmp_mo %>% relocate(date, .before = mo) #move date var to front
cmp_mo <- cmp_mo[order(cmp_mo$date),] #reorder dataset by date
cmp_mo <- cmp_mo %>% separate(location,c("lat","lon"),sep="-") #un-merge location column into latitude and longitude
cmp_mo$lat <- as.integer(cmp_mo$lat) #convert from char to int
cmp_mo$lon <- as.integer(cmp_mo$lon) #convert from char to int
head(cmp_mo)

# Visualize some data
setwd("/home/kareande/esci167-project")
sst = cmp_mo$sst
chl = cmp_mo$chl
time = cmp_mo$date
png(gsub(" ", "", paste("figs/SSTvsChl.png")))
plot(sst,chl,xlab="Sea Surface Temperature Anomaly (deg C)",ylab="Chlorophyl Concentration (chl-a/m^3)",
    main="Monthly Sea Surface Temperature Anomaly vs
    Chlorophyl Concentration")
dev.off()
png(gsub(" ", "", paste("figs/ChlvsTime.png")))
plot(time,chl,xlab="Time (Decimal Date)",ylab="Chlorophyl Concentration (chl-a/m^3)", type="l",
    main="Monthly Chlorophyl Concentration")
dev.off()
png(gsub(" ", "", paste("figs/SSTvsTime.png")))
plot(time,sst,xlab="Time (Decimal Date)",ylab="Sea Surface Temperature Anomaly (deg C)", type="l",
    main="Monthly Sea Surface Temperature Anomaly")
dev.off()

################## Fit OLS Model ##################
# all vars significant
# poor fit for extreme obs, autocorrelation

# Fit ordinary least squares (OLS) model
model = lm(sst ~ date + lat + lon +
             nit + sil + oxy + pho + npp + chl +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo)
res = model$residuals

# Plot time series and fit
png(gsub(" ", "", paste("figs/SST_OLS.png")))
plot(time,sst,xlab="Time (Years)",ylab="Sea Surface Temperature Anomaly", type="l",
    main="Monthly Sea Surface Temperature Anomaly with OLS Fit")
lines(time,model$fitted.values,col="red")
dev.off()

# Check independance of OLS model
png(gsub(" ", "", paste("figs/SST_OLS_res.png")))
plot(res, ylab="OLS Residuals",
     main="Monthly Sea Surface Temperature Anomaly OLS Residuals") #seasonal pattern
dev.off()
dwtest(model) #p-value < 2.2e-16, autocorrelation
summary(model)

# Check normality and constant variance
# for funsies
a = length(res)
b = round((length(res)-4998)/2)
c = a-b
res_sh = res[b:c]
shapiro.test(res_sh) #p-value = 5.235e-11, NOT normal

midres = length(res)/2
bartlett.test(res,c(rep(1,midres),
                     rep(2,midres))) #p-value < 2.2e-16, NOT constant

# Cochrane-transformation
# still not independent (p-value < 2.2e-16)
n = length(cmp_mo$sst)
rhohat = sum(model$residuals[-1]*model$residuals[-n])/sum(model$residuals^2)
sst.t = cmp_mo$sst[-1]-rhohat*cmp_mo$sst[-n]
model.t = lm(sst.t ~ date[-1] + lat[-1] + lon[-1] +
             nit[-1] + sil[-1] + oxy[-1] + pho[-1] + npp[-1] + chl[-1] +
             qnet[-1] + slp[-1] + sat[-1] + wndsp[-1] + sstRoC[-1], data=cmp_mo)

# Check independance of transformed model
res.t = model.t$residuals
png(gsub(" ", "", paste("figs/SST_OLS_Coch_res.png")))
plot(res.t, ylab="OLS Residuals",
     main="Monthly SST Anomaly Cochrane-Transformed
    OLS Residuals") #clustered around 0
dwtest(model.t) #p-value < 2.2e-16, autocorrelation
dev.off()

# Cochrane-transformation x100
# still not independent (p-value < 2.2e-16)
n = length(cmp_mo$sst)-100
rhohat = sum(model$residuals[-100]*model$residuals[-n])/sum(model$residuals^2)
sst.t2 = cmp_mo$sst[-100]-rhohat*cmp_mo$sst[-n]
model.t2 = lm(sst.t2 ~ date[-100] + lat[-100] + lon[-100] +
             nit[-100] + sil[-100] + oxy[-100] + pho[-100] + npp[-100] + chl[-100] +
             qnet[-100] + slp[-100] + sat[-100] + wndsp[-100] + sstRoC[-100], data=cmp_mo)

# Check independance of seriously transformed model
res.t2 = model.t2$residuals
plot(res.t2, ylab="OLS Residuals",
     main="Monthly SST Anomaly 100x Cochrane-Transformed
    OLS Residuals") #clustered around 0
dwtest(model.t) #p-value < 2.2e-16, autocorrelation

################## Fit GLS Model ##################
setwd("/home/kareande/esci167-project")
library("nlme") #acf, pacf
library("trend") #mk.test
library("dplyr") #tidymodels language
#install.packages()

# Determine order of autocorrelation
png(gsub(" ", "", paste("figs/SST_PACF.png")))
par(mfrow=c(1,2))
acf(res, main="Monthly SST Series ACF")
pacf(res, main="Monthly SST Series PACF") #2 or 28 orders of autocorrelation
dev.off()

# Fit GLS model to bypass independence assumption #DO OVERNIGHT
sst.gls = gls(sst ~ date + lat + lon +
             nit + sil + oxy + pho + chl +
             qnet  + sat + wndsp + sstRoC, data=cmp_mo, correlation=corARMA(p=28))
png(gsub(" ", "", paste("figs/SST_GLS.png")))
plot(time,sst,xlab="Time (Decimal Date)",ylab="Sea Surface Temperature Anomoly (C)", type="l",
    main="Monthly Sea Surface Temperature Anomaly with GLS Fit")
lines(time,sst.gls$fitted.values,col="red")
dev.off()
summary(sst.gls)

# Check independance of GLS model
res.gls = sst.gls$residuals
png(gsub(" ", "", paste("figs/SST_PACF.png")))
plot(res.gls, main="Monthly SST Anomaly GLS Residuals") #
dev.off()
dwtest(sst.gls) #p-value = 

# Check normality and constant variance of GLS model
shapiro.test(res.gls) #p-value = 
bartlett.test(res.gls,c(rep(1,length(res.gls)/2),
                    rep(2,length(res.gls)/2))) #p-value =

# Find the trend
png(gsub(" ", "", paste("figs/SST_GLS_trend.png")))
plot(time,sst,xlab="Time (Decimal Date)",ylab="Sea Surface Temperature Anomoly (C)", type="l",
    main="Monthly SST Anomaly GLS-Determined Trend")
lines(time, sst.gls$fitted, col="red")
dev.off()
summary(sst.gls) #estimate of the trend -0.2827
library(trend)
mk.test(res.gls) #p-value = 
sens.slope(res.gls) #


################## Fit QR model ##################
setwd("/home/kareande/esci167-project")
library("quantreg") #QR modeling
library("dplyr") #tidymodels language
#install.packages()

# Fit quantile regression model for 90th percentile (MHW event)
sst_q90 = rq(sst ~ date + lat + lon +
             nit + sil + oxy + pho + npp + chl +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.9)
summary(sst_q90)
res_q90 = sst_q90$residuals
dwtest(sst_q90) #p-value < 2.2e-16, not independent
ci_sstq9 = ci(sst_q90, type="parameter") # confidence interval

# Fit quantile regression model for median (SST)
sst_q50 = rq(sst ~ date + lat + lon +
             nit + sil + oxy + pho + npp + chl +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.5)
summary(sst_q50)
res_q50 = sst_q50$residuals
dwtest(sst_q50) #p-value < 2.2e-16, not independent
ci_sstq5 = ci(sst_q50, type="parameter") # confidence interval

# Fit quantile regression model for 10th percentile (SST)
sst_q10 = rq(sst ~ date + lat + lon +
             nit + sil + oxy + pho + npp + chl +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.1)
summary(sst_q10)
res_q10 = sst_q10$residuals
dwtest(sst_q10) #p-value < 2.2e-16, not independent


#plot SST QR residuals
png(gsub(" ", "", paste("figs/SST_QR_res.png")))
par(mfrow=c(2,1))
plot(res_q90, main="Monthly SST Anomaly QR Residuals
    90th Percentile") #U-shaped pattern
plot(res_q50, main="Monthly SST Anomaly QR Residuals
    Median Quantile") #U-shaped pattern
dev.off()

#Plot SST trends
png(gsub(" ", "", paste("figs/SST_QR_trends.png")))
plot(time,sst,xlab="Time (Decimal Date)",ylab="Sea Surface Temperature Anomoly (C)", type="l",
    main="Monthly SST Anomaly QR Trends")
lines(time,49.07826 + -0.02378*time,col="red",lwd=2) #90th
lines(time,43.47694 + -0.02122*time,col="gold",lwd=2) #Median
lines(time,31.66388 + -0.01565*time,col="blue",lwd=2) #10th
labels=c('90th Percentile',"Median Quantile", "10th Percentile")
colors=c("red","gold", "blue")
legend("bottomleft", inset=.05, labels,lty=c(1,1,1),lwd=c(2,2,1),col=colors)
dev.off()

# Fit quantile regression model for 10th percentile (LChl event)
chl_q10 = rq(chl ~ date + lat + lon +
             nit + sil + oxy + pho + npp + sst +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.1)
summary(chl_q10)
res_q10c = chl_q10$residuals
dwtest(chl_q10) #p-value < 2.2e-16, not independent

# Fit quantile regression model for median (Chl)
chl_q50 = rq(chl ~ date + lat + lon +
             nit + sil + oxy + pho + npp + sst +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.5)
summary(chl_q50)
res_q50c = chl_q50$residuals
dwtest(chl_q50) #p-value < 2.2e-16, not independent

# Fit quantile regression model for 90th percentile (Chl)
chl_q90 = rq(chl ~ date + lat + lon +
             nit + sil + oxy + pho + npp + sst +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.9)
summary(chl_q90)
res_q90c = chl_q90$residuals
dwtest(chl_q90) #p-value < 2.2e-16, not independent

#plot Chl QR residuals
png(gsub(" ", "", paste("figs/Chl_QR_res.png")))
par(mfrow=c(2,1))
plot(res_q10c, main="Monthly Chl Concentration QR Residuals
    10th Percentile") #
plot(res_q50, main="Monthly Chl Concentration QR Residuals
    Median Quantile") #
dev.off()

# Plot Chl trends
png(gsub(" ", "", paste("figs/Chl_QR_trends.png")))
plot(time,chl,type="l")
lines(time,-2.86775 + 0.00404*time,col="green",lwd=2) #10th
lines(time,-21.69045 + 0.01516*time,col="gold",lwd=2) #Median
lines(time,23.77714 + 0.00038*time,col="blue",lwd=2) #90th
labels=c('90th Percentile',"Median Quantile", "10th Percentile")
colors=c("green","gold", "blue")
legend("topleft", inset=.05, labels,lty=c(1,1,1),lwd=c(2,2,1),col=colors)
dev.off()


##################****** 2012-2019 ******##################
setwd("/home/kareande/lchl-mhw-events")
library("quantreg") #QR modeling
library("dplyr") #tidymodels
library("tidyr")
library("lubridate") #date and time manipulation
library("lmtest") #DW test for independence
#install.packages()

################## Data Processing ##################
# Get data from senior thesis project (OISST, NCEP, CMEMS)
cmp_df = read.csv("cmpndData/cmpnd_tst_0lag.csv",header=T) #get df
cmp_df <- subset(cmp_df, yr>2011) 
cmp_df <- cmp_df[ , ! names(cmp_df) %in% c("ds.y","location.y")] #remove duplicate columns
names(cmp_df)[names(cmp_df) == "ds.x"] <- "ds" #rename column
names(cmp_df)[names(cmp_df) == "location.x"] <- "location" #rename column
cmp_df$date <- as.Date(cmp_df$ds, origin = "1992-01-01") #convert days since to date
cmp_df <- cmp_df %>% dplyr::mutate(year = lubridate::year(cmp_df$date), #seperate date into year,
                         month = lubridate::month(cmp_df$date),         #month, and
                         day = lubridate::day(cmp_df$date))             #day columns
cmp_df <- cmp_df[ , ! names(cmp_df) %in% c("yr","lon","lat","ds","cmpCat","date")] #remove columns
cmp_df <- cmp_df %>% relocate(c(day, month, year), .before = doy) #move date vars to front
head(cmp_df)

# Calculate monthly means
cmp_mo = aggregate(cmp_df[,6:17],
                   by=list(cmp_df$month, cmp_df$year, cmp_df$location), #aggregate by month, year, and location
                   FUN=mean, data=cmp_df, na.rm=TRUE)
cmp_mo <- cmp_mo[order(cmp_mo$Group.2),] #reorder dataset by year
cmp_mo <- cmp_mo[order(cmp_mo$Group.1),] #reorder dataset by month
names(cmp_mo)[names(cmp_mo) == "Group.1"] <- "mo" #rename column
names(cmp_mo)[names(cmp_mo) == "Group.2"] <- "yr" #rename column
names(cmp_mo)[names(cmp_mo) == "Group.3"] <- "location" #rename column
cmp_mo$date <- gsub(" ", "", paste(cmp_mo$yr, "-", cmp_mo$mo, "-01")) #merge month and year
cmp_mo$date <- ymd(cmp_mo$date) #get date column
cmp_mo$date <- decimal_date(cmp_mo$date) #convert date to decimal date
cmp_mo <- cmp_mo %>% relocate(date, .before = mo) #move date var to front
cmp_mo <- cmp_mo[order(cmp_mo$date),] #reorder dataset by date
cmp_mo <- cmp_mo %>% separate(location,c("lat","lon"),sep="-") #un-merge location column into latitude and longitude
cmp_mo$lat <- as.integer(cmp_mo$lat) #convert from char to int
cmp_mo$lon <- as.integer(cmp_mo$lon) #convert from char to int
head(cmp_mo)

# Visualize some data
setwd("/home/kareande/esci167-project")
sst = cmp_mo$sst
chl = cmp_mo$chl
time = cmp_mo$date
png(gsub(" ", "", paste("figs/recentSSTvsChl.png")))
plot(sst,chl,xlab="Sea Surface Temperature Anomaly (deg C)",ylab="Chlorophyl Concentration (chl-a/m^3)",
    main="Monthly Sea Surface Temperature Anomaly vs
    Chlorophyl Concentration")
dev.off()
png(gsub(" ", "", paste("figs/recentChlvsTime.png")))
plot(time,chl,xlab="Time (Decimal Date)",ylab="Chlorophyl Concentration (chl-a/m^3)", type="l",
    main="Monthly Chlorophyl Concentration")
dev.off()
png(gsub(" ", "", paste("figs/recentSSTvsTime.png")))
plot(time,sst,xlab="Time (Decimal Date)",ylab="Sea Surface Temperature Anomaly (deg C)", type="l",
    main="Monthly Sea Surface Temperature Anomaly")
dev.off()

################## Fit QR model ##################
setwd("/home/kareande/esci167-project")
library("quantreg") #QR modeling
library("dplyr") #tidymodels language
#install.packages()

# Fit quantile regression model for 90th percentile (MHW event)
sst_q90 = rq(sst ~ date + lat + lon +
             nit + sil + oxy + pho + npp + chl +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.9)
summary(sst_q90) #-142.81806 0.07184
res_q90 = sst_q90$residuals
dwtest(sst_q90) #p-value < 2.2e-16, not independent

# Fit quantile regression model for median (SST)
sst_q50 = rq(sst ~ date + lat + lon +
             nit + sil + oxy + pho + npp + chl +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.5)
summary(sst_q50) #-218.20527 0.10924
res_q50 = sst_q50$residuals
dwtest(sst_q50) #p-value < 2.2e-16, not independent

# Fit quantile regression model for 10th percentile (SST)
sst_q10 = rq(sst ~ date + lat + lon +
             nit + sil + oxy + pho + npp + chl +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.1)
summary(sst_q10) #-266.77249 0.13292
res_q10 = sst_q10$residuals
dwtest(sst_q10) #p-value < 2.2e-16, not independent


#plot SST QR residuals
png(gsub(" ", "", paste("figs/recentSST_QR_res.png")))
par(mfrow=c(2,1))
plot(res_q90, main="Monthly SST Anomaly QR Residuals
    90th Percentile") #U-shaped pattern
plot(res_q50, main="Monthly SST Anomaly QR Residuals
    Median Quantile") #U-shaped pattern
dev.off()

#Plot SST trends
png(gsub(" ", "", paste("figs/recentSST_QR_trends.png")))
plot(time,sst,xlab="Time (Decimal Date)",ylab="Sea Surface Temperature Anomoly (C)", type="l",
    main="Monthly SST Anomaly QR Trends")
lines(time,-142.81806 + 0.07184*time,col="red",lwd=2) #90th
lines(time,-218.20527 + 0.10924*time,col="gold",lwd=2) #Median
lines(time,-266.77249 + 0.13292*time,col="blue",lwd=2) #10th
labels=c('90th Percentile',"Median Quantile", "10th Percentile")
colors=c("red","gold", "blue")
legend("bottomleft", inset=.05, labels,lty=c(1,1,1),lwd=c(2,2,1),col=colors)
dev.off()

# Fit quantile regression model for 10th percentile (LChl event)
chl_q10 = rq(chl ~ date + lat + lon +
             nit + sil + oxy + pho + npp + sst +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.1)
summary(chl_q10)
res_q10c = chl_q10$residuals
dwtest(chl_q10) #p-value < 2.2e-16, not independent

# Fit quantile regression model for median (Chl)
chl_q50 = rq(chl ~ date + lat + lon +
             nit + sil + oxy + pho + npp + sst +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.5)
summary(chl_q50)
res_q50c = chl_q50$residuals
dwtest(chl_q50) #p-value < 2.2e-16, not independent

# Fit quantile regression model for 90th percentile (Chl)
chl_q90 = rq(chl ~ date + lat + lon +
             nit + sil + oxy + pho + npp + sst +
             qnet + slp + sat + wndsp + sstRoC, data=cmp_mo, tau=0.9)
summary(chl_q90)
res_q90c = chl_q90$residuals
dwtest(chl_q90) #p-value < 2.2e-16, not independent

#plot Chl QR residuals
png(gsub(" ", "", paste("figs/recentChl_QR_res.png")))
par(mfrow=c(2,1))
plot(res_q10c, main="Monthly Chl Concentration QR Residuals
    10th Percentile") #
plot(res_q50, main="Monthly Chl Concentration QR Residuals
    Median Quantile") #
dev.off()

# Plot Chl trends
png(gsub(" ", "", paste("figs/recentChl_QR_trends.png")))
plot(time,chl,type="l")
lines(time,2.31648 + 0.00177*time,col="green",lwd=2) #10th
lines(time,125.86873 + -0.05613*time,col="gold",lwd=2) #Median
lines(time,86.85087 + -0.02826*time,col="blue",lwd=2) #90th
labels=c('90th Percentile',"Median Quantile", "10th Percentile")
colors=c("green","gold", "blue")
legend("topleft", inset=.05, labels,lty=c(1,1,1),lwd=c(2,2,1),col=colors)
dev.off()
