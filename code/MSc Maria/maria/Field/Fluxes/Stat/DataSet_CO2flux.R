## Info
# Smoothed data is the data frame with removed missing values for CO2 and DO. It contains also smoothed values for CO2 and DO daily fluxes.
# Flux_daily is the daily average of the raw data



library(car)
library(vegan)
library(nlme)
library(plyr)
library(ggplot2)


setwd("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/")
CO2_flux_data <- read.csv("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/CO2_flux_data.csv")


# Time conversion
# Check https://stat.ethz.ch/R-manual/R-devel/library/base/html/as.POSIXlt.html
CO2_flux_data$Average.Sampling.time <- as.POSIXct(CO2_flux_data$Average.Sampling.time,format = "%m/%d/%Y  %H:%M")
CO2_flux_data$Date <- as.Date(CO2_flux_data$Date, "%m/%d/%Y")

as.POSIXlt("2016-11-10") - as.POSIXlt("2016-05-18") 


head(CO2_flux_data)
tail(CO2_flux_data)
names(CO2_flux_data)


# check different variable by plotting them 
plot(CO2_flux_data[select.list(names(CO2_flux_data), multiple=TRUE, graphics=TRUE)])


# Make an daily average file
Flux_daily <- aggregate(. ~ Date, CO2_flux_data, mean)

plot(CO2_flux_data$Average.Sampling.time, CO2_flux_data$CHL.average, type = "l", col = "darkgreen", ylim = c(-1000, 800))
lines(CO2_flux_data$Average.Sampling.time, CO2_flux_data$Flux.umol.m2.h.1)

plot(CO2_flux_data$CHL.average, CO2_flux_data$Flux.umol.m2.h.1)
plot(CO2_flux_data$Wind.speed..10m..0.1m.s, CO2_flux_data$Flux.umol.m2.h.1)

summary(CO2_flux_data$Flux.umol.m2.h.1, by(CO2_flux_data$Date, FUN = mean, simplify = TRUE))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -823.80 -233.30 -140.20 -163.90  -64.12   59.38     627 





par(mfrow = c(2,1))
plot(CO2_flux_data$Average.Sampling.time, CO2_flux_data$CO2.water.umol.L, type = "l")
plot(CO2_flux_data$Average.Sampling.time, CO2_flux_data$DO.umol.L, type = "l")
par(mfrow = c(1,1))

plot(CO2_flux_data$CO2.water.umol.L, CO2_flux_data$DO.umol.L)
plot(CO2_flux_data[,c(7,11,25)])

plot(Flux_daily$Date, Flux_daily$CHL.average, type = "l", ylim = c(-1000, 800))
lines(Flux_daily$Date, Flux_daily$Flux.umol.m2.h.1)

summary(Flux_daily$Flux.umol.m2.h.1, by(Date))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -540.40 -272.50 -224.20 -203.70  -77.56    3.12 



#extract from the raw data the daily peaks.
Gas_min <- aggregate(. ~ Date, CO2_flux_data[,c(3,7,11)], min)
Gas_max <- aggregate(. ~ Date, CO2_flux_data[,c(3,7,11)], max)
Gas_median <- aggregate(. ~ Date, CO2_flux_data[,c(3,7,11)], median)

Gas_change <- data.frame(Gas_max, Gas_min[,c(2:3)], Gas_median[,c(2:3)])
colnames(Gas_change) [c(2:7)] <- c("CO2_max", "DO_max", "CO2_min", "DO_min", "CO2_median", "DO_median")

Gas_change$CO2_diff <- Gas_change$CO2_max - Gas_change$CO2_min
Gas_change$DO_diff <- Gas_change$DO_max - Gas_change$DO_min
summary(Gas_change)

par(mfrow = c(2,1))
plot(Gas_change$Date, Gas_change$DO_diff, ylab = "Daily absolute change in DO (DO max- DO min)", xlab = "Month (2016)")
plot(Gas_change$Date, Gas_change$CO2_diff, ylab = "Daily absolute change in C02 (C02 max- CO2 min)", xlab = "Month (2016)")
par(mfrow = c(1,1))

par(mfrow = c(2,1))
plot(Gas_change$DO_diff, Gas_change$CO2_diff, xlab = "Daily absolute change in DO (DO max- DO min)", ylab = "Daily absolute change in C02 (C02 max- CO2 min)")
plot(log(Gas_change$DO_diff), log(Gas_change$CO2_diff), xlab = "Daily absolute change in DO (DO max- DO min) log10", ylab = "Daily absolute change in C02 (C02 max- CO2 min) log 10")
par(mfrow = c(1,1))

cor(Gas_change$DO_diff, Gas_change$CO2_diff)






##Check the smooth function
lo <- loess(CO2_flux_data$CO2.water.umol.L ~ CO2_flux_data$point, span = 0.6)
plot(CO2_flux_data$point, CO2_flux_data$CO2.water.umol.L, type= 'l')
lines(predict(lo), col='red', lwd=2)
# Smoothing the whole dataset doesn't work with any method or span (hardness of the smooth)
# note that the smoothing doesn't match the end of the plot because of the missing values. 



###############################################
##Check the smooting function on daily basis 
smooth_CO2 <- data.frame(CO2_flux_data[c(30:53),])

par(mar = c(5,5,2,5))
CO2_smoothing_perDate_test <- ggplot(smooth_CO2, aes(Average.Sampling.time, CO2.water.umol.L)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit=fit<<-..y..), n=24) +
  theme(
    #panel.border = element_rect(fill = NA, colour = "black", size = 1), 
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0,5),
        axis.line.y = element_line(color="black", size = 0,5)) +
  labs(
    y = "Dissolved CO2 umol L-1",
    x = "Time in hours"
  )
CO2_smoothing_perDate_test
fit

ggplot_build(CO2_smoothing_perDate_test)$data[[2]]$y

#over the day the Loess smooting fit good
summary(lm(fit~CO2_flux_data$CO2.water.umol.L[c(30:53)]))








## To do 
# In order to match the points:
# Make new dataframe with no missing values
# And start at time 0.00 missing the first and the last date, so that the fit can be added to the dataframe.
# Do this for both CO2 and DO and create a graph for the data points and fitted values. 
# extract the min and max values and repeat as with the raw data-set


# There are missing values in the data set, so the dataframe was split up in several dataframes. These include high and low values for CO2.
CO2_smooth_data <- CO2_flux_data[complete.cases(CO2_flux_data[,7]),]
count_date <- count(CO2_smooth_data, 'Date')
count_date <- count_date[count_date$freq == 24,]
count_date <- count_date$Date

# from 19-05 till 14-06 
# from 24-06 till 28-06 
# from 30-06 till 23-08
# from 26-08 till 15-10
# from 27-10 till 31-10
# from 05-11 till 09-11



# create dataframe without missing daily values

CO2_smooth_data <- with(CO2_smooth_data, CO2_smooth_data[
                (Date >= "2016-05-19" & Date <= "2016-06-14") | 
                (Date >= "2016-06-24" & Date <= "2016-06-28") | 
                (Date >= "2016-06-30" & Date <= "2016-08-23") | 
                (Date >= "2016-08-26" & Date <= "2016-10-15") | 
                (Date >= "2016-10-27" & Date <= "2016-10-31") | 
                (Date >= "2016-11-05" & Date <= "2016-11-09") ,])


plot( CO2_smooth_data$Date, CO2_smooth_data$CO2.water.umol.L)

#apply the smoothing function on daily basis over the whole dataset
CO2_smoothing_perDate <- ggplot(CO2_smooth_data, aes(Average.Sampling.time, CO2.water.umol.L)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit=fit_CO2<<-..y..), n=24) +
  facet_wrap(~Date)
CO2_smoothing_perDate

ggplot_build(CO2_smoothing_perDate)$data[[2]]$y


#the fit looks good
plot(CO2_smooth_data$CO2.water.umol.L, type = 'l')
lines(fit_CO2, col='red')
summary(lm(CO2_smooth_data$CO2.water.umol.L ~ fit_CO2))

CO2_smooth_data$CO2_fit <- fit_CO2



#############################################################################################
# To compare CO2, HCO3, and DIC and DO, a dataframe was made with no missing values for both

names(CO2_flux_data)
Smoothed_data <- CO2_flux_data[complete.cases(CO2_flux_data[,c(7,11)]),]
count_date <- count(Smoothed_data, 'Date')
count_date <- count_date[count_date$freq == 24,]
count_date <- count_date$Date
count_date


# from 19-05 till 22-05 
# from 24-06 till 28-06 
# from 30-06 till 23-08
# from 26-08 till 15-10
# from 27-10 till 31-10
# from 05-11 till 09-11



# create seperate dataframes.

Smoothed_data <- with(Smoothed_data, Smoothed_data[
  (Date >= "2016-05-19" & Date <= "2016-05-22") | 
    (Date >= "2016-05-26" & Date <= "2016-05-27") | 
    (Date >= "2016-06-02" & Date <= "2016-06-09") | 
    (Date >= "2016-06-24" & Date <= "2016-06-28") | 
    (Date >= "2016-06-30" & Date <= "2016-07-02") |
    (Date >= "2016-07-07" & Date <= "2016-07-17") | 
    (Date >= "2016-07-21" & Date <= "2016-08-01") | 
    (Date >= "2016-08-04" & Date <= "2016-08-09") |
    (Date >= "2016-08-11" & Date <= "2016-08-22") |
    (Date >= "2016-09-01" & Date <= "2016-10-15") |
    (Date >= "2016-10-17" & Date <= "2016-10-18") |
    (Date >= "2016-11-05" & Date <= "2016-11-09")  ,])

# add DIC to the file
Smoothed_data$CO3umol <- Smoothed_data$CO3*1000000
Smoothed_data$HCO3umol <- Smoothed_data$HCO3*1000000
Smoothed_data$DICumol <- Smoothed_data$CO2.water.umol.L + Smoothed_data$HCO3umol + Smoothed_data$CO3umol

par(mfrow = c(4,1))
plot( Smoothed_data$Average.Sampling.time, Smoothed_data$CO2.water.umol.L)
plot( Smoothed_data$Average.Sampling.time, Smoothed_data$DO.umol.L)
plot( Smoothed_data$Average.Sampling.time, Smoothed_data$HCO3umol)
plot( Smoothed_data$Average.Sampling.time, Smoothed_data$DICumol)
par(mfrow = c(1,1))



###########################################################################
## check the fluxes for CO2 and HCO3

summary(Smoothed_data)
# There are no missing values anymore in HCO3


par(mfrow = c(2,1))
plot(Smoothed_data$Average.Sampling.time, Smoothed_data$HCO3umol, type = 'l') 
plot(Smoothed_data$Average.Sampling.time, Smoothed_data$CO2.water.umol.L, type = 'l') 
par(mfrow = c(1,1))

plot(Smoothed_data$HCO3umol, Smoothed_data$CO2.water.umol.L)
plot(Smoothed_data$CO2.water.umol.L, Smoothed_data$DO.umol.L)
plot(Smoothed_data$HCO3umol, Smoothed_data$DO.umol.L)

par(mfrow = c(2,1))
plot(Smoothed_data$Average.Sampling.time[c(1250:1350)], Smoothed_data$HCO3umol[c(1250:1350)], type = 'l') 
plot(Smoothed_data$Average.Sampling.time[c(1250:1350)], Smoothed_data$DO.umol.L[c(1250:1350)], type = 'l')
par(mfrow = c(1,1))

par(mfrow = c(2,1))
plot(Smoothed_data$Average.Sampling.time[c(1250:1350)], Smoothed_data$HCO3umol[c(1250:1350)], type = 'l') 
plot(Smoothed_data$Average.Sampling.time[c(1250:1350)], Smoothed_data$CO2.water.umol.L[c(1250:1350)], type = 'l')
par(mfrow = c(1,1))

par(mfrow = c(2,1))
plot(Smoothed_data$Average.Sampling.time[c(1:100)], Smoothed_data$HCO3umol[c(1:100)], type = 'l') 
plot(Smoothed_data$Average.Sampling.time[c(1:100)], Smoothed_data$DO.umol.L[c(1:100)], type = 'l')
par(mfrow = c(1,1))


par(mfrow = c(2,1))
plot(Flux_daily$Date, Flux_daily$HCO3, type = 'l') 
plot(Flux_daily$Date, Flux_daily$DO.umol.L, type = 'l')
par(mfrow = c(1,1))

summary(lm(Flux_daily$HCO3 ~ Flux_daily$DO.umol.L))





#########################################################################################
#apply the smoothing function on daily basis over the whole dataset 

#CO2
Smoothing_CO2_perDate <- ggplot(Smoothed_data, aes(Average.Sampling.time, CO2.water.umol.L)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit = Smooth_fit_CO2 <<-..y..), n=24) +
  facet_wrap(~Date)
Smoothing_CO2_perDate

ggplot_build(Smoothing_CO2_perDate)$data[[2]]$y

# DO umol/L #apply the smoothing function on daily basis over the whole dataset 
Smoothing_DO_perDate <- ggplot(Smoothed_data, aes(Average.Sampling.time, DO.umol.L)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit = Smooth_fit_DO <<-..y..), n=24) +
  facet_wrap(~Date)
Smoothing_DO_perDate

ggplot_build(Smoothing_DO_perDate)$data[[2]]$y


# HCO3 umol/L #apply the smoothing function on daily basis over the whole dataset 
Smoothing_HCO3_perDate <- ggplot(Smoothed_data, aes(Average.Sampling.time, HCO3*1000000)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit = Smooth_fit_HCO3 <<-..y..), n=24) +
  facet_wrap(~Date)
Smoothing_HCO3_perDate

ggplot_build(Smoothing_HCO3_perDate)$data[[2]]$y


# CO3 umol/L #apply the smoothing function on daily basis over the whole dataset 
Smoothing_CO3_perDate <- ggplot(Smoothed_data, aes(Average.Sampling.time, CO3*1000000)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit = Smooth_fit_CO3 <<-..y..), n=24) +
  facet_wrap(~Date)
Smoothing_CO3_perDate

ggplot_build(Smoothing_CO3_perDate)$data[[2]]$y


###########################################################
# Check the LOESS curve

#the fit looks good
plot(Smoothed_data$CO2.water.umol.L, type = 'l')
lines(Smooth_fit_CO2, col='red')
summary(lm(Smoothed_data$CO2.water.umol.L ~ Smooth_fit_CO2))


#the fit looks good
plot(Smoothed_data$DO.umol.L, type = 'l')
lines(Smooth_fit_DO, col='red')
summary(lm(Smoothed_data$DO.umol.L ~ Smooth_fit_DO))


#the fit looks good
plot(Smoothed_data$HCO3*1000000, type = 'l')
lines(Smooth_fit_HCO3, col='red')
summary(lm(Smoothed_data$HCO3*1000000 ~ Smooth_fit_HCO3))

#the fit looks good
plot(Smoothed_data$CO3*1000000, type = 'l')
lines(Smooth_fit_CO3, col='red')
summary(lm(Smoothed_data$CO3*1000000 ~ Smooth_fit_CO3))


Smoothed_data$CO2_fit <- Smooth_fit_CO2
Smoothed_data$DO_fit <- Smooth_fit_DO
Smoothed_data$HCO3_fit <- Smooth_fit_HCO3
Smoothed_data$CO3_fit <- Smooth_fit_CO3
Smoothed_data$DIC_fit <- Smooth_fit_CO2 + Smooth_fit_HCO3 + Smooth_fit_CO3





#####################################################################
## Daily changes in DIC and DO
## Build a new dataframe for the daily DO and CO2 range

Gas_min_smooth <- aggregate(. ~ Date, Smoothed_data[,c("Date", "DO_fit", "CO2_fit", "HCO3_fit", "DIC_fit")], min)
Gas_max_smooth <- aggregate(. ~ Date, Smoothed_data[,c("Date", "DO_fit", "CO2_fit", "HCO3_fit", "DIC_fit")], max)
Gas_median_smooth <- aggregate(. ~ Date, Smoothed_data[,c("Date", "DO_fit", "CO2_fit", "HCO3_fit", "DIC_fit")], median)


Gas_change_smooth <- data.frame(Gas_max_smooth, Gas_min_smooth[,c(2:ncol(Gas_min_smooth))], Gas_median_smooth[,c(2:ncol(Gas_median_smooth))])
colnames(Gas_change_smooth)[c(2:ncol(Gas_change_smooth))] <- c("DO_max", "CO2_max",  "HCO3_max", "DIC_max", "DO_min", "CO2_min",  "HCO3_min", "DIC_min", "DO_median", "CO2_median",  "HCO3_median", "DIC_median")
Gas_change_smooth


Gas_change_smooth$CO2_diff <- Gas_change_smooth$CO2_max - Gas_change_smooth$CO2_min
Gas_change_smooth$DO_diff <- Gas_change_smooth$DO_max - Gas_change_smooth$DO_min
Gas_change_smooth$HCO3_diff <- Gas_change_smooth$HCO3_max - Gas_change_smooth$HCO3_min
Gas_change_smooth$DIC_diff <- Gas_change_smooth$DIC_max - Gas_change_smooth$DIC_min


rm(Gas_min_smooth, Gas_max_smooth, Gas_median_smooth)



par(mfrow = c(4,1))
plot(Gas_change_smooth$Date, Gas_change_smooth$DO_diff)
plot(Gas_change_smooth$Date, Gas_change_smooth$CO2_diff)
plot(Gas_change_smooth$Date, Gas_change_smooth$DIC_diff)
plot(Gas_change_smooth$Date, Gas_change_smooth$HCO3_diff)
par(mfrow = c(2,2))

plot(Gas_change_smooth$DO_diff, Gas_change_smooth$CO2_diff)
plot(Gas_change_smooth$DO_diff, Gas_change_smooth$HCO3_diff)
plot(Gas_change_smooth$DO_diff, Gas_change_smooth$DIC_diff)
plot(Gas_change_smooth$CO2_diff, Gas_change_smooth$DIC_diff)

plot(lm(log(Gas_change_smooth$DO_diff) ~ log(Gas_change_smooth$CO2_diff)))
plot(lm(Gas_change_smooth$DO_diff ~ Gas_change_smooth$HCO3_diff))
plot(lm(Gas_change_smooth$DO_diff ~ Gas_change_smooth$DIC_diff))
plot(lm(log(Gas_change_smooth$CO2_diff) ~ log(Gas_change_smooth$DIC_diff)))

summary(lm(log(Gas_change_smooth$DO_diff) ~ log(Gas_change_smooth$CO2_diff)))
summary(lm(Gas_change_smooth$DO_diff ~ Gas_change_smooth$HCO3_diff))
summary(lm(Gas_change_smooth$DO_diff ~ Gas_change_smooth$DIC_diff))
summary(lm(log(Gas_change_smooth$CO2_diff) ~ log(Gas_change_smooth$DIC_diff)))



#########################################################################################################
## Check the relative change of the maximum fluctuation, compared to the average value of the daily median
Gas_change_smooth$CO2_diff_rel <- Gas_change_smooth$CO2_diff/Gas_change_smooth$CO2_median
Gas_change_smooth$DO_diff_rel <- Gas_change_smooth$DO_diff/Gas_change_smooth$DO_median
Gas_change_smooth$HCO3_diff_rel <- Gas_change_smooth$HCO3_diff/Gas_change_smooth$HCO3_median
Gas_change_smooth$DIC_diff_rel <- Gas_change_smooth$DIC_diff/Gas_change_smooth$DIC_median
summary(Gas_change_smooth)


par(mfrow = c(4,1))
plot(Gas_change_smooth$Date, Gas_change_smooth$DO_diff_rel)
plot(Gas_change_smooth$Date, Gas_change_smooth$CO2_diff_rel)
plot(Gas_change_smooth$Date, Gas_change_smooth$DIC_diff_rel)
plot(Gas_change_smooth$Date, Gas_change_smooth$HCO3_diff_rel)
par(mfrow = c(1,1))

## plot function
plot(Gas_change_smooth[select.list(names(Gas_change_smooth), multiple = TRUE, graphics = TRUE)])

par(mfrow = c(2,2))
plot(Gas_change_smooth$DO_diff_rel[-15], Gas_change_smooth$CO2_diff_rel[-15])
plot(Gas_change_smooth$DO_diff_rel, Gas_change_smooth$HCO3_diff_rel)
plot(Gas_change_smooth$DO_diff_rel, Gas_change_smooth$DIC_diff_rel)
plot(Gas_change_smooth$CO2_diff_rel[-15], Gas_change_smooth$DIC_diff_rel[-15])


plot(lm(Gas_change_smooth$DO_diff_rel[-15] ~ Gas_change_smooth$CO2_diff_rel[-15]))
#15 outlier high DO low CO2
plot(lm(Gas_change_smooth$DO_diff_rel ~ Gas_change_smooth$HCO3_diff_rel))
plot(lm(Gas_change_smooth$DO_diff_rel ~ Gas_change_smooth$DIC_diff_rel))
plot(lm(Gas_change_smooth$CO2_diff_rel[-15] ~ Gas_change_smooth$DIC_diff_rel[-15]))
#15 outlier

summary(lm(Gas_change_smooth$DO_diff_rel[-15] ~ Gas_change_smooth$CO2_diff_rel[-15]))
summary(lm(Gas_change_smooth$DO_diff_rel ~ Gas_change_smooth$HCO3_diff_rel))
summary(lm(Gas_change_smooth$DO_diff_rel ~ Gas_change_smooth$DIC_diff_rel))
summary(lm(Gas_change_smooth$CO2_diff_rel[-15] ~ Gas_change_smooth$DIC_diff_rel[-15]))







##############################################################################
#split up dataset, to off-season and high-season
#from the first of july untill october
Gas_change_smooth_Hseason <- Gas_change_smooth[c(16:93),]
Gas_change_smooth_Lseason <- Gas_change_smooth[-c(16:93),]



#When max values are plotted against each other there seems to be be similar linearity (on a different scale)
par(mfrow = c(2,1))
plot(Gas_change_smooth_Lseason$CO2_max, Gas_change_smooth_Lseason$DO_max)
plot(Gas_change_smooth_Hseason$CO2_max, Gas_change_smooth_Hseason$DO_max)

# When the data set is divided, there does not seem to be an obvision relationship between increasing or decreasing differecnces in dissolved gas.
plot(Gas_change_smooth_Lseason$CO2_diff, Gas_change_smooth_Lseason$DO_diff)
plot(Gas_change_smooth_Hseason$CO2_diff, Gas_change_smooth_Hseason$DO_diff)

# When the data set is divided, there does not seem to be an obvision relationship between increasing or decreasing differecnces in dissolved gas.
plot(Gas_change_smooth_Lseason$CO2_diff_rel, Gas_change_smooth_Lseason$DO_diff_rel)
plot(Gas_change_smooth_Hseason$CO2_diff_rel, Gas_change_smooth_Hseason$DO_diff_rel)

par(mfrow = c(1,1))

#check linear regression for devided data sets


lm_max_L <- lm(Gas_change_smooth_Lseason$DO_max ~ Gas_change_smooth_Lseason$CO2_max)
par(mfrow = c(2,2))
plot(lm_max_L)
par(mfrow = c(1,1))
summary(lm_max_L)

#hte high season has bad fit of the model
lm_max_H <- lm(Gas_change_smooth_Hseason$DO_max ~ Gas_change_smooth_Hseason$CO2_max)
par(mfrow = c(2,2))
plot(lm_max_H)
par(mfrow = c(1,1))
summary(lm_max_H)
























##################################################################################
# check the relation of DIC and DO

plot(Smoothed_data$Average.Sampling.time, Smoothed_data$DICumol/10, col = "darkorange", type = 'l')
points(Smoothed_data$Average.Sampling.time, Smoothed_data$DO.umol.L, col="darkblue", type = 'l')


par(mfrow = c(2,2))
plot(Smoothed_data$DO_fit, Smoothed_data$CO2_fit)
plot(Smoothed_data$DO_fit, Smoothed_data$HCO3_fit)
plot(Smoothed_data$DO_fit, Smoothed_data$DIC_fit)
plot(Smoothed_data$CO2_fit, Smoothed_data$DIC_fit)
par(mfrow = c(1,1))

plot(log(Smoothed_data$DO_fit), log(Smoothed_data$CO2_fit))
# There is no linear relation in the whole data set. Furthermore, the months are grouped in the plot.


## Check the slope of CO2 and DO over 24 hour periods

par(mfrow = c(2,2))
plot(Smoothed_data$DO_fit[c(1:24)], Smoothed_data$CO2_fit[1:24])
plot(Smoothed_data$DO_fit[c(50:73)], Smoothed_data$CO2_fit[50:73])
plot(Smoothed_data$DO_fit[c(1000:1024)], Smoothed_data$CO2_fit[1000:1024])
plot(Smoothed_data$DO_fit[c(2000:2024)], Smoothed_data$CO2_fit[2000:2024])
par(mfrow = c(1,1))

#extract the slope of the lm over a 24 hour range
ind_j <- 1:2736
ind_k <- 24:2760
ind_w <- 1:2760
lm_daily_slope <- data.frame()
lm_daily_slope[c(1:2736),] <- data.frame()

summary(lm(log(Smoothed_data$CO2_fit[c(1:24)]) ~ log(Smoothed_data$DO_fit[1:24])))
plot(log(Smoothed_data$CO2_fit[c(1:24)]) , log(Smoothed_data$DO_fit[1:24]))


for (i in seq_along(ind_j)){
  lm_daily_slope$slopeCO2[ind_w[i]] <- lm(Smoothed_data$CO2_fit[ind_j[i]:ind_k[i]] ~ Smoothed_data$DO_fit[ind_j[i]:ind_k[i]])$coefficients[[2]]
}

plot(lm_daily_slope$slopeCO2, ylab = "the change of CO2 relative to DO")

summary(lm_daily_slope$slopeCO2)



#check the relative change (max-min) compared to the median of each day
############## For the measured data
CO2_median <- aggregate(CO2.water.umol.L ~ Date, data = Smoothed_data, median)
CO2_median <- CO2_median[rep(row.names(CO2_median), 24), 1:2]
CO2_median <- CO2_median[order(CO2_median$Date),] 
Smoothed_data$CO2_median <- CO2_median$CO2.water.umol.L

DO_median <- aggregate(DO.umol.L ~ Date, data = Smoothed_data, median)
DO_median <- DO_median[rep(row.names(DO_median), 24), 1:2]
DO_median <- DO_median[order(DO_median$Date),] 
Smoothed_data$DO_median <- DO_median$DO.umol.L

HCO3_median <- aggregate(HCO3 ~ Date, data = Smoothed_data, median)
HCO3_median <- HCO3_median[rep(row.names(HCO3_median), 24), 1:2]
HCO3_median <- HCO3_median[order(HCO3_median$Date),] 
Smoothed_data$HCO3_median <- HCO3_median$HCO3


Smoothed_data$CO2_diff_rel <- (Smoothed_data$CO2.water.umol.L/Smoothed_data$CO2_median)
Smoothed_data$DO_diff_rel <- (Smoothed_data$DO.umol.L/Smoothed_data$DO_median)
Smoothed_data$HCO3_diff_rel <- (Smoothed_data$HCO3/Smoothed_data$HCO3_median)

plot(Smoothed_data$Average.Sampling.time, Smoothed_data$CO2_diff_rel)
plot(Smoothed_data$Average.Sampling.time, Smoothed_data$DO_diff_rel)
plot(Smoothed_data$Average.Sampling.time, Smoothed_data$HCO3_diff_rel)
plot(log(Smoothed_data$CO2_diff_rel), log(Smoothed_data$DO_diff_rel))
plot(log(Smoothed_data$CO2_diff_rel), log(Smoothed_data$HCO3_diff_rel))
plot(log(Smoothed_data$HCO3_diff_rel), log(Smoothed_data$DO_diff_rel))


# check the linear regression
# linearity in the model is lacking, which is still lacking with a log transformation
# There is linearity when the model is plit up different months

Smoothed_data_lm_adapted <- with(Smoothed_data, Smoothed_data[
  (Date >= "2016-05-01" & Date < "2016-11-07"),])
par(mfrow=c(2,2))
plot(lm(Smoothed_data_lm_adapted$CO2_diff_rel ~ Smoothed_data_lm_adapted$DO_diff_rel))
summary(lm(Smoothed_data_lm_adapted$CO2_diff_rel ~ Smoothed_data_lm_adapted$DO_diff_rel))
par(mfrow = c(1,1)) 

plot(Smoothed_data_lm_adapted$CO2_diff_rel, Smoothed_data_lm_adapted$DO_diff_rel)
# The relative changes in CO2 






#check the relative change (max-min) compared to the median of each day
############## For the fitted data 
CO2_median <- aggregate(CO2_fit ~ Date, data = Smoothed_data, median)
CO2_median <- CO2_median[rep(row.names(CO2_median), 24), 1:2]
CO2_median <- CO2_median[order(CO2_median$Date),] 
Smoothed_data$CO2_median <- CO2_median$CO2_fit

DO_median <- aggregate(DO_fit ~ Date, data = Smoothed_data, median)
DO_median <- DO_median[rep(row.names(DO_median), 24), 1:2]
DO_median <- DO_median[order(DO_median$Date),] 
Smoothed_data$DO_median <- DO_median$DO_fit

HCO3_median <- aggregate(HCO3_fit ~ Date, data = Smoothed_data, median)
HCO3_median <- HCO3_median[rep(row.names(HCO3_median), 24), 1:2]
HCO3_median <- HCO3_median[order(HCO3_median$Date),] 
Smoothed_data$HCO3_median <- HCO3_median$HCO3_fit

DIC_median <- aggregate(DIC_fit ~ Date, data = Smoothed_data, median)
DIC_median <- DIC_median[rep(row.names(DIC_median), 24), 1:2]
DIC_median <- DIC_median[order(DIC_median$Date),] 
Smoothed_data$DIC_median <- DIC_median$DIC_fit



# Check the median values
plot(Smoothed_data$Average.Sampling.time, Smoothed_data$CO2_fit)
points(Smoothed_data$Average.Sampling.time, Smoothed_data$CO2_median, col = "red")
lines(Smoothed_data$Average.Sampling.time, Smoothed_data$Wind.speed.m.s, col = "green")
Smoothed_data_median <- with(Smoothed_data, Smoothed_data[
  (Date >= "2016-06-24" & Date < "2016-06-25"),])
plot(Smoothed_data_median$Average.Sampling.time, Smoothed_data_median$CO2_fit)
points(Smoothed_data_median$Average.Sampling.time, Smoothed_data_median$CO2_median, col = "red")
# !!!!!!24-06-2016 has a low median and forms an outlier. Don't know yet why!!!!!!!!!!!!!!!!!
median(Smoothed_data_median$CO2_fit)
mean(max(Smoothed_data_median$CO2_fit):min(Smoothed_data_median$CO2_fit))



Smoothed_data$CO2_diff_rel <- (Smoothed_data$CO2_fit/Smoothed_data$CO2_median)-1
Smoothed_data$DO_diff_rel <- (Smoothed_data$DO_fit/Smoothed_data$DO_median)-1
Smoothed_data$HCO3_diff_rel <- (Smoothed_data$HCO3_fit/Smoothed_data$HCO3_median)-1
Smoothed_data$DIC_diff_rel <- (Smoothed_data$DIC_fit/Smoothed_data$DIC_median)-1





par(mfrow = c(4,1))
plot(Smoothed_data$Average.Sampling.time, Smoothed_data$CO2_diff_rel)
plot(Smoothed_data$Average.Sampling.time, Smoothed_data$DO_diff_rel)
plot(Smoothed_data$Average.Sampling.time, Smoothed_data$HCO3_diff_rel)
plot(Smoothed_data$Average.Sampling.time, Smoothed_data$DIC_diff_rel)
par(mfrow = c(2,2))
plot(Smoothed_data$CO2_diff_rel, Smoothed_data$DO_diff_rel)
plot(Smoothed_data$CO2_diff_rel, Smoothed_data$HCO3_diff_rel)
plot(Smoothed_data$HCO3_diff_rel, Smoothed_data$DO_diff_rel)
plot(Smoothed_data$DIC_diff_rel, Smoothed_data$DO_diff_rel)
par(mfrow = c(1,1))


plot(Smoothed_data_lm_adapted$CO2_diff_rel, Smoothed_data_lm_adapted$DO_diff_rel, xlim = c(0.5,1.3), ylim = c(0.5,1.3))


# check the linear regression
Smoothed_data_lm_adapted <- with(Smoothed_data, Smoothed_data[
  (CO2_diff_rel < 1 & DO_diff_rel < 1),]) 
# (log(CO2_diff_rel) < 0.5 & log(DO_diff_rel) < 0.5),]) 
plot(Smoothed_data_lm_adapted$CO2_diff_rel ~ Smoothed_data_lm_adapted$DO_diff_rel)
par(mfrow=c(2,2))
plot(lm(Smoothed_data_lm_adapted$CO2_diff_rel ~ Smoothed_data_lm_adapted$DO_diff_rel))
summary(lm(Smoothed_data_lm_adapted$DIC_diff_rel ~ Smoothed_data_lm_adapted$DO_diff_rel))
par(mfrow = c(1,1)) 
lm <- lm(Smoothed_data_lm_adapted$CO2_diff_rel ~ Smoothed_data_lm_adapted$DO_diff_rel)
hist(lm$residuals)
# The distribution of residual is slightly sqewed to the left. 
# By eliminating extreme vlaues the data becomes normally distrubuted, the slope and the intercept don't change more than 10%
# Because the data has so many points the slope is assumed to be correct. Therefore, only the extreme values of >= were removed.  





########################################
# Check the linear regression per month


months(Smoothed_data$Date)
plot(subset(Smoothed_data$CO2_diff_rel, months(Smoothed_data$Date) == "July") , subset(Smoothed_data$DO_diff_rel, months(Smoothed_data$Date) == "July"))
summary(lm(subset(Smoothed_data$CO2_diff_rel, months(Smoothed_data$Date) == "July") ~ subset(Smoothed_data$DO_diff_rel, months(Smoothed_data$Date) == "July")))
par(mfrow = c(2,2))
plot(lm(log(subset(Smoothed_data$CO2_diff_rel, months(Smoothed_data$Date) == "July")) ~ log(subset(Smoothed_data$DO_diff_rel, months(Smoothed_data$Date) == "July"))))
par(mfrow = c(1,1))


##########################
## Check the slope by 24 hours

ind_j <- 1:2736
ind_k <- 24:2760
ind_w <- 1:2760
lm_week <- data.frame()
lm_week[c(1:2736),] <- data.frame()

for (i in seq_along(ind_j)){
lm_week$slopeCO2[ind_w[i]] <- lm(Smoothed_data$CO2_diff_rel[ind_j[i]:ind_k[i]] ~
                                   Smoothed_data$DO_diff_rel[ind_j[i]:ind_k[i]])$coefficients[[2]]
}

for (i in seq_along(ind_j)){
  lm_week$slopeDIC[ind_w[i]] <- lm(Smoothed_data$DIC_diff_rel[ind_j[i]:ind_k[i]] ~
                                     Smoothed_data$DO_diff_rel[ind_j[i]:ind_k[i]])$coefficients[[2]]
}


par(mfrow = c(2,2))
plot(Smoothed_data$CO2_diff_rel[1:24], Smoothed_data$DO_diff_rel[1:24])
plot(Smoothed_data$CO2_diff_rel[4:28], Smoothed_data$DO_diff_rel[4:28])
plot(Smoothed_data$CO2_diff_rel[50:74], Smoothed_data$DO_diff_rel[50:74])
plot(Smoothed_data$CO2_diff_rel[100:124], Smoothed_data$DO_diff_rel[100:124])
par(mfrow = c(1,1))


par(mfrow = c(2,2))
plot(Smoothed_data$DIC_diff_rel[1:24], Smoothed_data$DO_diff_rel[1:24])
plot(Smoothed_data$DIC_diff_rel[4:28], Smoothed_data$DO_diff_rel[4:28])
plot(Smoothed_data$DIC_diff_rel[50:74], Smoothed_data$DO_diff_rel[50:74])
plot(Smoothed_data$DIC_diff_rel[100:124], Smoothed_data$DO_diff_rel[100:124])
par(mfrow = c(1,1))


par(mfrow = c(2,1))
plot(lm_week$slopeCO2, ylab = "The slope of the lm => relative_daily_change in CO2:DO")
mean(lm_week$slopeCO2)

plot(lm_week$slopeDIC)
mean(lm_week$slopeDIC)
par(mfrow = c(1,1))





########################################
# Check the fluctuations over daily time

plot(log(Smoothed_data$CO2_fit), log(Smoothed_data$DO_fit), type = 'l')
plot(Smoothed_data$CO2_fit, Smoothed_data$DO_fit, type = 'l', xlim = c(0,2))
lm_fit_day <- lm(Smoothed_data$CO2_fit ~ Smoothed_data$DO_fit)
par(mfrow = c(2,2))
plot(lm_fit_day)
par(mfrow = c(1,1))
summary(lm_fit_day)


# CO2 below 1 and above
Smoothed_data_HSeason <- with(Smoothed_data, Smoothed_data[
  (CO2_fit >= 0 & CO2_fit <= 1),])
Smoothed_data_LSeason <- with(Smoothed_data, Smoothed_data[
  (CO2_fit > 1),])

plot(log(Smoothed_data_HSeason$CO2_fit), log(Smoothed_data_HSeason$DO_fit))
plot(Smoothed_data_HSeason$CO2_fit, Smoothed_data_HSeason$DO_fit)
lm_fit_day_Low <- lm(log(Smoothed_data_HSeason$CO2_fit) ~ log(Smoothed_data_HSeason$DO_fit))
par(mfrow = c(2,2))
plot(lm_fit_day_Low)
par(mfrow = c(1,1))
summary(lm_fit_day_Low)

plot(log(Smoothed_data_LSeason$CO2_fit), log(Smoothed_data_LSeason$DO_fit))
plot(Smoothed_data_LSeason$CO2_fit, Smoothed_data_LSeason$DO_fit)
lm_fit_day_Low <- lm(log(Smoothed_data_LSeason$CO2_fit) ~ log(Smoothed_data_LSeason$DO_fit))
par(mfrow = c(2,2))
plot(lm_fit_day_Low)
par(mfrow = c(1,1))
summary(lm_fit_day_Low)



> plot(Gas_change_smooth$Date, Gas_change_smooth$CO2_max, col="red")
> points(Gas_change_smooth$Date, Gas_change_smooth$CO2_min, col = "green")
> lines(Gas_change_smooth$Date, Gas_change_smooth$CO2_median)









######################################
##plot graphs


plot(Smoothed_data$Average.Sampling.time, Smoothed_data$CO2.water.umol.L, col = "darkorange", ylim = c(0, 700))
points(Smoothed_data$Average.Sampling.time, Smoothed_data$DO.umol.L, col= "darkblue")


par(mar = c(5,5,2,5))
with(Smoothed_data, plot(Average.Sampling.time, CO2.water.umol.L, col="darkorange", pch=16, cex=0.8,
             ylab=expression(Dissolved_CO[2] (umol/L))))

par(new = T)
with(Smoothed_data, plot(Average.Sampling.time, DO.umol.L, pch=16, axes=F, xlab=NA, ylab=NA, cex=0.8, col="darkblue"))
axis(side = 4)
mtext(side = 4, line = 3, "Dissolved oxygen (umol/L)")
legend("topleft",
       legend=c(expression(CO[2] (umol/L)), expression(DO (umol/L))),
       lty=c(1,0), pch=c(NA, 16, cex = 1), col=c("darkorange", "darkblue"), cex = 1)
#Make this plot with excell


#install.packages("ggthemes")
library(ggthemes)
library(grid)

# color-blind-friendly palettes, to change gray into black replace "#999999" with -> "#000000"
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


CO2 <- ggplot(Smoothed_data, aes(Average.Sampling.time, CO2.water.umol.L)) + 
  geom_point(colour = "darkorange", size = 1) + #geom_rangeframe() +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line.x = element_line(color="black", size = 0,5),
        axis.line.y = element_line(color="black", size = 0,5)) +
  labs(
    x = "Time in months (2016)",
    y = "Dissolved CO2 in umol/L"
  )
 
DO <- ggplot(Smoothed_data, aes(Average.Sampling.time, DO.umol.L)) + 
  geom_point(colour = "darkblue", size = 1) + 
  theme_bw() +
theme(
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank()
) + scale_y_continuous(position = "right") +
  theme(axis.line.x = element_line(color="black", size = 0,5),
        axis.line.y = element_line(color="black", size = 0,5)) +
  labs(
    x = "Time in months (2016)",
    y = "Dissolved CO2 in umol/L"
  )
par(mar = c(5,5,2,5))
CO2
par(new=T)
DO


par(mar = c(5,5,2,5))  
CO2_DO_cor
par(new = T)
+ scale_y_continuous(position = "right")






CO2_DO_cor <- ggplot(Smoothed_data, aes(DO_diff_rel, DIC_diff_rel, colour = months(Smoothed_data$Date))) +
  geom_hline(yintercept = 0, colour = "grey" , linetype = 2 ) + geom_vline(xintercept = 0, colour = "grey" , linetype = 2 ) +
  geom_point() +  ylim(-1,1) + xlim(-1,1) + scale_colour_manual(values=cbPalette) +
  geom_abline(intercept = 0   , slope = 1) +
  theme(
    #panel.border = element_rect(fill = NA, colour = "black", size = 1), 
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
    theme(axis.line.x = element_line(color="black", size = 0,5),
        axis.line.y = element_line(color="black", size = 0,5)) +
  labs(
    x = "Relative change in dissolved DO umol L-1" ,
    y = "Relative change in dissolved DIC umol L-1",
    color = "Months (2016)"
  )
CO2_DO_cor



# color-blind-friendly palettes, to change gray into black replace "#999999" with -> "#000000"
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

CO2_DO_cor <- ggplot(Smoothed_data, aes(CO2_fit, DO_fit, colour = months(Smoothed_data$Date))) +
  geom_point() + scale_colour_manual(values=cbPalette) + 
  theme(
    #panel.border = element_rect(fill = NA, colour = "black", size = 1), 
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0,5),
        axis.line.y = element_line(color="black", size = 0,5)) +
  labs(
    x = "Dissolved CO2 umol L-1",
    y = "Dissolved DO umol L-1",
    color = "Months (2016)"
  )
CO2_DO_cor

















Smoothed_data_daily <- aggregate(. ~ Date, Smoothed_data, mean)
plot(Smoothed_data_daily$Date, Smoothed_data_daily$DO_diff_rel)

pca_Smoothed_data_daily <- prcomp(Smoothed_data_daily[,c(34,35,36)])
plot(pca_Smoothed_data_daily)
Anova(pca_Smoothed_data_daily)


ggbiplot(pca_Smoothed_data_daily, choices=c(1,2), obs.scale = 1, var.scale = 1,
         ellipse = TRUE, circle = F) +
  #xlim(-2, 2) + ylim(-2, 2) +
  geom_point(aes( color = Smoothed_data_daily$CHL.average), size=3) +
  geom_hline(aes(0), size=.2, col="grey", alpha=0.5, yintercept = 0) + 
  geom_vline(aes(0), size=.2, col="grey", alpha=0.5, xintercept = 0) +
  guides(shape=guide_legend("Season"), colour = guide_legend("Chl a (ug/L)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),  axis.line = element_line(colour = "black"), legend.direction= 'vertical', legend.position = 'right')






CO2_DO_cor <- ggplot(Smoothed_data_daily, aes(log(CO2_diff_rel), log(DO_diff_rel), colour = CHLV.Volts, size = pH)) +
  geom_point() + #scale_colour_manual(values=cbPalette) + 
  theme(
    #panel.border = element_rect(fill = NA, colour = "black", size = 1), 
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0,5),
        axis.line.y = element_line(color="black", size = 0,5)) +
  labs(
    x = "Dissolved CO2 umol L-1",
    y = "Dissolved DO umol L-1",
    color = "Months (2016)"
  )
CO2_DO_cor

