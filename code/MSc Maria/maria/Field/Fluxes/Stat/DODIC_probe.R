## Info
# Smoothed data is the data frame with removed missing values for CO2 and DO. It contains also smoothed values for CO2 and DO daily fluxes.
# Flux_daily is the daily average of the raw data


#If there are unstalled packages use: install.packages("name_package")
library(car)
library(vegan)
library(nlme)
library(plyr)
library(ggplot2)



#Chewbacca
setwd("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/")
CO2_flux_data <- read.csv("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/DIC_alk_data.csv")

#Przemek
setwd("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/")
CO2_flux_data <- read.csv("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/DIC_alk_data.csv", sep = ";")

#UvA
setwd("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/")
CO2_flux_data <- read.csv("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/DIC_alk_data.csv")


## DIC depents on the data set
CO2_flux_data$DICumol <- CO2_flux_data$DIC..Alk..pH.*1000000
CO2_flux_data$DICumol <- CO2_flux_data$DIC..pH.CO2.*1000000
CO2_flux_data$CO2umol <- CO2_flux_data$CO2.water.umol.L*1000000
CO2_flux_data$HCO3umol <- CO2_flux_data$HCO3*1000000
CO2_flux_data$CO3umol <- CO2_flux_data$CO3*1000000

  
# Time conversion
# Check https://stat.ethz.ch/R-manual/R-devel/library/base/html/as.POSIXlt.html
CO2_flux_data$Average.Sampling.time <- as.POSIXct(CO2_flux_data$Average.Sampling.time,format = "%m-%d-%Y  %H:%M")
CO2_flux_data$Date <- as.Date(CO2_flux_data$Date, "%m-%d-%Y")
CO2_flux_data$Time <- format(CO2_flux_data$Average.Sampling.time, "%H:%M")

as.POSIXlt("2016-11-10") - as.POSIXlt("2016-05-18") 

head(CO2_flux_data)
tail(CO2_flux_data)
names(CO2_flux_data)



################################################################################################
#                      check different variable by plotting them                               #  
     plot(CO2_flux_data[select.list(names(CO2_flux_data), multiple=TRUE, graphics=TRUE)])      #
                                                                                               #
################################################################################################


#################################################################################
#Check the data values

plot(CO2_flux_data$Average.Sampling.time, CO2_flux_data$CHL.aver)
plot(CO2_flux_data$Average.Sampling.time[CO2_flux_data$CHL.aver<1000 & CO2_flux_data$CHL.aver > 80], 
     CO2_flux_data$CHL.aver[CO2_flux_data$CHL.aver<1000 & CO2_flux_data$CHL.aver > 80], type = "l", col = "darkgreen")
summary(CO2_flux_data$DO.umol.L)

plot(CO2_flux_data$Average.Sampling.time, CO2_flux_data$pH)
text(CO2_flux_data$Average.Sampling.time, CO2_flux_data$pH, row.names(CO2_flux_data), cex=0.6, pos=4, col="red")
#rm points 6814, 7669, 8493, 9170, 10172, 11008, 12345, 12344, 14024, 14859, 14857, 15693, 16532, 17491, 18207, 18209, 19050

###Remove extreme values due to visitations

#Chlorophyll
CO2_flux_data$CHL.aver[CO2_flux_data$CHL.aver > 1000 | CO2_flux_data$CHL.aver < 80] <- NA

#pH
CO2_flux_data$pH[c(1799, 3478, 5144, 5145, 6814, 7669, 8493, 9170, 10172, 11008, 12345, 12344, 14024, 14859, 14857, 15693, 16532, 17491, 18207, 18209, 19050, 10174, 3800, 3802, 5983, 5981, 8495, 9336)] <- NA


###recheck data

#Chlorophyl
summary(CO2_flux_data$CHL.aver)
plot(CO2_flux_data$Average.Sampling.time, CO2_flux_data$CHL.aver)

#Chlorophyl
plot(CO2_flux_data$Average.Sampling.time, CO2_flux_data$pH)



###############################################################################################
##plot graphs


par(mfrow = c(4,1))
plot(CO2_flux_data$Average.Sampling.time, CO2_flux_data$DO.umol.L, type = "l")
plot(CO2_flux_data$Average.Sampling.time, CO2_flux_data$pH, type = "l")
plot(CO2_flux_data$Average.Sampling.time, CO2_flux_data$CO2umol, type = "l")
plot(CO2_flux_data$Average.Sampling.time, CO2_flux_data$DICumol, type = "l")
par(mfrow = c(1,1))

plot(CO2_flux_data$CO2.water.umol.L, CO2_flux_data$DO.umol.L)
names(CO2_flux_data)
plot(CO2_flux_data[,c(7:10,12)])

plot(CO2_flux_data$HCO3[CO2_flux_data$pH <= 9], CO2_flux_data$DO.umol.L[CO2_flux_data$pH <= 9])











################################################################################################
# Make an daily average file
Flux_daily <- aggregate(. ~ Date, CO2_flux_data, mean)



################################################################################################
#                      check different variable by plotting them                               #  
       plot(Flux_daily[select.list(names(Flux_daily), multiple=TRUE, graphics=TRUE)])          #
#
################################################################################################

summary(Flux_daily)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -823.80 -233.30 -140.20 -163.90  -64.12   59.38     627 




names(CO2_flux_data)
#extract from the raw data the daily peaks.
Gas_min <- aggregate(. ~ CO2_flux_data$Date, CO2_flux_data[,c("CO2umol", "HCO3umol", "CO3umol", "DICumol", "DO.umol.L")], min)
Gas_max <- aggregate(. ~ CO2_flux_data$Date, CO2_flux_data[,c("CO2umol", "HCO3umol", "CO3umol", "DICumol", "DO.umol.L")], max)
Gas_median <- aggregate(. ~ CO2_flux_data$Date, CO2_flux_data[,c("CO2umol", "HCO3umol", "CO3umol", "DICumol", "DO.umol.L")], median)

Gas_change <- data.frame(Gas_max, Gas_min[,c(2:ncol(Gas_min))], Gas_median[,c(2:ncol(Gas_median))])
names(Gas_change)
colnames(Gas_change) [c(1:ncol(Gas_change))] <- c("Date","CO2_max", "HCO3_max", "CO3_max", "DIC_max","DO_max", "CO2_min", "HCO3_min", "CO3_min", "DIC_min","DO_min","CO2_median", "HCO3_median", "CO3_median", "DIC_median","DO_median")

Gas_change$CO2_diff <- Gas_change$CO2_max - Gas_change$CO2_min
Gas_change$DO_diff <- Gas_change$DO_max - Gas_change$DO_min
Gas_change$DIC_diff <- Gas_change$DIC_max - Gas_change$DIC_min
summary(Gas_change)

par(mfrow = c(3,1))
plot(Gas_change$Date, Gas_change$DO_diff, ylab = "Daily absolute change in DO (DO max- DO min)", xlab = "Month (2016)", ylim = c(0,500))
plot(Gas_change$Date, Gas_change$DIC_diff, ylab = "Daily absolute change in DIC (DIC max- DIC min)", xlab = "Month (2016)", ylim = c(0,500))
plot(Gas_change$Date, Gas_change$CO2_diff, ylab = "Daily absolute change in C02 (CO2 max- CO2 min)", xlab = "Month (2016)")
par(mfrow = c(1,1))

par(mfrow = c(2,1))
plot(Gas_change$DO_diff, Gas_change$CO2_diff, xlab = "Daily absolute change in DO (DO max- DO min)", ylab = "Daily absolute change in C02 (C02 max- CO2 min)")
plot(Gas_change$DO_diff, Gas_change$DIC_diff, xlab = "Daily absolute change in DO (DO max- DO min)", ylab = "Daily absolute change in DIC (DIC max- DIC min)")
plot(log(Gas_change$DO_diff), log(Gas_change$CO2_diff), xlab = "Daily absolute change in DO (DO max- DO min) log10", ylab = "Daily absolute change in C02 (C02 max- CO2 min) log 10")
par(mfrow = c(1,1))

cor(Gas_change$DO_diff, Gas_change$DIC_diff)





###############################################
##Check the smooting function on daily basis 
smooth_CO2 <- data.frame(CO2_flux_data[c(30:125),])

par(mar = c(5,5,2,5))
CO2_smoothing_perDate_test <- ggplot(smooth_CO2, aes(Average.Sampling.time, CO2.water.umol.L)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit=fit<<-..y..), n=96) +
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
summary(lm(fit~CO2_flux_data$CO2.water.umol.L[c(30:125)]))







##########################################################################################################################
## Select for data points per day in order to extract the same number of points
## for now days with full data sets were used 

# There are missing values in the data set, so the dataframe was split up in several dataframes. These include high and low values for CO2.
Smoothed_data <- CO2_flux_data[complete.cases(CO2_flux_data[,c("CO2umol", "HCO3umol", "CO3umol", "DICumol", "DO.umol.L")]),]

count(count(Smoothed_data$Date)$freq)
counts <- count(count(Smoothed_data$Date)$freq)$x
Sub_count <- list()

for(i in counts){
  Sub_count[[i]] <- Smoothed_data[Smoothed_data$Date %in% as.Date(names(table(Smoothed_data$Date))[table(Smoothed_data$Date) == i]),]
}


# Try with more time the LOESS function per date instead of ggplot(LOESS)
Dates_inData <- as.Date(count(Smoothed_data$Date)$x)

Sub_Date <- list()
for(i in Dates_inData){
  Sub_Date[[i]] <- Smoothed_data[Smoothed_data$Date %in% as.Date(names(table(Smoothed_data$Date))) == "2016-11-06",]
}


counts

Sub_count_96 <- as.data.frame(Sub_count[[96]])

#Sub_count_20 <- as.data.frame(Sub_count[[45]])
#Sub_count_20 <- as.data.frame(Sub_count[[87]])
#Sub_count_20 <- as.data.frame(Sub_count[[20]])
#Sub_count_20 <- as.data.frame(Sub_count[[20]])
#Sub_count_20 <- as.data.frame(Sub_count[[20]])
#Sub_count_20 <- as.data.frame(Sub_count[[20]])




#apply the smoothing function on daily basis over the whole dataset
CO2_smoothing_96 <- ggplot(Sub_count_96, aes(Average.Sampling.time, CO2.water.umol.L)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit=fit_CO2<<-..y..), n=96) +
  facet_wrap(~Date)
CO2_smoothing_96

Sub_count_96$CO2_fit <- ggplot_build(CO2_smoothing_96)$data[[2]]$y

plot(Sub_count_96$Average.Sampling.time, Sub_count_96$CO2_fit, type = 'l')


#the fit looks good
plot(Sub_count_96$Average.Sampling.time, Sub_count_96$CO2.water.umol.L, type = 'l')
lines(Sub_count_96$Average.Sampling.time, Sub_count_96$CO2_fit, col='red')
summary(lm(Sub_count_96$CO2.water.umol.L ~ Sub_count_96$CO2_fit))





###########################################################################
## check the fluxes for CO2, HCO3, CO3, DIC, and DO




summary(Smoothed_data)
# There are no missing values anymore in HCO3



par(mfrow = c(2,1))
plot(Smoothed_data$Average.Sampling.time, Smoothed_data$HCO3umol, type = 'l') 
plot(Smoothed_data$Average.Sampling.time, Smoothed_data$CO2.water.umol.L, type = 'l') 
par(mfrow = c(1,1))

plot(Smoothed_data$HCO3umol, Smoothed_data$CO2.water.umol.L)
plot(Smoothed_data$CO2.water.umol.L, Smoothed_data$DO.umol.L)
plot(Smoothed_data$HCO3umol, Smoothed_data$DO.umol.L)
plot(Smoothed_data$DICumol, Smoothed_data$DO.umol.L)

par(mfrow = c(2,1))
plot(Smoothed_data$Average.Sampling.time[c(1250:1350)], Smoothed_data$DICumol[c(1250:1350)], type = 'l') 
plot(Smoothed_data$Average.Sampling.time[c(1250:1350)], Smoothed_data$DO.umol.L[c(1250:1350)], type = 'l')
par(mfrow = c(1,1))

par(mfrow = c(2,1))
plot(Smoothed_data$Average.Sampling.time[c(1250:1350)], Smoothed_data$DICumol[c(1250:1350)], type = 'l') 
plot(Smoothed_data$Average.Sampling.time[c(1250:1350)], Smoothed_data$CO2.water.umol.L[c(1250:1350)], type = 'l')
par(mfrow = c(1,1))

par(mfrow = c(2,1))
plot(Smoothed_data$Average.Sampling.time[c(1:100)], Smoothed_data$DICumol[c(1:100)], type = 'l') 
plot(Smoothed_data$Average.Sampling.time[c(1:100)], Smoothed_data$DO.umol.L[c(1:100)], type = 'l')
par(mfrow = c(1,1))



summary(lm(Smoothed_data$DICumol ~ Smoothed_data$DO.umol.L))




#########################################################################################
#apply the smoothing function on daily basis over the whole dataset 

##Check first the daily fluxes
Sub_count_96$Date
Sub_count_96_1 <-  with(Sub_count_96, Sub_count_96[(Date == "2016-06-08"),])


par(mar = c(5,5,2,5))
CO2_smoothing_perDate_test <- ggplot(Sub_count_96_1, aes(Average.Sampling.time, DICumol)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.2, formula = y~x,  aes(outfit=fit<<-..y..), n=96) +
  theme(
    #panel.border = element_rect(fill = NA, colour = "black", size = 1), 
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0,5),
        axis.line.y = element_line(color="black", size = 0,5)) +
  labs(
    y = "umol L-1",
    x = "Time in hours"
  )
CO2_smoothing_perDate_test
fit

ggplot_build(CO2_smoothing_perDate_test)$data[[2]]$y





#CO2
Smoothing_CO2_perDate <- ggplot(Sub_count_96, aes(Average.Sampling.time, CO2.water.umol.L)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit = Smooth_fit_CO2 <<-..y..), n=96) +
  facet_wrap(~Date)
Smoothing_CO2_perDate

Sub_count_96$CO2_fit <- ggplot_build(Smoothing_CO2_perDate)$data[[2]]$y

# DO umol/L #apply the smoothing function on daily basis over the whole dataset 
Smoothing_DO_perDate <- ggplot(Sub_count_96, aes(Average.Sampling.time, DO.umol.L)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit = Smooth_fit_DO <<-..y..), n=96) +
  facet_wrap(~Date)
Smoothing_DO_perDate

Sub_count_96$DO_fit <- ggplot_build(Smoothing_DO_perDate)$data[[2]]$y


# HCO3 umol/L #apply the smoothing function on daily basis over the whole dataset 
Smoothing_HCO3_perDate <- ggplot(Sub_count_96, aes(Average.Sampling.time, HCO3umol)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit = Smooth_fit_HCO3 <<-..y..), n=96) +
  facet_wrap(~Date)
Smoothing_HCO3_perDate

Sub_count_96$HCO3_fit <- ggplot_build(Smoothing_HCO3_perDate)$data[[2]]$y


# CO3 umol/L #apply the smoothing function on daily basis over the whole dataset 
Smoothing_CO3_perDate <- ggplot(Sub_count_96, aes(Average.Sampling.time, CO3umol)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit = Smooth_fit_CO3 <<-..y..), n=96) +
  facet_wrap(~Date)
Smoothing_CO3_perDate

Sub_count_96$CO3_fit <- ggplot_build(Smoothing_CO3_perDate)$data[[2]]$y


# DIC umol/L #apply the smoothing function on daily basis over the whole dataset 
Smoothing_DIC_perDate <- ggplot(Sub_count_96, aes(Average.Sampling.time, DICumol)) +
  geom_point() +
  stat_smooth(method = "loess", fill="lightblue", colour="darkblue", size=0.5, span = 0.4, aes(outfit = Smooth_fit_DIC <<-..y..), n=96) +
  facet_wrap(~Date)
Smoothing_DIC_perDate

Sub_count_96$DIC_fit <- ggplot_build(Smoothing_DIC_perDate)$data[[2]]$y




###########################################################
# Check the LOESS curve

#the fit looks good
plot(Sub_count_96$Average.Sampling.time, Sub_count_96$CO2.water.umol.L, type = 'l')
lines(Sub_count_96$Average.Sampling.time, Sub_count_96$CO2_fit, col='red')
summary(lm(Sub_count_96$CO2.water.umol.L ~ Sub_count_96$CO2_fit))


#the fit looks good
plot(Sub_count_96$Average.Sampling.time, Sub_count_96$DICumol, type = 'l')
lines(Sub_count_96$Average.Sampling.time, Sub_count_96$DIC_fit, col='red')
summary(lm(Sub_count_96$DICumol ~ Sub_count_96$DIC_fit))


#the fit looks good
plot(Sub_count_96$Average.Sampling.time, Sub_count_96$DO.umol.L, type = 'l')
lines(Sub_count_96$Average.Sampling.time, Sub_count_96$DO_fit, col='red')
summary(lm(Sub_count_96$DO.umol.L ~ Sub_count_96$DO_fit))

#the fit looks good
plot(Sub_count_96$Average.Sampling.time, Sub_count_96$HCO3umol, type = 'l')
lines(Sub_count_96$Average.Sampling.time, Sub_count_96$HCO3_fit, col='red')
summary(lm(Sub_count_96$CO2.water.umol.L ~ Sub_count_96$HCO3_fit))

#the fit looks good
plot(Sub_count_96$Average.Sampling.time, Sub_count_96$CO3umol, type = 'l')
lines(Sub_count_96$Average.Sampling.time, Sub_count_96$CO3_fit, col='red')
summary(lm(Sub_count_96$CO2.water.umol.L ~ Sub_count_96$CO3_fit))



#chech the relation of the fits. 
plot(Sub_count_96[,c("CO2_fit", "DIC_fit", "DO_fit")])
par(mfrow = c(2,2))
plot(lm(Sub_count_96$DO_fit ~ Sub_count_96$CO2_fit))
plot(lm(Sub_count_96$DO_fit ~ Sub_count_96$DIC_fit))
plot(lm(Sub_count_96$CO2_fit ~ Sub_count_96$DIC_fit))  # very bad fit




summary(lm(Sub_count_96$DO_fit ~ Sub_count_96$CO2_fit))
summary(lm(Sub_count_96$DO_fit ~ Sub_count_96$DIC_fit))
summary(lm(Sub_count_96$CO2_fit ~ Sub_count_96$DIC_fit))



#####################################################################
## Daily changes in DIC and DO
## Build a new dataframe for the daily DO and CO2 range

Gas_min_smooth <- aggregate(. ~ Date, Sub_count_96[,c("Date", "DO_fit", "CO2_fit", "HCO3_fit", "DIC_fit")], min)
Gas_max_smooth <- aggregate(. ~ Date, Sub_count_96[,c("Date", "DO_fit", "CO2_fit", "HCO3_fit", "DIC_fit")], max)
Gas_median_smooth <- aggregate(. ~ Date, Sub_count_96[,c("Date", "DO_fit", "CO2_fit", "HCO3_fit", "DIC_fit")], median)


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


plot(Sub_count_96$Average.Sampling.time, Sub_count_96$DIC_fit, col = "darkorange", type = 'l', ylim = c(0,5000))
points(Sub_count_96$Average.Sampling.time, Sub_count_96$DO.umol.L, col="darkblue", type = 'l')


par(mfrow = c(2,2))
plot(Sub_count_96$DO_fit, Sub_count_96$CO2_fit)
plot(Sub_count_96$DO_fit, Sub_count_96$HCO3_fit)
plot(Sub_count_96$DO_fit, Sub_count_96$DIC_fit)
plot(Sub_count_96$CO2_fit, Sub_count_96$DIC_fit)
par(mfrow = c(1,1))

plot(log(Smoothed_data$DO_fit), log(Smoothed_data$CO2_fit))
# There is no linear relation in the whole data set. Furthermore, the months are grouped in the plot.

###########################################
## Check the slope over a 24 hour periods


#DIC
names(table(Smoothed_data$Date))
par(mfrow = c(2,2))
plot(Sub_count_96$DO_fit[Sub_count_96$Date == "2016-05-20"], Sub_count_96$DIC_fit[Sub_count_96$Date == "2016-05-20"])
plot(Sub_count_96$DO_fit[Sub_count_96$Date == "2016-06-18"], Sub_count_96$DIC_fit[Sub_count_96$Date == "2016-06-18"])
plot(Sub_count_96$DO_fit[Sub_count_96$Date == "2016-08-09"], Sub_count_96$DIC_fit[Sub_count_96$Date == "2016-08-09"])
plot(Sub_count_96$DO_fit[Sub_count_96$Date == "2016-10-10"], Sub_count_96$DIC_fit[Sub_count_96$Date == "2016-10-10"])
par(mfrow = c(1,1))


#CO2
par(mfrow = c(2,2))
plot(Sub_count_96$DO_fit[Sub_count_96$Date == "2016-05-20"], Sub_count_96$CO2_fit[Sub_count_96$Date == "2016-05-20"])
plot(Sub_count_96$DO_fit[Sub_count_96$Date == "2016-06-18"], Sub_count_96$CO2_fit[Sub_count_96$Date == "2016-06-18"])
plot(Sub_count_96$DO_fit[Sub_count_96$Date == "2016-08-09"], Sub_count_96$CO2_fit[Sub_count_96$Date == "2016-08-09"])
plot(Sub_count_96$DO_fit[Sub_count_96$Date == "2016-10-10"], Sub_count_96$CO2_fit[Sub_count_96$Date == "2016-10-10"])
par(mfrow = c(1,1))



################################
#extract the slope over an hour period of time
ind_j <- 1:(nrow(Sub_count_96)-4)
ind_k <- 4:nrow(Sub_count_96)
ind_w <- 1:nrow(Sub_count_96)
lm_daily_slope <- data.frame()
lm_daily_slope[c(1:nrow(Sub_count_96)),] <- data.frame()

#CO2
CO2_slope_day <- summary(lm(CO2_fit ~ ., data=Sub_count_96))

plot(CO2_slope_day$coefficients[c(3:649)], ylim = c(-0.000001, 0.000001))

for (i in seq_along(ind_j)){
  lm_daily_slope$slopeCO2[ind_w[i]] <- lm(Sub_count_96$CO2_fit[ind_j[i]:ind_k[i]] ~ Sub_count_96$Average.Sampling.time[ind_j[i]:ind_k[i]])$coefficients[[2]]
}

summary(lm(Sub_count_96$CO2_fit[c(1:4)] ~ Sub_count_96$DO_fit[1:4]))
plot(Sub_count_96$CO2_fit[c(1:4)] , Sub_count_96$DO_fit[1:4])


#DIC
for (i in seq_along(ind_j)){
  lm_daily_slope$slopeCO2[ind_w[i]] <- lm(Sub_count_96$CO2_fit[ind_j[i]:ind_k[i]] ~ Sub_count_96$Average.Sampling.time[ind_j[i]:ind_k[i]])$coefficients[[2]]
}

summary(lm(Sub_count_96$DIC_fit[c(1:4)] ~ Sub_count_96$DO_fit[1:4]))
plot(Sub_count_96$DIC_fit[c(1:4)] , Sub_count_96$DO_fit[1:4])

for (i in seq_along(ind_j)){
  lm_daily_slope$slopeDIC[ind_w[i]] <- lm(Sub_count_96$DIC_fit[ind_j[i]:ind_k[i]] ~ Sub_count_96$DO_fit[ind_j[i]:ind_k[i]])$coefficients[[2]]
}



par(mfrow = c(2,2))
plot(lm_daily_slope$slopeCO2, ylab = "the change of CO2 relative to DO", xlab = "hours", ylim = c(-1, 1), type = 'l', xlim = c(1,48))
plot(lm_daily_slope$slopeDIC, ylab = "the change of DIC relative to DO", xlab = "hours", ylim = c(-100, 100), type = 'l', xlim = c(1,48))
plot(lm_daily_slope$slopeCO2, ylab = "the change of CO2 relative to DO", xlab = "hours", ylim = c(-0.5, 0.5), type = 'l', xlim = c(2000,2048))
plot(lm_daily_slope$slopeDIC, ylab = "the change of DIC relative to DO", xlab = "hours", ylim = c(-100, 100), type = 'l', xlim = c(2000,2048))
par(mfrow = c(1,1))


summary(lm_daily_slope$slopeCO2)



#check the relative change (max-min) compared to the median of each day
############## For the measured data
DIC_median <- aggregate(DIC_fit ~ Date, data = Sub_count_96, median)
DIC_median <- DIC_median[rep(row.names(DIC_median), 96), 1:2]
DIC_median <- DIC_median[order(DIC_median$Date),] 
Sub_count_96$DIC_median <- DIC_median$DIC_fit

CO2_median <- aggregate(CO2_fit ~ Date, data = Sub_count_96, median)
CO2_median <- CO2_median[rep(row.names(CO2_median), 96), 1:2]
CO2_median <- CO2_median[order(CO2_median$Date),] 
Sub_count_96$CO2_median <- CO2_median$CO2_fit

DO_median <- aggregate(DO_fit ~ Date, data = Sub_count_96, median)
DO_median <- DO_median[rep(row.names(DO_median), 96), 1:2]
DO_median <- DO_median[order(DO_median$Date),] 
Sub_count_96$DO_median <- DO_median$DO_fit

HCO3_median <- aggregate(HCO3 ~ Date, data = Sub_count_96, median)
HCO3_median <- HCO3_median[rep(row.names(HCO3_median), 24), 1:2]
HCO3_median <- HCO3_median[order(HCO3_median$Date),] 
Sub_count_96$HCO3_median <- HCO3_median$HCO3


Sub_count_96$CO2_diff_rel <- (Sub_count_96$CO2_fit/Sub_count_96$CO2_median)
Sub_count_96$DO_diff_rel <- (Sub_count_96$DO_fit/Sub_count_96$DO_median)
Sub_count_96$HCO3_diff_rel <- (Sub_count_96$HCO3_fit/Sub_count_96$HCO3_median)
Sub_count_96$DIC_diff_rel <- (Sub_count_96$DIC_fit/Sub_count_96$DIC_median)

plot(Sub_count_96$Average.Sampling.time, Sub_count_96$CO2_diff_rel)
plot(Sub_count_96$Average.Sampling.time, Sub_count_96$DIC_diff_rel)
plot(Sub_count_96$Average.Sampling.time, Sub_count_96$DO_diff_rel)
plot(Sub_count_96$Average.Sampling.time, Sub_count_96$HCO3_diff_rel)
plot(Sub_count_96$CO2_diff_rel-1, Sub_count_96$DO_diff_rel-1)
plot(Sub_count_96$DIC_diff_rel-1, Sub_count_96$DO_diff_rel-1)
plot(Sub_count_96$CO2_diff_rel-1, Sub_count_96$HCO3_diff_rel-1)
plot(Sub_count_96$HCO3_diff_rel-1, Sub_count_96$DO_diff_rel-1)


# check the linear regression
# linearity in the model is lacking, which is still lacking with a log transformation
# There is linearity when the model is plit up different months

Sub_count_96 <- with(Smoothed_data, Smoothed_data[
  (Date >= "2016-05-01" & Date < "2016-11-07"),])
par(mfrow=c(2,2))
plot(lm(Sub_count_96$DIC_diff_rel ~ Sub_count_96$DO_diff_rel))
summary(lm(Sub_count_96$DIC_diff_rel ~ Sub_count_96$DO_diff_rel))
par(mfrow = c(1,1)) 

plot(Smoothed_data_lm_adapted$CO2_diff_rel, Smoothed_data_lm_adapted$DO_diff_rel)
# The relative changes in CO2 




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








##########################################################################################
## Save document in excel 
library(xlsx)
#write.table(Sub_count_96, "C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/DO_DIC_fit_probe.txt", sep = ';')
write.table(Sub_count_96, "C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/DO_DIC_fit_alk.txt", sep = ';')




























########################################
# Check the linear regression per month


months(Sub_count_96$Date)
plot(subset(Sub_count_96$CO2_diff_rel, months(Sub_count_96$Date) == "mei") , subset(Sub_count_96$DO_diff_rel, months(Sub_count_96$Date) == "mei"))
summary(lm(subset(Sub_count_96$CO2_diff_rel, months(Sub_count_96$Date) == "July") ~ subset(Sub_count_96$DO_diff_rel, months(Sub_count_96$Date) == "July")))
par(mfrow = c(2,2))
plot(lm(log(subset(Sub_count_96$CO2_diff_rel, months(Sub_count_96$Date) == "July")) ~ log(subset(Sub_count_96$DO_diff_rel, months(Sub_count_96$Date) == "July"))))
par(mfrow = c(1,1))


##########################
## Check the slope by 24 hours

ind_j <- 1:2736
ind_k <- 24:2760
ind_w <- 1:2760
lm_week <- data.frame()
lm_week[c(1:2736),] <- data.frame()

for (i in seq_along(ind_j)){
  lm_week$slopeCO2[ind_w[i]] <- lm(Smoothed_data$CO2_diff_rel[ind_j[i]:ind_k[i]] ~ Smoothed_data$DO_diff_rel[ind_j[i]:ind_k[i]])$coefficients[[2]]
}

for (i in seq_along(ind_j)){
  lm_week$slopeDIC[ind_w[i]] <- lm(Smoothed_data$DIC_diff_rel[ind_j[i]:ind_k[i]] ~ Smoothed_data$DO_diff_rel[ind_j[i]:ind_k[i]])$coefficients[[2]]
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



plot(Gas_change_smooth$Date, Gas_change_smooth$CO2_max, col="red")
points(Gas_change_smooth$Date, Gas_change_smooth$CO2_min, col = "green")
lines(Gas_change_smooth$Date, Gas_change_smooth$CO2_median)









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






