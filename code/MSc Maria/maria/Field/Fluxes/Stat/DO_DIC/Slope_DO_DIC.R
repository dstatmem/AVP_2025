# Info this file is based on the output of the rscript DODIC_probe
# The max and min per day were calculated, plus the change in DO and DIC. 






#If there are unstalled packages use: install.packages("name_package")
library(car)
library(vegan)
library(nlme)
library(plyr)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)


library(readxl)
DO_DIC_fit <- read_excel("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/DO_DIC_fit_alk.xlsx")


#Chewbacca
setwd("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/DO_DIC/")
DO_DIC_fit <- read.csv("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/DO_DIC_fit_alk.csv")

#Przemek
setwd("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/DO_DIC/")
DO_DIC_fit <- read.csv("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/DO_DIC_fit_alk.csv")

#UvA
setwd("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/DO_DIC/")
DO_DIC_fit <- read.csv("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/DO_DIC_fit_alk.csv")


#time conversion
DO_DIC_fit$Average.Sampling.time <- as.POSIXct(DO_DIC_fit$Average.Sampling.time, format = "%m/%d/%Y %H:%M")  #In case of error check the seperator in the date i.e %m-%d-%Y or %m/%d/%Y
DO_DIC_fit$Date <- as.Date(DO_DIC_fit$Date, "%m/%d/%Y")
DO_DIC_fit$Time <- format(DO_DIC_fit$Average.Sampling.time, "%H:%M")

months(DO_DIC_fit$Date)

#DO_DIC_fit$Time1 <- hour(DO_DIC_fit$Average.Sampling.time) + minute(DO_DIC_fit$Average.Sampling.time)

#Remove the missing values 
#DO_DIC_fit <- DO_DIC_fit[complete.cases(DO_DIC_fit[,c("H")]),]




CO2_DO_cor <- ggplot(DO_DIC_fit, aes(DO_smooth_av,  DIC_smooth_av, colour = months(DO_DIC_fit$Date))) +
  geom_hline(yintercept = 0, colour = "grey" , linetype = 2 ) + geom_vline(xintercept = 0, colour = "grey" , linetype = 2 ) +
  geom_point() +
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


##########################################################################################################
# Extract the min, and max values

DO_DIC_fit %>%
  group_by(Date) %>%
  summarise(DO_min = min(DO_smooth_av),
            DO_max = max(DO_smooth_av),
            DIC_min = min(DIC_smooth_av),
            DIC_max = max(DIC_smooth_av),
            CO2_min = min(CO2_smooth_av),
            CO2_max = max(CO2_smooth_av)) -> Min_max

Min_max <- Min_max[-118,]




Min_max$CO2_diff <- Min_max$CO2_max - Min_max$CO2_min
Min_max$DO_diff <- Min_max$DO_max - Min_max$DO_min
Min_max$DIC_diff <- Min_max$DIC_max - Min_max$DIC_min
summary(Min_max)


###########################################################################################################
# find the location of min-max vlaues and add those to the min-max file


Min_max$Time_DOmax <- do.call(rbind, by(data = DO_DIC_fit,
                                        INDICES = list(DO_DIC_fit$Date),
                                        function(x) x[which.max(x$DO_smooth_av),
                                                      ]))$Average.Sampling.time
Min_max$Time_DICmax <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DIC_smooth_av), ]))$Average.Sampling.time
Min_max$Time_CO2max  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$CO2_smooth_av), ]))$Average.Sampling.time
Min_max$Time_DOmin <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DO_smooth_av), ]))$Average.Sampling.time
Min_max$Time_DICmin  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DIC_smooth_av), ]))$Average.Sampling.time
Min_max$Time_CO2min <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$CO2_smooth_av), ]))$Average.Sampling.time


#plot the max min values 
par(mfrow = c(3,1))
plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$DO_smooth_av, type = 'l')
points(Min_max$Time_DOmax, Min_max$DO_max, col = "red")
points(Min_max$Time_DOmin, Min_max$DO_min, col = "darkgreen")
plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$DIC_smooth_av, type = 'l')
points(Min_max$Time_DICmax, Min_max$DIC_max, col = "red")
points(Min_max$Time_DICmin, Min_max$DIC_min, col = "darkgreen")
plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$CO2_smooth_av, type = 'l')
points(Min_max$Time_CO2max, Min_max$CO2_max, col = "red")
points(Min_max$Time_CO2min, Min_max$CO2_min, col = "darkgreen")

# Some points do not fit, check between which time the most points lay
Min_max$Time_DOmax <- hour(do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DO_smooth_av), ]))$Average.Sampling.time)
Min_max$Time_DICmax <- hour(do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DIC_smooth_av), ]))$Average.Sampling.time)
Min_max$Time_CO2max  <- hour(do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$CO2_smooth_av), ]))$Average.Sampling.time)
Min_max$Time_DOmin <- hour(do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DO_smooth_av), ]))$Average.Sampling.time)
Min_max$Time_DICmin  <- hour(do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DIC_smooth_av), ]))$Average.Sampling.time)
Min_max$Time_CO2min <- hour(do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$CO2_smooth_av), ]))$Average.Sampling.time)

par(mfrow = c(2,3))
hist(Min_max$Time_DOmax, ylab = "DOmax count", xlab = "Time")
hist(Min_max$Time_DICmax, ylab = "DICmax count", xlab = "Time")
hist(Min_max$Time_CO2max, ylab = "CO2max count", xlab = "Time")
hist(Min_max$Time_DOmin, ylab = "DOmin count", xlab = "Time")
hist(Min_max$Time_DICmin, ylab = "DICmin count", xlab = "Time")
hist(Min_max$Time_CO2min, ylab = "CO2min count", xlab = "Time")

table(Min_max$Time_DOmax)
table(Min_max$Time_DICmax)
table(Min_max$Time_CO2max)
table(Min_max$Time_DOmin)
table(Min_max$Time_DICmin)
table(Min_max$Time_CO2min)

#redo the selection for min and max but then over a certain time frame within one day


Min_max_new <- data.frame()
Min_max_new <- Min_max_new[c(1:nrow(Min_max)),]

table(Min_max$Time_DOmax)
DO_DIC_fit_cropped_DOmax <- DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(11:23),]
Min_max_new$Time_DOmax <- do.call(rbind, by(data = DO_DIC_fit_cropped_DOmax, INDICES = list(DO_DIC_fit_cropped_DOmax$Date), function(x) x[which.max(x$DO_smooth_av), ]))$Average.Sampling.time

table(Min_max$Time_DICmax)
DO_DIC_fit_cropped_DICmax <- DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(0:16),]
Min_max_new$Time_DICmax <- do.call(rbind, by(data = DO_DIC_fit_cropped_DICmax, INDICES = list(DO_DIC_fit_cropped_DICmax$Date), function(x) x[which.max(x$DIC_smooth_av), ]))$Average.Sampling.time

table(Min_max$Time_CO2max)
DO_DIC_fit_cropped_CO2max <- DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(0:16),]
Min_max_new$Time_CO2max  <- do.call(rbind, by(data = DO_DIC_fit_cropped_CO2max, INDICES = list(DO_DIC_fit_cropped_CO2max$Date), function(x) x[which.max(x$CO2_smooth_av), ]))$Average.Sampling.time

table(Min_max$Time_DOmin)
DO_DIC_fit_cropped_DOmin <- DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(0:16),]
Min_max_new$Time_DOmin <- do.call(rbind, by(data = DO_DIC_fit_cropped_DOmin, INDICES = list(DO_DIC_fit_cropped_DOmin$Date), function(x) x[which.min(x$DO_smooth_av), ]))$Average.Sampling.time

table(Min_max$Time_DICmin)
DO_DIC_fit_cropped_DICmin <- DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(9:23),]
Min_max_new$Time_DICmin  <- do.call(rbind, by(data = DO_DIC_fit_cropped_DICmin, INDICES = list(DO_DIC_fit_cropped_DICmin$Date), function(x) x[which.min(x$DIC_smooth_av), ]))$Average.Sampling.time

table(Min_max$Time_CO2min)
DO_DIC_fit_cropped_CO2min <- DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(11:23),]
Min_max_new$Time_CO2min  <- do.call(rbind, by(data = DO_DIC_fit_cropped_DICmin, INDICES = list(DO_DIC_fit_cropped_DICmin$Date), function(x) x[which.min(x$CO2_smooth_av), ]))$Average.Sampling.time




DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(14:21),] %>%
  group_by(Date) %>%
  summarise(max(DO_smooth_av)) -> DO_max
DO_max <- DO_max[-118,]

DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(2:15),] %>%
  group_by(Date) %>%
  summarise(max(DIC_smooth_av)) -> DIC_max
DIC_max <- DIC_max[-118,]

DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(1:15),] %>%
  group_by(Date) %>%
  summarise(max(CO2_smooth_av)) -> CO2_max
CO2_max <- CO2_max[-118,]

DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(3:15),] %>%
  group_by(Date) %>%
  summarise(min(DO_smooth_av)) -> DO_min
DO_min <- DO_min[-118,]

DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(14:21),] %>%
  group_by(Date) %>%
  summarise(min(DIC_smooth_av)) -> DIC_min
DIC_min <- DIC_min[-118,]

DO_DIC_fit[hour(DO_DIC_fit$Average.Sampling.time) == c(13:22),] %>%
  group_by(Date) %>%
  summarise(min(CO2_smooth_av)) -> CO2_min
CO2_min <- CO2_min[-118,]






par(mfrow = c(3,1))
plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$DO_smooth_av, type = 'l', ylab = "Dissolved oxygen (umol/L)", xlab = "month (2016)")
points(Min_max_new$Time_DOmax, DO_max$`max(DO_smooth_av)`, col = "darkgreen")
points(Min_max_new$Time_DOmin, DO_min$`min(DO_smooth_av)`, col = "red")
plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$DIC_smooth_av, type = 'l', ylab = "Dissolved inorganic carbon (umol/L)", xlab = "month (2016)")
points(Min_max_new$Time_DICmax, DIC_max$`max(DIC_smooth_av)`, col = "darkgreen")
points(Min_max_new$Time_DICmin, DIC_min$`min(DIC_smooth_av)`, col = "red")
plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$CO2_smooth_av, type = 'l', ylab = "CO2 (umol/L)", xlab = "month (2016)")
points(Min_max_new$Time_CO2max, CO2_max$`max(CO2_smooth_av)`, col = "darkgreen")
points(Min_max_new$Time_CO2min, CO2_min$`min(CO2_smooth_av)`, col = "red")


# The peaks of DO predict the peak of DIC and CO2 well, therefore the 
par(mfrow = c(1,1))
plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$DO_smooth_av, type = 'l', ylab = "Dissolved oxygen (umol/L)", xlab = "month (2016)")
lines(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$DIC_smooth_av/2, col = "red")
lines(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$DIC_smooth_av/5, col = "darkgreen")

plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$DIC_smooth_av, type = 'l', ylab = "Dissolved inorganic carbon (umol/L)", xlab = "month (2016)")
points(Min_max_new$Time_DOmax, DIC_min$`min(DIC_smooth_av)`, col = "darkgreen")
points(Min_max_new$Time_DOmin, DIC_max$`max(DIC_smooth_av)`, col = "red")


plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$CO2_smooth_av, type = 'l', ylab = "CO2 (umol/L)", xlab = "month (2016)")
points(Min_max_new$Time_DOmax, CO2_min$`min(CO2_smooth_av)`, col = "darkgreen")
points(Min_max_new$Time_DOmin, CO2_max$`max(CO2_smooth_av)`, col = "red")



plot(day(Min_max_new$Time_DICmax), (Min_max_new$Time_DICmax - Min_max_new$Time_DOmin)/60)
(Min_max_new$Time_CO2max - Min_max_new$Time_DOmin)/60

(Min_max_new$Time_DICmin - Min_max_new$Time_DOmax)/60
(Min_max_new$Time_DICmax - Min_max_new$Time_DOmin)/60


# Check the min-max value if the time for DO is used
df_DODIC_min <- DO_DIC_fit[DO_DIC_fit$Average.Sampling.time %in% Min_max_new$Time_DOmin,]
df_DODIC_max <- DO_DIC_fit[DO_DIC_fit$Average.Sampling.time %in% Min_max_new$Time_DOmax,]

#plots based on the min max of DO
par(mfrow = c(1,1))
plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$DO_smooth_av, type = 'l', ylab = "Dissolved oxygen (umol/L)", xlab = "month (2016)")
points(df_DODIC_min$Average.Sampling.time, df_DODIC_min$DO_smooth_av, col = "darkgreen")
points(df_DODIC_max$Average.Sampling.time, df_DODIC_max$DO_smooth_av, col = "red")
plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$DIC_smooth_av, type = 'l', ylab = "Dissolved inorganic carbon (umol/L)", xlab = "month (2016)")
points(df_DODIC_min$Average.Sampling.time, df_DODIC_min$DIC_smooth_av, col = "darkgreen")
points(df_DODIC_max$Average.Sampling.time, df_DODIC_max$DIC_smooth_av, col = "red")
plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$CO2_smooth_av, type = 'l', ylab = "CO2 (umol/L)", xlab = "month (2016)")
points(df_DODIC_min$Average.Sampling.time, df_DODIC_min$CO2_smooth_av, col = "darkgreen")
points(df_DODIC_max$Average.Sampling.time, df_DODIC_max$CO2_smooth_av, col = "red")



#plots based on the min max of DO
par(mfrow = c(1,1))
plot(DO_DIC_fit$Average.Sampling.time[1:960], DO_DIC_fit$DO_smooth_av[1:960], type = 'l', ylab = "Dissolved oxygen (umol/L)", xlab = "month (2016)")
points(df_DODIC_min$Average.Sampling.time, df_DODIC_min$DO_smooth_av, col = "darkgreen")
points(df_DODIC_max$Average.Sampling.time, df_DODIC_max$DO_smooth_av, col = "red")
lines(DO_DIC_fit$Average.Sampling.time[1:960], DO_DIC_fit$DIC_smooth_av[1:960]/10, type = 'l', ylab = "Dissolved inorganic carbon (umol/L)", xlab = "month (2016)")
points(df_DODIC_min$Average.Sampling.time, df_DODIC_min$DIC_smooth_av/10, col = "darkgreen")
points(df_DODIC_max$Average.Sampling.time, df_DODIC_max$DIC_smooth_av/10, col = "red")
plot(DO_DIC_fit$Average.Sampling.time, DO_DIC_fit$CO2_smooth_av, type = 'l', ylab = "CO2 (umol/L)", xlab = "month (2016)")
points(df_DODIC_min$Average.Sampling.time, df_DODIC_min$CO2_smooth_av, col = "darkgreen")
points(df_DODIC_max$Average.Sampling.time, df_DODIC_max$CO2_smooth_av, col = "red")





#####################################################################################################################################
# Check the change per timestep

# The time difference is not always 15 minutes, so calculate and multiply the change per time step times the time difference to get the hourly rate. 

summary(DO_DIC_fit$TimeDif)

ind_j <- 1:(nrow(DO_DIC_fit)-1)
ind_k <- 2:nrow(DO_DIC_fit)
ind_w <- 1:nrow(DO_DIC_fit)


for (i in seq_along(ind_j)){
   DO_DIC_fit$TimeDif[ind_w[i]] <- 60/(as.double(DO_DIC_fit$Average.Sampling.time[ind_k[i]]-DO_DIC_fit$Average.Sampling.time[ind_j[i]]))
}

DO_DIC_fit$DOchange_hour <- DO_DIC_fit$DO_slope*DO_DIC_fit$TimeDif
DO_DIC_fit$DICchange_hour <- DO_DIC_fit$DIC_slope*DO_DIC_fit$TimeDif
DO_DIC_fit$CO2change_hour <- DO_DIC_fit$CO2_slope*DO_DIC_fit$TimeDif


## plot the change
par(mfrow = c(3,1))
plot(DO_DIC_fit$Average.Sampling.time,DO_DIC_fit$DOchange_hour, type = 'l')
plot(DO_DIC_fit$Average.Sampling.time,DO_DIC_fit$DICchange_hour, type = 'l')
plot(DO_DIC_fit$Average.Sampling.time,DO_DIC_fit$CO2change_hour, type = 'l')

#remove extreme values
DO_DIC_fit$DOchange_hour[DO_DIC_fit$DOchange_hour > 200] <- NA
DO_DIC_fit$DOchange_hour[DO_DIC_fit$DOchange_hour < -200] <- NA
DO_DIC_fit$DICchange_hour[DO_DIC_fit$DICchange_hour > 200] <- NA
DO_DIC_fit$DICchange_hour[DO_DIC_fit$DICchange_hour < -200] <- NA
DO_DIC_fit$CO2change_hour[DO_DIC_fit$CO2change_hour > 6] <- NA
DO_DIC_fit$CO2change_hour[DO_DIC_fit$CO2change_hour < -6] <- NA


par(mfrow = c(2,1))
plot(DO_DIC_fit$Average.Sampling.time[1000:1960],DO_DIC_fit$DO_smooth_av[1000:1960], type = 'l')
plot(DO_DIC_fit$Average.Sampling.time[1000:1960],DO_DIC_fit$DOchange_hour[1000:1960], type = 'l')
plot(DO_DIC_fit$Average.Sampling.time[1000:1960],DO_DIC_fit$DICchange_hour[1000:1960], type = 'l')
plot(DO_DIC_fit$Average.Sampling.time[1000:1960],DO_DIC_fit$DIC_smooth_av[1000:1960], type = 'l')
plot(DO_DIC_fit$DOchange_hour,DO_DIC_fit$CO2change_hour)
plot(DO_DIC_fit$CO2change_hour,DO_DIC_fit$DICchange_hour)


par(mfrow = c(2,2))
# there are some outlier series (which are probrably formed by a whole day, cause removing outlier points results in new outliers within the same range),
# plus the serie is not linear due to those outliers.  
plot(DO_DIC_fit$DOchange_hour,DO_DIC_fit$DICchange_hour)
plot(DO_DIC_fit$DOchange_hour,DO_DIC_fit$CO2change_hour)
plot(DO_DIC_fit$CO2change_hour,DO_DIC_fit$DICchange_hour)


plot(lm(DO_DIC_fit$DOchange_hour~DO_DIC_fit$DICchange_hour))
summary(lm(DO_DIC_fit$DOchange_hour~DO_DIC_fit$DICchange_hour))

DO_DIC_fit_red <- with(DO_DIC_fit, DO_DIC_fit[(DO_DIC_fit$DICchange_hour < 0.2 & DO_DIC_fit$DICchange_hour > -0.2) &
                                                (DO_DIC_fit$DOchange_hour < 0.5 & DO_DIC_fit$DOchange_hour > -0.5) &
                                                (DO_DIC_fit$CO2change_hour < 1.8 & DO_DIC_fit$CO2change_hour > -1.8),])

plot(DO_DIC_fit_red$DOchange_hour,DO_DIC_fit_red$DICchange_hour)
plot(DO_DIC_fit_red$DOchange_hour,DO_DIC_fit_red$CO2change_hour)
plot(DO_DIC_fit_red$CO2change_hour,DO_DIC_fit_red$DICchange_hour)

#Reducing the data by removing extreme values doesn't improve the model 
plot(lm(DO_DIC_fit_red$DOchange_hour~DO_DIC_fit_red$DICchange_hour))
summary(lm(DO_DIC_fit_red$DOchange_hour~DO_DIC_fit_red$DICchange_hour))


plot(lm(DO_DIC_fit$DOchange_hour[1:960]~DO_DIC_fit$DICchange_hour[1:960]))
hist(summary(lm(DO_DIC_fit$DOchange_hour[1:960]~DO_DIC_fit$DICchange_hour[1:960]))$residuals)
summary(lm(DO_DIC_fit$DOchange_hour[1000:1960]~DO_DIC_fit$DICchange_hour[1000:1960]))
####
#!!!!!!!!! Check the change of the slope per month !!!!!!!!!!

par(mfrow = c(3,1))
plot(DO_DIC_fit$Average.Sampling.time[1:960],DO_DIC_fit$DICchange_hour[1:960], type = 'l')
plot(DO_DIC_fit$Average.Sampling.time[1:960],DO_DIC_fit$DICchange_hour[1:960], type = 'l')
plot(DO_DIC_fit$Average.Sampling.time[1:960],DO_DIC_fit$CO2change_hour[1:960], type = 'l')




# Check the max change
par(mfrow =c(1,1))
plot(DO_DIC_fit$DOchange_hour, DO_DIC_fit$CHLV.Volt)
plot(DO_DIC_fit$DOchange_hour, DO_DIC_fit$I..umol.m2.s.)
summary(lm(DO_DIC_fit$DOchange_hour ~ DO_DIC_fit$I..umol.m2.s.))
plot(DO_DIC_fit$DICchange_hour, DO_DIC_fit$I..umol.m2.s.)
summary(lm(DO_DIC_fit$DICchange_hour ~ DO_DIC_fit$I..umol.m2.s.))

#check the change per mg chl/L
DO_DIC_fit$DOchange_hour_chl <- DO_DIC_fit$DOchange_hour/DO_DIC_fit$CHL.aver
plot(DO_DIC_fit$DOchange_hour_chl, DO_DIC_fit$I..umol.m2.s.)
summary(lm(DO_DIC_fit$DOchange_hour_chl ~ DO_DIC_fit$I..umol.m2.s.))


DO_DIC_fit %>%
  group_by(Date) %>%
  summarise(DO_min = min(DOchange_hour),
            DO_max = max(DOchange_hour),
            DIC_min = min(DICchange_hour),
            DIC_max = max(DICchange_hour),
            CO2_min = min(CO2change_hour),
            CO2_max = max(CO2change_hour)) -> Min_max_change
Min_max_change <- Min_max_change[-118,]


#extract the values from the dataset for the min and max change
Min_max_change$Time_DOmax <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DOchange_hour), ]))$Average.Sampling.time
Min_max_change$Time_DICmax <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DICchange_hour), ]))$Average.Sampling.time
Min_max_change$Time_CO2max  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$CO2change_hour), ]))$Average.Sampling.time
Min_max_change$Time_DOmin <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DOchange_hour), ]))$Average.Sampling.time
Min_max_change$Time_DICmin  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DICchange_hour), ]))$Average.Sampling.time
Min_max_change$Time_CO2min <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$CO2change_hour), ]))$Average.Sampling.time

Min_max_change$I_DOmax <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DOchange_hour), ]))$I..umol.m2.s.
Min_max_change$I_DICmax <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DICchange_hour), ]))$I..umol.m2.s.
Min_max_change$I_CO2max  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$CO2change_hour), ]))$I..umol.m2.s.
Min_max_change$I_DOmin <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DOchange_hour), ]))$I..umol.m2.s.
Min_max_change$I_DICmin  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DICchange_hour), ]))$I..umol.m2.s.
Min_max_change$I_CO2min <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$CO2change_hour), ]))$I..umol.m2.s.

Min_max_change$T.C_DOmax <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DOchange_hour), ]))$T.C
Min_max_change$T.C_DICmax <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DICchange_hour), ]))$T.C
Min_max_change$T.C_CO2max  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$CO2change_hour), ]))$T.C
Min_max_change$T.C_DOmin <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DOchange_hour), ]))$T.C
Min_max_change$T.C_DICmin  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DICchange_hour), ]))$T.C
Min_max_change$T.C_CO2min <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$CO2change_hour), ]))$T.C


Min_max$T.C_DOmax <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DO_smooth_av), ]))$T.C
Min_max$T.C_DICmax <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DIC_smooth_av), ]))$T.C
Min_max$T.C_CO2max  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$CO2_smooth_av), ]))$T.C
Min_max$T.C_DOmin <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DO_smooth_av), ]))$T.C
Min_max$T.C_DICmin  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DIC_smooth_av), ]))$T.C
Min_max$T.C_CO2min <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$CO2_smooth_av), ]))$T.C

Min_max$I_DOmax <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DO_smooth_av), ]))$I..umol.m2.s.
Min_max$I_DICmax <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$DIC_smooth_av), ]))$I..umol.m2.s.
Min_max$I_CO2max  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.max(x$CO2_smooth_av), ]))$I..umol.m2.s.
Min_max$I_DOmin <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DO_smooth_av), ]))$I..umol.m2.s.
Min_max$I_DICmin  <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$DIC_smooth_av), ]))$I..umol.m2.s.
Min_max$I_CO2min <- do.call(rbind, by(data = DO_DIC_fit, INDICES = list(DO_DIC_fit$Date), function(x) x[which.min(x$CO2_smooth_av), ]))$I..umol.m2.s.

write.table(Min_max, "C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/Min_max.txt", sep = ';')
write.table(Min_max_change, "C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/Min_max_change.txt", sep = ';')
