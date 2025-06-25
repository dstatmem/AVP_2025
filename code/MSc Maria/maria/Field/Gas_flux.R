setwd("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes")
Gas_flux_R <- read.csv("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes/Gas_flux_R.csv")
View(Gas_flux_R)
library(plyr)
library(car)




plot(Gas_flux_R)

lm_GasFlux <- lm(Gas_flux_R$Flux.umol.m2.h.1~Gas_flux_R$I..umol.m2.s. + Gas_flux_R$Chl.fluorescence + Gas_flux_R$Chl.Lab.ug.L)
summary(lm_GasFlux)


names(Gas_flux_R)
daily_flux <-ddply(Gas_flux_R, .(Date), summarize, Temp = mean(T..oC.), CO2_flux = mean(Flux.umol.m2.h.1), Chl_Lab = mean(Chl.Lab.ug.L))

########################################################################
## Check the distribution of the data

# check the distribution of the data and check for individual outliers.


##Temperature
hist(Gas_flux_R$T..oC.)
# looks fine




hist(Gas_flux_R$pH)
shapiro.test(log(Gas_flux_R$pH))
qqnorm(log(Gas_flux_R$pH));qqline(log(Gas_flux_R$pH))
plot(ecdf(Gas_flux_R$pH), do.points=FALSE, verticals=TRUE)

#the distrubution is far from normal
#therefore, check the distrubtion for high and low pH
High_pH <- Gas_flux_R$pH[Gas_flux_R$pH > 9.1]
Low_pH <- Gas_flux_R$pH[Gas_flux_R$pH <= 9.1]
#by splitting up the data in high and low pH the a normal distrubution can be aquired if that is needed for some tests.  
#The Shapiro test, however, still detects abnormallity which is not uncommen in large data-sets since the smallest diviation is detected. 
#The histogram and the qq-plot both show normal distrubtion
shapiro.test(High_pH)
hist(High_pH)
qqnorm(High_pH);qqline(High_pH)

shapiro.test(Low_pH)
hist(Low_pH)
qqnorm(Low_pH);qqline(Low_pH)

Gas_flux_R$pH-log <-
hist(Gas_flux_R$Wind.speed.m.s)
hist(Gas_flux_R$I.average)
hist(Gas_flux_R$Chl.fluorescence.ug.L)
hist(Gas_flux_R$Chl.Fluorescence.Volt)
hist(Gas_flux_R$Chl.Lab.ug.L)
hist(Gas_flux_R$X.CO2....405ppm)
hist(Gas_flux_R$CO2.water.mol.L)
hist(Gas_flux_R$Flux.umol.m2.h.1)



#####################################################################################
## RDA or CCA

# The data overall shows linear relations, plus the data are counts or measures. A CCA is based on unimodel data and on data relative data, such as species composition.
# Therefore, an RDA might be more appropiate for this data set. 
names(Gas_flux_R)
rda_co2_flux <- rda(Gas_flux_R$CO2.water.mol.L + Gas_flux_R$Flux.umol.m2.h.1~Gas_flux_R$T..oC. + Gas_flux_R$pH + Gas_flux_R$Wind.speed.m.s + Gas_flux_R$I.average + Gas_flux_R$Chl.Lab.ug.L + Gas_flux_R$Chl.Fluorescence.Volt)
