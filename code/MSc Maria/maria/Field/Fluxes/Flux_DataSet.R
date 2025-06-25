library(car)
library(vegan)
library(nlme)
install.packages("plyr")
library(plyr)
setwd("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes/")


CO2_flux_data <- read.csv("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/")
  

names(CO2_flux_data)
daily_flux_day <- ddply(Gas_flux_H, .(Date), summarize, Temp = mean(T.C), pH = mean(pH), CO2.Water = mean(CO2.water.umol.L), 
                   I = mean(I.in), Chl.Fluo.volt = mean(CHL.average), Chl.volt = mean(CHLV.Volts), Wind = mean(Wind.speed.m.s), 
                   Flux = mean(Flux.umol.m2.h.1), HCO3 = mean(HCO3), CO3 = mean(CO3))
