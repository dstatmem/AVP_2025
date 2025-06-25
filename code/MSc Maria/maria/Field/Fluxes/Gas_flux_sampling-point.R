setwd("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/")
Gas_flux_SP <- read.csv("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Fluxes/Stat/Gas_flux_SamplingPoint.csv")
View(Gas_flux_SP)
library(plyr)
library(car)
library(vegan)
library(lmtest)


plot(Gas_flux_SP)

names(Gas_flux_SP)
Gas_flux_SP <- data.frame(Gas_flux_SP[c(1,2,4:11),])


##########################################################################
# Testing the assumptions for linear regression

##linearity
plot(Gas_flux_SP)
plot(log10(Gas_flux_SP$pH), log10(Gas_flux_SP$CO2))

Gas_flux_SP$pH_log <- log10(Gas_flux_SP$pH)
Gas_flux_SP$CO2_log <- log10(Gas_flux_SP$CO2)
Gas_flux_SP$HCO3_log <- log10(Gas_flux_SP$HCO3)
Gas_flux_SP$CO3_log <- log10(Gas_flux_SP$CO3)

plot(Gas_flux_SP$HCO3, Gas_flux_SP$CO2)
plot(Gas_flux_SP$HCO3, Gas_flux_SP$CO3)
plot(Gas_flux_SP$CO2, Gas_flux_SP$CO3)

## Colinearity
lm_all1 <- lm(Gas_flux_SP$CO2Flux ~ Gas_flux_SP$NO2 + Gas_flux_SP$NO2NO3 + Gas_flux_SP$NH4 + Gas_flux_SP$PO4 + Gas_flux_SP$I + Gas_flux_SP$pH_log + Gas_flux_SP$Temp + Gas_flux_SP$Wind )
lm_all2 <- lm(Gas_flux_SP$CO2Flux ~  Gas_flux_SP$ChlLab + Gas_flux_SP$I + Gas_flux_SP$pH_log + Gas_flux_SP$Temp + Gas_flux_SP$Wind)
lm_all3 <- lm(Gas_flux_SP$CO2Flux ~  Gas_flux_SP$I + Gas_flux_SP$pH_log + Gas_flux_SP$Temp + Gas_flux_SP$Wind)
lm_all4 <- lm(Gas_flux_SP$CO2Flux ~  Gas_flux_SP$ChlFluo + Gas_flux_SP$I + Gas_flux_SP$Temp + Gas_flux_SP$Wind)

# For the fluxes build a model with Chl and pH seperatly, where there is corrected for pH and Chlorrophyll

bwplot(Gas_flux_SP$CO2Flux ~ Gas_flux_SP$ChlFluo * Gas_flux_SP$pH_log)

vif(lm_all1)
vif(lm_all2)
vif(lm_all3)
vif(lm_all4)
vif(lm_all5)





dwtest(lm_all4)


#ks.test(x, "pnorm", mean=1, sd=2)
ks.test(Gas_flux_SP$CO2Flux,Gas_flux_SP$ChlLab )


## homoscedasticity



## normal distribution of the error




###############################################################################################################
## models

lm_Flux_Chl <- lm(Gas_flux_SP$CO2Flux ~  Gas_flux_SP$ChlFluo * Gas_flux_SP$pH_log + Gas_flux_SP$I + Gas_flux_SP$Temp + Gas_flux_SP$Wind)
summary(lm_Flux_Chl)



######################################################################################################################


##############################################
#Just some orientation stuff

Gas <- data.frame(Gas_flux_SP$CO2, Gas_flux_SP$HCO3, Gas_flux_SP$CO3)
Nutrients <- data.frame(Gas_flux_SP[,c( 15:18)])
Chl <- data.frame(Gas_flux_SP[,c(10:11)])
Env <- data.frame(Gas_flux_SP[,c(2,8,9)])

rdaNut <- rda(Gas ,Nutrients)
rdaChl <- rda(Gas ,Chl)
rdaEnv <- rda(Gas ,Env)

rdaChl_Nut <- rda(Chl, Env)
anova(rdaChl_Nut)

anova(rdaNutrients)
anova(rdaChl)
anova(rdaEnv)


plot_var1 <- varpart(Gas_flux_SP$CO2 + Gas_flux_SP$HCO3 + Gas_flux_SP$CO3, ~  Gas_flux_SP$ChlFluo, scale = TRUE)


rdaGasNut <- rda(Gas_flux_SP$CO2 + Gas_flux_SP$HCO3 + Gas_flux_SP$CO3 ~ Gas_flux_SP$NO2NO3 + Gas_flux_SP$NO2 + Gas_flux_SP$NH4 + Gas_flux_SP$PO4, scale = TRUE)
rdaGasNut_Chl <- rda(Gas_flux_SP$CO2 + Gas_flux_SP$HCO3 + Gas_flux_SP$CO3 ~ Gas_flux_SP$NO2NO3 + Gas_flux_SP$NO2 + Gas_flux_SP$NH4 + Gas_flux_SP$PO4 + Condition(Gas_flux_SP$ChlLab + Gas_flux_SP$ChlFluo), scale = TRUE)
anova(rdaGasNut_Chl)
summary(rdaGasNut_Chl)
envfit(rdaGasNut_Chl ~ Gas_flux_SP$NO2NO3 + Gas_flux_SP$NO2 + Gas_flux_SP$NH4 + Gas_flux_SP$PO4, scale = TRUE)
# nutrients without condition for Chlorophyll are significant. 


rdaGasChl <- rda(Gas_flux_SP$CO2 + Gas_flux_SP$HCO3 + Gas_flux_SP$CO3 ~ Gas_flux_SP$ChlLab + Gas_flux_SP$ChlFluo, scale = TRUE)
rdaGasChl_Nut <- rda(Gas_flux_SP$CO2 + Gas_flux_SP$HCO3 + Gas_flux_SP$CO3 ~  Gas_flux_SP$ChlLab + Gas_flux_SP$ChlFluo + Condition(Gas_flux_SP$NO2NO3 + Gas_flux_SP$NO2 + Gas_flux_SP$NH4 + Gas_flux_SP$PO4 + Gas_flux_SP$pH), scale = TRUE)
anova(rdaGasChl_Nut)  #without ChlLab the model becomes significant. 
envfit(rdaGasChl_Nut ~ Gas_flux_SP$ChlLab + Gas_flux_SP$ChlFluo)


rdaGasEnv <- rda(Gas_flux_SP$CO2 + Gas_flux_SP$HCO3 + Gas_flux_SP$CO3 ~ Gas_flux_SP$Wind + Gas_flux_SP$Temp + Gas_flux_SP$I + Gas_flux_SP$pH)
anova(rdaGasEnv)
envfit(rdaGasEnv ~ Gas_flux_SP$Wind + Gas_flux_SP$Temp + Gas_flux_SP$I + Gas_flux_SP$pH)
ordistep(rdaGasEnv)

# Wind explains little variance in Chl
rdaChlEnv <- rda(Gas_flux_SP$ChlLab + Gas_flux_SP$ChlFluo ~ Gas_flux_SP$Temp + Gas_flux_SP$I + Gas_flux_SP$Wind)
anova(rdaChlEnv)
summary(rdaChlEnv)
envfit(rdaChlEnv~  Gas_flux_SP$Temp + Gas_flux_SP$I + Gas_flux_SP$Wind)
Var3 <- varpart(Gas_flux_SP$ChlLab + Gas_flux_SP$ChlFluo, ~ Gas_flux_SP$Temp,~ Gas_flux_SP$I,~ Gas_flux_SP$Wind)
plot(Var3)

#there is a difference in how much variance is explained by nutrient when looking at the fluorescence or the lab data
# Wind explains little variance in Chl
rdaChlNut <- rda(Gas_flux_SP$ChlFluo + Gas_flux_SP$ChlLab  ~ Gas_flux_SP$NO2 , scale= TRUE)
anova(rdaChlNut)
summary(rdaChlNut)
envfit(rdaChlNut~  Gas_flux_SP$NO2NO3 + Gas_flux_SP$NO2 + Gas_flux_SP$NH4 + Gas_flux_SP$PO4)
ordistep(rdaChlNut)

plot_var2 <- varpart( Gas_flux_SP$ChlFluo + Gas_flux_SP$ChlLab ,~Gas_flux_SP$NO2,~ Gas_flux_SP$NH4, scale=TRUE)
plot(plot_var2)
