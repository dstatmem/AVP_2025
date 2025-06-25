# respiration * Temperature
# from the field and the lab

library(car)
library(ggplot2)
library (vegan)
library (vegan3d)
library(devtools)
library(ggbiplot)
library(multcomp)
library(nlme)
library(pastecs)
library(lattice)
library(plyr)
library(lmtest)
library(Rmisc)

setwd("D:/UvA/Master/IBED stage/Jolanda/Stat/PI_curve")
PI_fit_R <- read.csv("D:/UvA/Master/IBED stage/Jolanda/Stat/PI_curve/coef_PI_model_R.csv")
Raw_resp_R <- read.csv("D:/UvA/Master/IBED stage/Jolanda/Stat/PI_curve/Raw_Resp_R.csv")

#View(PI_fit_R)
names(PI_fit_R)


##Check first for for outliers in the raw respiration data
#by fitting the expected value against the data

# check first the residuals against temperature
Raw_resp_R$Residuals <- (Raw_resp_R$Resp - Raw_resp_R$Fit)
hist(Raw_resp_R$Residuals)
#There is a normal distribution of the residuals
plot(Raw_resp_R$Residuals ~ Raw_resp_R$Temp)
#The residuals become larger with larger temperatures


#check for heteroskedasticity against temperature
bptest(Resp ~ Temp, data = Raw_resp_R)
#which is not significant but strong present




#check for the model fit, outlier in de model, and for heteroskedasticity

Raw_lm <- lm(Raw_resp_R$Resp*-1 ~ Raw_resp_R$Fit)
summary(Raw_lm)
#the model fits significantly but with a low R-squared, this becomes higher when the data is averaged.

#check for outliers
plot(Raw_lm)
Raw_cook <-cooks.distance(Raw_lm)
max(Raw_cook)
#no cook's distance is larger than 0.5, infact there small. So there is no worry that the outlier found by eye have a significant impact on the model









#####################################
#check the data which is averaged over the channels

plot(PI_fit_R$Resp ~ PI_fit_R$Temp)




#Include a positive respiaration vlaue in the data
PI_fit_R$Resp_neg <- PI_fit_R$Resp*-1

Resp <- Raw_resp_R$Resp*-1
Temp <- Raw_resp_R$Temp

r_fit <-(coef(nls(Resp ~  exp(const)* exp(Ea/ (0.000086173303 * (Temp + 273.15))),  start=list(Ea=-0.3986, const=13.437))))

# and the fits, three fits were included, that of the averaged data over the channels, that of the averaged data over the channels and days, and the R fit of the data averaged over the raw data
# the R file has been also checked for the averaged data but there is barely a difference.
PI_fit_R$fit_all <- exp(13.437)* exp(-0.3986/ (0.000086173303 * (PI_fit_R$Temp + 273.15)))
PI_fit_R$fit_av <- exp(13.928)* exp(-0.4106/ (0.000086173303 * (PI_fit_R$Temp + 273.15)))
PI_fit_R$fit_r <- exp(11.6166792)* exp(-0.3511325/ (0.000086173303 * (PI_fit_R$Temp + 273.15))) #averaged data
PI_fit_R$fit_r <- exp(11.9670563)* exp(-0.3596382/ (0.000086173303 * (PI_fit_R$Temp + 273.15))) #raw data

cor(PI_fit_R$Resp_neg, PI_fit_R$fit_av)

plot(PI_fit_R$Resp_neg ~ PI_fit_R$fit_av)
plot(PI_fit_R$Resp_neg ~ PI_fit_R$fit_r)

#Calcualte the residuals between the data and the fitted data
PI_fit_R$Residuals_all <- (PI_fit_R$Resp_neg - PI_fit_R$fit_all)
PI_fit_R$Residuals_av <- (PI_fit_R$Resp_neg - PI_fit_R$fit_av)
PI_fit_R$Residuals_r <- (PI_fit_R$Resp_neg - PI_fit_R$fit_r)

#plot to see if the distribution of the residuals shows a normal pattern
plot(PI_fit_R$Residuals_all ~ PI_fit_R$Temp)
plot(PI_fit_R$Residuals_av ~ PI_fit_R$Temp)
plot(PI_fit_R$Residuals_r ~ PI_fit_R$Temp)
# once the data i averaged over the channels there is no sign of heteroskedasticity over temperature!!!!!!!!!!



#chech the significance of the fit to the data, and chech the R-squared.  

# if you want to compare more models, to check which one best fits, an anova tells if there is -first significance difference between the models, and -second which model has lower std. 


lm_fit_all <- lm(Resp_neg ~ fit_all, data = PI_fit_R) 
summary(lm_fit_all)
plot(lm_fit_all)


lm_fit_av <- lm(Resp_neg ~ fit_av, data = PI_fit_R) 
summary(lm_fit_av)
plot(lm_fit_av)

lm_fit_r <- lm(Resp_neg ~ fit_r, data = PI_fit_R) 
summary(lm_fit_r)
plot(lm_fit_r)


anova(lm_fit_av, lm_fit_all, lm_fit_r)
#It doesn't matter which fit is used, the fit of the per date averaged data or the fit of the -only- per channel averaged data



#######################################################################
##Check data averaged over channels and date

#Create a new dataframe with the average values of Fit and resperation
PI_mean <- ddply(PI_fit_R, ~Date, summarise, Resp=mean(Resp_neg), Fit_all=mean(fit_all), Fit_av=mean(fit_av), Fit_r=mean(fit_r), Temp=mean(Temp))


#check the linear regresion of the of the for the average values
#lead the intercept trhought zero, since this is assumed in the model. 
lm_fit_mean_all <- lm(Resp ~ Fit_all, data = PI_mean) 
plot(PI_mean$Resp ~ PI_mean$Fit_all)
summary(lm_fit_mean_all)
plot(lm_fit_mean)
hist(lm_fit_mean$residuals)
bptest(Resp ~ Fit_all, data = PI_mean)
se <- sqrt(diag(vcov(lm_fit_mean_all)))

lm_fit_mean_av <- lm(Resp ~ Fit_av, data = PI_mean) 
summary(lm_fit_mean_av)


lm_fit_mean_r <- lm(Resp ~ Fit_r, data = PI_mean) 
summary(lm_fit_mean_r)


#Also check for the 
anova(lm_fit_mean_av, lm_fit_mean_all, lm_fit_mean_r)

#Chech the resididuals of the average values 
PI_mean$residuals <- (PI_mean$Resp - PI_mean$Fit_all)

plot(PI_mean$residuals ~ PI_mean$Temp)


#check for the effect of temperature on respiration, where there is corrected for the date (repeated measurements)
lme_fit <- lme(log10(Resp) ~ Temp, random= ~1|Date, data = PI_mean, method = "ML") 
plot(log10(PI_mean$Resp) ~ Temp)
summary(lme_fit)
plot(lme_fit)
bptest(log(PI_mean$Resp) ~ PI_mean$Temp)

## assing the fit 
Maximum_metabolic_rate = function(T) {
    exp(13.437)* exp(-0.3986/ (0.000086173303 * (T + 273.15)))
}





fit_line_Resp_temp <- Maximum_metabolic_rate(Temp)



#################!!!!!!!!!!!! check the figure !!!!!!!!!!!!!!!!!!!!####################
# plot the model fit

PI_mean$SE_fit_all_pos <- c(PI_mean$Fit_all- 2*0.27103593)
PI_mean$SE_fit_all_neg <- c(PI_mean$Fit_all + 2*0.27103593)

Data_fit_all <- ggplot(data = PI_mean) + 
  geom_point(aes(Temp, Resp)) + geom_smooth(aes(Temp, Fit_all), color = "darkblue") +
  theme(legend.direction = 'horizontal', legend.position = 'top', 
        panel.grid.major = element_line(color = "white"), panel.grid.minor = element_blank(),
        panel.background = element_blank(),  axis.line = element_line(color = "black")) +
        ylab("Respiration umol / ug-Chl / h") + xlab("Temperature C") 
Data_fit_all

Data_fit_all <-  Data_fit_all + geom_line(aes(Temp, SE_fit_all_pos), color = "lightblue")
Data_fit_all <-  Data_fit_all + geom_line(aes(Temp, SE_fit_all_neg), color = "lightblue")
Data_fit_all



geom_smooth(aes(Temp, SE_fit_all_pos), color = "lightblue") + 
  geom_smooth(aes(Temp, SE_fit_all_pos), color = "lightblue") +

summary(lm_fit_mean_all)


Data_fit_av <- ggplot(data = PI_mean) + 
  geom_point(aes(Temp, Resp)) + geom_smooth(aes(Temp, Fit_av), color = "darkblue") +
  theme(legend.direction = 'horizontal', legend.position = 'top', 
        panel.grid.major = element_line(color = "white"), panel.grid.minor = element_blank(),
        panel.background = element_blank(),  axis.line = element_line(color = "black")) +
        ylim(0.01,0.2) + ylab("Respiration umol / ug-Chl / h") + xlab("Temperature C") 
Data_fit_av

summary(lm_fit_mean_av)

#Add additional temperature for the exention of the fit
#for (i in 5:0) {
#PI_mean <- rbind(PI_mean, c(0,0,0,i))
#}

#PI_mean <- rbind(PI_mean, c("0","0","0",i))


#Add new fit parameters
#PI_mean$Fit <- exp(13.437)* exp(-0.3986/ (0.000086173303 * (PI_mean$Temp + 273.15)))


#Replot the plot with extended temperature
#Data_fit <- ggplot(data = PI_mean) + 
  #geom_point(aes(Temp, Resp)) + geom_smooth(aes(Temp, Fit), color = "darkblue") +
  #theme(legend.direction = 'vertical', legend.position = 'left', 
        #panel.grid.major = element_line(color = "white"), panel.grid.minor = element_blank(),
        #panel.background = element_blank(),  axis.line = element_line(color = "black")) +
        #ylim(0.01,0.2) + ylab("Respiration umol/L/h") + xlab("Temperature C") 
#Data_fit

#lm_fit <- predict(lm_fit_mean,se.fit=TRUE)
#plot(lm_fit$fit)

##graph info
# fit : exp(13.437)* exp(-0.3986/ (0.000086173303 * (Temp + 273.15)))
# R-square = 0.3818
# P-value < 0.01
# se residuals = 0.0291


### Check to do this from one data.frame PImean
#remove data.frame fit 
#