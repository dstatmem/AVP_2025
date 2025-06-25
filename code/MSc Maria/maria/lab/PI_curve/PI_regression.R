Result_PI_curve_lab <- read_csv("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/Result_PI_curve_lab.csv")







summary(Result_PI_curve_lab$P_max)
plot(Result_PI_curve_lab$P_max)
summary(Result_PI_curve_lab$h_i)
plot(Result_PI_curve_lab$h_i)
summary(Result_PI_curve_lab$Resp)
plot(Result_PI_curve_lab$Resp)



# Create a bootstrap function 
bootReg <- function(formula, data, indices) {
  d <- data[indices,]
  fit <- lm(formula, data = d)
  return(coef(fit))
}

library(boot)


Temp_range <- seq(steps = nrow(Result_PI_curve_lab),from = min(Result_PI_curve_lab$Temp)-1, to = max(Result_PI_curve_lab$Temp)+1 )


#######################################
# P_max
#######################################

# There does not seem to be a clear trend in P_max and temp



summary(Result_PI_curve_lab$P_max)



par(mfrow = c(2,2))
plot(lm(log(Result_PI_curve_lab$P_max) ~ Result_PI_curve_lab$Temp))
par(mfrow = c(1,1))
summary(lm(log(Result_PI_curve_lab$P_max) ~ Result_PI_curve_lab$Temp))


Pmax_lab_boot <- boot(statistic = bootReg, formula = log(P_max) ~ Temp, data = Result_PI_curve_lab, R = 2000)
boot.ci(Pmax_lab_boot, type = "bca", index = 1)
boot.ci(Pmax_lab_boot, type = "bca", index = 2)



plot(Result_PI_curve_lab$Temp, Result_PI_curve_lab$P_max, col = "orange", ylab = "Maximum photosynthesis [umol/L/h]", xlab = "Temperature [C]")
lines(Temp_range, exp(-1.39071) * exp(0.06772112 * Temp_range), col = "blue")                                                    # Based on I60
lines(Temp_range, exp(-1.743) * exp(0.0312 * Temp_range), col = "blue", lty = 3)  
lines(Temp_range, exp(-0.701) * exp(0.0955 * Temp_range), col = "blue", lty = 3)




####################################################33
# REspiration
######################################################



# Calculate the natural logarithm of r and the inverse of temp
Result_PI_curve_lab$r_ln <- log(Result_PI_curve_lab$Resp*-1)
Result_PI_curve_lab$invK <- 1/((Result_PI_curve_lab$Temp+273.15)*0.000086173303)

par(mfrow = c(2,2))
plot(lm(Result_PI_curve_lab$r_ln ~ Result_PI_curve_lab$invK))
par(mfrow = c(1,1))

summary(lm(Result_PI_curve_lab$r_ln ~ Result_PI_curve_lab$invK))



r_lab_boot <- boot(statistic = bootReg, formula = r_ln ~ invK, data = Result_PI_curve_lab, R = 2000)
boot.ci(r_lab_boot, type = "bca", index = 1)
boot.ci(r_lab_boot, type = "bca", index = 2)



plot(Result_PI_curve_lab$Temp, Result_PI_curve_lab$Resp*-1, col = "orange", ylab = "Respiration [umol/L/h]", xlab = "Temperature [C]", ylim = c(0,0.5))
lines(Temp_range, exp(15.15017)* exp(-0.44208/ (0.000086173303 * (Temp_range + 273.15))), col = "blue")                                                    # Based on I60
lines(Temp_range, exp(7.73)* exp(-0.5783/ (0.000086173303 * (Temp_range + 273.15))), col = "blue", lty = 3)  
lines(Temp_range, exp(20.63) * exp(-0.2546/ (0.000086173303 * (Temp_range + 273.15))), col = "blue", lty = 3)




##################################################
# h_i against chl votl
##################################################




par(mfrow = c(2,2))
plot(lm(Result_PI_curve_lab$h_i ~ Result_PI_curve_lab$Temp))                # This models shows non-linearity, however it fits in reality better than the inverse Chl model
par(mfrow = c(1,1))


summary(lm(Result_PI_curve_lab$h_i ~ Result_PI_curve_lab$Temp))


mean(Result_PI_curve_lab$h_i, na.rm = TRUE)
sd(Result_PI_curve_lab$h_i, na.rm = TRUE)/sqrt(nrow(Result_PI_curve_lab))
