library(readr)

# Calcium precipitation

CA_prec <- read_delim("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/CA_prec.csv", 
                           ";", escape_double = FALSE, trim_ws = TRUE)

par(mfrow = c(2,2))
plot(lm(CA_prec$Ca~CA_prec$`Alk ueq`))
par(mfrow = c(1,1))
summary(lm(CA_prec$Ca~CA_prec$`Alk ueq`))
summary(lm(CA_prec$`Alk ueq`~CA_prec$Ca))

par(mar = mar.default+1)

CA_prec$Ca_range <- seq(length=27, from = 600, to = 1700)

{
  
plot(CA_prec$`Alk ueq` ~ CA_prec$Ca, ylab = "Alkalinity ["~mu*" equivalent L"^"-1"*"]", xlab = "Calcium ["~mu*"mol L" ^"-1"*"]",  pch = 16, col = "darkblue")
lines(127.22657 + 2.17576 * CA_prec$Ca_range ~ CA_prec$Ca_range , col = "red")
lines(127.22657 + 80.77017*1.96 + (2.17576 + 0.1619156)* CA_prec$Ca_range ~ CA_prec$Ca_range , col = "red", lty = 3)
lines(127.22657 - 80.77017*1.96 + (2.17576 - 0.1619156)* CA_prec$Ca_range ~ CA_prec$Ca_range , col = "red", lty = 3)
legend("topleft", 
       inset = c(0.02,0.02), 
       cex = 1.1, 
       bty = "n", 
       legend = c(expression("Alkalinity "~mu*"eq"), "Regression line", "95% confidence interval"), 
       text.col = c("black", "black", "black"),
       col = c("darkblue", "red", "red"), 
       lwd=1, lty=c(NA,1,3), pch=c(16,NA,NA))
}



# active mixng layer
Zaml <- read_delim("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Zaml.csv", 
                        ";", escape_double = FALSE, trim_ws = TRUE)
par(mfrow = c(2,2))
plot(lm(Zaml$`Average Depth` ~ log(Zaml$`Average Temp`)))
par(mfrow = c(1,1))
summary(lm(log(Zaml$`Average Depth`) ~ Zaml$`Average Temp`))

Temp_range <- seq(steps = 24, from = min(Zaml$`Average Temp`)-5, to = max(Zaml$`Average Temp`)+5)

plot(Zaml$`Average Depth` ~ Zaml$`Average Temp`)
lines(exp(0.69741) * exp(-0.11611 * Zaml$Temp_range) ~ Zaml$Temp_range)


plot(Zaml$`Average Depth` ~ Zaml$`Average Temp`, ylab = "Actively mixing layer depth [m]", xlab = "Temperature [°C]", 
     cex.axis = 1.3, cex.lab = 1.3, pch = 16, col = "darkblue", ylim = c(0,1.8), xlim = c(4,25))
lines(exp(0.69741) * exp(-0.11611 * Temp_range) ~ Temp_range, col = "red")
lines(exp(0.69741 + 0.36711 * 1.96) * exp((-0.11611 + 0.02079*1.96) * Temp_range) ~ Temp_range, col = "red", lty = 3)
lines(exp(0.69741 - 0.36711 * 1.96) * exp((-0.11611 - 0.02079*1.96) * Temp_range) ~ Temp_range, col = "red", lty = 3)
legend("topright", 
       cex = 1.1, 
       bty = "n", 
       legend = c(expression("Z"[aml]), "Regression line", "95% confidence interval"), 
       text.col = c("black", "black", "black"),
       col = c("darkblue", "red", "red"), 
       lwd=1, lty=c(NA,1,3), pch=c(16,NA,NA))
     


# Extinction coefficient

Extinction_coef <- read_delim("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Extinction_coef.csv", 
                                   ";", escape_double = FALSE, trim_ws = TRUE)


par(mfrow = c(2,2))
plot(lm(Extinction_coef$Extinction_coef_m1 ~ Extinction_coef$Chl_ugL))
par(mfrow = c(1,1))
summary(lm(Extinction_coef$Extinction_coef_m1[-c(19,21)] ~ Extinction_coef$Chl_ugL[-c(19,21)]))
summary(lm(Extinction_coef$Extinction_coef_m1 ~ Extinction_coef$Chl_ugL))

Extinction_coef$Chl_range <- seq(length = nrow(Extinction_coef), from = min(Extinction_coef$Chl_ugL), to = max(Extinction_coef$Chl_ugL)+50)



plot(Extinction_coef$Extinction_coef_m1 ~ Extinction_coef$Chl_ugL, ylab = expression("Extinction coefficient [m"^-1*"]"), xlab = expression(paste("Chlorophyll ", italic(a), "["~mu~"g L"^"-1"*"]")),  pch = 16, col = "darkblue", ylim = c(1,12))
lines(1.385828 + 0.022218 * Extinction_coef$Chl_range ~ Extinction_coef$Chl_range, col = "red")
lines(1.385828 + 0.307098 + (0.022218 + 0.002304*1.96) * Extinction_coef$Chl_range ~ Extinction_coef$Chl_range, col = "red", lty = 3)
lines(1.385828 - 0.307098 + (0.022218 - 0.002304*1.96) * Extinction_coef$Chl_range ~ Extinction_coef$Chl_range, col = "red", lty = 3)
legend("topleft", 
       inset = c(0.0,0.0), 
       cex = 1.1, 
       bty = "n", 
       legend = c(expression("K"[tot]), "Regression line", "95% confidence interval"), 
       text.col = c("black", "black", "black"),
       col = c("darkblue", "red", "red"), 
       lwd=1, lty=c(NA,1,3), pch=c(16,NA,NA))




Profile_AVP <- read_delim("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Profile_AVP.csv", 
                              ";", escape_double = FALSE, col_types = cols(Date = col_date(format = "%m-%d-%Y")), 
                              trim_ws = TRUE)


Profile_AVP <- Profile_AVP[complete.cases(Profile_AVP[c("Temp", "DEP200 meters")]),]

Profile_AVP <- Profile_AVP[Profile_AVP$Date > "2016-05-17",]

length(unique(Profile_AVP$Nr))
summary(Profile_AVP$Temp)


layout(matrix(1:28, nrow = 4, byrow = TRUE))
par(mar = c(4,3,2,0.5))
for (i in unique(Profile_AVP$Date)){
  df <- subset(Profile_AVP[Profile_AVP$`DEP200 meters` > 0.05 & Profile_AVP$`DEP200 meters` < 1.5,], Date == i)
  plot(df$Temp, df$`DEP200 meters`, ylim = c(2,0), xlim = c(5,26), main = unique(df$Date, cex = 2),
        ylab = "", xlab = "", pch = 16, col = "darkblue", cex.lab = 1.2, cex.axis = 1.2)
  lines(df$Temp, df$`DEP200 meters`,  lty = 3)

} 

summary(df$Dommol_L)

layout(matrix(1:28, nrow = 4, byrow = TRUE))
par(mar = c(4,3,2,0.5))
for (i in unique(Profile_AVP$Date)){
  df <- subset(Profile_AVP[Profile_AVP$`DEP200 meters` > 0.05 & Profile_AVP$`DEP200 meters` < 1.5,], Date == i)
  plot(df$Dommol_L, df$`DEP200 meters`, ylim = c(2,0), xlim = c(0,1), main = unique(df$Date),
       ylab = "", xlab = "", pch = 16, col = "darkblue", cex.lab = 1.2, cex.axis = 1.2)
  lines(df$Dommol_L, df$`DEP200 meters`,  lty = 3)
} 


library(dplyr)
Profile_AVP[Profile_AVP$`DEP200 meters` > 0.6,] %>%
  group_by(Date) %>%
  summarise(diff = max(Temp)-min(Temp)) -> Temp_diff

Temp_diff
summary(Temp_diff)
sd(Temp_diff$diff)/length(sqrt(Temp_diff$diff))

library(dplyr)
Profile_AVP[Profile_AVP$`DEP200 meters` < 0.3 & Profile_AVP$`DEP200 meters` > 0.05,] %>%
  group_by(Date) %>%
  summarise(max_sur = max(Dommol_L),
            min_sur = min(Dommol_L),
            mean_sur = mean(Dommol_L)) -> DO_sur_diff
summary(DO_sur_diff)
sd(DO_sur_diff$mean_sur)/length(sqrt(DO_sur_diff$mean_sur))


library(dplyr)
Profile_AVP[Profile_AVP$`DEP200 meters` > 1.2 & Profile_AVP$`DEP200 meters` < 1.5,] %>%
  group_by(Date) %>%
  summarise(max_bent = max(Dommol_L),
            min_bent = min(Dommol_L),
            mean_bent = mean(Dommol_L)) -> DO_bent_diff
summary(DO_bent_diff)

