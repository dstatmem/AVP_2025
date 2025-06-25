#check for packages https://rdrr.io/all/cran/
library(AquaEnv)
library(LakeMetabolizer)
library(ggplot2)
library(readr)


#Chewbacca
DF_data <- read_delim("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/Merged_files_R.csv", 
                      ",", escape_double = FALSE, col_types = cols(`Average Sampling time` = col_datetime(format = "%d/%m/%Y %H:%M")), 
                      trim_ws = TRUE)

#Przemek
DF_data <- read_delim("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Merged_files_R.csv", 
                      ";", escape_double = FALSE, col_types = cols(`Average Sampling time` = col_datetime(format = "%d-%m-%Y %H:%M")), 
                      trim_ws = TRUE)
setwd("C:/Users/admin/Dropbox/MSc Maria/Maria/Lab/PI_curve/")




# Add date
DF_data$Date <- format(DF_data$`Average Sampling time`, "%Y-%m-%d")

#Time difference in hours
DF_data$Time/3600


#calculate the change in DO per hour
# The linear regression (LINEST in excel) hass been used over 5 points
DF_data$`d_DO_umol/L` <- DF_data$`d_DO_umol/L`/4                          # the average difference per time step 

# Convert to change per hour
ind_j <- 1:(nrow(DF_data)-1)
ind_k <- 2:nrow(DF_data)
ind_w <- 1:nrow(DF_data)


for (i in seq_along(ind_j)){
  DF_data$TimeDif[ind_w[i]] <- 60/(as.double(DF_data$`Average Sampling time`[ind_k[i]]-DF_data$`Average Sampling time`[ind_j[i]]))
}

DF_data$`d_DO_umol/L` <- DF_data$`d_DO_umol/L`*DF_data$TimeDif



plot(DF_data$`Average Sampling time`, DF_data$`d_DO_umol/L`)
plot(DF_data$`I (umol/m2/s)`[DF_data$Date == "2016-05-27"], DF_data$`d_DO_umol/L`[DF_data$Date == "2016-05-27"])




# create a new dataframe for the PI_curve
PI_field <- DF_data[, c("Average Sampling time", "Time", "Date", "d_DO_umol/L", "I (umol/m2/s)", "Chl_lab", "CHL aver mg/L", "T_C_water")]
PI_field <- PI_field[complete.cases(PI_field[,c("d_DO_umol/L","I (umol/m2/s)","CHL aver mg/L")]),]


# Calculate the DO change per ug chlorophyll
PI_field$d_DO_Chl <- PI_field$`d_DO_umol/L`/PI_field$`CHL aver mg/L`

plot(PI_field$`Average Sampling time`, PI_field$d_DO_Chl)
plot(DF_data$`I (umol/m2/s)`[DF_data$Date == "2016-05-27"], PI_field$d_DO_Chl[DF_data$Date == "2016-05-27"])
plot(PI_field$`Average Sampling time`, PI_field$`d_DO_umol/L`)

PI_field$Days <- as.Date(PI_field$Date) - as.Date(PI_field$Date[1])+1
PI_field$Days <- as.integer(PI_field$Days)
max(PI_field$Days)

PI_field_list <- list()
for (i in 1:max(PI_field$Days)) {
  PI_field_list[[i]] <- PI_field[which(PI_field$Days == i), ]
}

PI_field_list <- PI_field_list[lapply(PI_field_list,length)>0]




#######################################
#             PI curves
#######################################


## assign the function 
Michaelis_Menten = function(par, I) {
  par[1] + (par[2]* I) / (par[3] + I)
}


# create new lists
d_DO <- list()
for (i in 1:max(PI_field$Days)) {
  d_DO[[i]] <- PI_field$d_DO_Chl[PI_field$Days == i]
}

I <- list()
for (i in 1:max(PI_field$Days)) {
  I[[i]] <- PI_field$`I (umol/m2/s)`[PI_field$Days == i]
}



# Fit the non-linear model
PI_curve_coef <- list()

#I_v <- I[[179]]
#coef(nls(d_DO[[179]] ~ r + (P_max * I_v) / (h_i + I_v),  start=list(r=-0.4, P_max=0.3,h_i=400)))

for (i in 1:max(PI_field$Days)){
  tryCatch({
    I_v <- I[[i]]
    PI_curve_coef[[i]] <- coef(nls(d_DO[[i]] ~ r + (P_max * I_v) / (h_i + I_v),  start=list(r=-0.4, P_max=0.3,h_i=400)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

PI_curve_coef


PI_field_coef <-  data.frame(Reduce(rbind,PI_curve_coef))

plot(PI_field_coef$r)

write.table(PI_field_coef,"C:/Users/admin/Dropbox/MSc Maria/Maria/Lab/PI_curve/PI_field_coef_ChlField.txt", sep = ";")



## check for the daily minimum temperature, to link the temp to respiration 

Temp <- by(PI_field$T_C_water, PI_field$Days, function(x) min(x))



# output the fit from the model

fit_PI_field <- lapply(PI_curve_coef, function(x) Michaelis_Menten(x,I1))
lapply(fit_PI_field[[8]], function(x) plot(x ~ I[[8]],type="l"))
plot(fit_PI_field[[8]], I[[8]])
write.table(fit_PI_field,"fit_PI_field.txt", sep = "\t", row.names = FALSE)


PI_field_coef <-  data.frame(Reduce(rbind,PI_curve_coef))





#########################################################################################################
#
#                         Correct for the light 
#
#########################################################################################################



# Calculate the light at 60 cm depth
PI_field$I_in_water <- PI_field$`I (umol/m2/s)` * exp(-1 * (0.0222177372 * PI_field$Chl_lab + 1.3858280092) * 0.6)


plot(PI_field$`Average Sampling time`, PI_field$I_in_water)

plot(PI_field$I_in_water, PI_field$d_DO_Chl, ylim = c(-0.2,0.2))
plot(PI_field$I_in_water[DF_data$Date == "2016-05-27"], PI_field$d_DO_Chl[DF_data$Date == "2016-05-27"])
plot(PI_field$`Average Sampling time`, PI_field$`d_DO_umol/L`)



#######################################
#             PI curves
#######################################


## assign the function 
Michaelis_Menten = function(par, I) {
  par[1] + (par[2]* I) / (par[3] + I)
}


# create new lists
d_DO <- list()
for (i in 1:max(PI_field$Days)) {
  d_DO[[i]] <- PI_field$d_DO_Chl[PI_field$Days == i]
}

I <- list()
for (i in 1:max(PI_field$Days)) {
  I[[i]] <- PI_field$I_in_water [PI_field$Days == i]
}



# Fit the non-linear model
PI_curve_coef <- list()

#I_v <- I[[179]]
#coef(nls(d_DO[[179]] ~ r + (P_max * I_v) / (h_i + I_v),  start=list(r=-0.4, P_max=0.3,h_i=400)))

for (i in 1:max(PI_field$Days)){
  tryCatch({
    I_v <- I[[i]]
    PI_curve_coef[[i]] <- coef(nls(d_DO[[i]] ~ r + (P_max * I_v) / (h_i + I_v),  start=list(r=-0.4, P_max=0.3,h_i=400)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

PI_curve_coef


PI_field_coef <-  data.frame(Reduce(rbind,PI_curve_coef))

plot(PI_field_coef$r)

write.table(PI_field_coef,"C:/Users/admin/Dropbox/MSc Maria/Maria/Lab/PI_curve/PI_field_coef.txt", sep = ";")

# If the Chl_lab is used then the respiration decreases over the year. Even when there is corrected for light















## check for the daily average temperature

Temp <- list()
for(I in 1:max(PI_field$Days)){
  tryCatch({
    Temp[[i]] <- mean(PI_field_list[[i]]$T_C_water)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

with(PI_field, mean(PI_field$T_C_water ~ PI_field$Days))
Temp <- by(PI_field$T_C_water, PI_field$Days, function(x) mean(x))



# output the fit from the model
fit_PI_field <- lapply(PI_curve_coef, function(x) Michaelis_Menten(x,I1))
lapply(fit_PI_field, function(x) plot(x ~ I1,type="l"))
plot(fit_PI_field[[8]]~I1)
write.table(fit_PI_field,"fit_PI_field.txt", sep = "\t", row.names = FALSE)


PI_field_coef_lightWater <-  data.frame(Reduce(rbind,PI_curve_coef))

