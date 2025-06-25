#check for packages https://rdrr.io/all/cran/	
library(AquaEnv)	
library(LakeMetabolizer)	
library(ggplot2)	
library(readr)	
library(plyr)	
library(tidyverse)	
# install.packages("C:/Users/admin/Dropbox/MSc Maria/Maria/Misc/Msc2017Jolanda/Msc2017Jolanda.Rproj", repos = NULL)	
library("Msc2017J", lib.loc="~/R/win-library/3.2")	
#install.packages("REddyProc", repos="http://R-Forge.R-project.org")	

library(REddyProc)	
library(imputeTS)	

sEddyProc$sMDSGapFillAfterUstar	



# Calculate the ECS measured wind speed at 10m	
#Wind <- read_csv("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/Wind.csv",	
#                      col_types = cols(`Average Sampling time` = col_datetime(format = "%m/%d/%Y %H:%M")))	
#Wind$windECS_U10 <-wind.scale.base(wnd = Wind$wind, wnd.z = 2.2)	
#write.csv(x = Wind, file = "C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/Wind_Rout.csv")	


# setwd("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Basic_stat")	

# Change the setting of the graph position cause otherwise the graphs falls out the window 	
mar.default <- c(5,4,4,2) + 0.1	
par(mar = mar.default + c(0, 2, 0, 0))	


bootReg <- function(formula, data, indices) {	
  d <- data[indices,]	
  fit <- lm(formula, data = d)	
  return(coef(fit))	
}	




PI_curve <- function(ChangeGas_Chl, Sample, PAR, Gas) {	
  
  if(Gas == "DIC"){	
    ChangeGas_Chl <- ChangeGas_Chl*-1	
  }	
  
  else{ChangeGas_Chl <- ChangeGas_Chl}	
  
  DF_Metabolic_PI <-  data.frame(ChangeGas_Chl)	
  DF_Metabolic_PI$change_Chl <-  ChangeGas_Chl	
  DF_Metabolic_PI$Sample <-  Sample	
  DF_Metabolic_PI$PAR <-  PAR	
  
  DF_Metabolic_PI <- DF_Metabolic_PI[complete.cases(DF_Metabolic_PI["change_Chl"]),]	
  
  
  Nr_samples <- nrow(table(DF_Metabolic_PI$Sample))	
  
  # create new lists	
  NEP <- list()	
  for (i in 1:Nr_samples) {	
    NEP[[i]] <- DF_Metabolic_PI$change_Chl[DF_Metabolic_PI$Sample == i]	
  }	
  
  I <- list()	
  for (i in 1:Nr_samples) {	
    I[[i]] <- DF_Metabolic_PI$PAR[DF_Metabolic_PI$Sample == i]	
  }	
  
  
  
  # Fit the non-linear model	
  PI_field_coef_DIC <- list()	
  
  
  for (i in 1:Nr_samples){	
    tryCatch({	
      I_v <- I[[i]]	
      PI_field_coef_DIC[[i]] <- coef(nls(NEP[[i]] ~ r + (P_max * I_v) / (h_i + I_v),  start=list(r=-0.04, P_max=0.16, h_i=168)))	
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
  }	
  
  return(PI_field_coef_DIC)	
  
}	



Moving_average <- function(Data, Steps) {	
  
  Steps1 <- rep(1/Steps, Steps)	
  
  Out <- filter(Data, Steps1, sides=2)	
  
  return(Out)	
  
}	



Change_perHour <- function(Data, Date_Time, filter_extreem, Filter_time_above, Filter_time_below) {	
  
  
  # Calculate the difference in Data_vector	
  Diff_Data <- append(diff(Data), 0, 0)	
  
  
  # Calculate the time difference, per step.	
  ind_j <- 1:(length(Data)-1)	
  ind_k <- 2:length(Data)	
  ind_w <- 1:length(Data)	
  
  TimeDif <- as.array(Diff_Data)	
  
  for (i in seq_along(ind_j)){	
    TimeDif[ind_k[i]] <-	
      60/(as.double(Date_Time[ind_k[i]] - Date_Time[ind_j[i]]))	
  }	
  
  TimeDif[1] <- 0	
  
  
  # Filters out extreme time differences	
  if(filter_extreem == TRUE) {	
    
    
    if( Filter_time_above == 0 & Filter_time_below == 0 ) {	
      
    }	
    
    else{	
      
      if(Filter_time_above == 0 & Filter_time_below != 0) {	
        TimeDif[TimeDif > 60/Filter_time_below] = NA	
      }	
      
      else {	
        
        if(Filter_time_above != 0 & Filter_time_below == 0) {	
          TimeDif[TimeDif < 60/Filter_time_above] = NA	
        }	
        
        else {	
          TimeDif[TimeDif < 60/Filter_time_above | TimeDif > 60/Filter_time_below] = NA	
        }	
        
      }	
    }	
    
    Out <- Diff_Data *   TimeDif	
    
  }	
  
  else  {	
    
    Out <- Diff_Data * TimeDif	
    
  }	
  
  return(Out)	
  
}	





# WindRose.R	
require(ggplot2)	
require(RColorBrewer)	

plot.windrose <- function(data,	
                          spd,	
                          dir,	
                          spdres = 2,	
                          dirres = 30,	
                          spdmin = 2,	
                          spdmax = 20,	
                          spdseq = NULL,	
                          palette = "YlGnBu",	
                          countmax = NA,	
                          debug = 0){	
  
  
  # Look to see what data was passed in to the function	
  if (is.numeric(spd) & is.numeric(dir)){	
    # assume that we've been given vectors of the speed and direction vectors	
    data <- data.frame(spd = spd,	
                       dir = dir)	
    spd = "spd"	
    dir = "dir"	
  } else if (exists("data")){	
    # Assume that we've been given a data frame, and the name of the speed 	
    # and direction columns. This is the format we want for later use.    	
  }  	
  
  # Tidy up input data ----	
  n.in <- NROW(data)	
  dnu <- (is.na(data[[spd]]) | is.na(data[[dir]]))	
  data[[spd]][dnu] <- NA	
  data[[dir]][dnu] <- NA	
  
  # figure out the wind speed bins ----	
  if (missing(spdseq)){	
    spdseq <- seq(spdmin,spdmax,spdres)	
  } else {	
    if (debug >0){	
      cat("Using custom speed bins \n")	
    }	
  }	
  # get information about the number of bins, etc.	
  n.spd.seq <- length(spdseq)	
  n.colors.in.range <- n.spd.seq - 1	
  
  # create the color map	
  spd.colors <- colorRampPalette(brewer.pal(min(max(3,	
                                                    n.colors.in.range),	
                                                min(9,	
                                                    n.colors.in.range)),                                               	
                                            palette))(n.colors.in.range)	
  
  if (max(data[[spd]],na.rm = TRUE) > spdmax){    	
    spd.breaks <- c(spdseq,	
                    max(data[[spd]],na.rm = TRUE))	
    spd.labels <- c(paste(c(spdseq[1:n.spd.seq-1]),	
                          '-',	
                          c(spdseq[2:n.spd.seq])),	
                    paste(spdmax,	
                          "-",	
                          max(data[[spd]],na.rm = TRUE)))	
    spd.colors <- c(spd.colors, "grey50")	
  } else{	
    spd.breaks <- spdseq	
    spd.labels <- paste(c(spdseq[1:n.spd.seq-1]),	
                        '-',	
                        c(spdseq[2:n.spd.seq]))    	
  }	
  data$spd.binned <- cut(x = data[[spd]],	
                         breaks = spd.breaks,	
                         labels = spd.labels,	
                         ordered_result = TRUE)	
  # clean up the data	
  data. <- na.omit(data)	
  
  # figure out the wind direction bins	
  dir.breaks <- c(-dirres/2,	
                  seq(dirres/2, 360-dirres/2, by = dirres),	
                  360+dirres/2)  	
  dir.labels <- c(paste(360-dirres/2,"-",dirres/2),	
                  paste(seq(dirres/2, 360-3*dirres/2, by = dirres),	
                        "-",	
                        seq(3*dirres/2, 360-dirres/2, by = dirres)),	
                  paste(360-dirres/2,"-",dirres/2))	
  # assign each wind direction to a bin	
  dir.binned <- cut(data[[dir]],	
                    breaks = dir.breaks,	
                    ordered_result = TRUE)	
  levels(dir.binned) <- dir.labels	
  data$dir.binned <- dir.binned	
  
  # Run debug if required ----	
  if (debug>0){    	
    cat(dir.breaks,"\n")	
    cat(dir.labels,"\n")	
    cat(levels(dir.binned),"\n")       	
  }  	
  
  # deal with change in ordering introduced somewhere around version 2.2	
  if(packageVersion("ggplot2") > "2.2"){    	
    cat("\n")	
    data$spd.binned = with(data, factor(spd.binned, levels = rev(levels(spd.binned))))	
    spd.colors = rev(spd.colors)	
  }	
  
  # create the plot ----	
  p.windrose <- ggplot(data = data,	
                       aes(x = dir.binned,	
                           fill = spd.binned)) +	
    geom_bar() + 	
    scale_x_discrete(drop = FALSE,	
                     labels = waiver()) +	
    coord_polar(start = -((dirres/2)/360) * 2*pi) +	
    scale_fill_manual(name = "Wind Speed (m/s)", 	
                      values = spd.colors,	
                      drop = FALSE) +	
    theme(axis.title.x = element_blank())	
  
  # adjust axes if required	
  if (!is.na(countmax)){	
    p.windrose <- p.windrose +	
      ylim(c(0,countmax))	
  }	
  
  # print the plot	
  print(p.windrose)  	
  
  # return the handle to the wind rose	
  return(p.windrose)	
}	







Sampling_Dates <- as.Date(c("18-5-2016","25-5-2016", "1-6-2016", "15-6-2016", "23-6-2016", "29-6-2016", "6-7-2016", "13-7-2016", "20-7-2016", "27-7-2016", "3-8-2016", "10-8-2016","17-8-2016", "24-8-2016", "31-8-2016", "7-9-2016", "14-9-2016", "21-9-2016", "28-9-2016", "6-10-2016", "12-10-2016", "19-10-2016", "26-10-2016", "4-11-2016", "9-11-2016", "16-11-2016"), format = "%d-%m-%Y")	





#Data imput	

DF_data <- read_delim("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/Merged_files_R.csv",	
                      col_types = cols(`Average Sampling time` = col_datetime(format = "%d/%m/%Y %H:%M"),	
                                       Date = col_date(format = "%d/%m/%Y"),	
                                       Date_dusk = col_date(format = "%d/%m/%Y")),";")	


# Check and filter data.	

{	
  # The relative humidity in the Schiphol data does not fall below 40.	
  # However, the data should be reliable therefore fill in NA for values lower than the lowest humidity	
  #Therefore, filter out the low values (there are a lot of values of <0.01, which are considered errors)	
  DF_data$`RH [%]`[DF_data$`RH [%]` <= 34] = NA	
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ###Remove errors (extreme high or low values)	
  #Chlorophyll	
  DF_data$`CHL aver mg/L`[DF_data$`CHL aver mg/L` > 1000 | DF_data$`CHL aver mg/L` < 80] <- NA	
  DF_data$`CHLV Volt`[DF_data$`CHLV Volt` < 0.05 | DF_data$`CHLV Volt` > 0.65] <- NA	
  DF_data$ChlVoltNorm <- 509.93 * DF_data$`CHLV Volt`	
  
  plot(DF_data$fraction_day[DF_data$Date == "2016-05-22"], DF_data$`CHL aver mg/L`[DF_data$Date == "2016-05-22"])	
  
  plot(DF_data$ChlVoltNorm)	
  
  #pH	
  # Remove visiting days, by removing points that decrease fast over time 	
  # in the plot these points are at poitns the raft was visited 	
  
  
  summary(append(abs(diff(DF_data$pH)), 0, 0))	
  
  
  DF_data$pH[is.na(DF_data$pH)] <- 0	
  DF_data$pH[append(abs(diff(DF_data$pH)), 0, 0) > 0.8] <- NA	
  DF_data$pH[DF_data$pH == 0] <- NA	
  DF_data$pH[is.na(DF_data$pH)] <- 0	
  DF_data$pH[append(abs(diff(DF_data$pH)), 0, 0) > 0.8] <- NA	
  DF_data$pH[DF_data$pH == 0] <- NA	
  
  
  DF_data$pH[c(14176, 14177, 12837, 14850)] <- NA	
  
  
  plot(DF_data$pH)	
  
  # There are some single low atmospheric presure point, which are most likely due to failure of the equipment	
  # Filter out the 3 times (6points) where the atm. presure falls below 9900	
  DF_data$`Atm_pressure [PA]`[DF_data$`Atm_pressure [PA]` < 9900] <- NA	
  
  DF_data$T_C_air <- DF_data$T_K_air - 273.15	
}	


# Fill gaps from the CO2 probe 	

DF_data$`Original CO2 probe water mol/L` <- DF_data$`CO2 probe water mol/L`	
DF_data$`CO2 probe water mol/L`[is.na(DF_data$`CO2 probe water mol/L`)] <- DF_data$`CO2  water mol/L`[is.na(DF_data$`CO2 probe water mol/L`)]	
DF_data$`CO2 probe water mol/m3` <- DF_data$`CO2 probe water mol/L`*1000	


plot(DF_data$Date, DF_data$`CO2 probe water mol/L`)	
points(DF_data$Date, DF_data$`CO2  water mol/L`, col = "red")	
points(DF_data$Date, DF_data$`Original CO2 probe water mol/L`, col = "black")	

# Calculate the mixing depth	
{	
  Profile_AVP <- read_csv("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Profile_AVP.csv",	
                          col_types = cols(Date = col_date(format = "%m/%d/%Y")))	
  
  library(dplyr)	
  Profile_AVP %>%	
    group_by(Date) %>%	
    meta.depths(Profile_AVP$'Temp C'[c(302,329)], Profile_AVP$'Depth'[c(302,329)])	
  # The mixing depth prpduces Nan's: meaning no thermocline. The mixing depth = 1.8m	
  # Otherwise the output is the min depth and max depth of teh thermocline	
  detach("package:dplyr", unload=TRUE)	
}	


# Calculate the time difference between steps	
{	
  DF_data$irr = 1	
  DF_data$irr[DF_data$`I (umol/m2/s)` == 0] = 0	
  
  ind_a <- (1:nrow(DF_data))	
  ind_b <- (2:(nrow(DF_data)+1))	
  
  for (i in ind_a) {	
    DF_data$Time_step[ind_a[i]] <- 60/as.integer(DF_data$`Average Sampling time`[ind_b[i]] - DF_data$`Average Sampling time`[ind_a[i]]) 	
  }	
}	



#Check the daily light hours	
{library(dplyr)	
  DF_data %>%	
    group_by(Date) %>%	
    summarise(irr0 = sum(irr[irr == 0]+1)/mean(Time_step, rm.na = TRUE),	
              irr1 = sum(irr[irr == 1])/mean(Time_step, rm.na = TRUE)	
    ) -> light_hours_day	
}	
# write.table(light_hours_day, "../../../Maria/light_hours_day", sep = ";")	

detach("package:dplyr", unload=TRUE)	





###############################################################################################################3	
#         Add moving average smoothed data	
################################################################################################################	

#Add smoothed lines to the dataset using moving averages	



DF_data$DIC_Alk_smooth <- Moving_average(DF_data$`DIC Alk mol/L`*1000000, 3)	
DF_data$DO_mg_smooth <- Moving_average(DF_data$`DO mg/L`, 3)	
DF_data$DO_umol_smooth <- Moving_average(DF_data$`DO umol/L`, 3)	

# DIC	
plot(DF_data$`Average Sampling time`, DF_data$DIC_Alk_smooth, col = "grey")	
# DO	
plot(DF_data$`Average Sampling time`, DF_data$DO_mg_smooth, col = "grey")	
plot(DF_data$`Average Sampling time`, DF_data$DO_umol_smooth, col = "grey")	



# Wind Johnson 1998.	
Wind_U10_Johnson <- DF_data$`WindSpeed (ECS)` * (10/1.95)^(0.108)	










#Remove the missing values of one minute difference which causing missing values in the smoothed values	
#DF_data <- DF_data[complete.cases(DF_data$H[DF_data$`Average Sampling time` > "2016-08-18 00:29:00" & DF_data$`Average Sampling time` < "2016-08-23 15:15:00"])]	


# Calculate the difference in DIC and OD over time	
Change_perHour	

#Calculate the DIC changes per hour (taken from the difference between each time step)	
DF_data$delta.dic <- Change_perHour(DF_data$DIC_Alk_smooth, DF_data$`Average Sampling time`, filter_extreem = TRUE, Filter_time_above = 15, Filter_time_below = 10)	


#Calculate the DO changes	
DF_data$delta.do <- Change_perHour(DF_data$DO_umol_smooth, DF_data$`Average Sampling time`, filter_extreem = TRUE, Filter_time_above = 15, Filter_time_below = 10)	





plot(DF_data$`Average Sampling time`, DF_data$delta.dic)	
plot(DF_data$`Average Sampling time`, DF_data$delta.do)	








# Check the cumulative hourly change over a day	
DF_data$Date <- format(DF_data$`Average Sampling time`, "%Y-%m-%d")	
DF_data$Date_dusk <- as.Date(DF_data$Date_dusk, format= "%d-%m-%Y")	

library(dplyr)	
DF_data %>%	
  group_by(Date_dusk) %>%	
  summarise(DIC_change_mean = mean(delta.dic),	
            DO_change_mean = mean(delta.do)) -> Sum_GasChange	
detach("package:dplyr", unload=TRUE)	

library(dplyr)	
DF_data %>%	
  group_by(Date) %>%	
  summarise(DIC = mean(DIC_Alk_smooth, na.rm = TRUE)) -> mean_DIC	
detach("package:dplyr", unload=TRUE)	

# write.table(mean_DIC, "/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Routput/Daily_mean_DIC")	

plot(Sum_GasChange$DO_change_mean)	
plot(Sum_GasChange$DIC_change_mean)	








###########################################################################################################################	
#	
# Calculate the CO2 flux	
#	
###########################################################################################################################	



# The flux by Cole	
k600_cole <- k.cole.base(DF_data$windECS_U10)                  #m/d	
DF_data$kCO2_Cole_R <- k600.2.kGAS.base(k600 = k600_cole, temperature = DF_data$T_C_water, gas = 'CO2')  #m/d	
DF_data$CO2flux_ColeProbeR <- (DF_data$kCO2_Cole_R * DF_data$`Enhancement factor` * (DF_data$`CO2 probe water mol/m3` - DF_data$`CO2 atm (umol/m3)`)  * 1000000) / 24	


plot(DF_data$`Average Sampling time`, DF_data$CO2flux_ColeProbeR, type = 'l')	
lines(DF_data$`Average Sampling time`, DF_data$`CO2flux_probe (umol/m2/h)`, col = "red")	

plot(DF_data$`Average Sampling time`, k600_cole, type = 'l')	







#Calculate the K600 using the macIntyre model	


# The function is described as:	
# k.macIntyre(ts.data, wnd.z, Kd, atm.press, params=c(1.2,0.4872,1.4784))	
LakeMetabolizer::k.macIntyre.base	

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
#	
# Note that in the MacIntyre script the parameter are as followed:	
#  C1 = param3, C2 = param1, and C3 = param2	
# Where C1, C2, C3 are the parameters given in the articles	




# The following function is adapted to the paper of Winslow et al 2016 where they use the model from MacIntyre et al., 2010.	
# The model as in the script from Winslow et al 2016 appears to be different than the paper	
# However, the flux matches the classic flux of Cole et al., better. The flux from this function is ~ -1e+05 	

adapted_Macintyre_Winslow_K600 <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd, 	
                                            RH, sw, lwnet, params = c(1.2, 0.4872, 1.4784)) 	
{	
  S_B <- 5.67e-08	
  emiss <- 0.972	
  Kelvin = 273.15	
  albedo_SW <- 0.07	
  vonK <- 0.41	
  swRat <- 0.46	
  mnWnd <- 0.2	
  g <- 9.81	
  C_w <- 4186	
  rpcI <- wnd < mnWnd	
  wnd[rpcI] <- mnWnd	
  mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
  C_D <- mm$C_D	
  E <- mm$alh	
  H <- mm$ash	
  dUdt <- sw * 0.93 - E - H + lwnet	
  Qo <- sw * (1 - albedo_SW) * swRat	
  rho_w <- water.density(Ts)	
  if (wnd.z != 10) {	
    e1 <- sqrt(C_D)	
    u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
  }	
  else {	
    u10 = wnd	
  }	
  rhoAir <- 1.2	
  vonK <- 0.41	
  tau <- C_D * u10^2 * rhoAir	
  uSt <- sqrt(tau/rho_w)	
  q1 <- 2 - 2 * exp(z.aml * -Kd)	
  q2 <- z.aml * Kd	
  q3 <- exp(z.aml * -Kd)	
  H_star <- dUdt - Qo * (q1/q2 - q3)	
  thermalExpFromTemp <- function(Ts) {	
    V <- 1	
    dT <- 0.001	
    T1 <- water.density(Ts)	
    T2 <- water.density(Ts + dT)	
    V2 <- T1/T2	
    dv_dT <- (V2 - V)/dT	
    alpha <- dv_dT	
    return(alpha)	
  }	
  tExp <- thermalExpFromTemp(Ts)	
  B1 = H_star * tExp * g	
  B2 = rho_w * C_w	
  Bflx = B1/B2	
  getKinematicVis <- function(Ts) {	
    tempTable <- seq(0, 100, by = 5)	
    visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897, 	
                  0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477, 	
                  0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
    v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
    v <- v * 1e-06	
    return(v)	
  }	
  kinV <- getKinematicVis(Ts)	
  
  
  #Adapted part# 	
  #adapted to the paper of Winslow 2016	
  KeNm = uSt^3	
  SmE = params[2] * (params[3]*Bflx + 1.76 * KeNm/(vonK * z.aml))	
  SmE[SmE < 0] = 0	
  Sk = SmE * kinV	
  Sk = Sk * 100^4 * 3600^4	
  Sk600 = params[1] * Sk^(1/4)	
  k600 <- as.numeric(Sk600)	
  k600 <- k600 * 24/100	
  return(k600)	
}	







# Where the parameters are changeble. Heiskanen et al. (2014) changed the parameters to c(0.5,0.77,0.3) for a small (max z = 10) stratified boreal lake	
# The fit the parameter to the model, as done is below from line 2040	

# dateTime  vector of datetime in POSIXct format	
# Ts        numeric value of surface water temperature, degC	
# airT      numeric value of air temperature, degC	
# Uz        numeric value of wind speed, m/s	
# RH        numeric value of relative humidity, %	
# atm.press atmospheric pressure in mb	
# wnd.z     height of wind measurement, m	
# airT.z    height of air temperature measurement, m (optional)	
# RH.z      height of relative humidity measurement, m (optional)	
# Kd Light  attenuation coefficient (Units:m^-1)	


# Set time and date	
DF_data$`Average Sampling time` <- as.POSIXct(DF_data$`Average Sampling time`)	


# remove missing values	
DF_data_macIntyre <- DF_data[complete.cases(DF_data[, c("Atm_pressure [PA]", "T_C_water", "T_K_air", "windECS_U10", "RH [%]", "I (umol/m2/s)",  "Sal", "WindSpeed (ECS)")]),]	
summary(DF_data_macIntyre)	
# calculate the shortwave radiation [W/m2]	
DF_data_macIntyre$sw_PAR <- par.to.sw.base(DF_data_macIntyre$`I (umol/m2/s)`, coeff=0.473)     # Defaults to value from Britton and Dodd (1976).                Britton, C. M., and J. D. Dodd. Relationships of photosynthetically active radiation and shortwave irradiance. Agricultural Meteorology 17, no. 1 (1976): 1-7.	

# calculate the longwave radiation [W/m2]	
DF_data_macIntyre$lw.net <- calc.lw.net.base(dateTime = DF_data_macIntyre$`Average Sampling time`, sw = DF_data_macIntyre$sw_PAR,	
                                             Ts = DF_data_macIntyre$T_C_water, lat = 52.296214, atm.press = (DF_data_macIntyre$`Atm_pressure [PA]`/10),	
                                             airT = DF_data_macIntyre$T_K_air - 273.15,  RH = DF_data_macIntyre$`RH [%]`)	

DF_data_macIntyre$Zaml <- 2.0086 * exp(-0.116 * DF_data_macIntyre$T_C_water)	


#  C1 = param3, C2 = param1, and C3 = param2	

DF_data_macIntyre1 <- DF_data_macIntyre	

DF_data_macIntyre <-  DF_data_macIntyre1 	





# The k600 was check by hand and it matches as calculated by hand, note that this K600 is given in m/d	

DF_data_macIntyre$k600.macIntyre <- k.macIntyre.base(dateTime = DF_data_macIntyre$`Average Sampling time`, wnd.z = 10,	
                                                     atm.press = DF_data_macIntyre$`Atm_pressure [PA]`/10, Ts = DF_data_macIntyre$T_C_water,	
                                                     z.aml = DF_data_macIntyre$Zaml,	
                                                     # z.aml = 0.15,	
                                                     airT = DF_data_macIntyre$T_K_air-273.15, wnd = DF_data_macIntyre$windECS_U10,	
                                                     RH = DF_data_macIntyre$`RH [%]`, sw =  DF_data_macIntyre$sw_PAR, lwnet =  DF_data_macIntyre$lw.net,	
                                                     Kd = 0.0222177372 * DF_data_macIntyre$Chl_lab + 1.3858280092, params = c(1.5600460, 0.2579, 0.5930818))   # m/d	




MAcIntyre_Bfl <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                           RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
{	
  S_B <- 5.67e-08	
  emiss <- 0.972	
  Kelvin = 273.15	
  albedo_SW <- 0.07	
  vonK <- 0.41	
  swRat <- 0.46	
  mnWnd <- 0.2	
  g <- 9.81	
  C_w <- 4186	
  rpcI <- wnd < mnWnd	
  wnd[rpcI] <- mnWnd	
  mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
  C_D <- mm$C_D	
  E <- mm$alh	
  H <- mm$ash	
  dUdt <- sw * 0.93 - E - H + lwnet	
  Qo <- sw * (1 - albedo_SW) * swRat	
  rho_w <- water.density(Ts)	
  if (wnd.z != 10) {	
    e1 <- sqrt(C_D)	
    u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
  }	
  else {	
    u10 = wnd	
  }	
  rhoAir <- 1.2	
  vonK <- 0.41	
  tau <- C_D * u10^2 * rhoAir	
  uSt <- sqrt(tau/rho_w)	
  q1 <- 2 - 2 * exp(z.aml * -Kd)	
  q2 <- z.aml * Kd	
  q3 <- exp(z.aml * -Kd)	
  H_star <- dUdt - Qo * (q1/q2 - q3)	
  thermalExpFromTemp <- function(Ts) {	
    V <- 1	
    dT <- 0.001	
    T1 <- water.density(Ts)	
    T2 <- water.density(Ts + dT)	
    V2 <- T1/T2	
    dv_dT <- (V2 - V)/dT	
    alpha <- dv_dT	
    return(alpha)	
  }	
  tExp <- thermalExpFromTemp(Ts)	
  B1 = H_star * tExp * g	
  B2 = rho_w * C_w	
  Bflx = B1/B2	
  getKinematicVis <- function(Ts) {	
    tempTable <- seq(0, 100, by = 5)	
    visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                  0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                  0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
    v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
    v <- v * 1e-06	
    return(v)	
  }	
  kinV <- getKinematicVis(Ts)	
  KeNm = uSt^3	
  SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
  SmE[SmE < 0] = 0	
  Sk = SmE * kinV	
  Sk = Sk * 100^4 * 3600^4	
  Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
  k600 <- as.numeric(Sk600)	
  k600 <- k600 * 24/100	
  return(Bflx)	
}	


MAcIntyre_KeNm <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                            RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
{	
  S_B <- 5.67e-08	
  emiss <- 0.972	
  Kelvin = 273.15	
  albedo_SW <- 0.07	
  vonK <- 0.41	
  swRat <- 0.46	
  mnWnd <- 0.2	
  g <- 9.81	
  C_w <- 4186	
  rpcI <- wnd < mnWnd	
  wnd[rpcI] <- mnWnd	
  mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
  C_D <- mm$C_D	
  E <- mm$alh	
  H <- mm$ash	
  dUdt <- sw * 0.93 - E - H + lwnet	
  Qo <- sw * (1 - albedo_SW) * swRat	
  rho_w <- water.density(Ts)	
  if (wnd.z != 10) {	
    e1 <- sqrt(C_D)	
    u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
  }	
  else {	
    u10 = wnd	
  }	
  rhoAir <- 1.2	
  vonK <- 0.41	
  tau <- C_D * u10^2 * rhoAir	
  uSt <- sqrt(tau/rho_w)	
  q1 <- 2 - 2 * exp(z.aml * -Kd)	
  q2 <- z.aml * Kd	
  q3 <- exp(z.aml * -Kd)	
  H_star <- dUdt - Qo * (q1/q2 - q3)	
  thermalExpFromTemp <- function(Ts) {	
    V <- 1	
    dT <- 0.001	
    T1 <- water.density(Ts)	
    T2 <- water.density(Ts + dT)	
    V2 <- T1/T2	
    dv_dT <- (V2 - V)/dT	
    alpha <- dv_dT	
    return(alpha)	
  }	
  tExp <- thermalExpFromTemp(Ts)	
  B1 = H_star * tExp * g	
  B2 = rho_w * C_w	
  Bflx = B1/B2	
  getKinematicVis <- function(Ts) {	
    tempTable <- seq(0, 100, by = 5)	
    visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                  0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                  0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
    v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
    v <- v * 1e-06	
    return(v)	
  }	
  kinV <- getKinematicVis(Ts)	
  KeNm = uSt^3	
  SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
  SmE[SmE < 0] = 0	
  Sk = SmE * kinV	
  Sk = Sk * 100^4 * 3600^4	
  Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
  k600 <- as.numeric(Sk600)	
  k600 <- k600 * 24/100	
  return(KeNm)	
}	




MAcIntyre_Sk <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                          RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
{	
  S_B <- 5.67e-08	
  emiss <- 0.972	
  Kelvin = 273.15	
  albedo_SW <- 0.07	
  vonK <- 0.41	
  swRat <- 0.46	
  mnWnd <- 0.2	
  g <- 9.81	
  C_w <- 4186	
  rpcI <- wnd < mnWnd	
  wnd[rpcI] <- mnWnd	
  mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
  C_D <- mm$C_D	
  E <- mm$alh	
  H <- mm$ash	
  dUdt <- sw * 0.93 - E - H + lwnet	
  Qo <- sw * (1 - albedo_SW) * swRat	
  rho_w <- water.density(Ts)	
  if (wnd.z != 10) {	
    e1 <- sqrt(C_D)	
    u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
  }	
  else {	
    u10 = wnd	
  }	
  rhoAir <- 1.2	
  vonK <- 0.41	
  tau <- C_D * u10^2 * rhoAir	
  uSt <- sqrt(tau/rho_w)	
  q1 <- 2 - 2 * exp(z.aml * -Kd)	
  q2 <- z.aml * Kd	
  q3 <- exp(z.aml * -Kd)	
  H_star <- dUdt - Qo * (q1/q2 - q3)	
  thermalExpFromTemp <- function(Ts) {	
    V <- 1	
    dT <- 0.001	
    T1 <- water.density(Ts)	
    T2 <- water.density(Ts + dT)	
    V2 <- T1/T2	
    dv_dT <- (V2 - V)/dT	
    alpha <- dv_dT	
    return(alpha)	
  }	
  tExp <- thermalExpFromTemp(Ts)	
  B1 = H_star * tExp * g	
  B2 = rho_w * C_w	
  Bflx = B1/B2	
  getKinematicVis <- function(Ts) {	
    tempTable <- seq(0, 100, by = 5)	
    visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                  0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                  0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
    v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
    v <- v * 1e-06	
    return(v)	
  }	
  kinV <- getKinematicVis(Ts)	
  KeNm = uSt^3	
  SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
  SmE[SmE < 0] = 0	
  Sk = SmE * kinV	
  Sk = Sk * 100^4 * 3600^4	
  Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
  k600 <- as.numeric(Sk600)	
  k600 <- k600 * 24/100	
  return(Sk)	
}	



MAcIntyre_Sme <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                           RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
{	
  S_B <- 5.67e-08	
  emiss <- 0.972	
  Kelvin = 273.15	
  albedo_SW <- 0.07	
  vonK <- 0.41	
  swRat <- 0.46	
  mnWnd <- 0.2	
  g <- 9.81	
  C_w <- 4186	
  rpcI <- wnd < mnWnd	
  wnd[rpcI] <- mnWnd	
  mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
  C_D <- mm$C_D	
  E <- mm$alh	
  H <- mm$ash	
  dUdt <- sw * 0.93 - E - H + lwnet	
  Qo <- sw * (1 - albedo_SW) * swRat	
  rho_w <- water.density(Ts)	
  if (wnd.z != 10) {	
    e1 <- sqrt(C_D)	
    u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
  }	
  else {	
    u10 = wnd	
  }	
  rhoAir <- 1.2	
  vonK <- 0.41	
  tau <- C_D * u10^2 * rhoAir	
  uSt <- sqrt(tau/rho_w)	
  q1 <- 2 - 2 * exp(z.aml * -Kd)	
  q2 <- z.aml * Kd	
  q3 <- exp(z.aml * -Kd)	
  H_star <- dUdt - Qo * (q1/q2 - q3)	
  thermalExpFromTemp <- function(Ts) {	
    V <- 1	
    dT <- 0.001	
    T1 <- water.density(Ts)	
    T2 <- water.density(Ts + dT)	
    V2 <- T1/T2	
    dv_dT <- (V2 - V)/dT	
    alpha <- dv_dT	
    return(alpha)	
  }	
  tExp <- thermalExpFromTemp(Ts)	
  B1 = H_star * tExp * g	
  B2 = rho_w * C_w	
  Bflx = B1/B2	
  getKinematicVis <- function(Ts) {	
    tempTable <- seq(0, 100, by = 5)	
    visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                  0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                  0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
    v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
    v <- v * 1e-06	
    return(v)	
  }	
  kinV <- getKinematicVis(Ts)	
  KeNm = uSt^3	
  SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
  SmE[SmE < 0] = 0	
  Sk = SmE * kinV	
  Sk = Sk * 100^4 * 3600^4	
  Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
  k600 <- as.numeric(Sk600)	
  k600 <- k600 * 24/100	
  return(SmE)	
}	



MAcIntyre_Ust <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                           RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
{	
  S_B <- 5.67e-08	
  emiss <- 0.972	
  Kelvin = 273.15	
  albedo_SW <- 0.07	
  vonK <- 0.41	
  swRat <- 0.46	
  mnWnd <- 0.2	
  g <- 9.81	
  C_w <- 4186	
  rpcI <- wnd < mnWnd	
  wnd[rpcI] <- mnWnd	
  mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
  C_D <- mm$C_D	
  E <- mm$alh	
  H <- mm$ash	
  dUdt <- sw * 0.93 - E - H + lwnet	
  Qo <- sw * (1 - albedo_SW) * swRat	
  rho_w <- water.density(Ts)	
  if (wnd.z != 10) {	
    e1 <- sqrt(C_D)	
    u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
  }	
  else {	
    u10 = wnd	
  }	
  rhoAir <- 1.2	
  vonK <- 0.41	
  tau <- C_D * u10^2 * rhoAir	
  uSt <- sqrt(tau/rho_w)	
  q1 <- 2 - 2 * exp(z.aml * -Kd)	
  q2 <- z.aml * Kd	
  q3 <- exp(z.aml * -Kd)	
  H_star <- dUdt - Qo * (q1/q2 - q3)	
  thermalExpFromTemp <- function(Ts) {	
    V <- 1	
    dT <- 0.001	
    T1 <- water.density(Ts)	
    T2 <- water.density(Ts + dT)	
    V2 <- T1/T2	
    dv_dT <- (V2 - V)/dT	
    alpha <- dv_dT	
    return(alpha)	
  }	
  tExp <- thermalExpFromTemp(Ts)	
  B1 = H_star * tExp * g	
  B2 = rho_w * C_w	
  Bflx = B1/B2	
  getKinematicVis <- function(Ts) {	
    tempTable <- seq(0, 100, by = 5)	
    visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                  0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                  0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
    v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
    v <- v * 1e-06	
    return(v)	
  }	
  kinV <- getKinematicVis(Ts)	
  KeNm = uSt^3	
  SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
  SmE[SmE < 0] = 0	
  Sk = SmE * kinV	
  Sk = Sk * 100^4 * 3600^4	
  Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
  k600 <- as.numeric(Sk600)	
  k600 <- k600 * 24/100	
  return(uSt)	
}	




MAcIntyre_kinV <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                            RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
{	
  S_B <- 5.67e-08	
  emiss <- 0.972	
  Kelvin = 273.15	
  albedo_SW <- 0.07	
  vonK <- 0.41	
  swRat <- 0.46	
  mnWnd <- 0.2	
  g <- 9.81	
  C_w <- 4186	
  rpcI <- wnd < mnWnd	
  wnd[rpcI] <- mnWnd	
  mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
  C_D <- mm$C_D	
  E <- mm$alh	
  H <- mm$ash	
  dUdt <- sw * 0.93 - E - H + lwnet	
  Qo <- sw * (1 - albedo_SW) * swRat	
  rho_w <- water.density(Ts)	
  if (wnd.z != 10) {	
    e1 <- sqrt(C_D)	
    u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
  }	
  else {	
    u10 = wnd	
  }	
  rhoAir <- 1.2	
  vonK <- 0.41	
  tau <- C_D * u10^2 * rhoAir	
  uSt <- sqrt(tau/rho_w)	
  q1 <- 2 - 2 * exp(z.aml * -Kd)	
  q2 <- z.aml * Kd	
  q3 <- exp(z.aml * -Kd)	
  H_star <- dUdt - Qo * (q1/q2 - q3)	
  thermalExpFromTemp <- function(Ts) {	
    V <- 1	
    dT <- 0.001	
    T1 <- water.density(Ts)	
    T2 <- water.density(Ts + dT)	
    V2 <- T1/T2	
    dv_dT <- (V2 - V)/dT	
    alpha <- dv_dT	
    return(alpha)	
  }	
  tExp <- thermalExpFromTemp(Ts)	
  B1 = H_star * tExp * g	
  B2 = rho_w * C_w	
  Bflx = B1/B2	
  getKinematicVis <- function(Ts) {	
    tempTable <- seq(0, 100, by = 5)	
    visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                  0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                  0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
    v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
    v <- v * 1e-06	
    return(v)	
  }	
  kinV <- getKinematicVis(Ts)	
  KeNm = uSt^3	
  SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
  SmE[SmE < 0] = 0	
  Sk = SmE * kinV	
  Sk = Sk * 100^4 * 3600^4	
  Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
  k600 <- as.numeric(Sk600)	
  k600 <- k600 * 24/100	
  return(kinV)	
}	


# write.table(DF_data_macIntyre, "../Data_files/Routput/DF_data_macIntyre")	



DF_data_macIntyre$kCO2_macIntyre_R <- k600.2.kGAS.base(k600 = DF_data_macIntyre$k600.macIntyre , temperature = DF_data_macIntyre$T_C_water, gas = 'CO2') #m/d	

# The CO2 flux in umol/m2/h	

# !!!!!! warning the umol/m3 was an early typo in the data set. It is mol/m3!!!!!!!!!!!!!!	
# The concenrtations are both converted to mol/m3. 	
# The gas transfer velocity is in m/d see the help of the function	
# The outcome is mutiplied by 1000000 to get umol/m2/d	
DF_data_macIntyre$CO2flux_macIntyreR <- (DF_data_macIntyre$kCO2_macIntyre_R * DF_data_macIntyre$`Enhancement factor` *	
                                           (DF_data_macIntyre$`CO2 probe water mol/m3` - DF_data_macIntyre$`CO2 atm (umol/m3)`)*1000000)/24	


DF_data_macIntyre$CO2flux_probe_macIntyreR <- (DF_data_macIntyre$kCO2_macIntyre_R * DF_data_macIntyre$`Enhancement factor` *	
                                                 (DF_data_macIntyre$`Original CO2 probe water mol/L`*1000 - DF_data_macIntyre$`CO2 atm (umol/m3)`)*1000000)/24	

DF_data_macIntyre$CO2_ALK_flux_macIntyreR <- (DF_data_macIntyre$kCO2_macIntyre_R * DF_data_macIntyre$`Enhancement factor` *	
                                                (DF_data_macIntyre$`CO2  water mol/L`*1000 - DF_data_macIntyre$`CO2 atm (umol/m3)`)*1000000)/24	


# Check the fluxes of both the CO2 flux based on the probe (CO2flux_macIntyreR) and the CO2 based on the calculated flux based on the alkalinity (CO2_ALK_flux_macIntyreR)	
plot(DF_data_macIntyre$CO2flux_probe_macIntyreR, DF_data_macIntyre$CO2_ALK_flux_macIntyreR)	
summary(lm(DF_data_macIntyre$CO2_ALK_flux_macIntyreR~DF_data_macIntyre$CO2flux_probe_macIntyreR))	
# There is no great difference between the CO2 probe data and the calculated CO2	
# CO2 probe flux = -4.040e+01 + 1 * CO2 Alk flux, with an R2 of 0.995	

# What is the max influx for average values	
(mean(DF_data_macIntyre$kCO2_macIntyre_R, na.rm = TRUE) * mean(DF_data_macIntyre$`Enhancement factor`, na.rm = TRUE) * (0-mean(DF_data_macIntyre$`CO2 atm (umol/m3)`, na.rm = TRUE))*1000000)/1000/1.8/24	
(max(DF_data_macIntyre$kCO2_macIntyre_R, na.rm = TRUE) * max(DF_data_macIntyre$`Enhancement factor`, na.rm = TRUE) * (0-max(DF_data_macIntyre$`CO2 atm (umol/m3)`, na.rm = TRUE))*1000000)/1000/1.8/24	




plot(DF_data$`Average Sampling time`, DF_data$CO2flux_ColeProbeR, col = "darkgrey")	
lines(DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$CO2flux_macIntyreR != 0], DF_data_macIntyre$CO2flux_macIntyreR[DF_data_macIntyre$CO2flux_macIntyreR != 0], type = 'l', col = "blue")	

plot(DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$CO2flux_macIntyreR != 0], DF_data_macIntyre$CO2flux_macIntyreR[DF_data_macIntyre$CO2flux_macIntyreR != 0], type = 'l', col = "blue")	
points(DF_data$`Average Sampling time`, DF_data$CO2flux_ColeProbeR, col = "darkgrey")	


plot(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$`DIC Alk mol/L`*1000000, ylim = c(-11000, 4000))	
points(DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$CO2flux_macIntyreR != 0], DF_data_macIntyre$CO2flux_macIntyreR[DF_data_macIntyre$CO2flux_macIntyreR != 0], col="blue")	


plot(DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$k600.macIntyre != 0], DF_data_macIntyre$k600.macIntyre[DF_data_macIntyre$k600.macIntyre != 0], col = "darkgrey")	













##################################################################################################################3	
#	
# O2 flux	
#	
###################################################################################################################	

# Calculate the O2 flux m-1/d-1	
DF_data_macIntyre$kO2_macIntyre_R <- k600.2.kGAS.base(k600 = DF_data_macIntyre$k600.macIntyre , temperature = DF_data_macIntyre$T_C_water, gas = 'O2')	




# create a data file without missing values	
DF_Metabolic <- DF_data_macIntyre[complete.cases(DF_data_macIntyre[, c("DO mg/L")]),]	


#ts.data, Object of class data.frame with two named columns "datetime" and "wtr" (water temp in deg C).	
ts.data <- data.frame(1:nrow(DF_data_macIntyre))	
ts.data$datetime <- DF_data_macIntyre$`Average Sampling time`	
ts.data$wtr <- DF_data_macIntyre$T_C_water	
ts.data <- ts.data[,c(2:3)]	





# Calculate the O2 equilibrium based on the articel Cole et al., 2000 Persistence of net heterotrophy in lakes in mg/L	
DF_data_macIntyre$`DO_sat_mg/L` <- o2.at.sat(ts.data,baro = DF_data_macIntyre$`Atm_pressure [PA]`/10,  salinity = DF_data_macIntyre$Sal, model = "garcia-benson")[,2]	

# Calculate the DO saturation in mg/L	
DF_data_macIntyre$`DO_sat_umol/L` <- o2.at.sat(ts.data,baro = DF_data_macIntyre$`Atm_pressure [PA]`/10,  salinity = DF_data_macIntyre$Sal, model = "garcia-benson")[,2] * 0.0312511719189*1000	






# Calculate the DO flux in umol/m2/h	
DF_data_macIntyre$DOflux_umol <- DF_data_macIntyre$kO2_macIntyre_R *	
  (((DF_data_macIntyre$DO_mg_smooth * 1000 - DF_data_macIntyre$`DO_sat_mg/L` * 1000) *	
      0.0312511719189 * 1000 ) )/ 24	

plot(DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$DOflux_umol != 0], DF_data_macIntyre$DOflux_umol[DF_data_macIntyre$DOflux_umol != 0], col = "darkblue", type = 'l', ylab = "Oxygen flux ["~mu~"mol L"^"-1"*"h"^"-1"*"]", xlab = "Months [2016]", ylim = c(-10000, 50000))	
abline(h = 0 , col = "red")	




# Remove the zero's from the file	
DF_data_macIntyre$DOflux_umol[DF_data_macIntyre$DOflux_umol == 0] = NA	


# check the relation with the CO2 flux	
plot(DF_data_macIntyre$CO2flux_macIntyreR, DF_data_macIntyre$DOflux_umol)	
summary(lm(DF_data_macIntyre$CO2flux_macIntyreR~DF_data_macIntyre$DOflux_umol))	
p.adjust(2.2e-16, method = "bonferroni", n = 11660)	
# The flux of CO2 is about 20% of that of O2	


# Estimate Std. Error t value Pr(>|t|)	
# (Intercept)                   -7.580e+02  9.382e+00  -80.79   <2e-16 ***	
#  DF_data_macIntyre$DOflux_umol -2.014e-01  1.281e-03 -157.30   <2e-16 ***	
#  ---	
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1	

# Residual standard error: 690 on 11093 degrees of freedom	
# (4910 observations deleted due to missingness)	
# Multiple R-squared:  0.6905,	Adjusted R-squared:  0.6904
# F-statistic: 1.716e+04 on 1 and 10912 DF,  p-value: < 2.2e-16	




#################################################################################################################################	
#	
#  Check and plot the water and environment characteristics 	
#	
##################################################################################################################################	



summary(DF_data$`Average Sampling time`)	

DF_data$Date <- as.POSIXct(DF_data$Date)	


# Chlorophyll	

summary(DF_data$Chl_lab)	
sd(DF_data$Chl_lab)/sqrt(length(DF_data$Chl_lab))	
which.min(DF_data$Chl_lab)	
DF_data$Date[which.max(DF_data$Chl_lab)]	


summary(DF_data$ChlVoltNorm)	

DF_data$`Average Sampling time`[format(DF_data$Date, "%m") == "07"][[which.max(DF_data$ChlVoltNorm[format(DF_data$Date, "%m") == "07"])]] 	
DF_data$`Average Sampling time`[format(DF_data$Date, "%m") == "09"][[which.max(DF_data$ChlVoltNorm[format(DF_data$Date, "%m") == "09"])]] 	
DF_data$`Average Sampling time`[format(DF_data$Date, "%m") == "07"][[which.min(DF_data$ChlVoltNorm[format(DF_data$Date, "%m") == "07"])]] 	




mar.2ysides <- c(5,5,3,5) + 0.1	
par(mar = mar.2ysides)	
Sampling_Dates <- as.Date(Sampling_Dates)	

{	
  plot(DF_data$`Average Sampling time`, DF_data$ChlVoltNorm, col = "yellow3", type = 'l', ylim = c(0,450), 	
       ylab = "Chlorophyll "~paste(italic(a)) *"concentration ["~mu*"g L"^"-1"*"]", xlab = "Time [Months]", cex.lab = 1.2, cex.axis = 1.2)	
  lines(DF_data$Date[as.Date(DF_data$Date) == Sampling_Dates], DF_data$Chl_lab[as.Date(DF_data$Date) == Sampling_Dates], col = "darkgreen", lwd = 3)	
  points(DF_data$Date[as.Date(DF_data$Date) == Sampling_Dates], DF_data$Chl_lab[as.Date(DF_data$Date) == Sampling_Dates], cex = 1.5, bg = "black", pch = 16)	
  axis(side = 4, cex.axis = 1.2)	
  mtext(side = "4", "Chlorophyll"~paste(italic(a)) *"fluorescence [Arb.U]", line = 3, cex = 1.2)	
  legend(	
    "topright",	
    inset = 0.02,	
    bty = "n",	
    cex = 1.1,	
    lty = 1,	
    legend = c(expression(paste("Chl ", italic(a))*" concentration"), expression(paste("Chl ", italic(a)," fluorescence"))),	
    col = c("darkgreen", "yellow3")	
    
    
  )	
  
  
  
}	


par(mar = mar.default + c(0, 1, 0, 0))	

# Temperature 	
{	
  summary(DF_data$T_C_water)	
  sd(DF_data$T_C_water, na.rm = TRUE)/nrow(DF_data[complete.cases(DF_data["T_C_water"]),])	
  
  summary(DF_data$T_C_air)	
  sd(DF_data$T_C_air, na.rm = TRUE)/nrow(DF_data[complete.cases(DF_data["T_C_air"]),])	
  
  
  t.test(DF_data$T_C_air[DF_data$irr == 0 & months(DF_data$Date) == "oktober" | months(DF_data$Date) == "november" ], DF_data$T_C_water[DF_data$irr == 0 & months(DF_data$Date) == "oktober" | months(DF_data$Date) == "november" ])	
  
  t.test(DF_data$T_C_air, DF_data$T_C_water)	
  p.adjust(2.2e-16, method = "bonferroni", n = length(DF_data$T_C_water))	
  
  plot(DF_data$`Average Sampling time`, DF_data$T_C_air, type = 'l', col = "darkorange", ylim = c(0,45), ylab = "Temperature [ºC degrees]", xlab = "Months [2016]")	
  lines(DF_data$`Average Sampling time`, DF_data$T_C_water, type = 'l', col = "darkblue")	
  legend(	
    "topright",	
    inset = c(0.0, 0.0),	
    cex = 1.1,	
    bty = "n",	
    lty = 1,	
    col = c("darkorange", "darkblue"),	
    legend = c("Air temperature ", "Water temperature")	
  )	
}	

# Wind 	

{	
  summary(DF_data$windECS_U10)	
  DF_data$`Average Sampling time`[which.max(DF_data$windECS_U10)]	
  DF_data$`Average Sampling time`[which.min(DF_data$windECS_U10)]	
  sd(DF_data$windECS_U10, na.rm = TRUE)/nrow(DF_data[complete.cases(DF_data["windECS_U10"]),])	
  
  
  wind <- DF_data[complete.cases(DF_data$windECS_U10, DF_data$`wind_dir [deg_from_north]`),]	
  wind$date <- as.POSIXct(wind$Date)	
  
  t.test(wind$windECS_U10[wind$`Average Sampling time` < "2016-08-30"], 	
         wind$windECS_U10[wind$`Average Sampling time` >= "2016-08-30"])	
  
  p.adjust(2.2e-16, method = "bonferroni", n = length(wind$windECS_U10[wind$`Average Sampling time` < "2016-08-30"]))	
  
  sd(wind$windECS_U10[wind$`Average Sampling time` < "2016-08-30"], na.rm = TRUE)/length(wind$windECS_U10[wind$`Average Sampling time` < "2016-08-30"])	
  sd(wind$windECS_U10[wind$`Average Sampling time` >= "2016-08-30"], na.rm = TRUE)/length(wind$windECS_U10[wind$`Average Sampling time` >= "2016-08-30"])	
  
  plot(DF_data$`Average Sampling time`, DF_data$windECS_U10)	
  
  plot(DF_data$`Average Sampling time`, DF_data$`wind_dir [deg_from_north]`)	
  points(DF_data$`Average Sampling time`[DF_data$windECS_U10 < 6], DF_data$`wind_dir [deg_from_north]`[DF_data$windECS_U10 < 6], col = "darkgreen")	
  points(DF_data$`Average Sampling time`[DF_data$windECS_U10 > 6], DF_data$`wind_dir [deg_from_north]`[DF_data$windECS_U10 > 6], col = "darkred")	
  
  
  
  
  
  Windrose <- plot.windrose(data = wind, spd = wind$windECS_U10, dir = wind$`wind_dir [deg_from_north]`, spdmax = 12.04,spdmin = 1, spdseq = c(0,3,6,9,12))	
  
  
  library(openair)	
  
  windRose(mydata = wind, ws = "windECS_U10", wd = "wind_dir [deg_from_north]", 	
           auto.text= FALSE, paddle = FALSE, annotate = FALSE, breaks=c(0,3,6,9,13), type = "month", smooth  =TRUE)	
  
  detach("package:openair", unload = TRUE)	
  rm(wind)	
}	


#pH  and alkalinity	

{	
  #pH	
  summary(DF_data$pH)	
  sd(DF_data$T_C_water, na.rm = TRUE)/nrow(DF_data[complete.cases(DF_data["pH"]),])	
  
  # Alkalinity	
  plot(DF_data$Alk_Alk)	
  summary(DF_data$Alk_Alk)	
  sd(DF_data$Alk_Alk, na.rm = TRUE)/nrow(DF_data[complete.cases(DF_data["Alk_Alk"]),])	
  
  #for distance/margin the mar c (down, left, up, right)	
  
  layout(matrix(c(1,2),c(2)),	
         heights = c(2,2), respect = FALSE)	
  
  par(mar=c(1.5, 5, 3, 2))	
  plot(DF_data$`Average Sampling time`, DF_data$pH, type = 'l', col = "darkred", 	
       ylab = "pH", xlab = "", xaxt='n', ylim = c(8,10.5), cex.lab = 1.2, cex.axis = 1.2)	
  legend(	
    "topright",	
    inset = c(0.0, 0.0),	
    cex = 1.1,	
    bty = "n",	
    legend = "A ")	
  
  
  Sampling_Dates <- as.POSIXct(Sampling_Dates)	
  Alk_titr <- c(0.003653578	,
                0.00362152	,
                0.003597191	,
                0.002516293	,
                0.001776456	,
                0.001642561	,
                0.001698642	,
                0.001729858	,
                0.001606521	,
                0.001780728	,
                0.001655168	,
                0.001581944	,
                0.001685932	,
                0.001784072	,
                0.001828363	,
                0.001795344	,
                0.001768224	,
                0.001948843	,
                0.002017333	,
                0.002095987	,
                0.00215288	,
                0.002264173	,
                0.002359147	,
                0.00255	,
                0.002605	,
                0.002617773	
  )	
  
  
  
  Sampling_Dates <- as.POSIXct(Sampling_Dates)	
  par(mar=c(4.5, 5, 0, 2))	
  plot(DF_data$`Average Sampling time`, DF_data$Alk_Alk*1000, type = 'l', col = "white", lwd = 2	
       , xlab = "Time [Months]",  ylab = expression("Alk [meq L"^"-1"*"]"), cex.lab = 1.2, cex.axis = 1.2)	
  legend(	
    "topright",	
    inset = c(0.0, 0.0),	
    cex = 1.1,	
    bty = "n",	
    legend = "B ")	
  points(Sampling_Dates, Alk_titr*1000, pch = 20, col = "black", cex = 1.5)	
  lines(Sampling_Dates, Alk_titr*1000)	
  
  
  par(mfrow=c(1,1))	
  par(mar=c(5,4,4,2) + 0.1)	
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  DF_data$`Average Sampling time`[which.max(DF_data$pH)]	
  DF_data$`Average Sampling time`[which.min(DF_data$pH)]	
  
  
  Sampling_Dates <- format(Sampling_Dates, "%Y-%m-%d")	
  
  
  
  
  
  
  DF_data$date <- as.POSIXct(DF_data$`Average Sampling time`)	
  
  
  
  ##################################	
  # Check the evaporation and precipitation	
  
  {	
    
    
    
    # plot evaporation and precipitation	
    Neerslag_verdamping <- read_delim("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Neerslag-verdamping.csv", 	
                                      ";", escape_double = FALSE, col_types = cols(Day = col_date(format = "%d-%m-%Y")), 	
                                      trim_ws = TRUE)	
    
    Neerslag_verdamping <- Neerslag_verdamping[c(1:366),]	
    
    
    
    
    summary(Neerslag_verdamping$`neerslag [mm/d]`[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"])	
    sd(Neerslag_verdamping$`neerslag [mm/d]`[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"])/sqrt(length(Neerslag_verdamping$`neerslag [mm/d]`[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"]))	
    Neerslag_verdamping$Day[which.max(Neerslag_verdamping$`neerslag [mm/d]`)]	
    
    summary(Neerslag_verdamping$`verdamping penman [mm/d]`[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"])	
    sd(Neerslag_verdamping$`verdamping penman [mm/d]`[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"])/sqrt(length(Neerslag_verdamping$`verdamping penman [mm/d]`[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"]))	
    Neerslag_verdamping$Day[which.max(Neerslag_verdamping$`verdamping penman [mm/d]`)]	
    
    layout(matrix(c(1,2),c(2)),	
           heights = c(2,2), respect = FALSE)	
    
    par(mar=c(1, 5, 3.5, 2))	
    plot(Neerslag_verdamping$Day[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"], 	
         Neerslag_verdamping$`neerslag [mm/d]`[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"], 	
         type = 'l', col = "firebrick3", ylab = "Precipitation [mm/d]", xlab = "", xaxt='n', ylim = c(0,50))	
    legend(	
      "topright",	
      inset = c(0.0, 0.0),	
      cex = 1.1,	
      bty = "n",	
      legend = "A ")	
    
    
    
    par(mar=c(4.5, 5, 0, 2))	
    plot(Neerslag_verdamping$Day[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"],	
         Neerslag_verdamping$`verdamping penman [mm/d]`[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"],	
         type = 'l', col = "dodgerblue3", lwd = 1	
         , xlab = "Months [2016]",  ylab = expression("Evaporation [mm/d]"))	
    legend(	
      "topright",	
      inset = c(0.0, 0.0),	
      cex = 1.1,	
      bty = "n",	
      legend = "B ")	
    
  }	
  
  
  
  par(mfrow=c(1,1))	
  par(mar=c(5,4,4,2) + 0.1)	
  
  
  
}	

##################################################################################################################################################3	
#	
# Calculate metabolism	
#	
###################################################################################################################################################	


{	
  
  # add column with night (0) and daytime (1)	
  
  DF_data_macIntyre$irr = 1	
  DF_data_macIntyre$irr[DF_data_macIntyre$`I (umol/m2/s)` == 0] = 0	
  
  datetime <- DF_data_macIntyre$`Average Sampling time`	
  
  
  # Check the primary production in umol/L/h	
  metab.bookkeep(datetime = DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$Date == "2016-05-19"],	
                 do.obs = DF_data_macIntyre$DO_mg_smooth[DF_data_macIntyre$Date == "2016-05-19"],	
                 do.sat = DF_data_macIntyre$`DO_sat_mg/L`[DF_data_macIntyre$Date == "2016-05-19"],	
                 k.gas = DF_data_macIntyre$kO2_macIntyre_R[DF_data_macIntyre$Date == "2016-05-19"],	
                 z.mix = 1.8, irr = DF_data_macIntyre$irr[DF_data_macIntyre$Date == "2016-05-19"]) * 0.0312511719189 * 1000 / 24	
  
  
  
  
  # Create a new dataframe	
  DF_Metabolic <- DF_data_macIntyre[complete.cases(DF_data_macIntyre[, c("kO2_macIntyre_R")]),]	
  
  
  
  
  # Chech which days are complete	
  
  table(table(DF_Metabolic$Date_dusk) == 96)	
  
  # Give days a number	
  DF_Metabolic$Days <- as.Date(DF_Metabolic$Date_dusk) - as.Date(DF_Metabolic$Date_dusk[1])+1	
  DF_Metabolic$Days <- as.integer(DF_Metabolic$Days)	
  max(DF_Metabolic$Days)	
  
  
  # There is a bug in plyr and dplyr. If the following script function gives an errror, then detached the package dplyr	
  detach("package:dplyr", unload=TRUE)	
  
  library(plyr)	
  # Select for only complete days	
  Complete_dates <- data.frame(count(DF_Metabolic$Days)$x[count(DF_Metabolic$Days)$freq  == 96])	
  
  
  # calculate the metabolism in umol/L/h for the whole lake	
  DF_Metabolic_coef <- list()	
  for (i in count(DF_Metabolic$Days)$x[count(DF_Metabolic$Days)$freq  == 96]) {	
    DF_Metabolic_sub <- DF_Metabolic[c(DF_Metabolic$Days == i),]	
    tryCatch({	
      DF_Metabolic_coef[[i]] <- metab.bookkeep(datetime = DF_Metabolic_sub$`Average Sampling time`,	
                                               do.obs = DF_Metabolic_sub$DO_mg_smooth,	
                                               do.sat = DF_Metabolic_sub$`DO_sat_mg/L`,	
                                               k.gas = DF_Metabolic_sub$kO2_macIntyre_R,	
                                               z.mix = 1.8,	
                                               irr = DF_Metabolic_sub$irr)*0.0312511719189*1000/24 }	
      , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
  }	
  
  
  DF_Metabolic_coef	
  
  DF_Metabolic_coef <-  Reduce(rbind,mapply(c, DF_Metabolic_coef, SIMPLIFY=FALSE))	
  DF_Metabolic_coef <- as.data.frame(DF_Metabolic_coef)	
  
  mean(Reduce(rbind,mapply(c,DF_Metabolic_coef$GPP, SIMPLIFY=FALSE)))	
  median(Reduce(rbind,mapply(c,DF_Metabolic_coef$GPP, SIMPLIFY=FALSE)))	
  
  mean(Reduce(rbind,mapply(c,DF_Metabolic_coef$R, SIMPLIFY=FALSE)))	
  median(Reduce(rbind,mapply(c,DF_Metabolic_coef$R, SIMPLIFY=FALSE)))	
  
  mean(Reduce(rbind,mapply(c,DF_Metabolic_coef$NEP, SIMPLIFY=FALSE)))	
  median(Reduce(rbind,mapply(c,DF_Metabolic_coef$NEP, SIMPLIFY=FALSE)))	
  
  
  plot(Reduce(rbind,mapply(c,DF_Metabolic_coef$NEP, SIMPLIFY=FALSE)))	
  
  #GPP numeric estimate of Gross Primary Production, mgO2L???1d???1	
  #R numeric estimate of Respiration, mgO2L???1d???1	
  #NEP numeric estimate of Net Ecosystem production, mgO2L???1d???1	
  
  #        GPP         R     NEP	
  #     5.497985 -3.591273 0.9061053        Umol/l/h    19-05-2016	
  
  
  # The temperature dependend respiration from the PI_curve is:	
  DF_Metabolic$resp_lab <- exp(13.437)* exp(-0.3986/ (0.000086173303 * (DF_Metabolic$T_C_water + 273.15)))*DF_Metabolic$`CHL aver mg/L`	
  
  
  ########################################################3	
  # Used for the Report!!!!!!!!!	
  ##########################################################	
  
  
  # Check the model of  with a Kalman filter application (which is a smoothing function that is applied on the DO data)	
  
  DF_Metabolic_coef_Kal1 <- list()	
  for (i in count(DF_Metabolic$Days)$x[count(DF_Metabolic$Days)$freq  == 96]) {	
    DF_Metabolic_sub <- DF_Metabolic[c(DF_Metabolic$Days == i),]	
    tryCatch({	
      DF_Metabolic_coef_Kal1[[i]] <-   metab.kalman(do.obs = DF_Metabolic_sub$`DO mg/L`,	
                                                    do.sat = DF_Metabolic_sub$`DO_sat_mg/L`,	
                                                    k.gas = DF_Metabolic_sub$kO2_macIntyre_R,	
                                                    z.mix = 1.8,	
                                                    irr = DF_Metabolic_sub$`I (umol/m2/s)`,	
                                                    wtr = DF_Metabolic$T_C_water)}	
      , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
  }	
  
  for (i in count(DF_Metabolic$Days)$x[count(DF_Metabolic$Days)$freq  == 96]) {	
    DF_Metabolic_sub <- DF_Metabolic[c(DF_Metabolic$Days == i),]	
    tryCatch({ DF_Metabolic_coef_Kal1[[i]]$Date <- DF_Metabolic_sub$Date}	
             , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
  }	
  
  
  DF_Metabolic_coef_Kal1	
  
  
  DF_Metabolic_coef_Kal1 <-  Reduce(rbind,mapply(c, DF_Metabolic_coef_Kal1, SIMPLIFY=FALSE))	
  DF_Metabolic_coef_Kal1 <- as.data.frame(DF_Metabolic_coef_Kal1)	
  DF_Metabolic_coef_Kal1 <-  Reduce(rbind,mapply(c, DF_Metabolic_coef_Kal1$metab, DF_Metabolic_coef_Kal1$Date, SIMPLIFY=FALSE))	
  DF_Metabolic_coef_Kal1 <- DF_Metabolic_coef_Kal1[,c(1:4)]	
  DF_Metabolic_coef_Kal <- data.frame(as.numeric(DF_Metabolic_coef_Kal1[,1]))	
  DF_Metabolic_coef_Kal[2] <- data.frame(as.numeric(DF_Metabolic_coef_Kal1[,2]))	
  DF_Metabolic_coef_Kal[3] <- data.frame(as.numeric(DF_Metabolic_coef_Kal1[,3]))	
  DF_Metabolic_coef_Kal[4] <- data.frame(as.POSIXct(DF_Metabolic_coef_Kal1[,4]))	
  
  names(DF_Metabolic_coef_Kal) <- c("GPP","R", "NEP", "Date")	
  
  
  
  par(mfrow = c(3,1))	
  plot(DF_Metabolic_coef_Kal$Date, DF_Metabolic_coef_Kal$GPP)	
  plot(DF_Metabolic_coef_Kal$Date, DF_Metabolic_coef_Kal$R)	
  plot(DF_Metabolic_coef_Kal$Date, DF_Metabolic_coef_Kal$NEP)	
  abline(h = 0, col = "red")	
}	

# the Avergae metabolic values in mmol/L/d  	

{	
  mean(DF_Metabolic_coef_Kal[,1], na.rm = TRUE)*0.0312511719189 	
  median(DF_Metabolic_coef_Kal[,1], na.rm = TRUE)*0.0312511719189 	
  sd(DF_Metabolic_coef_Kal[,1], na.rm = TRUE)/sqrt(length(DF_Metabolic_coef_Kal[,1][is.na(DF_Metabolic_coef_Kal[,1]) == FALSE]))*0.0312511719189 	
  mean(DF_Metabolic_coef_Kal[,2], na.rm = TRUE)*0.0312511719189 	
  median(DF_Metabolic_coef_Kal[,2], na.rm = TRUE)*0.0312511719189 	
  sd(DF_Metabolic_coef_Kal[,2], na.rm = TRUE)/sqrt(length(DF_Metabolic_coef_Kal[,2][is.na(DF_Metabolic_coef_Kal[,2]) == FALSE]))*0.0312511719189  	
  mean(DF_Metabolic_coef_Kal[,3], na.rm = TRUE)*0.0312511719189 	
  median(DF_Metabolic_coef_Kal[,3], na.rm = TRUE)*0.0312511719189 	
  sd(DF_Metabolic_coef_Kal[,3], na.rm = TRUE)/sqrt(length(DF_Metabolic_coef_Kal[,3][is.na(DF_Metabolic_coef_Kal[,3]) == FALSE]))*0.0312511719189 	
}	

# The GPP is twice as high than R	
# The Kalma n filter results in simmilar NEP results as that fro Cole. However, GPP estimates and R are lower, probaly due to the smoothing function 	
# Because the unsmoothed data (as used here) is niosy the Kalman filter is used. 	













###################################################################################################################	
#	
# Benthic fluxes	
#	
###################################################################################################################	


# For now it is assumed that the mineralisation has a dDO:dCO2 of 1:1	
# Li, J., Crowe, S. A., Miklesh, D., Kistner, M., Canfield, D. E., & Katsev, S. (2012). Carbon mineralization and oxygen dynamics in sediments with deep oxygen penetration, Lake Superior. Limnology and Oceanography, 57(6), 1634-1650.	

# Verslag Tom de Ruyter	
DF_Metabolic$`CO2flux_benthic_umol/L/h` <- 10 ^(0.0515 * DF_Metabolic$T_C_water + 1.9949) * 1000 * -1/1800/24        # conversion from mg/m2/d to umol/L/h multiplied by -1 because it relfects a negative flux, C into the the lake	
DF_Metabolic$`O2flux_benthic_umol/L/h` <- 10 ^(0.0515 * DF_Metabolic$T_C_water + 1.9949) * 1000/1800/24              # conversion from mg/m2/d to umol/L/h O2 leaves the lake into the sediment, therefore the flux is positive	


# Gudasz	
#DF_Metabolic$`CO2flux_benthic_umol/L/h` <- ((10 ^(0.0362 * DF_Metabolic$T_C_water - 1.635)/44.0095) * 1000 * -1)/1800/24        # conversion from mg/m2/d to umol/L/h multiplied by -1 because it relfects a negative flux, C into the the lake	
#DF_Metabolic$`O2flux_benthic_umol/L/h` <- ((10 ^(0.0362 * DF_Metabolic$T_C_water - 1.635)/44.0095) * 1000)/1800/24              # conversion from mg/m2/d to umol/L/h O2 leaves the lake into the sediment, therefore the flux is positive	

# Gudasz, C., Bastviken, D., Steger, K., Premke, K., Sobek, S., & Tranvik, L. J. (2010). Temperature-controlled organic carbon mineralization in lake sediments. Nature, 466(7305), 478-481.	
# the estimated flux from the respicon is 0.2777778 umol/L/h, which is aprox a 4 fold higher than in articles	




# den Heyer, C., & Kalff, J. (1998). Organic matter mineralization rates in sediments: A within-and among-lake study. Limnology and oceanography, 43(4), 695-705.	
#DF_Metabolic$`CO2flux_benthic_umol/L/h` <- (0.98 * DF_Metabolic$T_C_water + 1.62 * 1000 * -1)/1800/24	









#################################################################################################################################################	
#	
# Estimate precipitation	
#	
#################################################################################################################################################	



# Zie excel file en paper :	

# Striegl, R. G., Reddy, M. M., Schuster, P. F., & Carter, V. (1994). Carbon budget for a groundwater-fed lake: Calcification supports summer photosynthesis. Limnol. Oceanogr, 39(6), 1319-1332.	


























#########################################################################################################################	
#	
# Set out PI-curves for the field, flux + d_DO	
#	
#########################################################################################################################	


# Used in report	














































































































{	
  # Calculate the irradiance at 60 cm	
  DF_Metabolic$`I(60)` <- DF_Metabolic$`I (umol/m2/s)` * exp(-(0.0222 * DF_Metabolic$Chl_lab + 1.3858) * 0.6)	
  
  
  
  
  # Calculate the DO flux in umol/L/h	
  DF_Metabolic$`DOflux_umol/L/h` <- DF_data_macIntyre$DOflux_umol/1800	
  plot(DF_data_macIntyre$`Average Sampling time`, DF_Metabolic$`DOflux_umol/L/h`, type = 'l')	
  
  
  plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$delta.do, type = 'l')	
  lines(DF_Metabolic$`Average Sampling time`,DF_Metabolic$`DOflux_umol/L/h`, col = "red")	
  
  
  
  # Calculate the total DO change in umol/L/h	
  DF_Metabolic$Tot_DO_change <- DF_Metabolic$delta.do + DF_Metabolic$`DOflux_umol/L/h` + DF_Metabolic$`O2flux_benthic_umol/L/h`	
  
  
  
  
  
  # plot the NEP per hour	
  plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$Tot_DO_change)	
  plot(DF_Metabolic$`I(60)`, DF_Metabolic$Tot_DO_change)	
  plot(DF_Metabolic$`I(60)`[DF_Metabolic$Date == "2016-05-19"], DF_Metabolic$Tot_DO_change[DF_Metabolic$Date == "2016-05-19"])	
  plot(DF_Metabolic$`I(60)`[DF_Metabolic$Days == 2], DF_Metabolic$Tot_DO_change[DF_Metabolic$Days == 2])	
  
  # The mean NEP is compareble to that of Bade	
  mean(DF_Metabolic$Tot_DO_change[DF_Metabolic$Date == "2016-10-24"])	
  # -0.3664637 the NEP from BAde is  -0.475337	
  
  
  
  
  
  
  
  # Check the PI_curve on the derived NEP in mg/h/L	
  ## assign the function	
  Michaelis_Menten = function(par, I) {	
    par[1] + (par[2]* I) / (par[3] + I)	
  }	
  
  DF_Metabolic_PI <- DF_Metabolic[complete.cases(DF_Metabolic[, c("Tot_DO_change")]),]	
  
  
  # Singular gradients are caused by too low error estimations of the residuals. This can be prevented by introducing some scatter. In this case it doesn't resolve anything, and produces more non-sence values	
  # add a bit of noise to avoid zero-residual data	
  # jitter_x <- jitter(DF_Metabolic_PI$`I(60)`,amount=0.02)	
  # jitter_y <- jitter(DF_Metabolic_PI$Tot_DO_change,amount=0.1)	
  # DF_Metabolic_PI$jitter_x <- jitter_x	
  # DF_Metabolic_PI$jitter_y <- jitter_y	
  
  # nls is also sensitive for initial parameter which are off.	
  # Therfore, another way to work around the errors of nls is to use nlsLM. This however also produces a lot of non-sence value such as a negative p_max	
  # library(minpack.lm)	
  
  
  
  # TO DO check to remove the parameters.	
  # DF_Metabolic_PI$Tot_DO_change_smooth <- filter(DF_Metabolic_PI$Tot_DO_change, half_hour, sides=2)	
  # DF_Metabolic_filter <- DF_Metabolic_PI$Tot_DO_change[DF_Metabolic_PI$Tot_DO_change > -0.3 & DF_Metabolic_PI$Tot_DO_change < 0.2]	
  # DF_Metabolic_filter$I <- DF_Metabolic_PI$`I(60)`[DF_Metabolic_PI$Tot_DO_change > -0.3 & DF_Metabolic_PI$Tot_DO_change < 0.2]	
  
  
  # The rate is given in umol/L/h	
  
  PI_curve_coef <- PI_curve(DF_Metabolic_PI$Tot_DO_change, Sample = DF_Metabolic_PI$Days, PAR = DF_Metabolic_PI$`I(60)`, Gas = "DO")	
  
  PI_curve_coef	
  
  
  
  
  # Create a list for the minimum temperature	
  Temp_min <- list()	
  for (i in 1:max(DF_Metabolic$Days)) {	
    Temp_min[[i]] <- min(DF_Metabolic$T_C_water[DF_Metabolic$Days == i])	
  }	
  
  
  Temp_max <- list()	
  for (i in 1:max(DF_Metabolic$Days)) {	
    Temp_max[[i]] <- max(DF_Metabolic$T_C_water[DF_Metabolic$Days == i])	
  }	
  
  
  Temp_av <- list()	
  for (i in 1:max(DF_Metabolic$Days)) {	
    Temp_av[[i]] <- mean(DF_Metabolic$T_C_water[DF_Metabolic$Days == i])	
  }	
  
  Light <- list()	
  for (i in 1:max(DF_Metabolic$Days)) {	
    Light[[i]] <- mean(DF_Metabolic$`I(60)`[DF_Metabolic$Days == i], na.rm = TRUE)	
  }	
  
  
  PI_field_coef <-  data.frame(Reduce(rbind,mapply(c, PI_curve_coef, Temp_min, Temp_max, Temp_av, Light, SIMPLIFY=FALSE)))	
  PI_field_coef$P_max[PI_field_coef$r > 0] = NA	
  PI_field_coef$h_i[PI_field_coef$r > 0] = NA	
  PI_field_coef$r[PI_field_coef$r > 0] = NA	
  
  
  
  # Check the metabolic rates	
  {	
    summary(PI_field_coef$P_max)	
    summary(PI_field_coef$h_i)	
    summary(PI_field_coef$r)	
    
    
    # The mean GPP in mmol/L/d	
    PI_field_coef$GPP <- (PI_field_coef$P_max * PI_field_coef$V7) / (PI_field_coef$h_i + PI_field_coef$V7)/1000*24	
    mean((PI_field_coef$P_max * PI_field_coef$V7) / (PI_field_coef$h_i + PI_field_coef$V7)/1000*24, na.rm = TRUE)	
    sd((PI_field_coef$P_max * PI_field_coef$V7) / (PI_field_coef$h_i + PI_field_coef$V7)/1000*24, na.rm = TRUE)/sqrt(nrow(PI_field_coef))	
    
    # The respiration in mmol/L/d	
    mean(PI_field_coef$r, na.rm = TRUE)/1000*24	
    sd(PI_field_coef$r/1000*24, na.rm = TRUE)/sqrt(nrow(PI_field_coef))	
    
    PI_field_coef$NEP <-PI_field_coef$GPP - PI_field_coef$r	
    
    mean(PI_field_coef$NEP, na.rm = TRUE)/1000*24	
    sd(PI_field_coef$NEP/1000*24, na.rm = TRUE)/sqrt(nrow(PI_field_coef))	
    
    
    # the rate of respiration (0.09624 mmol/L/d) is near to the estimated R by kalmans filter. 	
    
    
    
    PI_field_coef$P_R <- abs(PI_field_coef$P_max/PI_field_coef$r)	
    plot(PI_field_coef$P_R)	
    summary(PI_field_coef$P_R)	
    
    
    par(mfrow = c(3,1))	
    plot(PI_field_coef[,5], PI_field_coef$r * -1, xlab = "Temperature", ylab = "Respiration")	
    lines(PI_field_coef[,5], exp(13.437)* exp(-0.3986/ (0.000086173303 * (PI_field_coef[,4] + 273.15))), col = "red")	
    plot(PI_field_coef[,5], PI_field_coef$P_max, xlab = "Temperature", ylab = "P_max")	
    plot(PI_field_coef[,6], PI_field_coef$h_i, xlab = "Temperature", ylab = "h_i")	
    
    plot(PI_field_coef$r)	
    plot(PI_field_coef$P_max)	
    plot(PI_field_coef$h_i)	
    par(mfrow = c(1,1))	
    
    plot(PI_field_coef$r[PI_field_coef$P_max < 200], PI_field_coef$P_max[PI_field_coef$P_max < 200])	
    # plot(PI_field_coef$r, PI_field_coef$P_max, ylim = c(0,1), xlim = c(-0.3,0))	
    summary(lm(PI_field_coef$P_max[PI_field_coef$P_max < 200] ~ PI_field_coef$r[PI_field_coef$P_max < 200]))	
    
    library(ggplot2)	
    r_pmax <- ggplot(PI_field_coef, aes(PI_field_coef$r, PI_field_coef$P_max), colour = PI_field_coef$V6) +	
      geom_point()	
    
    ggplot(PI_field_coef, aes(x= r, y= P_max, colour=V6)) + geom_point(size=2) +	
      scale_colour_gradientn(colours=rainbow(4))	
    
  }	
  
  
  
}	


# Bentic flux = 0.5 -0.25 mmol/m2/h	
# min sediment = 1.5 - 2.5 mmol/m2/h	

# Herczeg, 1987 % gas exchange	


# CO2 =  flux atm + bent + min + R - phot	









##############################################################################################################################################	
#	
# Check the changes in DO and CO2	
#	
##############################################################################################################################################	


# Calculate the CO2 flux per umol/L/h	
DF_Metabolic$`CO2flux_umol/L/h` <- DF_Metabolic$CO2flux_macIntyreR/1800	

DF_Metabolic	


# Calculate the total CO2 change	
DF_Metabolic$Tot_DIC_change <- DF_Metabolic$delta.dic + DF_Metabolic$`CO2flux_umol/L/h` + DF_Metabolic$`CO2flux_benthic_umol/L/h`  + DF_Metabolic$Ca_prec_umolL - DF_Metabolic$`DIC_IN-OUT_umol/h`	
DF_Metabolic$Tot_DO_change <- DF_Metabolic$delta.do + DF_Metabolic$`DOflux_umol/L/h`  + DF_Metabolic$`O2flux_benthic_umol/L/h`	


plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$Tot_DIC_change, type = 'l', ylim = c(-200, 200), col = "red")	
lines(DF_Metabolic$`Average Sampling time`, DF_Metabolic$Tot_DO_change*-1)	















#################################################################################################################3	
#	
# Plot and summurise DO, DIC, light and metabolic	
#	
######################################################################################################################	

{	
  
  mar.default <- c(5,4,4,2) + 0.1	
  par(mar = mar.default + c(0, 2, 0, 0))	
  
  DF_data$`I(60)` <- DF_data$`I (umol/m2/s)` * exp(-(0.0222 * DF_data$Chl_lab + 1.3858) * 0.6)	
  
  
  library(dplyr)	
  DF_data %>%	
    group_by(Date) %>%	
    summarise(I = max(`I (umol/m2/s)`),	
              I60 = max(`I(60)`)) -> PAR_max	
  
  library(dplyr)	
  DF_data %>%	
    group_by(Date) %>%	
    summarise(I = mean(`I (umol/m2/s)`)*60*60*24/1000000,	
              I60 = mean(`I(60)`)*60*60*24/1000000) -> PAR_dailyDose	
  
  detach("package:dplyr", unload = TRUE)	
  
  
  
  # plot the max light irradiance	
  {	
    plot(PAR_max$Date[-1], PAR_max$I[-1],  col = "darkblue", pch = 20,	
         ylab = "PAR ["~mu~"E m"^"-1 "*"s"^"-1"*"]", xlab = "Month [2016]", cex = 1.1, cex.axis = 1.1, cex.lab =1.1, ylim = c(0,2100))	
    lines(PAR_max$Date[-1], PAR_max$I[-1],  col = "darkblue", lty = 3)	
    points(PAR_max$Date[-1], PAR_max$I60[-1],  col = "darkgreen", pch = 20)	
    lines(PAR_max$Date[-1], PAR_max$I60[-1],  col = "darkgreen", lty = 3)	
    legend(	
      "topright",	
      legend = c("PAR surface", "PAR 60 cm depth"),	
      pch = 20,	
      col = c("darkblue", "darkgreen"),	
      cex = 1.1,	
      bty  = "n")	
    
    # plot the daily dose	
    
    # plot the light irradiance	
    {	
      plot(PAR_dailyDose$Date[-1], PAR_dailyDose$I[-1],  col = "darkblue", pch = 20,	
           ylab = expression("PAR [mol m"^"-2 "*"d"^"-1"*"]"), xlab = "Time [Month]", cex = 1.1, cex.axis = 1.1, cex.lab =1.1, ylim = c(0,70))	
      lines(PAR_dailyDose$Date[-1], PAR_dailyDose$I[-1],  col = "darkblue", lty = 3)	
      points(PAR_dailyDose$Date[-1], PAR_dailyDose$I60[-1],  col = "darkgreen", pch = 20)	
      lines(PAR_dailyDose$Date[-1], PAR_dailyDose$I60[-1],  col = "darkgreen", lty = 3)	
      legend(	
        "topright",	
        legend = c("PAR surface", "PAR 60 cm depth"),	
        pch = 20,	
        col = c("darkblue", "darkgreen"),	
        cex = 1.1,	
        bty  = "n")	
    }	
    
    #DF_Metabolic$`I(60)`[DF_Metabolic$`I(60)` > 250] = NA	
    
    summary(PAR_dailyDose)	
    sd(PAR_dailyDose$I, na.rm = TRUE)/sqrt(length(PAR_dailyDose))	
    sd(PAR_dailyDose$I60, na.rm = TRUE)/sqrt(length(PAR_dailyDose))	
    
    DF_Metabolic$`I (umol/m2/s)`	
    mean(ddply(DF_Metabolic, .(Date), summarise, max=max(`I (umol/m2/s)`), mean = mean(`I (umol/m2/s)`))$max)	
    sd(ddply(DF_Metabolic, .(Date), summarise, max=max(`I (umol/m2/s)`), mean = mean(`I (umol/m2/s)`))$max)/length(ddply(DF_Metabolic, .(Date), summarise, max=max(`I (umol/m2/s)`), mean = mean(`I (umol/m2/s)`))$max)	
    which.max(ddply(DF_Metabolic, .(Date), summarise, max=max(`I (umol/m2/s)`), mean = mean(`I (umol/m2/s)`))$max)	
    which.min(ddply(DF_Metabolic, .(Date), summarise, max=max(`I (umol/m2/s)`), mean = mean(`I (umol/m2/s)`))$max[-c(1,34, 176, 157, 86, 131)])	
    
    ddply(DF_Metabolic, .(Date), summarise, max=max(`I (umol/m2/s)`), mean = mean(`I (umol/m2/s)`))$max==105.8	
    
    mean(ddply(DF_Metabolic, .(Date), summarise, max=max(`I(60)`), mean = mean(`I(60)`))$max, na.rm = TRUE)	
    sd(ddply(DF_Metabolic, .(Date), summarise, max=max(`I(60)`), mean = mean(`I(60)`))$max, na.rm = TRUE)/length(ddply(DF_Metabolic, .(Date), summarise, max=max(`I(60)`), mean = mean(`I(60)`))$max)	
    
    ddply(DF_Metabolic, .(Date), summarise, max=max(`I(60)`, na.rm = TRUE), mean = mean(`I(60)`, na.rm = TRUE))	
    which.max(ddply(DF_Metabolic, .(Date), summarise, max=max(`I(60)`, na.rm = TRUE), mean = mean(`I(60)`, na.rm = TRUE))$max)	
    which.min(ddply(DF_Metabolic, .(Date), summarise, max=max(`I(60)`, na.rm = TRUE), mean = mean(`I(60)`, na.rm = TRUE))$max)	
  }	
  
  
  
  
  # DIC and O2	
  
  {	
    
    library(dplyr)	
    DF_data_macIntyre %>%	
      group_by(Date_dusk) %>%	
      summarise(mean_DO_change = mean(DO_umol_smooth - `DO_sat_umol/L`, na.rm =TRUE),	
                mean_CO2_change = mean(`CO2 probe water mol/L`*1000000-`CO2 atm (umol/L)`*1000000, na.rm =TRUE)) -> DO_CO2	
    detach("package:dplyr", unload = TRUE)	
    
    plot(DO_CO2$Date_dusk, DO_CO2$mean_DO_change)	
    abline(h=0, v= c(DO_CO2$Date_dusk[DO_CO2$Date_dusk == "2016-05-24"], DO_CO2$Date_dusk[DO_CO2$Date_dusk == "2016-09-29"]), col = "red")	
    plot(DO_CO2$Date_dusk, DO_CO2$mean_CO2_change)	
    abline(h=0, v=  DO_CO2$Date_dusk[DO_CO2$Date_dusk == "2016-11-07"], col = "red", lty = 3)	
    
    
    
    # Summary for DO 	
    summary(DF_data_macIntyre$DO_umol_smooth) 	
    DF_data_macIntyre$`Average Sampling time`[which.max(DF_data_macIntyre$DO_umol_smooth)]	
    DF_data_macIntyre$`Average Sampling time`[which.min(DF_data_macIntyre$DO_umol_smooth)]	
    sd(DF_data_macIntyre$DO_umol_smooth, na.rm = TRUE)/length(DF_data_macIntyre$DO_umol_smooth)	
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # Summary for CO2	
    summary(DF_data_macIntyre$`CO2 probe water mol/L`*1000000) 	
    DF_data_macIntyre$`Average Sampling time`[which.max(DF_data_macIntyre$`CO2 probe water mol/L`*1000000)]	
    DF_data_macIntyre$`Average Sampling time`[which.min(DF_data_macIntyre$`CO2 probe water mol/L`*1000000)]	
    sd(DF_data_macIntyre$`CO2 probe water mol/L`*1000000, na.rm = TRUE)/length(DF_data_macIntyre$`CO2 probe water mol/L`*1000000)	
    
    
    
    
    
    
    
    
    
    
    # Summary for DIC	
    summary(DF_data_macIntyre$`DIC Alk mol/L`*1000000) 	
    DF_data_macIntyre$`Average Sampling time`[which.max(DF_data_macIntyre$`DIC Alk mol/L`*1000000)]	
    DF_data_macIntyre$`Average Sampling time`[which.min(DF_data_macIntyre$`DIC Alk mol/L`*1000000)]	
    sd(DF_data_macIntyre$`DIC Alk mol/L`*1000000, na.rm = TRUE)/length(DF_data_macIntyre$`DIC Alk mol/L`*1000000)	
    
    
    # Plot the graph	
    {       	
      layout(matrix(c(1,2,3),c(3)),	
             heights = c(1,1,1), respect = FALSE)	
      
      #layout(matrix(c(1,2),c(2)),	
      #       heights = c(1,1), respect = FALSE)	
      
      # plot the DO 	
      par(mar=c(0, 5, 5, 2), xpd = FALSE)	
      
      plot(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$DO_mg_smooth * 0.0312511719189*1000, col = "blue", type = 'l', 	
           ylab = "DO ["~mu~"mol L"^"-1"*"]", xlab = "", xaxt ='n', cex.axis=1.3, cex.lab = 1.5, ylim = c(00, 900))	
      
      rect(xleft = DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-05-23"], 	
           xright = DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-09-29"],	
           ybottom = -100, ytop = 1000, col = rgb(0,0,1, alpha = 0.03), border = NA)	
      abline(v= c(DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-05-23"], 	
                  DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-09-29"]), col = "grey", lty = 3)	
      lines(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$`DO_sat_umol/L`, col = "red" )	
      lines(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$DO_mg_smooth * 0.0312511719189*1000, col = "blue")	
      legend(	
        inset = c(0.05, -0.0),	
        "topleft",	
        horiz=TRUE,	
        cex = 1.3,	
        bty = "n",	
        lty = 1,	
        col = c("blue", "red"),	
        legend = c("DO measured ", "DO at equilibrium")	
        
        
      )	
      
      
      legend(	
        "topright",	
        inset = c(0.0, 0.0),	
        cex = 1.5,	
        bty = 'n',	
        box.col = "white",	
        bg = "white",	
        legend = "A ")	
      
      
      #plot CO2	
      par(mar=c(2.5, 5, 2.5, 2 ))	
      plot(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$`Original CO2 probe water mol/L`*1000000, col = "orange", 	
           type = 'l', ylab = "CO"[2]*" ["~mu~"mol L"^"-1"*"]", ylim = c(0,50), xlab = "", xaxt = 'n', cex.axis=1.3, cex.lab = 1.5)	
      rect(xleft = DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-05-19"], 	
           xright = DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-11-07"],	
           ybottom = -10, ytop = 850, col = rgb(0,0,1, alpha = 0.03), border = NA)	
      abline(v= c(DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-05-19"], 	
                  DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-11-07"]), col = "grey", lty = 3)	
      lines(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$`CO2  water mol/L`*1000000, col = "black" )	
      lines(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$`Original CO2 probe water mol/L`*1000000, col = "orange" )	
      lines(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$`CO2 atm (umol/L)`*1000000, col = "red" )	
      abline(v=  DO_CO2$Date_dusk[DO_CO2$Date_dusk == "2016-11-07"], col = "red", lty = 3)	
      legend(	
        inset = c(0.05, -0.0),	
        "topleft",	
        horiz=TRUE,	
        cex = 1.3,	
        bty = "n",	
        box.col = "white",	
        lty = 1,	
        col = c("orange", "black", "red"),	
        bg = "white",	
        legend = c(expression("CO"[2]*" measured"),expression("CO"[2]*" calculated"), expression("CO"[2]*" at equilibrium"))	
      )	
      legend(	
        "topright",	
        inset = c(0.0, 0.0),	
        cex = 1.5,	
        bty = 'n',	
        box.col = "white",	
        bg = "white",	
        legend = "B ")	
      
      
      
      #plot DIC	
      par(mar=c(5, 5, 0, 2))	
      plot(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$DIC_Alk_smooth, col = "darkorange4", 	
           type = 'l', ylab = "DIC ["~mu~"mol L"^"-1"*"]", xlab = ("Time [Months]"),  cex.axis=1.3, cex.lab = 1.5,	
           ylim = c(100,5000))	
      legend(	
        inset = c(0.05, -0.0),	
        "topleft",	
        horiz=TRUE,	
        cex = 1.3,	
        bty = "n",	
        lty = 1,	
        col = c("darkorange4"),	
        legend = c(expression("DIC measured  "))	
      )	
      legend(	
        "topright",	
        inset = c(0.0, 0.0),	
        cex = 1.5,	
        bty = "n",	
        box.col = "white",	
        bg = "white",	
        legend = "C ")	
      
      
      par(mfrow=c(1,1))	
      par(mar=c(5,4,4,2) + 0.1)	
    }	
    
  }	
  
  
  
  # DO and CO2	
  
  
  layout(matrix(c(1,2),c(2)),	
         heights = c(1,1), respect = FALSE)	
  
  # plot the DO 	
  par(mar=c(1, 5, 4, 2), xpd = FALSE)	
  
  plot(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$DO_mg_smooth * 0.0312511719189*1000, col = "blue", type = 'l', 	
       ylab = "DO ["~mu~"mol L"^"-1"*"]", xlab = "", xaxt ='n', cex.axis=1.2, cex.lab = 1.5, ylim = c(0, 1000))	
  
  rect(xleft = DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-05-23"], 	
       xright = DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-09-29"],	
       ybottom = -50, ytop = 1050, col = rgb(0,0,1, alpha = 0.03), border = NA)	
  abline(v= c(DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-05-23"], 	
              DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-09-29"]), col = "grey", lty = 3)	
  lines(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$`DO_sat_umol/L`, col = "red" )	
  lines(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$DO_mg_smooth * 0.0312511719189*1000, col = "blue")	
  legend(	
    inset = c(0.05, -0.0),	
    "topleft",	
    horiz=FALSE,	
    cex = 1.3,	
    bty = "n",	
    lty = 1,	
    col = c("blue", "red"),	
    legend = c("DO measured ", "DO at equilibrium")	
    
    
  )	
  
  
  
  #plot CO2	
  par(mar=c(5, 5, 0, 2 ))	
  plot(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$`CO2 probe water mol/L`*1000000, col = "orange", 	
       type = 'l', ylab = "CO"[2]*" ["~mu~"mol L"^"-1"*"]", ylim = c(0,50), xlab = "Months [2016]", cex.axis=1.2, cex.lab = 1.5)	
  rect(xleft = DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-05-19"], 	
       xright = DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-11-07"],	
       ybottom = -10, ytop = 850, col = rgb(0,0,1, alpha = 0.03), border = NA)	
  abline(v= c(DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-05-19"], 	
              DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-11-07"]), col = "grey", lty = 3)	
  lines(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$`CO2 atm (umol/L)`*1000000, col = "red" )	
  abline(v=  DO_CO2$Date_dusk[DO_CO2$Date_dusk == "2016-11-07"], col = "red", lty = 3)	
  legend(	
    inset = c(0.05, -0.0),	
    "topleft",	
    horiz=FALSE,	
    cex = 1.3,	
    bty = "n",	
    box.col = "white",	
    lty = 1,	
    col = c("orange",  "red"),	
    bg = "white",	
    legend = c(expression("CO"[2]*" measured"), expression("CO"[2]*" at equilibrium"))	
  )	
  
  
  
  
  
  
  
  
  
  
  # Plot Chl, DO, and pH	
  
  
  # Plot the graph	
  
  {	
    layout(matrix(c(1,2),c(2)),	
           heights = c(1.5,1), respect = FALSE)	
    
    par(mar=c(0.5, 5, 3, 2), xpd = FALSE)	
    
    plot(DF_data$Date[as.Date(DF_data$Date) == Sampling_Dates], DF_data$Chl_lab[as.Date(DF_data$Date) == Sampling_Dates], col = "Darkgreen", type = 'l', lwd = 3, ylim = c(0,450), 	
         ylab = expression("Chlorophyll" ~paste(italic( a))*"["~mu~"g L"^"-1"*"]"), xlab = "", xaxt = 'n', cex.lab = 1.5, cex.axis = 1.5)	
    points(DF_data$Date[as.Date(DF_data$Date) == Sampling_Dates], DF_data$Chl_lab[as.Date(DF_data$Date) == Sampling_Dates], cex = 1.7, bg = "black", pch = 16)	
    
    
    
    
    par(mar=c(5, 5, 0, 2))	
    plot(DF_data$`Average Sampling time`, DF_data$pH, type = 'l', col = "darkred", ylab = "pH", xlab = "Month [2016]", ylim = c(8,10.5),	
         cex.axis = 1.5, cex.lab = 1.5)	
    
    
    
  }	
  
  
  # Plot DO with Chlorophyll	
  
  
  
  
  par(mar=c(mar.2ysides) +0.5, xpd = FALSE)	
  
  plot(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$DO_mg_smooth * 0.0312511719189*1000, col = "blue", type = 'l', 	
       ylab = "DO ["~mu~"mol L"^"-1"*"]", xlab = "Months [2016]", cex.axis=1, cex.lab = 1, ylim = c(0, 1000))	
  
  rect(xleft = DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-05-23"], 	
       xright = DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-09-29"],	
       ybottom = -50, ytop = 1050, col = rgb(0,0,1, alpha = 0.03), border = NA)	
  abline(v= c(DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-05-23"], 	
              DF_data_macIntyre$`Average Sampling time`[DF_data_macIntyre$`Average Sampling time` == "2016-09-29"]), col = "grey", lty = 3)	
  lines(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$`DO_sat_umol/L`, col = "red" )	
  lines(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$DO_mg_smooth * 0.0312511719189*1000, col = "blue")	
  
  legend(	
    inset = c(0.05, -0.0),	
    "topleft",	
    cex = 1,	
    bty = "n",	
    lty = 1,	
    col = c("blue", "red", "darkgreen"),	
    legend = c("DO measured ", "DO at equilibrium", expression("Chl "~alpha~"")))	
  
  par(new = T)	
  with(DF_data, plot(DF_data$date, DF_data$Chl_lab, type = 'l', lwd = 2, axes=F, xlab=NA, ylab=NA, cex=1, col = "darkgreen",	
                     ylim = c(0,500)))	
  axis(side = 4)	
  mtext(side = 4, line = 3, "Chlorophyll "~alpha~ " ["~mu~"g L"^"-1"*"]")	
  
  
  {	
    # Plot DIC with Chlorophyll	
    
    layout(matrix(c(1,2),c(2)),	
           heights = c(2,1), respect = FALSE)	
    
    par(mar=c(1.5, 5, 3, 5), xpd = FALSE)	
    
    plot(DF_data_macIntyre$`Average Sampling time`, DF_data_macIntyre$DIC_Alk_smooth, col = "darkorange", type = 'l', 	
         ylab = "DIC ["~mu~"mol L"^"-1"*"]", xlab = "", xaxt='n', cex.axis=1.2, cex.lab = 1.3)	
    
    legend(	
      inset = c(0.05, -0.0),	
      "topright",	
      cex = 1.2,	
      bty = "n",	
      pch = c(NA,16),	
      lty = c(1,1),	
      col = c("darkorange", "darkgreen", "darkgreen"),	
      legend = c("DIC", expression("Chl "~alpha~"")))	
    
    par(new = T)	
    with(DF_data, plot(DF_data$date, DF_data$Chl_lab, type = 'l', lwd = 2, axes=F, xlab=NA, ylab=NA, cex=1, col = "darkgreen",	
                       ylim = c(0,500)))	
    points(DF_data$date[as.Date(DF_data$date) == as.Date(Sampling_Dates)], DF_data$Chl_lab[as.Date(DF_data$date) == as.Date(Sampling_Dates)], pch = 16, col = "darkgreen")	
    axis(side = 4, cex.axis = 1.2)	
    mtext(side = 4, line = 3, "Chlorophyll "~alpha~ " ["~mu~"g L"^"-1"*"]", cex = 1.3)	
    
    
    par(mar=c(4.5, 5, 0, 5))	
    plot(DF_data$`Average Sampling time`, DF_data$pH, type = 'l', col = "darkred", ylab = "pH", xlab = "Months [2016]",  ylim = c(8,10.5),	
         cex.axis=1.2, cex.lab = 1.3)	
    
  }	
  
  
  
}	

par(mfrow=c(1,1))	
mar.default <- c(5,4,4,2) + 0.1	
par(mar = mar.default + c(0, 2, 0, 0))	

















######################################################################################################################################################################3	
#	
# Use the ECS flux data	
#	
######################################################################################################################################################################	


plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$CO2FluxECS_umol_m2_h, type = 'l', ylab ="Rough ECS CO2 flux (umol/m2/h)", xlab = "Months (2016)")	
lines(DF_Metabolic$`Average Sampling time`, 0*DF_Metabolic$CO2FluxECS_umol_m2_h, col = "red")	

# add the flux from the vu	
DF_Metabolic$CO2flux_ECS_MA_day <-  Moving_average(as.double(DF_Metabolic$CO2FluxECS_umol_m2_h), Steps = 48)	

plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$CO2flux_ECS_MA_day, type = 'l',  ylab ="Half day moving average ECS CO2 flux (umol/m2/h)", xlab = "Months (2016)")	
lines(DF_Metabolic$`Average Sampling time`, DF_Metabolic$CO2flux_macIntyreR, col= "blue")	
abline(h=0, col = "red")	

plot(DF_Metabolic[c("CO2flux_macIntyreR", "CO2FluxECS_umol_m2_h", "CO2flux_ECS_MA_day")])	

plot(DF_Metabolic[c("wind_dir [deg_from_north]", "windECS_U10", "CO2FluxECS_umol_m2_h")])	



# Plot the flux in relation to the wind direction 	
{	
  plot(DF_Metabolic$`wind_dir [deg_from_north]`, DF_Metabolic$CO2FluxECS_umol_m2_h/1000,	
       ylab = expression("EC CO"[2]*" flux [mmol m"^"-2"*" h"^"-1"*"]"), xlab = "Wind direction in degrees from North", ylim = c(-3e+02, 3e+02), cex.axis  = 1.1, cex.lab = 1.1)	
  #lines(loess.smooth(DF_Metabolic$`wind_dir [deg_from_north]`, DF_Metabolic$CO2FluxECS_umol_m2_h)$x, loess.smooth(DF_Metabolic$`wind_dir [deg_from_north]`, DF_Metabolic$CO2FluxECS_umol_m2_h)$y, col = "red")	
  points(DF_Metabolic$`wind_dir [deg_from_north]`[DF_Metabolic$`wind_dir [deg_from_north]` > 260 & DF_Metabolic$`wind_dir [deg_from_north]` < 315 |	
                                                    DF_Metabolic$`wind_dir [deg_from_north]` > 173 & DF_Metabolic$`wind_dir [deg_from_north]` < 186],	
         DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$`wind_dir [deg_from_north]` > 260 & DF_Metabolic$`wind_dir [deg_from_north]` < 315 |	
                                             DF_Metabolic$`wind_dir [deg_from_north]` > 173 & DF_Metabolic$`wind_dir [deg_from_north]` < 186]/1000, col = "red")	
  points(DF_Metabolic$`wind_dir [deg_from_north]`[DF_Metabolic$`wind_dir [deg_from_north]` > 347 | DF_Metabolic$`wind_dir [deg_from_north]` < 27],	
         DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$`wind_dir [deg_from_north]` > 347 | DF_Metabolic$`wind_dir [deg_from_north]` < 27]/1000, col = "darkgreen")	
  abline(h=0, col = "red")	
  legend(	
    "bottomleft",	
    bty = 'n',	
    pch = 1,	
    col = c("darkgreen", "red"),	
    legend = c("Long fetch", "Short fetch"),	
    cex = 1.1	
  )	
  
  
  t.test_wind <- t.test(DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$`wind_dir [deg_from_north]` > 347 | DF_Metabolic$`wind_dir [deg_from_north]` < 27],	
                        DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$`wind_dir [deg_from_north]` > 260 & DF_Metabolic$`wind_dir [deg_from_north]` < 315 |	
                                                            DF_Metabolic$`wind_dir [deg_from_north]` > 173 & DF_Metabolic$`wind_dir [deg_from_north]` < 186])	
  
  
  p.adjust(t.test_wind[3], method = "bonferroni", n = length(DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$`wind_dir [deg_from_north]` > 347 | DF_Metabolic$`wind_dir [deg_from_north]` < 27]))	
  # there is a relation between wind direction and co2 flux in to the water.	
  # when the wind comes from the North, the open side from the lake the influx of CO2 is heighest.	
}	








###############################	
#	
# Estimate the k600	
#	
###############################	

{	
  # Convert the mol/L to umol/m2  -> Fc = ka(Cw - Ca)  so k=Fc/(a*(Cw-Ca))	
  DF_Metabolic$k_CO2flux_ECS <- DF_Metabolic$CO2FluxECS_umol_m2_h/(DF_Metabolic$`Enhancement factor`*(DF_Metabolic$`CO2 probe water mol/L`*1000000*1000*1.8 - DF_Metabolic$`CO2 atm (umol/L)`*1000000*1000*1.8))	
  
  # Multiply by 100 to convert from m/h to cm/h	
  DF_Metabolic$k600_CO2flux_ECS <- (DF_Metabolic$k_CO2flux_ECS)/((getSchmidt(DF_Metabolic$T_C_water, gas = "CO2")/600)^(-0.5))*100	
  
  hist(DF_Metabolic$k600_CO2flux_ECS, xlim = c(-100 , 100), breaks=10000)	
  
  plot(DF_Metabolic$`wind_dir [deg_from_north]`[DF_Metabolic$k600_CO2flux_ECS < 40 | DF_Metabolic$k600_CO2flux_ECS > -40], DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$k600_CO2flux_ECS < 40 | DF_Metabolic$k600_CO2flux_ECS > -40])	
  lines(loess.smooth(DF_Metabolic$`wind_dir [deg_from_north]`[DF_Metabolic$k600_CO2flux_ECS < 40 | DF_Metabolic$k600_CO2flux_ECS > -40],	
                     DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$k600_CO2flux_ECS < 40 | DF_Metabolic$k600_CO2flux_ECS > -40])$x,	
        loess.smooth(DF_Metabolic$`wind_dir [deg_from_north]`[DF_Metabolic$k600_CO2flux_ECS < 40 | DF_Metabolic$k600_CO2flux_ECS > -40],	
                     DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$k600_CO2flux_ECS < 40 | DF_Metabolic$k600_CO2flux_ECS > -40])$y, col = "red")	
  
  
  # Plot the filtered data	
  plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$CO2FluxECS_umol_m2_h/1000, 	
       ylab = "Atmospheric CO"[2]*" flux ["~mu~"mol m"^"-2"*" h"^"-1"*"]", xlab = "Time [Months 2016]", ylim = c(-3e+02, 3e+02), cex.axis = 1.1, cex.lab = 1.1)	
  points(DF_Metabolic$`Average Sampling time`[DF_Metabolic$k600_CO2flux_ECS > 40 | DF_Metabolic$k600_CO2flux_ECS < -40], DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$k600_CO2flux_ECS > 40 | DF_Metabolic$k600_CO2flux_ECS < -40]/1000, col = "red")	
  legend("topleft", 	
         inset = c(0.0,0.0), 	
         cex = 1.1, 	
         bty = "n", 	
         legend = c(expression("CO"[2]*" flux [k600 < 40 cm h"^"-1"*"]"), expression("CO"[2]*" flux [k600 > 40 cm h"^"-1"*"]")), 	
         text.col = c("black", "black"),	
         col = c("black", "red"),	
         pch = c(1,1)	
  )	
  
  
  
}	






length(DF_Metabolic$CO2FluxECS_filter_umol_m2_h[complete.cases(DF_Metabolic$CO2FluxECS_filter_umol_m2_h)])/length(DF_Metabolic$CO2FluxECS_umol_m2_h[complete.cases(DF_Metabolic$CO2FluxECS_umol_m2_h)])	






####################################################################################################################################3	
#	
# Calculate atmospheric stability	
#	
#########################################################################################################################################	









































































































































































































































































































{	
  
  MAcIntyre_Bfl <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                             RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
  {	
    S_B <- 5.67e-08	
    emiss <- 0.972	
    Kelvin = 273.15	
    albedo_SW <- 0.07	
    vonK <- 0.41	
    swRat <- 0.46	
    mnWnd <- 0.2	
    g <- 9.81	
    C_w <- 4186	
    rpcI <- wnd < mnWnd	
    wnd[rpcI] <- mnWnd	
    mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
    C_D <- mm$C_D	
    E <- mm$alh	
    H <- mm$ash	
    dUdt <- sw * 0.93 - E - H + lwnet	
    Qo <- sw * (1 - albedo_SW) * swRat	
    rho_w <- water.density(Ts)	
    if (wnd.z != 10) {	
      e1 <- sqrt(C_D)	
      u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
    }	
    else {	
      u10 = wnd	
    }	
    rhoAir <- 1.2	
    vonK <- 0.41	
    tau <- C_D * u10^2 * rhoAir	
    uSt <- sqrt(tau/rho_w)	
    q1 <- 2 - 2 * exp(z.aml * -Kd)	
    q2 <- z.aml * Kd	
    q3 <- exp(z.aml * -Kd)	
    H_star <- dUdt - Qo * (q1/q2 - q3)	
    thermalExpFromTemp <- function(Ts) {	
      V <- 1	
      dT <- 0.001	
      T1 <- water.density(Ts)	
      T2 <- water.density(Ts + dT)	
      V2 <- T1/T2	
      dv_dT <- (V2 - V)/dT	
      alpha <- dv_dT	
      return(alpha)	
    }	
    tExp <- thermalExpFromTemp(Ts)	
    B1 = H_star * tExp * g	
    B2 = rho_w * C_w	
    Bflx = B1/B2	
    getKinematicVis <- function(Ts) {	
      tempTable <- seq(0, 100, by = 5)	
      visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                    0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                    0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
      v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
      v <- v * 1e-06	
      return(v)	
    }	
    kinV <- getKinematicVis(Ts)	
    KeNm = uSt^3	
    SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
    SmE[SmE < 0] = 0	
    Sk = SmE * kinV	
    Sk = Sk * 100^4 * 3600^4	
    Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
    k600 <- as.numeric(Sk600)	
    k600 <- k600 * 24/100	
    return(Bflx)	
  }	
  
  
  MAcIntyre_KeNm <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                              RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
  {	
    S_B <- 5.67e-08	
    emiss <- 0.972	
    Kelvin = 273.15	
    albedo_SW <- 0.07	
    vonK <- 0.41	
    swRat <- 0.46	
    mnWnd <- 0.2	
    g <- 9.81	
    C_w <- 4186	
    rpcI <- wnd < mnWnd	
    wnd[rpcI] <- mnWnd	
    mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
    C_D <- mm$C_D	
    E <- mm$alh	
    H <- mm$ash	
    dUdt <- sw * 0.93 - E - H + lwnet	
    Qo <- sw * (1 - albedo_SW) * swRat	
    rho_w <- water.density(Ts)	
    if (wnd.z != 10) {	
      e1 <- sqrt(C_D)	
      u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
    }	
    else {	
      u10 = wnd	
    }	
    rhoAir <- 1.2	
    vonK <- 0.41	
    tau <- C_D * u10^2 * rhoAir	
    uSt <- sqrt(tau/rho_w)	
    q1 <- 2 - 2 * exp(z.aml * -Kd)	
    q2 <- z.aml * Kd	
    q3 <- exp(z.aml * -Kd)	
    H_star <- dUdt - Qo * (q1/q2 - q3)	
    thermalExpFromTemp <- function(Ts) {	
      V <- 1	
      dT <- 0.001	
      T1 <- water.density(Ts)	
      T2 <- water.density(Ts + dT)	
      V2 <- T1/T2	
      dv_dT <- (V2 - V)/dT	
      alpha <- dv_dT	
      return(alpha)	
    }	
    tExp <- thermalExpFromTemp(Ts)	
    B1 = H_star * tExp * g	
    B2 = rho_w * C_w	
    Bflx = B1/B2	
    getKinematicVis <- function(Ts) {	
      tempTable <- seq(0, 100, by = 5)	
      visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                    0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                    0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
      v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
      v <- v * 1e-06	
      return(v)	
    }	
    kinV <- getKinematicVis(Ts)	
    KeNm = uSt^3	
    SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
    SmE[SmE < 0] = 0	
    Sk = SmE * kinV	
    Sk = Sk * 100^4 * 3600^4	
    Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
    k600 <- as.numeric(Sk600)	
    k600 <- k600 * 24/100	
    return(KeNm)	
  }	
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  MAcIntyre_Sk <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                            RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
  {	
    S_B <- 5.67e-08	
    emiss <- 0.972	
    Kelvin = 273.15	
    albedo_SW <- 0.07	
    vonK <- 0.41	
    swRat <- 0.46	
    mnWnd <- 0.2	
    g <- 9.81	
    C_w <- 4186	
    rpcI <- wnd < mnWnd	
    wnd[rpcI] <- mnWnd	
    mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
    C_D <- mm$C_D	
    E <- mm$alh	
    H <- mm$ash	
    dUdt <- sw * 0.93 - E - H + lwnet	
    Qo <- sw * (1 - albedo_SW) * swRat	
    rho_w <- water.density(Ts)	
    if (wnd.z != 10) {	
      e1 <- sqrt(C_D)	
      u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
    }	
    else {	
      u10 = wnd	
    }	
    rhoAir <- 1.2	
    vonK <- 0.41	
    tau <- C_D * u10^2 * rhoAir	
    uSt <- sqrt(tau/rho_w)	
    q1 <- 2 - 2 * exp(z.aml * -Kd)	
    q2 <- z.aml * Kd	
    q3 <- exp(z.aml * -Kd)	
    H_star <- dUdt - Qo * (q1/q2 - q3)	
    thermalExpFromTemp <- function(Ts) {	
      V <- 1	
      dT <- 0.001	
      T1 <- water.density(Ts)	
      T2 <- water.density(Ts + dT)	
      V2 <- T1/T2	
      dv_dT <- (V2 - V)/dT	
      alpha <- dv_dT	
      return(alpha)	
    }	
    tExp <- thermalExpFromTemp(Ts)	
    B1 = H_star * tExp * g	
    B2 = rho_w * C_w	
    Bflx = B1/B2	
    getKinematicVis <- function(Ts) {	
      tempTable <- seq(0, 100, by = 5)	
      visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                    0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                    0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
      v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
      v <- v * 1e-06	
      return(v)	
    }	
    kinV <- getKinematicVis(Ts)	
    KeNm = uSt^3	
    SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
    SmE[SmE < 0] = 0	
    Sk = SmE * kinV	
    Sk = Sk * 100^4 * 3600^4	
    Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
    k600 <- as.numeric(Sk600)	
    k600 <- k600 * 24/100	
    return(Sk)	
  }	
  
  
  
  MAcIntyre_Sme <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                             RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
  {	
    S_B <- 5.67e-08	
    emiss <- 0.972	
    Kelvin = 273.15	
    albedo_SW <- 0.07	
    vonK <- 0.41	
    swRat <- 0.46	
    mnWnd <- 0.2	
    g <- 9.81	
    C_w <- 4186	
    rpcI <- wnd < mnWnd	
    wnd[rpcI] <- mnWnd	
    mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
    C_D <- mm$C_D	
    E <- mm$alh	
    H <- mm$ash	
    dUdt <- sw * 0.93 - E - H + lwnet	
    Qo <- sw * (1 - albedo_SW) * swRat	
    rho_w <- water.density(Ts)	
    if (wnd.z != 10) {	
      e1 <- sqrt(C_D)	
      u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
    }	
    else {	
      u10 = wnd	
    }	
    rhoAir <- 1.2	
    vonK <- 0.41	
    tau <- C_D * u10^2 * rhoAir	
    uSt <- sqrt(tau/rho_w)	
    q1 <- 2 - 2 * exp(z.aml * -Kd)	
    q2 <- z.aml * Kd	
    q3 <- exp(z.aml * -Kd)	
    H_star <- dUdt - Qo * (q1/q2 - q3)	
    thermalExpFromTemp <- function(Ts) {	
      V <- 1	
      dT <- 0.001	
      T1 <- water.density(Ts)	
      T2 <- water.density(Ts + dT)	
      V2 <- T1/T2	
      dv_dT <- (V2 - V)/dT	
      alpha <- dv_dT	
      return(alpha)	
    }	
    tExp <- thermalExpFromTemp(Ts)	
    B1 = H_star * tExp * g	
    B2 = rho_w * C_w	
    Bflx = B1/B2	
    getKinematicVis <- function(Ts) {	
      tempTable <- seq(0, 100, by = 5)	
      visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                    0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                    0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
      v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
      v <- v * 1e-06	
      return(v)	
    }	
    kinV <- getKinematicVis(Ts)	
    KeNm = uSt^3	
    SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
    SmE[SmE < 0] = 0	
    Sk = SmE * kinV	
    Sk = Sk * 100^4 * 3600^4	
    Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
    k600 <- as.numeric(Sk600)	
    k600 <- k600 * 24/100	
    return(SmE)	
  }	
  
  
  
  MAcIntyre_Ust <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                             RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
  {	
    S_B <- 5.67e-08	
    emiss <- 0.972	
    Kelvin = 273.15	
    albedo_SW <- 0.07	
    vonK <- 0.41	
    swRat <- 0.46	
    mnWnd <- 0.2	
    g <- 9.81	
    C_w <- 4186	
    rpcI <- wnd < mnWnd	
    wnd[rpcI] <- mnWnd	
    mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
    C_D <- mm$C_D	
    E <- mm$alh	
    H <- mm$ash	
    dUdt <- sw * 0.93 - E - H + lwnet	
    Qo <- sw * (1 - albedo_SW) * swRat	
    rho_w <- water.density(Ts)	
    if (wnd.z != 10) {	
      e1 <- sqrt(C_D)	
      u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
    }	
    else {	
      u10 = wnd	
    }	
    rhoAir <- 1.2	
    vonK <- 0.41	
    tau <- C_D * u10^2 * rhoAir	
    uSt <- sqrt(tau/rho_w)	
    q1 <- 2 - 2 * exp(z.aml * -Kd)	
    q2 <- z.aml * Kd	
    q3 <- exp(z.aml * -Kd)	
    H_star <- dUdt - Qo * (q1/q2 - q3)	
    thermalExpFromTemp <- function(Ts) {	
      V <- 1	
      dT <- 0.001	
      T1 <- water.density(Ts)	
      T2 <- water.density(Ts + dT)	
      V2 <- T1/T2	
      dv_dT <- (V2 - V)/dT	
      alpha <- dv_dT	
      return(alpha)	
    }	
    tExp <- thermalExpFromTemp(Ts)	
    B1 = H_star * tExp * g	
    B2 = rho_w * C_w	
    Bflx = B1/B2	
    getKinematicVis <- function(Ts) {	
      tempTable <- seq(0, 100, by = 5)	
      visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                    0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                    0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
      v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
      v <- v * 1e-06	
      return(v)	
    }	
    kinV <- getKinematicVis(Ts)	
    KeNm = uSt^3	
    SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
    SmE[SmE < 0] = 0	
    Sk = SmE * kinV	
    Sk = Sk * 100^4 * 3600^4	
    Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
    k600 <- as.numeric(Sk600)	
    k600 <- k600 * 24/100	
    return(uSt)	
  }	
  
  
  
  
  MAcIntyre_kinV <- function (wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd,	
                              RH, sw, lwnet, params = c(1.5600460, 0.2579, 0.5930818))	
  {	
    S_B <- 5.67e-08	
    emiss <- 0.972	
    Kelvin = 273.15	
    albedo_SW <- 0.07	
    vonK <- 0.41	
    swRat <- 0.46	
    mnWnd <- 0.2	
    g <- 9.81	
    C_w <- 4186	
    rpcI <- wnd < mnWnd	
    wnd[rpcI] <- mnWnd	
    mm <- calc.zeng(dateTime, Ts, airT, wnd, RH, atm.press, wnd.z)	
    C_D <- mm$C_D	
    E <- mm$alh	
    H <- mm$ash	
    dUdt <- sw * 0.93 - E - H + lwnet	
    Qo <- sw * (1 - albedo_SW) * swRat	
    rho_w <- water.density(Ts)	
    if (wnd.z != 10) {	
      e1 <- sqrt(C_D)	
      u10 <- wnd/(1 - e1/vonK * log(10/wnd.z))	
    }	
    else {	
      u10 = wnd	
    }	
    rhoAir <- 1.2	
    vonK <- 0.41	
    tau <- C_D * u10^2 * rhoAir	
    uSt <- sqrt(tau/rho_w)	
    q1 <- 2 - 2 * exp(z.aml * -Kd)	
    q2 <- z.aml * Kd	
    q3 <- exp(z.aml * -Kd)	
    H_star <- dUdt - Qo * (q1/q2 - q3)	
    thermalExpFromTemp <- function(Ts) {	
      V <- 1	
      dT <- 0.001	
      T1 <- water.density(Ts)	
      T2 <- water.density(Ts + dT)	
      V2 <- T1/T2	
      dv_dT <- (V2 - V)/dT	
      alpha <- dv_dT	
      return(alpha)	
    }	
    tExp <- thermalExpFromTemp(Ts)	
    B1 = H_star * tExp * g	
    B2 = rho_w * C_w	
    Bflx = B1/B2	
    getKinematicVis <- function(Ts) {	
      tempTable <- seq(0, 100, by = 5)	
      visTable <- c(1.792, 1.519, 1.308, 1.141, 1.007, 0.897,	
                    0.804, 0.727, 0.661, 0.605, 0.556, 0.513, 0.477,	
                    0.444, 0.415, 0.39, 0.367, 0.347, 0.328, 0.311, 0.296)	
      v <- data.frame(approx(tempTable, visTable, xout = Ts))[2]	
      v <- v * 1e-06	
      return(v)	
    }	
    kinV <- getKinematicVis(Ts)	
    KeNm = uSt^3	
    SmE = params[1] * -Bflx + params[2] * KeNm/(vonK * z.aml)	
    SmE[SmE < 0] = 0	
    Sk = SmE * kinV	
    Sk = Sk * 100^4 * 3600^4	
    Sk600 = params[3] * 600^(-0.5) * Sk^(1/4)	
    k600 <- as.numeric(Sk600)	
    k600 <- k600 * 24/100	
    return(kinV)	
  }	
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  DF_Metabolic$Mac_Intyre_Bfl <- MAcIntyre_Bfl(dateTime = DF_Metabolic$`Average Sampling time`, wnd.z = 10,	
                                               atm.press = DF_Metabolic$`Atm_pressure [PA]`/10, Ts = DF_Metabolic$T_C_water,	
                                               z.aml = DF_Metabolic$Zaml,	
                                               # z.aml = 0.15,	
                                               airT = DF_Metabolic$T_C_air, wnd = DF_Metabolic$windECS_U10,	
                                               RH = DF_Metabolic$`RH [%]`, sw =  DF_Metabolic$sw_PAR, lwnet =  DF_Metabolic$lw.net,	
                                               Kd = 0.0222177372 * DF_Metabolic$Chl_lab + 1.3858280092)	
  
  
  DF_Metabolic$Mac_Intyre_Kenm <- MAcIntyre_KeNm(dateTime = DF_Metabolic$`Average Sampling time`, wnd.z = 10,	
                                                 atm.press = DF_Metabolic$`Atm_pressure [PA]`/10, Ts = DF_Metabolic$T_C_water,	
                                                 z.aml = DF_Metabolic$Zaml,	
                                                 # z.aml = 0.15,	
                                                 airT = DF_Metabolic$T_C_air, wnd = DF_Metabolic$windECS_U10,	
                                                 RH = DF_Metabolic$`RH [%]`, sw =  DF_Metabolic$sw_PAR, lwnet =  DF_Metabolic$lw.net,	
                                                 Kd = 0.0222177372 * DF_Metabolic$Chl_lab + 1.3858280092)	
  
  
  
  DF_Metabolic$Mac_Intyre_uST <- MAcIntyre_Ust(dateTime = DF_Metabolic$`Average Sampling time`, wnd.z = 10,	
                                               atm.press = DF_Metabolic$`Atm_pressure [PA]`/10, Ts = DF_Metabolic$T_C_water,	
                                               z.aml = DF_Metabolic$Zaml,	
                                               # z.aml = 0.15,	
                                               airT = DF_Metabolic$T_C_air, wnd = DF_Metabolic$windECS_U10,	
                                               RH = DF_Metabolic$`RH [%]`, sw =  DF_Metabolic$sw_PAR, lwnet =  DF_Metabolic$lw.net,	
                                               Kd = 0.0222177372 * DF_Metabolic$Chl_lab + 1.3858280092)	
  
  
  
  plot(DF_Metabolic$`Average Sampling time`[DF_Metabolic$Date == "2016-06-01"], DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$Date == "2016-06-01"])	
  plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$Mac_Intyre_Bfl)	
  # The bouyancy flux is over all negative.	
  summary(DF_Metabolic$Mac_Intyre_Bfl)	
  sd(DF_Metabolic$Mac_Intyre_Bfl, na.rm = TRUE)/sqrt(length(DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$Mac_Intyre_Bfl != NA]))	
  
  
  
  sd(DF_Metabolic$Mac_Intyre_Bfl, na.rm = TRUE)/sqrt(length(DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$Mac_Intyre_Bfl != NA]))	
  
  
  
  DF_Metabolic$US <- 1.95/(DF_Metabolic$Mac_Intyre_Kenm/(0.41*DF_Metabolic$Mac_Intyre_Bfl))	
  DF_Metabolic$L <- DF_Metabolic$Mac_Intyre_Kenm/(0.41*DF_Metabolic$Mac_Intyre_Bfl)	
  
  plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$L)	
  plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$US)	
  
  
  
  
  
  
  plot(DF_Metabolic$`Average Sampling time`[DF_Metabolic$`Average Sampling time` > "2016-06-01" & DF_Metabolic$`Average Sampling time` < "2016-06-30"], DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$`Average Sampling time` > "2016-06-01" & DF_Metabolic$`Average Sampling time` < "2016-06-30"])	
  points(DF_Metabolic$`Average Sampling time`[DF_Metabolic$US < -0.55 | DF_Metabolic$US > 0.22], DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$US < -0.55 | DF_Metabolic$US > 0.22], col = "red")	
  
  
  
  # Remove data points for stability and k600	
  
  
  # Atmospgeric stability	
  quantile(DF_Metabolic$US, c(0.1,.15, 0.2, 0.25, 0.3, 0.35, 0.4, .50,  .75, 0.8, .85, 0.9,  .99))	
  
  
  DF_Metabolic$K600_filter_US <- DF_Metabolic$k600_CO2flux_ECS	
  
  # Removing stable conditions removes mostly winddirections from the residential area	
  
  DF_Metabolic$K600_filter_US[DF_Metabolic$US < -15 | DF_Metabolic$US > 0] = NA	
  DF_Metabolic$CO2FluxECS_filter_umol_m2_h <- DF_Metabolic$CO2FluxECS_umol_m2_h	
  DF_Metabolic$CO2FluxECS_filter_umol_m2_h[DF_Metabolic$US < -15 | DF_Metabolic$US > 0] = NA	
  DF_Metabolic$CO2FluxECS_filter_umol_L_h <- DF_Metabolic$CO2FluxECS_umol_m2_h/1800	
  DF_Metabolic$CO2FluxECS_filter_umol_L_h[DF_Metabolic$US < -15 | DF_Metabolic$US > 0] = NA	
  #plot(DF_Metabolic$`wind_dir [deg_from_north]`,DF_Metabolic$K600_filter_US)	
  
  
  
  # Remove extreme K600 values	
  
  DF_Metabolic$K600_filter_US[DF_Metabolic$k600_CO2flux_ECS > 40 | DF_Metabolic$k600_CO2flux_ECS < -40] = NA	
  DF_Metabolic$CO2FluxECS_filter_umol_m2_h[DF_Metabolic$k600_CO2flux_ECS > 40 | DF_Metabolic$k600_CO2flux_ECS < -40] = NA	
  DF_Metabolic$CO2FluxECS_filter_umol_L_h[DF_Metabolic$k600_CO2flux_ECS > 40 | DF_Metabolic$k600_CO2flux_ECS < -40] = NA	
  
  DF_Metabolic$CO2flux_macIntyreR[DF_Metabolic$CO2flux_macIntyreR == 0] = NA	
  
  
  # only 3906 points lay outside the range of k600 < -50 and > 50, this is 27.03 %	
  
  DF_Metabolic$CO2flux_ECS_MA_day <-  Moving_average(DF_Metabolic$CO2FluxECS_filter_umol_m2_h, Steps = 2)	
  
  plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$CO2FluxECS_filter_umol_m2_h, type = 'l' )	
  lines(DF_Metabolic$`Average Sampling time`, DF_Metabolic$CO2flux_macIntyreR, col= "red")	
  
  plot(DF_Metabolic$CO2FluxECS_filter_umol_m2_h~DF_Metabolic$CO2flux_macIntyreR, ylim = c(-10000, 10000))	
  summary(lm(DF_Metabolic$CO2flux_ECS_MA_day~DF_Metabolic$CO2flux_macIntyreR))	
  
  hist(DF_Metabolic$CO2FluxECS_filter_umol_m2_h, breaks = 100)	
  
  plot(DF_Metabolic$`wind_dir [deg_from_north]`, DF_Metabolic$K600_filter_US)	
  
  sum(count(DF_Metabolic$k600_CO2flux_ECS[DF_Metabolic$k600_CO2flux_ECS > 40 | DF_Metabolic$k600_CO2flux_ECS < -40])$freq)/nrow(DF_Metabolic)	
  sum(count(DF_Metabolic$k600_CO2flux_ECS[DF_Metabolic$US < -15 | DF_Metabolic$US > 0])$freq)/nrow(DF_Metabolic$k600_CO2flux_ECS)	
  
  sum(count(DF_Metabolic$CO2FluxECS_filter_umol_m2_h)$freq)/sum(count(DF_data$CO2FluxECS_umol_m2_h)$freq)	
  
  
  
  
  
  
  plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$CO2FluxECS_umol_m2_h/1000, 	
       ylab = expression("Atmospheric CO"[2]*" flux [mmol m"^"-2"*" h"^"-1"*"]"), xlab = "Time [Months 2016]", ylim = c(-3e+02,3e+02), cex.lab = 1.1, cex.axis = 1.1)	
  points(DF_Metabolic$`Average Sampling time`[DF_Metabolic$US < -15 | DF_Metabolic$US > 0], DF_Metabolic$CO2FluxECS_umol_m2_h[DF_Metabolic$US < -15 | DF_Metabolic$US > 0]/1000, col = "red")	
  legend("topleft", 	
         inset = c(0.0,0.0), 	
         cex = 1.1, 	
         bty = "n", 	
         legend = c(expression("CO"[2]*" flux"), expression("CO"[2]*" flux [AS < -15, AS > 0]")), 	
         text.col = c("black", "black"),	
         col = c("black", "red"),	
         pch = c(1,1))	
  
}	


library(i)	
plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$CO2FluxECS_filter_umol_m2_h, col = rgb(0,0,0), pch = 149,	
     ylab = "Atmospheric CO"[2]*" flux ["~mu~"mol m"^"-2"*" h"^"-1"*"]", xlab = "Months [2016]")	
points(DF_Metabolic$`Average Sampling time`, DF_Metabolic$CO2flux_macIntyreR, col = "blue", pch = 149)	
abline(h = 0, col = rgb(1,0,0,0.7))	
legend("topleft", 	
       inset = c(0.0,0.0), 	
       cex = 1.1, 	
       bty = "n", 	
       legend = c(expression("CO"[2]*" flux ECS"), expression("CO"[2]*" flux model")), 	
       text.col = c("black", "black"),	
       col = c("black", "blue"),	
       pch = c(149))	




# The total CO2 in flux into the lake	
{	
  
  # The average flux per mmol/m2/d	
  summary(DF_Metabolic$CO2FluxECS_filter_umol_m2_h)/1000*24	
  mean(DF_Metabolic$CO2FluxECS_filter_umol_m2_h, na.rm = TRUE)/1000*24	
  sd(DF_Metabolic$CO2FluxECS_filter_umol_m2_h, na.rm = TRUE)/sqrt(nrow(DF_Metabolic))/1000*24	
  # The total flux umol/period sampling	
  mean(DF_Metabolic$CO2FluxECS_filter_umol_m2_h/1000000, na.rm = TRUE) * as.numeric(DF_Metabolic$`Average Sampling time`[nrow(DF_Metabolic)] - DF_Metabolic$`Average Sampling time`[1])*24*541120 	
  sd(DF_Metabolic$CO2FluxECS_filter_umol_m2_h/1000000, na.rm = TRUE )/sqrt(nrow(DF_Metabolic)) * as.numeric(DF_Metabolic$`Average Sampling time`[nrow(DF_Metabolic)] - DF_Metabolic$`Average Sampling time`[1])*24*541120	
  
  
  
  # The average flux per umol/m2/d	
  summary(DF_Metabolic$CO2flux_macIntyreR)/1000*24	
  mean(DF_Metabolic$CO2flux_macIntyreR, na.rm = TRUE)/1000*24	
  sd(DF_Metabolic$CO2flux_macIntyreR, na.rm = TRUE)/sqrt(nrow(DF_Metabolic))/1000*24	
  
  # The total flux umol/period sampling	
  mean(DF_Metabolic$CO2flux_macIntyreR, na.rm = TRUE) * as.numeric(DF_Metabolic$`Average Sampling time`[nrow(DF_Metabolic)] - DF_Metabolic$`Average Sampling time`[1])*24/1000000*541120 	
  sd(DF_Metabolic$CO2flux_macIntyreR/1000000, na.rm = TRUE )/sqrt(nrow(DF_Metabolic)) * as.numeric(DF_Metabolic$`Average Sampling time`[nrow(DF_Metabolic)] - DF_Metabolic$`Average Sampling time`[1])*24*541120	
  
  
  # O2	
  
  # The average flux per umol/m2/d	
  summary(DF_Metabolic$DOflux_umol)/1000*24	
  mean(DF_Metabolic$DOflux_umol, na.rm = TRUE)/1000*24	
  sd(DF_Metabolic$DOflux_umol, na.rm = TRUE)/sqrt(nrow(DF_Metabolic))/1000*24	
  
  # The total flux umol/period sampling	
  mean(DF_Metabolic$DOflux_umol, na.rm = TRUE) * as.numeric(DF_Metabolic$`Average Sampling time`[nrow(DF_Metabolic)] - DF_Metabolic$`Average Sampling time`[1])*24/1000000000*541120 	
  sd(DF_Metabolic$DOflux_umol/1000000, na.rm = TRUE )/sqrt(nrow(DF_Metabolic)) * as.numeric(DF_Metabolic$`Average Sampling time`[nrow(DF_Metabolic)] - DF_Metabolic$`Average Sampling time`[1])*24*541120	
  
  
  
}	


#############################################################################################################################	
#	
# Estimate the constants from the Mac Intyre model based on the ECS data	
#	
#############################################################################################################################	

{	
  
  DF_Metabolic_ECS <- DF_Metabolic[complete.cases(DF_Metabolic["CO2FluxECS_filter_umol_m2_h"]),]	
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  DF_Metabolic_ECS$Mac_Intyre_SK <- MAcIntyre_Sk(dateTime = DF_Metabolic_ECS$`Average Sampling time`, wnd.z = 10,	
                                                 atm.press = DF_Metabolic_ECS$`Atm_pressure [PA]`/10, Ts = DF_Metabolic_ECS$T_C_water,	
                                                 z.aml = DF_Metabolic_ECS$Zaml,	
                                                 # z.aml = 0.15,	
                                                 airT = DF_Metabolic_ECS$T_K_air-273.15, wnd = DF_Metabolic_ECS$windECS_U10,	
                                                 RH = DF_Metabolic_ECS$`RH [%]`, sw =  DF_Metabolic_ECS$sw_PAR, lwnet =  DF_Metabolic_ECS$lw.net,	
                                                 Kd = 0.0222177372 * DF_Metabolic_ECS$Chl_lab + 1.3858280092)	
  
  DF_Metabolic_ECS$Mac_Intyre_SmE <- MAcIntyre_Sme(dateTime = DF_Metabolic_ECS$`Average Sampling time`, wnd.z = 10,	
                                                   atm.press = DF_Metabolic_ECS$`Atm_pressure [PA]`/10, Ts = DF_Metabolic_ECS$T_C_water,	
                                                   z.aml = DF_Metabolic_ECS$Zaml,	
                                                   # z.aml = 0.15,	
                                                   airT = DF_Metabolic_ECS$T_K_air-273.15, wnd = DF_Metabolic_ECS$windECS_U10,	
                                                   RH = DF_Metabolic_ECS$`RH [%]`, sw =  DF_Metabolic_ECS$sw_PAR, lwnet =  DF_Metabolic_ECS$lw.net,	
                                                   Kd = 0.0222177372 * DF_Metabolic_ECS$Chl_lab + 1.3858280092, params = c(1.5600460, 0.2579, 0.5930818))	
  
  
  DF_Metabolic_ECS$Mac_Intyre_kinV <- MAcIntyre_kinV(dateTime = DF_Metabolic_ECS$`Average Sampling time`, wnd.z = 10,	
                                                     atm.press = DF_Metabolic_ECS$`Atm_pressure [PA]`/10, Ts = DF_Metabolic_ECS$T_C_water,	
                                                     z.aml = DF_Metabolic_ECS$Zaml,	
                                                     # z.aml = 0.15,	
                                                     airT = DF_Metabolic_ECS$T_K_air-273.15, wnd = DF_Metabolic_ECS$windECS_U10,	
                                                     RH = DF_Metabolic_ECS$`RH [%]`, sw =  DF_Metabolic_ECS$sw_PAR, lwnet =  DF_Metabolic_ECS$lw.net,	
                                                     Kd = 0.0222177372 * DF_Metabolic_ECS$Chl_lab + 1.3858280092)	
  
  
  
  Flux <- DF_Metabolic_ECS$CO2FluxECS_filter_umol_m2_h	
  Bfl <- DF_Metabolic_ECS$Mac_Intyre_Bfl	
  Kenm <- DF_Metabolic_ECS$Mac_Intyre_Kenm	
  Zaml <- DF_Metabolic_ECS$Zaml	
  SmE <- DF_Metabolic_ECS$Mac_Intyre_SmE	
  kinV <- DF_Metabolic_ECS$Mac_Intyre_kinV	
  k600 <- DF_Metabolic_ECS$K600_filter_US	
  
  
  
  
  coef(nls(SmE ~ C2 * -Bfl + C3 * Kenm/(0.41 * Zaml),  start=list(C2=0.7, C3=0.3)))	
  
  #        C2        C3	
  #  0.5930818    1.5600460	
  
  
  # estimate C1	
  
  SmE[SmE < 0] = 0	
  Sk = SmE * kinV	
  Sk = Sk * 100^4 * 3600^4	
  
  
  # The K600 is given in cm/h	
  
  library(minpack.lm)	
  
  
  
  
  nlsLM(k600 ~ C1 * 600^(-0.5) * ((0.5930818 * -Bfl + 1.5600460 * Kenm/(0.41 * Zaml)) * kinV * 100^4 * 3600^4)^(1/4),  start=list(C1 = 1.5))	
  
  # C1	
  # 0.2579	
  
  
  mean(0.259 * 600^(-0.5) * ((0.5930818 * -Bfl + 1.5600460 * Kenm/(0.41 * Zaml)) * kinV * 100^4 * 3600^4)^(1/4))	
  
  # C1 0.2579	
  # C2 0.5930818	
  # C3 1.5600460	
  
  
  
  # with the mistake in the code of Winslow the vector will be:	
  # c(1.5600460, 0.2579, 0.5930818)	
  
}	



#########################################################################################################################3	
#	
# buoyancy flux, wind and k600 to compare with MacIntyre 2010	
#	
############################################################################################################################	
{	
  DF_wind <- data.frame(DF_Metabolic$`Average Sampling time`)	
  DF_wind$windECS_U10 <- DF_Metabolic$windECS_U10	
  DF_wind$Mac_Intyre_Bfl <- DF_Metabolic$Mac_Intyre_Bfl	
  DF_wind$K600_filter_US <- DF_Metabolic$K600_filter_US	
  DF_wind$K600_MI <- DF_Metabolic$k600.macIntyre/24*100	
  
  mean(DF_wind$K600_filter_US, na.rm = TRUE)	
  sd(DF_wind$K600_filter_US, na.rm = TRUE)/sqrt(length(DF_wind$windECS_U10))	
  
  k600_wind_B <- data.frame(c(1:21))	
  k600_wind_B$wind <- seq(length=21, from=0, to=10)	
  
  
  ind_a <- seq(length=19, from=0.5, to=9.5)	
  ind_b <- seq(length=19, from=1, to=10)	
  ind_c <- seq(length=21, from=1, to=21)	
  
  for( i in seq_along(ind_c)) {	
    k600_wind_B$k600Bneg[ind_c[i]] <- mean(DF_wind$K600_filter_US[DF_wind$windECS_U10 > ind_a[i] & DF_wind$windECS_U10 < ind_b[i] &	
                                                                    DF_wind$Mac_Intyre_Bfl < 0], na.rm = TRUE)	
  }	
  
  for( i in seq_along(ind_c)) {	
    k600_wind_B$k600BnegSD[ind_c[i]] <- sd(DF_wind$K600_filter_US[DF_wind$windECS_U10 > ind_a[i] & DF_wind$windECS_U10 < ind_b[i] &	
                                                                    DF_wind$Mac_Intyre_Bfl < 0], na.rm = TRUE)/sqrt(length(DF_wind$K600_filter_US[DF_wind$windECS_U10 > ind_a[i] & DF_wind$windECS_U10 < ind_b[i] &	
                                                                                                                                                    DF_wind$Mac_Intyre_Bfl < 0]))	
  }	
  
  for( i in seq_along(ind_c)) {	
    k600_wind_B$k600Bpos[ind_c[i]] <- mean(DF_wind$K600_filter_US[DF_wind$windECS_U10 > ind_a[i] & DF_wind$windECS_U10 < ind_b[i] &	
                                                                    DF_wind$Mac_Intyre_Bfl > 0], na.rm = TRUE)	
  }	
  
  
  for( i in seq_along(ind_c)) {	
    k600_wind_B$k600BnegMI[ind_c[i]] <- mean(DF_wind$K600_MI[DF_wind$windECS_U10 > ind_a[i] & DF_wind$windECS_U10 < ind_b[i] &	
                                                               DF_wind$Mac_Intyre_Bfl < 0], na.rm = TRUE)	
  }	
  
  
  for( i in seq_along(ind_c)) {	
    k600_wind_B$k600BnegMISD[ind_c[i]] <- sd(DF_wind$K600_MI[DF_wind$windECS_U10 > ind_a[i] & DF_wind$windECS_U10 < ind_b[i] &	
                                                               DF_wind$Mac_Intyre_Bfl < 0], na.rm = TRUE)/sqrt(length(DF_wind$K600_MI[DF_wind$windECS_U10 > ind_a[i] & DF_wind$windECS_U10 < ind_b[i] &	
                                                                                                                                        DF_wind$Mac_Intyre_Bfl < 0]))	
  }	
  
  for( i in seq_along(ind_c)) {	
    k600_wind_B$k600BposMI[ind_c[i]] <- mean(DF_wind$K600_MI[DF_wind$windECS_U10 > ind_a[i] & DF_wind$windECS_U10 < ind_b[i] &	
                                                               DF_wind$Mac_Intyre_Bfl > 0])	
  }	
  
  for( i in seq_along(ind_c)) {	
    k600_wind_B$k600BposMISD[ind_c[i]] <- sd(DF_wind$K600_MI[DF_wind$windECS_U10 > ind_a[i] & DF_wind$windECS_U10 < ind_b[i] &	
                                                               DF_wind$Mac_Intyre_Bfl > 0], na.rm = TRUE)/sqrt(length(DF_wind$K600_MI[DF_wind$windECS_U10 > ind_a[i] & DF_wind$windECS_U10 < ind_b[i] &	
                                                                                                                                        DF_wind$Mac_Intyre_Bfl > 0]))	
  }	
  
  
  Jonsson2008 <- data.frame(c(0:5))	
  Jonsson2008$K600 <- c(-4.36, 1.744186, 5.5232558, 6.8313953, 7.5581395, 10.319767)	
  Jonsson2008$wind <- c(0.658307, 1.523511, 2.501567, 3.517241, 4.401254, 5.435737)	
  
  MacIntyre2010 <- data.frame(c(0:5))	
  MacIntyre2010$wind <- c(1,2,3,4,5,6)	
  MacIntyre2010$windNegB <- c(2.583732, 6.685321, 8.82253, 11.71886, 9.838231, 14.52018)	
  MacIntyre2010$windPosB <- c(NA, 5.167463, 7.036909, 6.584857, 9.079303, 8.761249)	
  
  par(mar = mar.default + 1)	
  # Plot the k600 against the wind	
  plot(k600_wind_B$wind, k600_wind_B$k600Bneg, col = "darkblue", xlim = c(0.5,6), ylim = c(-5,17), pch = 16,	
       ylab = expression("K"[600]*" [cm h"^"-1"*"]"), xlab = expression("Wind speed [m s"^"-1"*"]"), cex.axis = 1.4, cex.lab = 1.5)	
  lines(k600_wind_B$wind[2:13], k600_wind_B$k600Bneg[2:13], col = "darkblue")	
  points(k600_wind_B$wind, k600_wind_B$k600Bpos, col = "darkred", pch = 16)	
  lines(k600_wind_B$wind[1:13], k600_wind_B$k600Bpos[1:13], col = "darkred")	
  points(Jonsson2008$wind, Jonsson2008$K600, col = "black", pch = 7)	
  lines(Jonsson2008$wind, Jonsson2008$K600, col = "black", lty = 1)	
  points(MacIntyre2010$wind, MacIntyre2010$windNegB, col = "darkred", pch =16)	
  lines(MacIntyre2010$wind, MacIntyre2010$windNegB, col = "darkred", lty = 1)	
  points(k600_wind_B$wind, k600_wind_B$k600BnegMI, col = "darkblue")	
  lines(k600_wind_B$wind[2:13], k600_wind_B$k600BnegMI[2:13], col = "darkblue", lty = 3)	
  points(k600_wind_B$wind, k600_wind_B$k600BposMI, col = "darkblue", pch = 3, cex = 1.7)	
  lines(k600_wind_B$wind[2:13], k600_wind_B$k600BposMI[2:13], col = "darkblue", lty = 3)	
  points(MacIntyre2010$wind, MacIntyre2010$windPosB, col = "darkred", pch = 3, cex = 1.7)	
  lines(MacIntyre2010$wind, MacIntyre2010$windPosB, col = "darkred", lty = 1)	
  legend("topleft", 	
         inset = c(0.02,0.0), 	
         cex = 1.1, 	
         bty = "n", 	
         legend = c(expression("Measured K"[600]*" ["~beta~"flux < 0]"), expression("Modeled K"[600]*" ["~beta~"flux < 0]"), expression("Modeled K"[600]*" ["~beta~"flux > 0]"), "Jonsson et al., 2008","MacIntyre et al., 2010 ["~beta~"flux < 0]", "MacIntyre et al., 2010 ["~beta~"flux > 0]"), 	
         text.col = "black",	
         col = c("darkblue", "darkblue", "darkblue", "black", "darkred", "darkred"),	
         pch = c(16,1,3,7,16,3), lty = c(1,3,3,1,1,1))	
  
  
  t.test(DF_wind$K600_MI[DF_wind$Mac_Intyre_Bfl < 0 & DF_wind$windECS_U10 <= 6 & DF_wind$windECS_U10 >= 1], DF_wind$K600_MI[DF_wind$Mac_Intyre_Bfl > 0 & DF_wind$windECS_U10 <= 6 & DF_wind$windECS_U10 >= 1])	
  
  
  summary(k600_wind_B$k600Bneg[k600_wind_B$wind < 7])	
  sd(k600_wind_B$k600Bneg[k600_wind_B$wind < 7])/sqrt(6)	
  summary(Jonsson2008$K600)	
  summary(MacIntyre2010$windNegB)	
  
  hist(DF_Metabolic$`WindSpeed (ECS)`, labels = TRUE)	
  
}	




########################################################################################################################################3	
#	
# Check the median and mean values for the ECS flux	
#	
#########################################################################################################################################	

# Mean and median values	

library(imputeTS)	
DF_Metabolic$CO2FluxECS_filter_umol_m2_h_kalman <- na.kalman(DF_Metabolic$CO2FluxECS_filter_umol_m2_h, model = "StructTS", smooth = TRUE)	


{	
  
  DF_Metabolic_ECS_median <- DF_Metabolic[c("Average Sampling time", "Date", "CO2flux_ColeProbeR","CO2flux_macIntyreR", "CO2FluxECS_filter_umol_m2_h", "CO2FluxECS_filter_umol_m2_h_kalman",	
                                            "ChlVoltNorm", "Chl_lab", "WindSpeed (ECS)", "T_C_water", "Mac_Intyre_Bfl", "pH", "DIC_Alk_smooth", "DOflux_umol")]	
  DF_Metabolic_ECS_median <- DF_Metabolic_ECS_median[complete.cases(DF_Metabolic_ECS_median[c("CO2flux_macIntyreR")]),]	
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  library(dplyr)	
  DF_Metabolic_ECS_median %>%	
    group_by(Date) %>%	
    summarise(median_flux_cole = median(CO2flux_ColeProbeR),	
              median_flux_MacIntyre = median(CO2flux_macIntyreR, na.rm = TRUE),	
              median_flux_MacIntyre_Alk = median(CO2flux_macIntyreR, na.rm = TRUE),	
              median_flux_ECS = median(CO2FluxECS_filter_umol_m2_h, na.rm = TRUE),	
              mean_flux_MacIntyre = mean(CO2flux_macIntyreR, na.rm = TRUE),	
              mean_flux_MacIntyre_Alk = mean(CO2flux_macIntyreR, na.rm = TRUE),	
              mean_flux_ECS = mean(CO2FluxECS_filter_umol_m2_h, na.rm = TRUE),	
              mean_flux_ECS_kal = mean(CO2FluxECS_filter_umol_m2_h_kalman, na.rm = TRUE),	
              mean_flux_DO = mean(DOflux_umol, na.rm = TRUE),	
              median_DIC_alk = median(DIC_Alk_smooth, na.rm = TRUE),	
              median_Chl = median(ChlVoltNorm, na.rm = TRUE),	
              median_Chl_lab = median(Chl_lab),	
              mean_wind = mean(`WindSpeed (ECS)`, na.rm = TRUE),	
              mean_BF = mean(Mac_Intyre_Bfl, na.rm = TRUE),	
              median_pH = median(pH, na.rm = TRUE),	
              median_temp = median(T_C_water)) -> Flux	
  detach("package:dplyr", unload = TRUE)	
  
  
  
  Flux$Date <- as.POSIXct(Flux$Date)	
  
  write.table(Flux, "C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Routput/Flux.txt")	
  
  
  par(mar = mar.default +1)	
  plot(Flux$Date, Flux$mean_flux_ECS/1000*24, ylab = expression("Daily averaged flux [mmol m"^"-2"*" d"^"-1"*"]"), xlab = "Months [2016]", pch = 19, 	
       col = rgb(0,0,0,0.9), type = 'l', lwd = 2, ylim = c(-500, 500), cex.axis = 1.4, cex.lab = 1.4)	
  #points(Flux$Date_dusk, Flux$median_flux_cole, col = "red")	
  points(Flux$Date, Flux$mean_flux_DO/1000*24, col = "dodgerblue1", pch = 15)	
  points(Flux$Date[Flux$mean_flux_MacIntyre != 0], Flux$mean_flux_MacIntyre[Flux$mean_flux_MacIntyre != 0]/1000*24, col = rgb(0,0,1,0.7), pch = 19)	
  lines(Flux$Date, Flux$mean_flux_ECS/1000*24, ylab = expression("mean CO"[2]*" flux [mmol m"^"-2"*" d"^"-1"*"]"), xlab = "Months [2016]", pch = 19, 	
        col = rgb(0,0,0,0.9), type = 'l', lwd = 2, ylim = c(-500, 500), cex.axis = 1.4, cex.lab = 1.4)	
  abline(h = 0, col = "red", lwd = 2, lty = 3)	
  legend(	
    "bottomright",	
    bty = "n",	
    legend = c(expression("CO"[2]*" Flux ECS"), expression("CO"[2]*" Flux model"), expression("O"[2]*" Flux model")),	
    cex = 1.4,	
    pch = c(NA,19,15), lty = c(1,NA,NA), lwd(NA, 2,2),	
    col = c(rgb(0,0,0,0.7), rgb(0,0,1,0.7), "steelblue3")	
    
  )	
  
  plot(Flux$Date, Flux$mean_flux_ECS_kal*24/1000, ylab = expression("CO"[2]*" flux [mmol m"^"-2"*" d"^"-1"*"]"), xlab = "Months [2016]", 	
       pch = 19, col = rgb(0,0,0,0.9), type = 'l', lwd = 2, ylim = c(-500, 500), cex.lab = 1.2, cex.axis = 1.2)	
  #points(Flux$Date_dusk, Flux$median_flux_cole, col = "red")	
  points(Flux$Date[Flux$mean_flux_MacIntyre != 0], Flux$mean_flux_MacIntyre[Flux$mean_flux_MacIntyre != 0]*24/1000, col = rgb(0,0,1,0.7), pch = 19)	
  abline(h = 0, col = "red", lwd = 2, lty = 3)	
  legend(	
    "topleft",	
    bty = "n",	
    cex = 1.1,	
    legend = c(expression("CO"[2]*" Flux ECS"), expression("CO"[2]*" Flux model")),	
    pch = c(NA,19), lty = c(1,NA), lwd(NA, 2),	
    col = c(rgb(0,0,0,0.7), rgb(0,0,1,0.7))	
    
  )	
  
  # Check which months are highest and lowest	
  
  library(dplyr)	
  DF_Metabolic_ECS_median %>%	
    group_by(months(as.POSIXct(Date))) %>%	
    summarise(median_flux_cole = median(CO2flux_ColeProbeR),	
              median_flux_MacIntyre = median(CO2flux_macIntyreR, na.rm = TRUE),	
              median_flux_ECS = median(CO2FluxECS_filter_umol_m2_h, na.rm = TRUE),	
              mean_flux_ECS = mean(CO2FluxECS_filter_umol_m2_h, na.rm = TRUE),	
              mean_flux_ECS_kal = mean(CO2FluxECS_filter_umol_m2_h_kalman, na.rm = TRUE),	
              mean_flux_MacIntyre = mean(CO2flux_macIntyreR, na.rm = TRUE),	
              s.e.flux_MacIntyre = sd(CO2flux_macIntyreR, na.rm = TRUE)/sqrt(length(CO2flux_macIntyreR)),	
              mean_O2flux_MacIntyre = mean(DOflux_umol, na.rm = TRUE),	
              s.e.O2flux_MacIntyre = sd(DOflux_umol, na.rm = TRUE)/sqrt(length(DOflux_umol)),	
              mean_flux_ECS = mean(CO2FluxECS_filter_umol_m2_h, na.rm = TRUE),	
              s.e.flux_ECS = sd(CO2FluxECS_filter_umol_m2_h, na.rm = TRUE)/sqrt(length(CO2FluxECS_filter_umol_m2_h)),	
              median_temp = median(T_C_water),	
              mean_Chl = mean(Chl_lab)) -> Flux_month	
  
  Flux_month[c("months(as.POSIXct(Date))","mean_flux_ECS", "s.e.flux_ECS")]	
  Flux_month[c("months(as.POSIXct(Date))","mean_flux_MacIntyre", "s.e.flux_MacIntyre")]	
  plot(Flux_month$mean_flux_ECS_kal)	
  
  Flux_month <-  Flux_month[order(match(Flux_month$`months(as.POSIXct(Date))`,month.name)),]	
  
  Flux_month$`months(as.POSIXct(Date))`	
  Flux_month[,c(5:10)]/1000*24	
  
  
  length(Flux$Date[Flux$mean_flux_ECS > 0])	
  length(Flux$Date)	
}	



##################################################################################################################################	
#	
#  Check the regression models for the ECS flux	
#	
###################################################################################################################################	

{	
  ##################################################################################################	
  
  
  
  
  # Check the relation between both fluxes	
  
  ## Median	
  # Check for the modeled flux based on the probe data	
  par(mfrow=c(2,2))	
  plot(lm(Flux$median_flux_ECS ~ Flux$median_flux_cole))	
  par(mfrow=c(1,1))	
  
  summary(lm(Flux$median_flux_ECS ~ Flux$median_flux_cole))	
  
  plot(Flux$mean_flux_ECS ~ Flux$mean_flux_MacIntyre)	
  
  # Mean	
  
  # Check for the modeled flux based on the probe data	
  
  par(mfrow=c(2,2))	
  plot(lm(Flux$mean_flux_ECS ~  Flux$mean_flux_MacIntyre))	
  par(mfrow=c(1,1))	
  
  
  par(mfrow=c(2,2))	
  plot(lm(Flux$mean_flux_ECS ~ sign(Flux$mean_flux_MacIntyre)*Flux$mean_flux_MacIntyre^2))	
  par(mfrow=c(1,1))	
  
  par(mfrow=c(2,2))	
  plot(lm(sign(Flux$mean_flux_ECS)*log(abs(Flux$mean_flux_ECS)) ~ Flux$mean_flux_MacIntyre))	
  par(mfrow=c(1,1))	
  
  summary(lm(Flux$mean_flux_ECS ~ Flux$mean_flux_MacIntyre))	
  summary(lm(Flux$mean_flux_ECS ~ 0 + Flux$mean_flux_MacIntyre))	
  
  # Check for the modeled flux based on the probe data	
  par(mfrow=c(2,2))	
  plot(lm(Flux$mean_flux_ECS ~ Flux$mean_flux_MacIntyre_Alk))	
  par(mfrow=c(1,1))	
  
  summary(lm(Flux$mean_flux_ECS ~ Flux$mean_flux_MacIntyre_Alk))	
  summary(lm(Flux$mean_flux_ECS ~ 0 + Flux$mean_flux_MacIntyre_Alk))	
  
  
  # The mean modeled fluxes are identical 	
  summary(lm(Flux$mean_flux_MacIntyre ~ Flux$mean_flux_MacIntyre_Alk))	
  
  # The linear regression between the median values for the ECS and the model MacIntyre flux	
  
  
  #Coefficients:	
  #  Estimate Std. Error t value Pr(>|t|)    	
  # Flux$median_flux_MacIntyre[-c(127, 23, 22, 21)]  1.19679    0.07369   16.24   <2e-16 ***	
  #  ---	
  #  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1	
  
  # Residual standard error: 1976 on 124 degrees of freedom	
  # (48 observations deleted due to missingness)	
  # Multiple R-squared:  0.6802,	Adjusted R-squared:  0.6777 
  # F-statistic: 263.8 on 1 and 124 DF,  p-value: < 2.2e-16	
  library(boot)	
  
  Flux$mean_flux_ECS_mmold <- Flux$mean_flux_ECS/1000*24	
  Flux$mean_flux_MacIntyre_mmold  <- Flux$mean_flux_MacIntyre/1000*24	
  
  plot(Flux$mean_flux_ECS, Flux$mean_flux_MacIntyre)	
  
  
  lm_flux <- lm(Flux$mean_flux_ECS_mmold ~ Flux$mean_flux_MacIntyre_mmold)	
  lm_flux <- lm(Flux$mean_flux_MacIntyre_mmold ~ Flux$mean_flux_ECS_mmold)	
  
  plot(lm_flux)	
  summary(lm_flux)	
  
  Flux_ECSMacIntyre_boot <- boot(statistic = bootReg, formula =  mean_flux_MacIntyre_mmold ~ mean_flux_ECS_mmold, data = Flux, R = 10000)	
  boot.ci(Flux_ECSMacIntyre_boot, type = "bca", index = 1)	
  boot.ci(Flux_ECSMacIntyre_boot, type = "bca", index = 2)	
  
  
  range_ECflux <- seq(steps = nrow(Flux+400), from = min(Flux$mean_flux_ECS_mmold-200, na.rm = TRUE), to = max(Flux$mean_flux_ECS_mmold+200, na.rm = TRUE))	
  
  # Plot the regression of both modeled CO2 fluxes	
  {	
    par(mfrow = c(1,1))	
    par(mar= mar.default + 0.5)	
    
    plot(Flux$mean_flux_ECS_mmold, Flux$mean_flux_MacIntyre_mmold, ylab = expression("Daily averaged modeled flux [mmol m"^"-2"*" d"^"-1"*"]" ), xlab = expression("Daily averaged measured flux [mmol m"^"-2"*" d"^"-1"*"]"), 	
         ylim = c(-300,200), xlim =c(-300,200 ), pch = 16, col = "darkblue", cex.axis = 1.1, cex.lab = 1.1)	
    abline(lm_flux, col="red")	
    lines(range_ECflux, -34.33 + 0.0676 * range_ECflux, col="red", lty =3)	
    lines(range_ECflux, -25.08 + 0.1757 * range_ECflux, col="red", lty =3)	
    legend(	
      "bottomleft",	
      inset = c(0.0,0.0),	
      cex = 1.1, 	
      bty = "n", 	
      legend = c(expression("CO"[2]*" Flux"), expression("Regression line"), expression("95% confidence interval")), 	
      text.col = "black",	
      col = c("darkblue", "red", "red"), 	
      pch = c(16,NA, NA), lty = c(NA,1,3))	
    
    
  }	
  
  # plot modeled O2 against modeled CO2	
  
  lm_flux <- lm(Flux$mean_flux_MacIntyre_mmold ~ Flux$mean_flux_DO)	
  
  par(mfrow = c(2,2))	
  plot(lm_flux)	
  summary(lm_flux)	
  
  Flux_ECSMacIntyre_boot <- boot(statistic = bootReg, formula =  mean_flux_MacIntyre_mmold ~ mean_flux_ECS_mmold, data = Flux, R = 10000)	
  boot.ci(Flux_ECSMacIntyre_boot, type = "bca", index = 1)	
  boot.ci(Flux_ECSMacIntyre_boot, type = "bca", index = 2)	
  
  
  {	
    par(mfrow = c(1,1))	
    par(mar= mar.default + 0.5)	
    
    plot(Flux$mean_flux_DO, Flux$mean_flux_MacIntyre_mmold, ylab = expression("Daily averaged modeled CO"[2]*" flux [mmol m"^"-2"*" d"^"-1"*"]" ), xlab = expression("Daily averaged modeled O"[2]*" flux [mmol m"^"-2"*" d"^"-1"*"]"), 	
         ylim = c(-300,200), pch = 16, col = "darkblue", cex.axis = 1.1, cex.lab = 1.1)	
    abline(lm_flux, col="red")	
    abline(summary(lm_flux)$coefficients[,1] + summary(lm_flux)$coefficients[,2]*1.96, col = "red", lty = 3)	
    abline(summary(lm_flux)$coefficients[,1] - summary(lm_flux)$coefficients[,2]*1.96, col = "red", lty = 3)	
    legend(	
      "bottomleft",	
      inset = c(0.0,0.0),	
      cex = 1.1, 	
      bty = "n", 	
      legend = c(expression("Gas Flux"), expression("Regression line"), expression("95% confidence interval")), 	
      text.col = "black",	
      col = c("darkblue", "red", "red"), 	
      pch = c(16,NA, NA), lty = c(NA,1,3))	
    
    
  }	
  
  
  
  DF_Metabolic_ECS$Tot_DIC_change_ECS <- (DF_Metabolic_ECS$delta.dic + DF_Metabolic_ECS$CO2FluxECS_filter_umol_L_h + DF_Metabolic_ECS$`CO2flux_benthic_umol/L/h` + DF_Metabolic_ECS$Ca_prec_umolL - DF_Metabolic_ECS$`DIC_IN-OUT_umol/h`)	
  DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS <- (DF_Metabolic_ECS$delta.dic + DF_Metabolic_ECS$CO2FluxECS_filter_umol_L_h + DF_Metabolic_ECS$`CO2flux_benthic_umol/L/h` + DF_Metabolic_ECS$Ca_prec_umolL - DF_Metabolic_ECS$`DIC_IN-OUT_umol/h`)/DF_Metabolic_ECS$ChlVoltNorm	
  
  
  DF_Metabolic_ECS$Tot_DIC_change <- (DF_Metabolic_ECS$delta.dic + DF_Metabolic_ECS$CO2flux_macIntyreR/1800 + DF_Metabolic_ECS$`CO2flux_benthic_umol/L/h` + DF_Metabolic_ECS$Ca_prec_umolL - DF_Metabolic_ECS$`DIC_IN-OUT_umol/h`)	
  DF_Metabolic_ECS$Tot_DIC_change_Chl <- (DF_Metabolic_ECS$delta.dic + DF_Metabolic_ECS$CO2flux_macIntyreR/1800 + DF_Metabolic_ECS$`CO2flux_benthic_umol/L/h` + DF_Metabolic_ECS$Ca_prec_umolL - DF_Metabolic_ECS$`DIC_IN-OUT_umol/h`)/DF_Metabolic_ECS$ChlVoltNorm	
  
  
  DF_Metabolic_ECS$Tot_DO_change <- (DF_Metabolic_ECS$delta.do + DF_Metabolic_ECS$DOflux_umol/1800 - DF_Metabolic_ECS$`CO2flux_benthic_umol/L/h`)	
  DF_Metabolic_ECS$Tot_DO_change_Chl <- (DF_Metabolic_ECS$delta.do + DF_Metabolic_ECS$DOflux_umol/1800 - DF_Metabolic_ECS$`CO2flux_benthic_umol/L/h` )/DF_Metabolic_ECS$ChlVoltNorm	
  
  
  
  plot(DF_Metabolic_ECS$Tot_DIC_change, DF_Metabolic_ECS$Tot_DIC_change_ECS)	
  abline(a = 0, b = 1, col = "red")	
  
  plot(DF_Metabolic_ECS$CO2FluxECS_filter_umol_L_h/DF_Metabolic_ECS$delta.dic, ylim = c(-15, 15))	
  
  
  
  # Conpare the change to the change in DO	
  plot(DF_Metabolic_ECS$`Average Sampling time`, DF_Metabolic_ECS$Tot_DIC_change_Chl, type = 'l')	
  lines(DF_Metabolic_ECS$`Average Sampling time`, DF_Metabolic_ECS$Tot_DO_change_Chl*-1, col="blue")	
  
  plot(DF_Metabolic_ECS$`Average Sampling time`, DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS*-1, type = 'l')	
  lines(DF_Metabolic_ECS$`Average Sampling time`, DF_Metabolic_ECS$Tot_DO_change_Chl, col="blue")	
  
  
  plot(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS_log*-1~DF_Metabolic_ECS$Tot_DO_change_Chl)	
  summary(lm(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS*-1 ~ DF_Metabolic_ECS$Tot_DO_change_Chl))	
  
  
  
  par(mfrow = c(2,2))	
  plot(lm(sqrt(abs(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2552:2555)]))*sign(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2552:2555)]) ~ DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2552:2555)]))	
  par(mfrow = c(1,1))	
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  summary(lm(sqrt(abs(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2552:2555)]))*sign(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2552:2555)]) ~ DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2552:2555)]))	
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  lm_d_DIC_DO <- lm(sqrt(abs(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2552:2555)]))*sign(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2552:2555)]) ~ DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2552:2555)])	
  p.adjust(2.2e-16, method = "bonferroni", n = 6560)	
  
  plot(sqrt(abs(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2552:2555)]))*sign(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2552:2555)]) ~ DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2552:2555)])	
  
  # plot the changes in DIC and DO 	
  {	
    plot(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2551:2555)]~DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2551:2555)], pch = 16, col = "darkblue", ylim = c(-2,2), xlim = c(-2,2),	
         ylab = "Rate of change in DIC ["~mu*"mol L"^"-1"*" h"^"-1"*"]", xlab = "Rate of change in DO ["~mu*"mol L"^"-1"*" h"^"-1"*"]", 	
         cex.lab = 1.2, cex.axis = 1.2)	
    abline(0,-1, lty = 3, col = "darkgrey")	
    points(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2551:2555)]~DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2551:2555)], pch = 16, col = "darkblue")	
    abline(sign(summary(lm_d_DIC_DO)$coef[,1])*summary(lm_d_DIC_DO)$coef[,1]^2, col = "red")	
    abline(sign(summary(lm_d_DIC_DO)$coef[,1])*summary(lm_d_DIC_DO)$coef[,1]^2 + summary(lm_d_DIC_DO)$coef[,2]*1.96, lty = 3, col = "red")	
    abline(sign(summary(lm_d_DIC_DO)$coef[,1])*summary(lm_d_DIC_DO)$coef[,1]^2 - summary(lm_d_DIC_DO)$coef[,2]*1.96, lty = 3, col = "red")	
    legend(	
      "bottomleft",	
      bty = 'n',	
      pch = c(16, NA, NA),	
      lty = c(NA, 1, 3),	
      col = c("darkblue", "red", "red"),	
      legend = c("Rate of change in gas", "Regression line", "95% confidence interval")	
    )	
  } 	
  
  
  
  
  par(mfrow = c(2,2))	
  plot(lm(sqrt(abs(DF_Metabolic_ECS$Tot_DIC_change_Chl[-c(2552:2555)]))*sign(DF_Metabolic_ECS$Tot_DIC_change_Chl[-c(2552:2555)]) ~ DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2552:2555)]))	
  par(mfrow = c(1,1))	
  
  
  summary(lm(sqrt(abs(DF_Metabolic_ECS$Tot_DIC_change_Chl[-c(2552:2555)]))*sign(DF_Metabolic_ECS$Tot_DIC_change_Chl[-c(2552:2555)]) ~ DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2552:2555)]))	
  
  lm_d_DIC_DO_MI <- lm(sqrt(abs(DF_Metabolic_ECS$Tot_DIC_change_Chl[-c(2552:2555)]))*sign(DF_Metabolic_ECS$Tot_DIC_change_Chl[-c(2552:2555)]) ~ DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2552:2555)])	
  p.adjust(2.2e-16, method = "bonferroni", n = 6561)	
  
  plot(sqrt(abs(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2552:2555)]))*sign(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2552:2555)]) ~ DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2552:2555)])	
  
  {	
    par(mfrow = c(2,1))	
    par(mar = c(1.5,5,3,2))	
    
    {	
      plot(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2551:2555)]~DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2551:2555)], pch = 16, col = "darkblue", ylim = c(-2,2), xlim = c(-2,2),	
           ylab = "Rate of change in DIC ["~mu*"mol L"^"-1"*" h"^"-1"*"]", xlab = "", xaxt = 'n',	
           cex.lab = 1.2, cex.axis = 1.2)	
      abline(0,-1, lty = 3, col = "darkgrey")	
      points(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2551:2555)]~DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2551:2555)], pch = 16, col = "darkblue")	
      abline(sign(summary(lm_d_DIC_DO)$coef[,1])*summary(lm_d_DIC_DO)$coef[,1]^2, col = "red")	
      abline(sign(summary(lm_d_DIC_DO)$coef[,1])*summary(lm_d_DIC_DO)$coef[,1]^2 + summary(lm_d_DIC_DO)$coef[,2]*1.96, lty = 3, col = "red")	
      abline(sign(summary(lm_d_DIC_DO)$coef[,1])*summary(lm_d_DIC_DO)$coef[,1]^2 - summary(lm_d_DIC_DO)$coef[,2]*1.96, lty = 3, col = "red")	
      legend(	
        "bottomleft",	
        bty = 'n',	
        pch = c(16, NA, NA),	
        lty = c(NA, 1, 3),	
        col = c("darkblue", "red", "red"),	
        legend = c("Rate of change in gas", "Regression line", "95% confidence interval")	
      )	
      legend(	
        "topright",	
        bty = 'n',	
        legend = "A",	
        cex = 1.2	
      )	
    } 	
    
    par(mar = c(4,5,0.5,2))	
    
    {	
      plot(DF_Metabolic_ECS$Tot_DIC_change_Chl[-c(2551:2555)]~DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2551:2555)], pch = 16, col = "darkblue", ylim = c(-2,2), xlim = c(-2,2),	
           ylab = "Rate of change in DIC ["~mu*"mol L"^"-1"*" h"^"-1"*"]", xlab = "Rate of change in DO ["~mu*"mol L"^"-1"*" h"^"-1"*"]", 	
           cex.lab = 1.2, cex.axis = 1.2)	
      abline(0,-1, lty = 3, col = "darkgrey")	
      points(DF_Metabolic_ECS$Tot_DIC_change_Chl_ECS[-c(2551:2555)]~DF_Metabolic_ECS$Tot_DO_change_Chl[-c(2551:2555)], pch = 16, col = "darkblue")	
      abline(sign(summary(lm_d_DIC_DO_MI)$coef[,1])*summary(lm_d_DIC_DO_MI)$coef[,1]^2, col = "red")	
      abline(sign(summary(lm_d_DIC_DO_MI)$coef[,1])*summary(lm_d_DIC_DO_MI)$coef[,1]^2 + summary(lm_d_DIC_DO_MI)$coef[,2]*1.96, lty = 3, col = "red")	
      abline(sign(summary(lm_d_DIC_DO_MI)$coef[,1])*summary(lm_d_DIC_DO_MI)$coef[,1]^2 - summary(lm_d_DIC_DO_MI)$coef[,2]*1.96, lty = 3, col = "red")	
      legend(	
        "bottomleft",	
        bty = 'n',	
        pch = c(16, NA, NA),	
        lty = c(NA, 1, 3),	
        col = c("darkblue", "red", "red"),	
        legend = c("Rate of change in gas", "Regression line", "95% confidence interval")	
      )	
      legend(	
        "topright",	
        bty = 'n',	
        legend = "B",	
        cex = 1.2	
      )	
    } 	
    
  }	
  
  
}	









##################################################################################################################################	
#	
#  Plot ECS data, BF etc	
#	
###################################################################################################################################	



# incoming and outgoing DIC in AVP	
{	
  
  summary(DF_data$`DIC_IN-OUT_umol/h`*-1)	
  sd(DF_data$`DIC_IN-OUT_umol/h`*-1, na.rm = TRUE)/sqrt(length(DF_data$`DIC_IN-OUT_umol/h`*-1))	
  
  par(mfrow = c(1,1))	
  par(mar = mar.default + c(0,1,0,3) )	
  
  plot(DF_data$Date, DF_data$`DIC_IN-OUT_umol/h`*-1, ylim = c(-3, 3), 	
       ylab = "DIC ["~mu~"mol L"^"-1"*" h"^"-1"*"]", xlab = "Months [2016]", pch = 16, col = "darkblue", type = 'l')	
  abline(h = 0, col ="grey")	
  #rect(xleft = DF_data$Date[DF_data$Date == as.POSIXct("2016-05-18")], 	
  #     xright = DF_data$Date[DF_data$Date == as.POSIXct("2016-11-17")],	
  #     ybottom = 0, ytop = -2.15, col = rgb(0.45,1,0.05, alpha = 0.002), border = NA, mar = -5)	
  
  
  #rect(xleft = DF_data$Date[DF_data$Date == as.POSIXct("2016-05-18")], 	
  #     xright = DF_data$Date[DF_data$Date == as.POSIXct("2016-11-17")],	
  #     ybottom = 0, ytop = 2.15, col = rgb(1,0.45,0.05, alpha = 0.002), border = NA, mar = -0.1)	
  
  mtext(side = 4, "Uptake DIC                                       . ", col = "darkgreen")	
  mtext(side = 4, " .                                        Ouput DIC", col = "darkred")	
  points(DF_data$Date, DF_data$`DIC_IN-OUT_umol/h`*-1, pch = 16, col = "darkblue", cex = 0.5)	
}	




# Bouyancy flux and temp	
{	
  
  DF_data$`Average Sampling time`[which.max(DF_data$T_C_water - DF_data$T_C_air)]	
  
  
  layout(matrix(c(1,2),c(2)),	
         heights = c(2,2), respect = FALSE)	
  par(mar=c(0.25, 5, 1, 2))	
  plot(DF_data$`Average Sampling time`, DF_data$T_C_air, type = 'l', col = "darkorange", ylim = c(0,50), 	
       ylab = "Temp [ºC degrees]", xlab = "", xaxt = 'n')	
  lines(DF_data$`Average Sampling time`, DF_data$T_C_water, type = 'l', col = "darkblue")	
  legend(	
    "topright",	
    inset = c(0.0, 0.0),	
    cex = 1.1,	
    bty = "n",	
    lty = 1,	
    col = c("darkorange", "darkblue"),	
    legend = c("Air temperature ", "Water temperature"))	
  
  legend(	
    "topleft",	
    inset = c(0.0, 0.0),	
    cex = 1.1,	
    bty = "n",	
    legend = "A ")	
  
  par(mar=c(4.5, 5, 0, 2))	
  plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$Mac_Intyre_Bfl, type = 'l', col = "darkblue", 	
       ylim = c(-4e-07, 3e-07), ylab = expression("BF [m"^"2"*" s"^"-3"*"]"), xlab = "Months [2016]")	
  abline(h = 0, col = "red")	
  
  legend(	
    "topleft",	
    inset = c(0.0, 0.0),	
    cex = 1.1,	
    bty = "n",	
    legend = "B ")	
  
  
  par(mfrow=c(1,1))	
  par(mar=c(5,4,4,2) + 0.1)	
  
  
  plot(DF_Metabolic$sw_PAR[DF_Metabolic$irr == 1], DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$irr == 1])	
  
  summary(DF_Metabolic$Mac_Intyre_Bfl)	
  
  DF_Metabolic$`Average Sampling time`[which.min(DF_Metabolic$Mac_Intyre_Bfl)]	
  mean(DF_Metabolic$Mac_Intyre_Bfl)	
  sd(DF_Metabolic$Mac_Intyre_Bfl)/length(DF_Metabolic$Mac_Intyre_Bfl)	
  
  mean(DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$irr == 1])	
  sd(DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$irr == 1])/length(DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$irr == 1])	
  
  mean(DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$irr == 0])	
  sd(DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$irr == 0])/length(DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$irr == 0])	
  
  t.test(DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$irr == 1], DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$irr == 0])	
  p.adjust(2.2e-16, method = "bonferroni", n = length(DF_Metabolic$Mac_Intyre_Bfl[DF_Metabolic$irr == 1]))	
  
  DF_Metabolic_ECS$`Average Sampling time`[which.min(DF_Metabolic_ECS$Mac_Intyre_Bfl)]	
  DF_Metabolic_ECS$`Average Sampling time`[which.min(DF_Metabolic_ECS$Mac_Intyre_Bfl[DF_Metabolic_ECS$`Average Sampling time` != as.POSIXct("2016-07-22")])]	
  
  
  tail(DF_Metabolic$Mac_Intyre_Bfl[order(DF_Metabolic$Mac_Intyre_Bfl)])	
  tail(DF_Metabolic$`Average Sampling time`[order(DF_Metabolic$Mac_Intyre_Bfl)])	
  
  library(dplyr)	
  DF_Metabolic %>%	
    group_by(Date_dusk) %>%	
    summarise(BF_mean = mean(Mac_Intyre_Bfl),	
              BF_min = min(Mac_Intyre_Bfl),	
              BF_max = max(Mac_Intyre_Bfl),	
              flux_mean = mean(CO2FluxECS_filter_umol_m2_h, na.rm = TRUE)) -> BF	
  detach("package:dplyr", unload = TRUE)	
  
  
  plot(BF)	
  summary(lm(BF$flux_mean ~ BF$BF_mean))	
  
  summary(DF_Metabolic$Mac_Intyre_Bfl)	
  
  
  
  
  
  names(BF)[1] <- "Day"	
  BF_Eva <- merge(BF, Neerslag_verdamping)	
  
  
  
  
}	


# Temp BF, wind, Evap, and precip	

plot(DF_data$Time, DF_data$windECS_U10)	

library(dplyr)	

DF_data %>%	
  group_by(Time) %>%	
  summarise(mean(windECS_U10)) -> Wind_day	
Wind_day$Time/3600	
detach("package:dplyr", unload = TRUE)	

#save the plot as 800 by 950	

{	
  
  layout(matrix(c(1:5),c(5)),	
         heights = c(1.4,1,2.6,0.9,1.2), respect = FALSE)	
  par(mar=c(1, 5, 0.5, 2))	
  plot(DF_data$`Average Sampling time`, DF_data$T_C_air, type = 'l', col = "darkorange", ylim = c(0,50), 	
       ylab = "Temp [ºC degrees]", xlab = "", xaxt = 'n', cex.lab = 1.4, cex.axis = 1.1)	
  lines(DF_data$`Average Sampling time`, DF_data$T_C_water, type = 'l', col = "darkblue")	
  legend(	
    "topright",	
    inset = c(0.02, 0.0),	
    cex = 1.3,	
    bty = "n",	
    lty = 1,	
    col = c("darkorange", "darkblue"),	
    legend = c("Air temperature", "Water temperature"))	
  
  legend(	
    "topleft",	
    inset = c(0.0, 0.0),	
    cex = 1.3,	
    bty = "n",	
    legend = "A ")	
  
  par(mar=c(1, 5, 0.5, 2))	
  plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$Mac_Intyre_Bfl, type = 'l', col = "darkblue", 	
       ylim = c(-4e-07, 3e-07), ylab = expression(~beta~"Flux [m"^"2"*" s"^"-3"*"]"), xlab = "", xaxt = 'n', 	
       cex.lab = 1.4, cex.axis = 1.1)	
  abline(h = 0, col = "red")	
  
  legend(	
    "topleft",	
    inset = c(0.0, 0.0),	
    cex = 1.3,	
    bty = "n",	
    legend = "B")	
  
  
  par(mar=c(1, 5, 0.5, 2))	
  plot(DF_data$`Average Sampling time`, DF_data$windECS_U10, ylab = expression("Windspeed [m s"^"-1"*"]", 	
                                                                               cex.lab = 1.4, cex.axis = 1.1), 	
       ylim = c(0,16), xlab = "", xaxt = 'n', type = 'l', col = "gray17", cex.lab = 1.3)	
  legend(	
    "topleft",	
    cex = 1.3,	
    bty = 'n',	
    legend = "C")	
  
  par(mar=c(1.5, 5, 0.5, 2))	
  plot(Neerslag_verdamping$Day[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"], 	
       Neerslag_verdamping$`neerslag [mm/d]`[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"], 	
       type = 'l', col = "firebrick3", ylab = expression("Precip [mm d"^"-1"*"]"), xlab = "", xaxt='n', ylim = c(0,50), 	
       cex.lab = 1.4, cex.axis = 1.1)	
  legend(	
    "topleft",	
    inset = c(0.0, 0.0),	
    cex = 1.3,	
    bty = "n",	
    legend = "D")	
  
  
  par(mar=c(5, 5, 0, 2))	
  plot(Neerslag_verdamping$Day[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"],	
       Neerslag_verdamping$`verdamping penman [mm/d]`[Neerslag_verdamping$Day > "2016-5-18" & Neerslag_verdamping$Day < "2016-11-17"],	
       type = 'l', col = "dodgerblue3", lwd = 1	
       , xlab = "Time [Months]",  ylab = expression("Evap [mm d"^"-1"*"]"), cex.lab = 1.4, cex.axis = 1.1)	
  legend(	
    "topleft",	
    inset = c(0.0, 0.0),	
    cex = 1.3,	
    bty = "n",	
    legend = "E")	
  
  
  
  
  par(new = TRUE)	
  par(mfrow=c(1,1))	
  par(mar=c(27.7,35,17.7,1.7))	
  
  par(new = TRUE)	
  par(mfrow=c(1,1))	
  par(mar=c(27,33,17.3,1.7))	
  plot(Wind_day$Time/3600, Wind_day$`mean(windECS_U10)`, ylab = "", xlab = "", cex.lab = 0.7, cex.axis = 0.7, las = 1)	
  text(12.5, -2, "Hours", cex = 0.7, xpd=NA)	
  
  
}	

par(mfrow=c(1,1))	
par(mar=c(5,4,4,2) + 0.1)	

# Plot the mean and normal CO2 flux	
{	
  
  
  
  
  summary(Flux$mean_flux_ECS/1000*24)	
  Flux$Date[which.max(Flux$mean_flux_ECS)]	
  Flux$Date[which.min(Flux$mean_flux_ECS)]	
  
  
  summary(Flux$mean_flux_ECS)	
  Flux$Date[which.max(Flux$mean_flux_ECS)]	
  Flux$Date[which.min(Flux$mean_flux_ECS)]	
  
  
  
  layout(matrix(c(1,2),c(2)),	
         heights = c(2,1), respect = FALSE)	
  
  par(mar=c(1.5, 5, 3, 2))	
  #axis(2, at = x <- seq(-1e4, 1e4, by = 1e4),  labels = paste(x/10000, "e+04", sep = "" ))	
  #points(Flux$Date_dusk, Flux$median_flux_cole, col = "red")	
  plot(Flux$Date, Flux$mean_flux_ECS/1000*24, ylab = expression("Daily averaged flux [mmol m"^"-2"*" d"^"-1"*"]"),xaxt = 'n', xlab = "", pch = 19, 	
       col = rgb(0,0,0,0.9), type = 'l', lwd = 2, ylim = c(-500, 500), cex.axis = 1.1, cex.lab = 1.1)	
  #points(Flux$Date_dusk, Flux$median_flux_cole, col = "red")	
  points(Flux$Date, Flux$mean_flux_DO/1000*24, col = "dodgerblue1", pch = 15)	
  points(Flux$Date[Flux$mean_flux_MacIntyre != 0], Flux$mean_flux_MacIntyre[Flux$mean_flux_MacIntyre != 0]/1000*24, col = rgb(0,0,1,0.7), pch = 19)	
  lines(Flux$Date, Flux$mean_flux_ECS/1000*24, ylab = expression("mean CO"[2]*" flux [mmol m"^"-2"*" d"^"-1"*"]"), xlab = "Months [2016]", pch = 19, 	
        col = rgb(0,0,0,0.9), type = 'l', lwd = 2, ylim = c(-500, 500), cex.axis = 1.1, cex.lab = 1.1)	
  abline(h = 0, col = "red", lwd = 2, lty = 3)	
  legend(	
    "bottomleft",	
    bty = "n",	
    legend = c(expression("Measured CO"[2]*" flux"), expression("Modeled CO"[2]*" flux"), expression("Modeled O"[2]*" flux")),	
    cex = 1.1,	
    pch = c(NA,19,15), lty = c(1,NA,NA), lwd(NA, 2,2),	
    col = c(rgb(0,0,0,0.7), rgb(0,0,1,0.7), "steelblue3"))	
  legend(	
    "topright",	
    inset = c(0.0, 0.0),	
    cex = 1.1,	
    bty = "n",	
    legend = "A ")	
  
  
  
  par(mar=c(4.5, 5, 0, 2))	
  plot(DF_Metabolic$`Average Sampling time`, DF_Metabolic$CO2FluxECS_filter_umol_m2_h/1000*24, col = rgb(0,0,0), pch = 149,	
       ylab = expression("CO"[2]*" flux [mmol m"^"-2"*" d"^"-1"*"]"), xlab = "Time [Months]",  ylim= c(-5500, 5500), cex.axis = 1.1, cex.lab = 1.1)	
  points(DF_Metabolic$`Average Sampling time`, DF_Metabolic$CO2flux_macIntyreR/1000, col = rgb(0,0,1), pch = 149)	
  abline(h = 0, col = rgb(1,0,0,0.7))	
  legend("bottomleft", 	
         inset = c(0.02,0.0), 	
         cex = 1.1, 	
         bty = "n", 	
         legend = c(expression("Measured CO"[2]*" flux"), expression("Modeled CO"[2]*" flux")), 	
         text.col = c("black", "black"),	
         col = c("black", "blue"),	
         pch = c(149))	
  legend(	
    "topright",	
    inset = c(0.0, 0.0),	
    cex = 1.1,	
    bty = "n",	
    legend = "B ")	
  
  
  
  par(mfrow=c(1,1))	
  par(mar=c(5,4,4,2) + 0.1)	
  
}	

..	





#################################################################################################################################################################	
#################################################################################################################################################################	
##########################################################################      NEP       #######################################################################	
#################################################################################################################################################################	
#################################################################################################################################################################	






####################################################################################################################	
#	
# Check the metabolism by the max min change per day based on the different fluxes	
#	
####################################################################################################################	

DF_Metabolic$Tot_DIC_change_ECS_kal <- (DF_Metabolic$delta.dic + DF_Metabolic$CO2FluxECS_filter_umol_m2_h_kalman/1800 + DF_Metabolic$`CO2flux_benthic_umol/L/h` + DF_Metabolic$Ca_prec_umolL - DF_Metabolic$`DIC_IN-OUT_umol/h`)	


DF_Metabolic$Tot_DIC_change_ECS <- (DF_Metabolic$delta.dic + DF_Metabolic$CO2FluxECS_filter_umol_L_h + DF_Metabolic$`CO2flux_benthic_umol/L/h` + DF_Metabolic$Ca_prec_umolL - DF_Metabolic$`DIC_IN-OUT_umol/h`)	

DF_Metabolic$Tot_DO_change <- (DF_Metabolic$delta.do + DF_Metabolic$`DOflux_umol/L/h` + DF_Metabolic$`CO2flux_benthic_umol/L/h`*-1)	


DF_Metabolic_ECS$Tot_DIC_change_ECS <- (DF_Metabolic_ECS$delta.dic + DF_Metabolic_ECS$CO2FluxECS_filter_umol_L_h + DF_Metabolic_ECS$`CO2flux_benthic_umol/L/h` + DF_Metabolic_ECS$Ca_prec_umolL - DF_Metabolic_ECS$`DIC_IN-OUT_umol/h`)	

DF_Metabolic_ECS$Tot_DIC_change_ECS_Chl <- (DF_Metabolic_ECS$delta.dic + DF_Metabolic_ECS$CO2FluxECS_filter_umol_L_h + DF_Metabolic_ECS$`CO2flux_benthic_umol/L/h` + DF_Metabolic_ECS$Ca_prec_umolL - DF_Metabolic_ECS$`DIC_IN-OUT_umol/h`)/DF_Metabolic_ECS$Chl_lab	






# the ECS flux was used but could not be used because of the filtering of the data. 	
# This resulted in a positive NEP, therefore use the raw data and the MI flux	


DF_data$`CO2flux_benthic_umol/L/h` <- ((10 * exp(0.0362 * DF_data$T_C_water + 1.635)/44.0095) * 1000 * -1)/1800/24        # conversion from mg/m2/d to umol/L/h multiplied by -1 because it relfects a negative flux, C into the the lake	
DF_data$`O2flux_benthic_umol/L/h` <- ((10 * exp(0.0362 * DF_data$T_C_water + 1.635)/44.0095) * 1000 )/1800/24        # conversion from mg/m2/d to umol/L/h multiplied by -1 because it relfects a negative flux, C into the the lake	


# DF_data$Tot_DO_change <- DF_data$delta.do + DF_data$`DOflux_umol/L/h` + DF_Metabolic$`O2flux_benthic_umol/L/h`	



Metabolic_Flux <- DF_Metabolic[complete.cases(DF_Metabolic[,c("delta.dic", "CO2flux_macIntyreR", "CO2flux_benthic_umol/L/h", "Ca_prec_umolL", "DIC_IN-OUT_umol/h")]),]	
Metabolic_Flux$Tot_DIC_change <- (Metabolic_Flux$delta.dic + Metabolic_Flux$CO2flux_macIntyreR/1800 + Metabolic_Flux$`CO2flux_benthic_umol/L/h` + Metabolic_Flux$Ca_prec_umolL - Metabolic_Flux$`DIC_IN-OUT_umol/h`)	



# check without the model, just the max value and min value on a day	


MetabolicRate_perHour <- function(date, time, light, gasChange_perHour) {	
  
  if(length(light != gas)){ warning(paste("length differs"))}	
  DF <- data.frame(date = date, time = time, light = light, gas = gasChange_perHour)	
  
  DF$irr <- 1	
  if(DF$light == 0){ DF$irr = 0}	
  
  require(dplyr)	
  DF %>%	
    group_by(date, irr) %>%	
    summarise(mean_change = mean(gas),	
              sun_up = min(time),	
              sun_down = max(time)) -> DF_Mean	
  
  DF_mean$Sunhours <- DF_mean$sun_down - DF_mean$sun_up	
  
}	

#perhaps later	

#mean_DO_change = mean(Tot_DO_change, na.rm = TRUE),	
#min_DO_change = min(Tot_DO_change, na.rm = TRUE),	
#max_DO_change = max(Tot_DO_change, na.rm = TRUE),	


library(dplyr)	
Metabolic_Flux %>%	
  group_by(Date_dusk, irr) %>%	
  summarise(sun_up = min(as.numeric(Time)),	
            sun_down = max(as.numeric(Time)),	
            max_I_surface = mean(`I (umol/m2/s)`)*60*60*24,	
            max_I_60 = mean(`I(60)`)*60*60*24,	
            mean_DIC_change = mean(Tot_DIC_change, na.rm = TRUE),	
            min_DIC_change = min(Tot_DIC_change, na.rm = TRUE),	
            max_DIC_change = max(Tot_DIC_change, na.rm = TRUE),	
            sum_DIC_change = sum(Tot_DIC_change, na.rm = TRUE),	
            mean_DO_change = mean(Tot_DO_change, na.rm = TRUE),	
            mean_T_C_water = mean(T_C_water, na.rm = TRUE),	
            min_T_C_water = min(T_C_water, na.rm = TRUE),	
            max_T_C_water = max(T_C_water, na.rm = TRUE),	
            mean_DW_mg_L = mean(DW_mg_L, na.rm = TRUE),	
            mean_flux_MI = mean(CO2flux_macIntyreR, na.rm = TRUE),	
            mean_ChlVolt_N = mean(ChlVoltNorm, na.rm = TRUE),	
            mean_DIC = mean(DIC_Alk_smooth, na.rm = TRUE),	
            Tot_change_DIC = (max(DIC_Alk_smooth, na.rm = TRUE) - min(DIC_Alk_smooth, na.rm = TRUE)),	
            Tot_change_DO = (max(DO_umol_smooth, na.rm = TRUE) - min(DO_umol_smooth, na.rm = TRUE)),	
            mean_chl_lab =  mean(Chl_lab, na.rm = TRUE),	
            mean_chl_volt = mean(ChlVoltNorm, na.rm = TRUE)) -> FLux_model	




library(dplyr)	
DF_Metabolic %>%	
  group_by(Date_dusk, irr) %>%	
  summarise(sun_up = min(as.numeric(Time)),	
            sun_down = max(as.numeric(Time)),	
            mean_DIC_change = mean(Tot_DIC_change_ECS_kal, na.rm = TRUE),	
            min_DIC_change = min(Tot_DIC_change_ECS_kal, na.rm = TRUE),	
            max_DIC_change = max(Tot_DIC_change_ECS_kal, na.rm = TRUE),	
            sum_DIC_change = sum(Tot_DIC_change_ECS_kal, na.rm = TRUE),	
            mean_T_C_water = mean(T_C_water, na.rm = TRUE),	
            min_T_C_water = min(T_C_water, na.rm = TRUE),	
            max_T_C_water = max(T_C_water, na.rm = TRUE),	
            mean_DW_mg_L = mean(DW_mg_L, na.rm = TRUE),	
            mean_flux_MI = mean(CO2FluxECS_filter_umol_L_h, na.rm = TRUE),	
            mean_ChlVolt_N = mean(ChlVoltNorm, na.rm = TRUE),	
            mean_DIC = mean(DIC_Alk_smooth, na.rm = TRUE),	
            Tot_change_DIC = (max(DIC_Alk_smooth, na.rm = TRUE) - min(DIC_Alk_smooth, na.rm = TRUE)),	
            Tot_change_DO = (max(DO_umol_smooth, na.rm = TRUE) - min(DO_umol_smooth, na.rm = TRUE)),	
            mean_chl_lab =  mean(Chl_lab, na.rm = TRUE),	
            mean_chl_volt = mean(ChlVoltNorm, na.rm = TRUE)) -> FLux_ECS	




library(dplyr)	
DF_Metabolic %>%	
  group_by(Date_dusk) %>%	
  summarise(sun_up = min(as.numeric(Time)),	
            sun_down = max(as.numeric(Time)),	
            mean_DIC_change = mean(Tot_DIC_change_ECS_kal, na.rm = TRUE),	
            min_DIC_change = min(Tot_DIC_change_ECS_kal, na.rm = TRUE),	
            max_DIC_change = max(Tot_DIC_change_ECS_kal, na.rm = TRUE),	
            sum_DIC_change = sum(Tot_DIC_change_ECS_kal, na.rm = TRUE),	
            mean_T_C_water = mean(T_C_water, na.rm = TRUE),	
            min_T_C_water = min(T_C_water, na.rm = TRUE),	
            max_T_C_water = max(T_C_water, na.rm = TRUE),	
            mean_DW_mg_L = mean(DW_mg_L, na.rm = TRUE),	
            mean_flux_ECS = mean(CO2FluxECS_filter_umol_L_h, na.rm = TRUE),	
            mean_ChlVolt_N = mean(ChlVoltNorm, na.rm = TRUE),	
            mean_DIC = mean(DIC_Alk_smooth, na.rm = TRUE),	
            Tot_change_DIC = (max(DIC_Alk_smooth, na.rm = TRUE) - min(DIC_Alk_smooth, na.rm = TRUE)),	
            Tot_change_DO = (max(DO_umol_smooth, na.rm = TRUE) - min(DO_umol_smooth, na.rm = TRUE)),	
            mean_chl_lab =  mean(Chl_lab, na.rm = TRUE),	
            mean_chl_volt = mean(ChlVoltNorm, na.rm = TRUE)) -> FLux_ECS_NEPDay	


plot(FLux_ECS_NEPDay$Date_dusk, FLux_ECS_NEPDay$mean_DIC_change)	


abline(h=0)	



FLux_model$SunHours <- (FLux_model$sun_down -  FLux_model$sun_up)/3600	
plot(FLux_model$Date_dusk[FLux_model$irr == 1], FLux_model$SunHours[FLux_model$irr == 1])	

FLux_ECS$SunHours <- (FLux_ECS$sun_down -  FLux_ECS$sun_up)/3600	
plot(FLux_ECS$Date_dusk[FLux_ECS$irr == 1], FLux_ECS$SunHours[FLux_ECS$irr == 1])	


Metabolic_model <- ungroup(FLux_model)	
Metabolic_ECS <- ungroup(FLux_ECS)	
detach("package:dplyr", unload = TRUE)	

write.table(Metabolic_model, "C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Metabolic_model.txt", sep = ";")	
#write.table(Metabolic_ECS, "C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/Metabolic_ECS.txt", sep = ";")	

# The file below was addapted in excell to include GPP 	
#Metabolic <- read_delim("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Metabolic.csv", 	
#                       ";", escape_double = FALSE, col_types = cols(Date = col_date(format = "%m/%d/%Y")), 	
#                      trim_ws = TRUE)	

Metabolic_model <- read_delim("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/Metabolic_model.csv", 	
                              
                              
                              
                              
                              ";", escape_double = FALSE, col_types = cols(Date = col_date(format = "%m/%d/%Y")), 	
                              trim_ws = TRUE)	

Metabolic_ECS <- read_delim("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/Metabolic_ECS.csv", 	
                            ";", escape_double = FALSE, col_types = cols(Date = col_date(format = "%m/%d/%Y")), 	
                            trim_ws = TRUE)	

Metabolic_model_DO <- read_delim("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/Metabolic_model_DO.csv", 	
                                 ";", escape_double = FALSE, col_types = cols(Date = col_date(format = "%d-%m-%Y")), 	
                                 trim_ws = TRUE)	

# Metabolism in umol/L/h based on DIC, the values R_day etc are in umol/L/d	

Metabolic_model <- Metabolic_model[complete.cases(Metabolic_model[c("R", "NEP")]),]	

# Metabolism in umol/L/h based on DIC, the values R_day etc are in umol/L/d	

Metabolic_model_DO <- Metabolic_model_DO[complete.cases(Metabolic_model_DO[c("R", "NEP")]),]	

#Metabolic_model <- Metabolic_model[complete.cases(Metabolic_model[c("R", "NEP")]),]	
Metabolic_ECS <- Metabolic_ECS[complete.cases(Metabolic_ECS[c("R", "NEP")]),]	

#Plot the raw per hour calculated metabolic rates	
par(mfrow = c(3,1))	
plot(Metabolic_model$Date, Metabolic_model$GPP)	
plot(Metabolic_model$Date, Metabolic_model$R)	
plot(Metabolic_model$Date, Metabolic_model$NEP)	
abline(h = 0, col = "red")	
par(mfrow = c(1,1))	


par(mfrow = c(3,1))	
plot(Metabolic_ECS$Date, Metabolic_ECS$GPP)	
plot(Metabolic_ECS$Date, Metabolic_ECS$R)	
plot(Metabolic_ECS$Date, Metabolic_ECS$NEP)	
abline(h = 0, col = "red")	
par(mfrow = c(1,1))	



#Calculate the metabolic per day, with in acount the sun hours. 	
par(mfrow = c(3,1))	
plot(Metabolic_model_DO$Date, Metabolic_model_DO$GPP_day)	
plot(Metabolic_model_DO$Date, Metabolic_model_DO$R_day)	
plot(Metabolic_model_DO$Date, Metabolic_model_DO$NEP_day)	
abline(h = 0, col = "red")	
par(mfrow = c(1,1))	

par(mfrow = c(3,1))	
plot(Metabolic_model$Date, Metabolic_model$GPP_day)	
plot(Metabolic_model$Date, Metabolic_model$R_day)	
plot(Metabolic_model$Date, Metabolic_model$NEP_day)	
abline(h = 0, col = "red")	
par(mfrow = c(1,1))	

par(mfrow = c(3,1))	
plot(Metabolic_ECS$Date, Metabolic_ECS$GPP_day)	
plot(Metabolic_ECS$Date, Metabolic_ECS$R_day)	
plot(Metabolic_ECS$Date, Metabolic_ECS$NEP_day)	
abline(h = 0, col = "red")	
par(mfrow = c(1,1))	




plot(Metabolic_ECS$Date, Metabolic_ECS$GPP_day/Metabolic_ECS$Chl_lab) 	
plot(Metabolic_ECS$Chl_lab, Metabolic_ECS$GPP_day/Metabolic_ECS$Chl_lab)	
plot(Metabolic_ECS$Date, Metabolic_ECS$R_day/Metabolic_ECS$Chl_lab) 	
plot(Metabolic_ECS$Date, Metabolic_ECS$NEP_day/Metabolic_ECS$Chl_lab) 	
plot(Metabolic_ECS$Chl_lab, Metabolic_ECS$NEP_day) 	
abline(h=0, col = "red")	


plot(Metabolic_model$Date, Metabolic_model$GPP_day/Metabolic_model$Chl_lab) 	
plot(Metabolic_model$Chl_lab, Metabolic_model$GPP_day/Metabolic_model$Chl_lab)	
plot(Metabolic_model$Date, Metabolic_model$R_day/Metabolic_model$Chl_lab) 	
plot(Metabolic_model$Date, Metabolic_model$NEP_day/Metabolic_model$Chl_lab) 	
plot(Metabolic_model$Chl_lab, Metabolic_model$NEP_day) 	
abline(h=0, col = "red")	

# this was adapted las minute there the name hasn´t been changed	
# Daily DIC in out consists of the daily DIC in in umol/L/d	

Daily_DIC_in_out <- read_delim("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/AVP_balans/Carbon_balans/CarbonBudget_InOut.csv", 	
                               ";", escape_double = FALSE, 	
                               trim_ws = TRUE)	

Daily_DIC_in_out$Date <- as.Date(Daily_DIC_in_out$Date, format = "%m/%d/%Y")	
View(Daily_DIC_in_out)	
names(Daily_DIC_in_out)	
Daily_DIC_in_out <- Daily_DIC_in_out[complete.cases(Daily_DIC_in_out[,c(1:6)]),]	


library(plyr)	
met <- join(Metabolic_model, Metabolic_ECS,  by = "Date")	
met <- join(met, Metabolic_model_DO,  by = "Date")	
met <- join(met, Daily_DIC_in_out[,c(1:6)],  by = "Date")	

names(met)	
names(met) <-  c( "Date", "Sun_hours", "NEP_model", "R_model", "GPP_model", "I_surface",  "Chl_lab", "T_C", "I",  "R_day__model" ,   	
                  "NEP_day_model",   "GPP_day_model"  , "diff_chl",    "change_chl", "Sun_hours", "NEP_ECS",  "R_ECS",  "GPP_ECS", "I_surface", "Chl_lab",   "T_C"  , "I_60",     	
                  "R_day_ECS" ,   "NEP_day_ECS" , "GPP_day_ECS", "diff_chl",  "change_chl",  "Sun_hours",  "NEP_DO" ,"R_DO", "GPP_DO", "I_surface", "Chl_lab",     	
                  "T_C", "I_60", "R_day_DO" ,   "NEP_day_DO" , "GPP_day_DO", "diff_chl", "change_chl", names(Daily_DIC_in_out[c(2:6)]))	

plot(met$Date, met$Inlet*-1)	
#The inlet is higher during spring and decreases during the year	

# DIC-In-Out is in mmol L/h and NEP in mmol/L/d	

plot(met$Inlet*-24/1000,met$NEP_day_ECS/-1000/met$Chl_lab)	
plot(met$Inlet*-24/1000,met$NEP_day_model/-1000/met$Chl_lab)	
plot(met$Inlet*-24/1000,met$NEP_day_DO/-1000/met$Chl_lab)	

# The NEP increases with increasing inlet 	
# The NEP got corrected for chlorophyll to eliminate effect of increased chlorophyll concetration or seasonal pattern of inlet	

plot(met$Inlet*-24/1000,met$change_chl)	
plot(met$Inlet*-24/1000,met$Chl_lab)	
plot(met$Incoming_MI*-24/1000,met$change_chl)	

plot(met$Date, met$Inlet*-1, ylim = c(-0.5, 2))	
points(met$Date, met$change_chl, col = "red")	

# there is however, no pattern visible with chlorophyll concentrations	
# but the change in chlorophyll does increase with increasing inlet of DIC. 	
# However, this might not be an effect of DIC inlet, it could be due to inlet of nutrients and more likely just a colinearity due to other factors	



par(mfrow = c(2,2))	
plot(lm(log(abs(met$NEP_day_ECS/met$Chl_lab))*sign(met$NEP_day_ECS) ~ met$Inlet))	
plot(lm(log(abs(met$NEP_day_model/met$Chl_lab))*sign(met$NEP_day_ECS) ~ met$Inlet))	
plot(lm(log(abs(met$NEP_day_DO/met$Chl_lab))*sign(met$NEP_day_ECS) ~ met$Inlet))	
par(mfrow = c(1,1))	


par(mfrow = c(2,2))	
plot(lm(met$NEP_day_ECS ~ met$Inlet))	
plot(lm(met$NEP_day_model/met$Chl_lab ~ met$Inlet))	
plot(lm(met$NEP_day_DO/met$Chl_lab ~ met$Inlet))	
par(mfrow = c(1,1))	



#convert to mmol/m2/d	
met$NEP_day_ECS_m2 <- met$NEP_day_ECS*-1.8	
met$NEP_day_model_m2 <- met$NEP_day_model*-1.8 	
met$NEP_day_DO_m2 <- met$NEP_day_DO*-1.8	

met$Inlet_m2 <- met$Inlet*-24*1.8	

plot(log(abs(met$NEP_day_ECS/met$Chl_lab))*sign(met$NEP_day_ECS) ~ met$Inlet)	


lm_DIC_EC <- lm(log(abs(met$NEP_day_ECS_m2/met$Chl_lab))*sign(met$NEP_day_ECS) ~ met$Inlet_m2)	
lm_DIC_model <- lm(log(abs(met$NEP_day_model_m2/met$Chl_lab))*sign(met$NEP_day_ECS) ~ met$Inlet_m2)	
lm_DIC_DO <- lm(log(abs(met$NEP_day_DO_m2/met$Chl_lab))*sign(met$NEP_day_ECS) ~ met$Inlet_m2)	


lm_DIC_EC <- lm(met$NEP_day_ECS_m2 ~ met$Inlet_m2)	
lm_DIC_model <- lm(met$NEP_day_model_m2 ~ met$Inlet_m2)	
lm_DIC_DO <- lm(met$NEP_day_DO_m2 ~ met$Inlet_m2)	


summary(lm_DIC_EC)	
summary(lm_DIC_model)	
summary(lm_DIC_DO)	


plot(met$NEP_day_ECS, met$NEP_day_model)	
summary(lm(met$NEP_day_ECS ~ met$NEP_day_model))	
summary(lm(met$GPP_day_ECS ~ met$GPP_day_model))	
summary(lm(met$R_day_ECS ~ met$R_day__model))	

summary(lm(met$NEP_day_ECS ~ met$NEP_day_DO))	
summary(lm(met$GPP_day_ECS ~ met$GPP_day_DO))	
summary(lm(met$R_day_ECS ~ met$R_day_DO))	

summary(lm(met$NEP_day_DO ~ met$NEP_day_model))	
summary(lm(met$GPP_day_DO ~ met$GPP_day_model))	
summary(lm(met$R_day_DO ~ met$R_day__model))	





summary(Metabolic_ECS$GPP_day)/1000	
summary(Metabolic_ECS$R_day)/1000	
summary(Metabolic_ECS$NEP_day)/1000	

summary(Metabolic_model$GPP_day)/1000	
summary(Metabolic_model$R_day)/1000	
summary(Metabolic_model$NEP_day)/1000	

summary(Metabolic_model_DO$GPP_day)/1000	
summary(Metabolic_model_DO$R_day)/1000	
summary(Metabolic_model_DO$NEP_day)/1000	


t.test(Metabolic_ECS$GPP_day, Metabolic_model$GPP_day)	
t.test(Metabolic_ECS$R_day, Metabolic_model$R_day)	
t.test(Metabolic_ECS$NEP_day, Metabolic_model$NEP_day)	

t.test(Metabolic_model_DO$GPP_day, Metabolic_model$GPP_day)	
t.test(Metabolic_model_DO$R_day, Metabolic_model$R_day)	
t.test(Metabolic_model_DO$NEP_day, Metabolic_model$NEP_day)	

t.test(Metabolic_ECS$GPP_day, Metabolic_model_DO$GPP_day)	
t.test(Metabolic_ECS$R_day, Metabolic_model_DO$R_day)	
t.test(Metabolic_ECS$NEP_day, Metabolic_model_DO$NEP_day)	

sd(Metabolic_ECS$GPP_day/1000, na.rm = TRUE)/length(Metabolic_ECS$GPP_day)	
sd(Metabolic_ECS$R_day/1000, na.rm = TRUE)/length(Metabolic_ECS$R_day)	
sd(Metabolic_ECS$NEP_day/1000, na.rm = TRUE)/length(Metabolic_ECS$NEP_day)	

sd(Metabolic_model$GPP_day/1000, na.rm = TRUE)/length(Metabolic_model$GPP_day)	
sd(Metabolic_model$R_day/1000, na.rm = TRUE)/length(Metabolic_model$R_day)	
sd(Metabolic_model$NEP_day/1000, na.rm = TRUE)/length(Metabolic_model$NEP_day)	


summary(Metabolic_ECS$GPP_day/Metabolic_ECS$Chl_lab)/24	
summary(Metabolic_ECS$R_day/Metabolic_ECS$Chl_lab)/24	
summary(Metabolic_ECS$NEP_day/Metabolic_ECS$Chl_lab)/24	




# plot the ecosytem production 	

{	
  layout(matrix(c(1,2,3,4,5,6),2,3),	
         heights = c(1,1), respect = FALSE)	
  par(mar = c(3, 5, 2, 0))	
  plot(Metabolic_ECS$Date[Metabolic_ECS$GPP_day < 0], Metabolic_ECS$GPP_day[Metabolic_ECS$GPP_day < 0]/-1000, col = "darkgreen", pch = 21, ylim = c(-1.4,1.4),	
       ylab = expression("Metabolic rate [mmol L"^"-1"*" d"^"-1"*"]"), xlab = "", xaxt = 'n',cex.axis = 1.3, cex.lab = 1.5)	
  abline(h = -0.0, col = "red", lwd = 2, lty = 3)	
  abline(v = as.Date("2016-08-10"), col = "darkgrey", lwd = 2.5, lty = 3)	
  points(Metabolic_ECS$Date[Metabolic_ECS$R_day > 0], Metabolic_ECS$R_day[Metabolic_ECS$R_day > 0]/-1000, col = "darkred", pch = 22)	
  points(Metabolic_ECS$Date[Metabolic_ECS$GPP_day < 0], Metabolic_ECS$GPP_day[Metabolic_ECS$GPP_day < 0]/-1000, col = "darkgreen", pch = 21)	
  points(Metabolic_ECS$Date, Metabolic_ECS$NEP_day/-1000, col = "darkblue", pch = 16)	
  legend(	
    "topright",	
    bty = 'n',	
    pch = c(22,21,16),	
    col = c("darkred", "darkgreen", "darkblue"),	
    legend = c("ER", "GPP", "NEP"),	
    cex = 1.3	
  )	
  legend(	
    "topleft",	
    bty = 'n',	
    legend = "A",	
    cex = 1.4	
  )	
  
  
  # plot the ecosytem production corrected for chlorophyll	
  par(mar = c(5, 5, 0, 0))	
  plot(Metabolic_ECS$Date[Metabolic_ECS$GPP_day < 0], Metabolic_ECS$GPP_day[Metabolic_ECS$GPP_day < 0]/Metabolic_ECS$Chl_lab[Metabolic_ECS$GPP_day < 0]/-1000, ylim = c(-0.012,0.012), col = "darkgreen", pch = 21,	
       ylab = expression("Metabolic rate [mmol "~mu*"g Chl"^"-1"*" d"^"-1"*"]"), xlab = "Time [Months]", cex.axis = 1.2, cex.lab = 1.5)	
  abline(v = as.Date("2016-08-10"), col = "darkgrey", lwd = 2.5, lty = 3)	
  abline(h = -0.0, col = "red", lwd = 2, lty = 3)	
  points(Metabolic_ECS$Date[Metabolic_ECS$R_day > 0], Metabolic_ECS$R_day[Metabolic_ECS$R_day > 0]/Metabolic_ECS$Chl_lab[Metabolic_ECS$R_day > 0]/-1000, col = "darkred", pch = 22)	
  points(Metabolic_ECS$Date[Metabolic_ECS$GPP_day < 0], Metabolic_ECS$GPP_day[Metabolic_ECS$GPP_day < 0]/Metabolic_ECS$Chl_lab[Metabolic_ECS$GPP_day < 0]/-1000, col = "darkgreen", pch = 21)	
  
  points(Metabolic_ECS$Date, Metabolic_ECS$NEP_day/Metabolic_ECS$Chl_lab/-1000, col = "darkblue", pch = 16)	
  legend(	
    "topright",	
    bty = 'n',	
    pch = c(22,21,16),	
    col = c("darkred", "darkgreen", "darkblue"),	
    legend = c("ER", "GPP", "NEP"),	
    cex = 1.3	
  )	
  
  legend(	
    "topleft",	
    bty = 'n',	
    legend = "D",	
    cex = 1.4	
  )	
  text(as.Date("2016-06-25"), -0.011, "Population increase", col = "black", cex = 1.1)	
  text(as.Date("2016-09-27"), -0.011, "Population decrease", col = "black", cex = 1.1)	
  
  
  
  # Plot the same for the modeled NEP	
  par(mar = c(3, 2, 2, 2.5))	
  plot(Metabolic_model$Date[Metabolic_model$GPP_day < 0], Metabolic_model$GPP_day[Metabolic_model$GPP_day < 0]/-1000, col = "darkgreen", pch = 21, ylim = c(-1.4,1.4),	
       ylab = expression(""), xlab = "", yaxt = 'n', xaxt = 'n',cex.axis = 1.2, cex.lab = 1.5)	
  abline(h = -0.0, col = "red", lwd = 2, lty = 3)	
  abline(v = as.Date("2016-08-10"), col = "darkgrey", lwd = 2.5, lty = 3)	
  points(Metabolic_model$Date[Metabolic_model$R_day > 0], Metabolic_model$R_day[Metabolic_model$R_day > 0]/-1000, col = "darkred", pch = 22)	
  points(Metabolic_model$Date[Metabolic_model$GPP_day < 0], Metabolic_model$GPP_day[Metabolic_model$GPP_day < 0]/-1000, col = "darkgreen", pch = 21)	
  points(Metabolic_model$Date, Metabolic_model$NEP_day/-1000, col = "darkblue", pch = 16)	
  legend(	
    "topright",	
    bty = 'n',	
    pch = c(22,21,16),	
    col = c("darkred", "darkgreen", "darkblue"),	
    legend = c("ER", "GPP", "NEP"),	
    cex = 1.3	
  )	
  legend(	
    "topleft",	
    bty = 'n',	
    legend = "B",	
    cex = 1.4	
  )	
  
  
  # plot the ecosytem production corrected for chlorophyll	
  par(mar = c(5, 2, 0, 2.5))	
  plot(Metabolic_model$Date[Metabolic_model$GPP_day < 0], Metabolic_model$GPP_day[Metabolic_model$GPP_day < 0]/Metabolic_model$Chl_lab[Metabolic_model$GPP_day < 0]/-1000, ylim = c(-0.012,0.012), col = "darkgreen", pch = 21,	
       ylab = expression(""), yaxt = 'n', xlab = "Time [Months]", cex.axis = 1.2, cex.lab = 1.5)	
  abline(v = as.Date("2016-08-10"), col = "darkgrey", lwd = 2.5, lty = 3)	
  abline(h = -0.0, col = "red", lwd = 2, lty = 3)	
  points(Metabolic_model$Date[Metabolic_model$R_day > 0], Metabolic_model$R_day[Metabolic_model$R_day > 0]/Metabolic_model$Chl_lab[Metabolic_model$R_day > 0]/-1000, col = "darkred", pch = 22)	
  points(Metabolic_model$Date[Metabolic_model$GPP_day < 0], Metabolic_model$GPP_day[Metabolic_model$GPP_day < 0]/Metabolic_model$Chl_lab[Metabolic_model$GPP_day < 0]/-1000, col = "darkgreen", pch = 21)	
  
  points(Metabolic_model$Date, Metabolic_model$NEP_day/Metabolic_model$Chl_lab/-1000, col = "darkblue", pch = 16)	
  legend(	
    "topright",	
    bty = 'n',	
    pch = c(22,21,16),	
    col = c("darkred", "darkgreen", "darkblue"),	
    legend = c("ER", "GPP", "NEP"),	
    cex = 1.3	
  )	
  legend(	
    "topleft",	
    bty = 'n',	
    legend = "E",	
    cex = 1.4	
  )	
  text(as.Date("2016-06-25"), -0.011, "Population increase", col = "black", cex = 1.1)	
  text(as.Date("2016-09-27"), -0.011, "Population decrease", col = "black", cex = 1.1)	
  
  
  
  
  
  # Plot the same for the DO modeled NEP	
  par(mar = c(3, 0, 2, 5))	
  plot(Metabolic_model_DO$Date[Metabolic_model_DO$GPP_day < 0], Metabolic_model_DO$GPP_day[Metabolic_model_DO$GPP_day < 0]/-1000, col = "darkgreen", pch = 21, ylim = c(-1.2,1.2),	
       ylab = expression(""), xlab = "", yaxt = 'n', xaxt = 'n',cex.axis = 1.2, cex.lab = 1.5)	
  abline(h = -0.0, col = "red", lwd = 2, lty = 3)	
  abline(v = as.Date("2016-08-10"), col = "darkgrey", lwd = 2.5, lty = 3)	
  points(Metabolic_model_DO$Date[Metabolic_model_DO$R_day > 0], Metabolic_model_DO$R_day[Metabolic_model_DO$R_day > 0]/-1000, col = "darkred", pch = 22)	
  points(Metabolic_model_DO$Date[Metabolic_model_DO$GPP_day < 0], Metabolic_model_DO$GPP_day[Metabolic_model_DO$GPP_day < 0]/-1000, col = "darkgreen", pch = 21)	
  points(Metabolic_model_DO$Date, Metabolic_model_DO$NEP_day/-1000, col = "darkblue", pch = 16)	
  legend(	
    "topright",	
    bty = 'n',	
    pch = c(22,21,16),	
    col = c("darkred", "darkgreen", "darkblue"),	
    legend = c("ER", "GPP", "NEP"),	
    cex = 1.3	
  )	
  legend(	
    "topleft",	
    bty = 'n',	
    legend = "C",	
    cex = 1.4	
  )	
  
  
  # plot the ecosytem production corrected for chlorophyll	
  par(mar = c(5, 0, 0, 5))	
  plot(Metabolic_model_DO$Date[Metabolic_model_DO$GPP_day < 0], Metabolic_model_DO$GPP_day[Metabolic_model_DO$GPP_day < 0]/Metabolic_model_DO$Chl_lab[Metabolic_model_DO$GPP_day < 0]/-1000, ylim = c(-0.012,0.012), col = "darkgreen", pch = 21,	
       ylab = expression(""), yaxt = 'n', xlab = "Time [Months]", cex.axis = 1.2, cex.lab = 1.5)	
  abline(v = as.Date("2016-08-10"), col = "darkgrey", lwd = 2.5, lty = 3)	
  abline(h = -0.0, col = "red", lwd = 2, lty = 3)	
  points(Metabolic_model_DO$Date[Metabolic_model_DO$R_day > 0], Metabolic_model_DO$R_day[Metabolic_model_DO$R_day > 0]/Metabolic_model_DO$Chl_lab[Metabolic_model_DO$R_day > 0]/-1000, col = "darkred", pch = 22)	
  points(Metabolic_model_DO$Date[Metabolic_model_DO$GPP_day < 0], Metabolic_model_DO$GPP_day[Metabolic_model_DO$GPP_day < 0]/Metabolic_model_DO$Chl_lab[Metabolic_model_DO$GPP_day < 0]/-1000, col = "darkgreen", pch = 21)	
  
  points(Metabolic_model_DO$Date, Metabolic_model_DO$NEP_day/Metabolic_model_DO$Chl_lab/-1000, col = "darkblue", pch = 16)	
  legend(	
    "topright",	
    bty = 'n',	
    pch = c(22,21,16),	
    col = c("darkred", "darkgreen", "darkblue"),	
    legend = c("ER", "GPP", "NEP"),	
    cex = 1.3	
  )	
  legend(	
    "topleft",	
    bty = 'n',	
    legend = "F",	
    cex = 1.3	
  )	
  text(as.Date("2016-06-25"), -0.011, "Population increase", col = "black", cex = 1.1)	
  text(as.Date("2016-09-27"), -0.011, "Population decrease", col = "black", cex = 1.1)	
}	


par(mfrow = c(1,1), mar = mar.default +1)	










library(lubridate)	
library(dplyr)	

# Link NEP to the fluxes	
Flux$Date <- as.Date(Flux$Date)	

library(plyr)	
NEP_Flux_ECS <- join(Metabolic_ECS, Flux, by = "Date")	
NEP_Flux_model <- join(Metabolic_model, Flux, by = "Date")	
NEP_Flux_model_DO <- join(Metabolic_model_DO, Flux, by = "Date")	


Flux$Date <- as.POSIXct(Flux$Date)	


# The fluxes are given in umol/m2/h	
# so convert them to mmol/L/d or NEP to umol/m2/h	

par(mfrow = c(3,1), xpd = FALSE)	
par(mar = mar.default +1)	
plot(NEP_Flux_ECS$NEP_day[-c(45,49)]*-1.8, NEP_Flux_ECS$mean_flux_ECS[-c(45,49)]/1000*24)	
abline(h=0, v=0)	
plot(NEP_Flux_model$NEP_day*-1.8, NEP_Flux_model$mean_flux_MacIntyre/1000*24)	
abline(h=0, v=0)	
plot(NEP_Flux_model_DO$NEP_day*-1.8, NEP_Flux_model_DO$mean_flux_DO/1000*24)	
abline(h=0, v=0)	



plot(NEP_Flux_model_DO$mean_flux_DO[NEP_Flux_model_DO$NEP_day < 0]/1000, NEP_Flux_model_DO$mean_flux_MacIntyre[NEP_Flux_model_DO$NEP_day < 0]/1000, col = "darkgreen")	
abline(h=0, v=0)	
points(NEP_Flux_model_DO$mean_flux_DO[NEP_Flux_model_DO$NEP_day > 0]/1000, NEP_Flux_model_DO$mean_flux_MacIntyre[NEP_Flux_model_DO$NEP_day > 0]/1000, col = "darkred", pch = 19)	


# the daily average flux and NEP in mmol/m2/d	

NEP_Flux_ECS$mean_flux_ECS_mmold	
NEP_Flux_ECS$NEP_day_mmold <- NEP_Flux_ECS$NEP_day*-1.8	

par(mfrow = c(2,2))	
plot(lm(NEP_Flux_ECS$mean_flux_ECS_mmold[-49] ~  (NEP_Flux_ECS$NEP_day_mmold[-49])))	
summary(lm(NEP_Flux_ECS$mean_flux_ECS_mmold[-49] ~  NEP_Flux_ECS$NEP_day_mmold[-49]))	

plot(lm(NEP_Flux_model$mean_flux_MacIntyre/1000*24 ~  NEP_Flux_model$NEP_day*1.8))	
summary(lm(NEP_Flux_model$mean_flux_MacIntyre/1000*24 ~  NEP_Flux_model$NEP_day*1.8))	

NEP_Flux_model_DO$mean_flux_DO_mmold <- NEP_Flux_model_DO$mean_flux_DO/1000*24	
NEP_Flux_model_DO$NEP_day_mmold <- NEP_Flux_model_DO$NEP_day*-1.8	

plot(lm(NEP_Flux_model_DO$mean_flux_DO_mmold ~  NEP_Flux_model_DO$NEP_day_mmold))	
summary(lm(NEP_Flux_model_DO$mean_flux_DO_mmold ~  NEP_Flux_model_DO$NEP_day_mmold))	



# Check linear relation NEP with the temperature and wind	


summary(lm(NEP_Flux_ECS$GPP_day ~ NEP_Flux_ECS$mean_flux_ECS))	

par(mfrow = c(2,2))	
plot(lm(NEP_Flux_ECS$mean_flux_ECS ~  NEP_Flux_ECS$mean_wind * NEP_Flux_ECS$mean_BF))	
summary(lm(NEP_Flux_ECS$mean_flux_ECS ~ NEP_Flux_ECS$mean_wind * NEP_Flux_ECS$mean_BF))	

plot(lm(NEP_Flux_ECS$NEP_day ~ NEP_Flux_ECS$mean_wind * NEP_Flux_ECS$mean_BF))	
summary(lm(NEP_Flux_ECS$NEP_day ~ NEP_Flux_ECS$mean_wind * NEP_Flux_ECS$mean_BF))	

NEP_Flux_ECS$changeChl[2:132] <- diff(log(NEP_Flux_ECS$median_Chl_lab))	

plot(lm(NEP_Flux_ECS$NEP_day ~ NEP_Flux_ECS$mean_wind * NEP_Flux_ECS$mean_BF * NEP_Flux_ECS$median_Chl_lab ))	


library(nlme)	


NEP_Flux_ECS <- NEP_Flux_ECS[complete.cases(NEP_Flux_ECS[c("median_Chl_lab", "changeChl", "mean_BF", "mean_wind", "NEP_day")]),]	

summary(lme(NEP_day ~ mean_wind * mean_BF, random= ~1|median_Chl_lab/changeChl, data = NEP_Flux_ECS,  method = "ML"))	
summary(lm(NEP_Flux_ECS$mean_flux_ECS ~ NEP_Flux_ECS$mean_wind * NEP_Flux_ECS$mean_BF ))	



summary(lm(NEP_Flux_ECS$NEP_day ~ NEP_Flux_ECS$mean_wind * NEP_Flux_ECS$mean_BF * NEP_Flux_ECS$median_Chl_lab * NEP_Flux_ECS$changeChl))	
summary(lm(NEP_Flux_ECS$mean_flux_ECS ~ NEP_Flux_ECS$mean_wind * NEP_Flux_ECS$mean_BF * NEP_Flux_ECS$median_Chl_lab * NEP_Flux_ECS$changeChl))	


plot(lm(sign(NEP_Flux_ECS$NEP_day)*log(abs(NEP_Flux_ECS$NEP_day/NEP_Flux_ECS$median_Chl_lab)) ~ NEP_Flux_ECS$mean_BF))	
summary(lm(sign(NEP_Flux_ECS$NEP_day)*log(abs(NEP_Flux_ECS$NEP_day/NEP_Flux_ECS$median_Chl_lab)) ~ NEP_Flux_ECS$mean_BF * NEP_Flux_ECS$mean_wind))	

plot(lm(sign(NEP_Flux_model$NEP_day)*log(abs(NEP_Flux_model$NEP_day/NEP_Flux_model$median_Chl_lab)) ~ NEP_Flux_model$mean_BF))	
summary(lm(sign(NEP_Flux_model$NEP_day)*log(abs(NEP_Flux_model$NEP_day/NEP_Flux_model$median_Chl_lab)) ~ NEP_Flux_model$mean_BF * NEP_Flux_model$mean_wind))	








# in mmol/m2/d	

library(dplyr)	
NEP_Flux_model %>%	
  group_by(week(Date)) %>%	
  summarise(	
    NEP_day = mean(NEP_day*-1*1.8, na.rm = TRUE),	
    GPP_day = mean(GPP_day*-1*1.8, na.rm = TRUE),	
    R_day = mean(GPP_day*-1*1.8, na.rm = TRUE),	
    mean_flux_MacIntyre = mean(mean_flux_MacIntyre/1000*24, na.rm = TRUE),	
    mean_flux_ECS = mean(mean_flux_ECS/1000*24, na.rm = TRUE),	
    mean_flux_DO = mean(mean_flux_DO/1000*24, na.rm = TRUE)) -> NEP_Flux_model_week	


library(dplyr)	
NEP_Flux_model_DO %>%	
  group_by(week(Date)) %>%	
  summarise(	
    NEP_day = mean(NEP_day*-1*1.8, na.rm = TRUE),	
    GPP_day = mean(GPP_day*-1*1.8, na.rm = TRUE),	
    R_day = mean(GPP_day*-1*1.8, na.rm = TRUE),	
    mean_flux_MacIntyre = mean(mean_flux_MacIntyre/1000*24, na.rm = TRUE),	
    mean_flux_ECS = mean(mean_flux_ECS/1000*24, na.rm = TRUE),	
    mean_flux_DO = mean(mean_flux_DO/1000*24, na.rm = TRUE)) -> NEP_Flux_model_DO_week	


library(dplyr)	
NEP_Flux_ECS %>%	
  group_by(week(Date)) %>%	
  summarise(	
    NEP_day = mean(NEP_day*-1*1.8, na.rm = TRUE),	
    GPP_day = mean(GPP_day*-1*1.8, na.rm = TRUE),	
    R_day = mean(GPP_day*-1*1.8, na.rm = TRUE),	
    mean_flux_MacIntyre = mean(mean_flux_MacIntyre/1000*24, na.rm = TRUE),	
    mean_flux_ECS = mean(mean_flux_ECS/1000*24, na.rm = TRUE),	
    mean_flux_DO = mean(mean_flux_DO/1000*24, na.rm = TRUE)) -> NEP_Flux_ECS_week	



# plot the model 	
par(mfrow = c(2,2)) 	
plot(NEP_Flux_model_week$GPP_day, NEP_Flux_model_week$mean_flux_MacIntyre)	
abline(h=0, v=0)	
plot(NEP_Flux_model_week$NEP_day, NEP_Flux_model_week$mean_flux_MacIntyre)	
abline(h=0, v=0)	
plot(NEP_Flux_model_DO_week$NEP_day, NEP_Flux_model_DO_week$mean_flux_DO)	
abline(h=0, v=0)	
plot(NEP_Flux_model_week$mean_flux_MacIntyre, NEP_Flux_model_week$mean_flux_DO)	
abline(h=0, v=0)	

#plot the ECS	
par(mfrow = c(2,2)) 	
plot(NEP_Flux_ECS_week$GPP_day, NEP_Flux_ECS_week$mean_flux_ECS)	
abline(h=0, v=0)	
plot(NEP_Flux_ECS_week$NEP_day, NEP_Flux_ECS_week$mean_flux_ECS)	
abline(h=0, v=0)	
plot(NEP_Flux_ECS_week$NEP_day, NEP_Flux_ECS_week$mean_flux_DO)	
abline(h=0, v=0)	
plot(NEP_Flux_ECS_week$mean_flux_MacIntyre, NEP_Flux_ECS_week$mean_flux_ECS)	
abline(h=0, v=0)	




par(mfrow = c(2,2))	
plot(lm(NEP_Flux_ECS_week$mean_flux_ECS ~ NEP_Flux_ECS_week$NEP_day))	
summary(lm(NEP_Flux_ECS_week$mean_flux_ECS[-c(15,4)] ~ NEP_Flux_ECS_week$NEP_day[-c(15,4)]))	
plot(lm(NEP_Flux_ECS_week$mean_flux_ECS ~ NEP_Flux_ECS_week$GPP_day))	
summary(lm(NEP_Flux_ECS_week$mean_flux_ECS ~ NEP_Flux_ECS_week$GPP_day))	
plot(lm(NEP_Flux_ECS_week$mean_flux_ECS ~ NEP_Flux_ECS_week$R_day))	
summary(lm(NEP_Flux_ECS_week$mean_flux_ECS ~ NEP_Flux_ECS_week$R_day))	

# The models do not fit as good as with the ECS data	
plot(lm(NEP_Flux_model_week$mean_flux_MacIntyre ~ NEP_Flux_model_week$NEP_day))	
summary(lm(NEP_Flux_model_week$mean_flux_MacIntyre ~ NEP_Flux_model_week$NEP_day))	
plot(lm(NEP_Flux_model_week$mean_flux_MacIntyre ~ NEP_Flux_model_week$GPP_day))	
summary(lm(NEP_Flux_model_week$mean_flux_MacIntyre ~ NEP_Flux_model_week$GPP_day))	
plot(lm(NEP_Flux_model_week$mean_flux_MacIntyre ~ NEP_Flux_model_week$R_day))	
summary(lm(NEP_Flux_model_week$mean_flux_MacIntyre ~ NEP_Flux_model_week$R_day))	



plot(lm(NEP_Flux_model_DO_week$mean_flux_MacIntyre ~ NEP_Flux_model_DO_week$NEP_day))	
summary(lm(NEP_Flux_model_DO_week$mean_flux_MacIntyre[-5] ~ NEP_Flux_model_DO_week$NEP_day[-5]))	
summary(lm(NEP_Flux_model_DO_week$mean_flux_MacIntyre[-5] ~ NEP_Flux_model_DO_week$GPP_day[-5]))	


lm_NEP_ECS <- lm(NEP_Flux_ECS_week$mean_flux_ECS[-c(15,4)] ~ NEP_Flux_ECS_week$NEP_day[-c(15,4)])	
lm_NEP_mod <- lm(NEP_Flux_model_week$mean_flux_MacIntyre[-16] ~ NEP_Flux_model_week$NEP_day[-16])	
lm_NEP_DO <- lm(NEP_Flux_model_DO_week$mean_flux_DO[-c(16,10,5)] ~ NEP_Flux_model_DO_week$NEP_day[-c(16,10,5)])	

lm_GPP_ECS <- lm(NEP_Flux_ECS_week$mean_flux_ECS ~ NEP_Flux_ECS_week$GPP_day)	
lm_GPP_mod <- lm(NEP_Flux_model_week$mean_flux_ECS ~ NEP_Flux_model_week$GPP_day)	
lm_GPP_DO <- lm(NEP_Flux_model_DO_week$mean_flux_DO[-5] ~ NEP_Flux_model_DO_week$GPP_day[-5])	

lm_R_ECS <- lm(NEP_Flux_ECS_week$mean_flux_ECS ~ NEP_Flux_ECS_week$R_day)	
lm_R_mod <- lm(NEP_Flux_model_week$mean_flux_ECS ~ NEP_Flux_model_week$R_day)	
lm_R_DO <- lm(NEP_Flux_model_DO_week$mean_flux_DO[-5] ~ NEP_Flux_model_DO_week$R_day[-5])	


plot(lm_NEP_ECS)	
plot(lm_NEP_mod)	
plot(lm_NEP_DO)	

plot(lm_GPP_ECS)	
plot(lm_GPP_mod)	
plot(lm_GPP_DO)	

plot(lm_R_ECS)	
plot(lm_R_mod)	
plot(lm_R_DO)	


lmtest::bptest(lm_NEP_ECS)	
lmtest::bptest(lm_NEP_mod)	
lmtest::bptest(lm_NEP_DO)	

lmtest::bptest(lm_GPP_ECS)	
lmtest::bptest(lm_GPP_mod)	
lmtest::bptest(lm_GPP_DO)	

lmtest::bptest(lm_R_ECS)	
lmtest::bptest(lm_R_mod)	
lmtest::bptest(lm_R_DO)	

# No heteroscedasticity	

shapiro.test(residuals(lm_NEP_ECS))	
shapiro.test(residuals(lm_NEP_mod))	
shapiro.test(residuals(lm_NEP_DO))	

shapiro.test(residuals(lm_GPP_ECS))	
shapiro.test(residuals(lm_GPP_mod))	
shapiro.test(residuals(lm_GPP_DO))	

shapiro.test(residuals(lm_R_ECS))	
shapiro.test(residuals(lm_R_mod))	
shapiro.test(residuals(lm_R_DO))	





summary(lm_NEP_ECS)	
summary(lm_NEP_mod)	
summary(lm_NEP_DO)	

summary(lm_GPP_ECS)	
summary(lm_GPP_mod)	
summary(lm_GPP_DO)	

summary(lm_R_ECS)	
summary(lm_R_mod)	
summary(lm_R_DO)	


plot(NEP_Flux_ECS_week$mean_flux_ECS ~ NEP_Flux_ECS_week$NEP_day)	
plot(NEP_Flux_model_week$mean_flux_ECS ~ NEP_Flux_model_week$NEP_day)	
plot(NEP_Flux_ECS_week$mean_flux_ECS ~ NEP_Flux_ECS_week$GPP_day)	
plot(NEP_Flux_model_week$mean_flux_ECS[-15] ~ NEP_Flux_model_week$GPP_day[-15])	
plot(NEP_Flux_ECS_week$mean_flux_ECS ~ NEP_Flux_ECS_week$R_day)	
plot(NEP_Flux_model_week$mean_flux_ECS[-15] ~ NEP_Flux_model_week$R_day[-15])	

met_sub <- met[c("Date",names(Daily_DIC_in_out)[c(2:6)],"NEP_day_ECS","NEP_day_model","NEP_day_DO",	
                 "GPP_day_ECS","GPP_day_model","GPP_day_DO",	
                 "R_day_ECS","R_day__model","R_day_DO")]	

# convert units to mmol m2 d	
library(dplyr)	
met_sub %>%	
  group_by(week(Date)) %>%	
  summarise(	
    Inlet = mean(Inlet*24*-1.8, na.rm = TRUE),	
    NEP_day_ECS = mean(NEP_day_ECS*1.8, na.rm = TRUE),	
    NEP_day_model = mean(NEP_day_model*1.8, na.rm = TRUE),	
    NEP_day_DO = mean(NEP_day_DO*1.8, na.rm = TRUE),	
    GPP_day_ECS = mean(GPP_day_ECS*1.8, na.rm = TRUE),	
    GPP_day_model = mean(GPP_day_model*1.8, na.rm = TRUE),	
    GPP_day_DO = mean(GPP_day_DO*1.8, na.rm = TRUE),	
    R_day_ECS = mean(R_day_ECS*1.8, na.rm = TRUE),	
    R_day_model = mean(R_day__model*1.8, na.rm = TRUE),	
    R_day_DO = mean(R_day_DO*1.8, na.rm = TRUE)) -> met_week	




plot(met_sub$Inlet,met_sub$NEP_day_ECS)	
plot(met_sub$Inlet,met_sub$NEP_day_model)	
plot(met_sub$Inlet,met_sub$NEP_day_DO)	


par(mfrow = c(2,2))	
plot(lm(met_sub$NEP_day_ECS~ met_sub$Inlet))	
plot(lm(log(met_sub$NEP_day_model) ~ met_sub$Daily_DIC_in_out))	
plot(lm(met_sub$NEP_day_DO ~ met_sub$Daily_DIC_in_out))	
par(mfrow = c(1,1))	

lm_DIC_EC <- (lm(met$NEP_day_ECS_m2 ~ met$Daily_DIC_in_out_m2))	
lm_DIC_model <- (lm(met$NEP_day_model_m2 ~ met$Daily_DIC_in_out_m2))	
lm_DIC_DO <- (lm(met$NEP_day_DO_m2[-5] ~ met$Daily_DIC_in_out_m2[-5]))	


summary(lm_DIC_EC)	
summary(lm_DIC_model)	
summary(lm_DIC_DO)	

met$GPP_day_ECS_m2 <- met$GPP_day_ECS*-1.8	
met$GPP_day_model_m2 <- met$GPP_day_model*-1.8	
met$GPP_day_DO_m2 <- met$GPP_day_DO*-1.8	


met$R_day_ECS_m2 <- met$R_day_ECS*1.8	
met$R_day_model_m2 <- met$R_day__model*1.8	
met$R_day_DO_m2 <- met$R_day_DO*1.8	


plot(lm(met$GPP_day_ECS_m2 ~ met$Inlet_m2))	
plot(lm(met$GPP_day_model_m2 ~ met$Inlet_m2))	
plot(lm(met$GPP_day_DO_m2 ~ met$Inlet_m2))	

summary(lm(met$GPP_day_ECS_m2 ~ met$Inlet_m2))	
summary(lm(met$GPP_day_model_m2 ~ met$Inlet_m2))	
summary(lm(met$GPP_day_DO_m2 ~ met$Inlet_m2))	



plot(lm(log(met$R_day_ECS_m2) ~ met$Inlet_m2))	
plot(lm(log(met$R_day_model_m2) ~ met$Inlet_m2))	
plot(lm(log(met$R_day_DO_m2) ~ met$Inlet_m2))	

summary(lm(log(met$R_day_ECS_m2) ~ met$Inlet_m2))	
summary(lm(log(met$R_day_model_m2) ~ met$Inlet_m2))	
summary(lm(log(met$R_day_DO_m2) ~ met$Inlet_m2))	




summary(lm(log(met_week$GPP_day_ECS*-1) ~ met_week$Inlet))	
summary(lm(log(met_week$GPP_day_model*-1) ~ met_week$Inlet))	
summary(lm(met_week$GPP_day_DO*-1 ~ met_week$Inlet))	


summary(lm(log(met_week$R_day_ECS[-c(3,25)]) ~ met_week$Inlet[-c(3,25)]))	
summary(lm(met_week$R_day_model ~ met_week$Inlet))	
summary(lm(met_week$R_day_DO ~ met_week$Inlet))	


# Check the relation between DIC inlfow and atmospheric flux	
# in mmol m2 d	

# Creates a datafiles by day 	
library(dplyr)	
DF_Metabolic %>%	
  group_by(Date) %>%	
  summarise(	
    Flux_CO2_ECS = mean(CO2FluxECS_umol_m2_h/1000*24, na.rm = TRUE),	
    Flux_CO2_model = mean(CO2flux_macIntyreR/1000*24, na.rm = TRUE),	
    Flux_O2_model = mean(DOflux_umol/1000*24, na.rm = TRUE),	
    DIC_in_out = mean(`DIC_IN-OUT_umol/h`*1800/1000*24, na.rm = TRUE) 	
  )-> DIC_InOut_Flux	

DIC_InOut_Flux$Date <- as.Date(DIC_InOut_Flux$Date)	

DICIn_Flux <- join(met_sub, DIC_InOut_Flux, by = "Date")	



# Creates a datafiles by week	
library(dplyr)	
DICIn_Flux %>%	
  group_by(lubridate::week(Date)) %>%	
  summarise(	
    Flux_CO2_ECS = mean(CO2FluxECS_umol_m2_h/1000*24, na.rm = TRUE),	
    Flux_CO2_model = mean(CO2flux_macIntyreR/1000*24, na.rm = TRUE),	
    Flux_O2_model = mean(DOflux_umol/1000*24, na.rm = TRUE),	
    Inlet_m2 = mean(`DIC_IN-OUT_umol/h`*1800/1000*24, na.rm = TRUE) 	
  )-> DIC_InOut_Flux	



plot(DICIn_Flux$Flux_CO2_ECS~DICIn_Flux$Inlet_m2)	
plot(DICIn_Flux$Flux_CO2_model~DICIn_Flux$Inlet_m2)	
plot(DICIn_Flux$Flux_O2_model~DICIn_Flux$Inlet_m2)	


par(mfrow = c(2,2))	
plot(lm(DICIn_Flux$Flux_CO2_ECS~DICIn_Flux$Inlet_m2))   # Transformation doesn't make the model better	
plot(lm(DICIn_Flux$Flux_CO2_model~DICIn_Flux$Inlet_m2))	
plot(lm(DICIn_Flux$Flux_O2_model~DICIn_Flux$Inlet_m2))	
par(mfrow = c(1,1))	

DICIn_Flux$Inlet_m2 <- DICIn_Flux$Inlet*-1.8	

lm_ECS_DICInOut <- (lm(DICIn_Flux$Flux_CO2_ECS/1000~DICIn_Flux$Inlet_m2))	
lm_model_DICInOut <- (lm(DICIn_Flux$Flux_CO2_model/1000~DICIn_Flux$Inlet_m2))	
lm_O2_DICInOut <- (lm(DICIn_Flux$Flux_O2_model/1000~DICIn_Flux$Inlet_m2))	

summary(lm_ECS_DICInOut)	
summary(lm_model_DICInOut)	
summary(lm_O2_DICInOut)	

##############################################	


# plot NEP against the DIC in-out	


{	
  par(mfrow = c(1,2))	
  par(mar = c(5, 5, 2, 2))	
  plot(met$Daily_DIC_in_out_m2, met$NEP_day_ECS_m2, ylim = c(-700,800), xlim = c(-90,90),	
       ylab = expression("NEP-C"[EC]*" [mmol m"^"-2"*"d"^"-1"*"]"), xlab = expression("DIC flow [mmol m"^"-2"*"d"^"-1"*"]"), cex.lab = 1.3, cex.axis = 1.3, pch = 19, col = "blue")	
  rect(xleft = 0, xright = 350, ytop = 1000, ybottom = 0, col = rgb(0.1,0.7,0.0, alpha = 0.1) , border = NA)	
  rect(xleft = -200, xright = 0, ytop = 0, ybottom = -1000, col = rgb(0.7,0.1,0.0, alpha = 0.1), border = NA)	
  abline(h=0, v=0, lty = 3)	
  abline(lm_DIC_EC, col = "red")	
  abline(summary(lm_DIC_EC)$coefficients[,1] + summary(lm_DIC_EC)$coefficients[,2]*1.96, col = "red", lty = 3)	
  abline(summary(lm_DIC_EC)$coefficients[,1] - summary(lm_DIC_EC)$coefficients[,2]*1.96, col = "red", lty = 3)	
  text(x = 55, y = 750, "Autotrophy", cex = 1.1)	
  text(x = -55, y = -650, "Heterotrophy", cex = 1.1)	
  legend(	
    "topleft",	
    legend = "A",	
    bty = 'n', 	
    cex = 1.4	
  )	
  
  par(mar = c(5, 5, 2, 2))	
  plot(met$Daily_DIC_in_out_m2, met$NEP_day_model_m2, ylim = c(-700,800), xlim = c(-90,90),	
       ylab = expression("NEP-C"[model]*" [mmol m"^"-2"*"d"^"-1"*"]"), xlab = expression("DIC flow [mmol m"^"-2"*"d"^"-1"*"]"), cex.lab = 1.3, cex.axis = 1.3, pch = 1, col = "black")	
  rect(xleft = 0, xright = 350, ytop = 1000, ybottom = 0, col = rgb(0.1,0.7,0.0, alpha = 0.1) , border = NA)	
  rect(xleft = -200, xright = 0, ytop = 0, ybottom = -1000, col = rgb(0.7,0.1,0.0, alpha = 0.1), border = NA)	
  abline(h=0, v=0, lty = 3)	
  abline(lm_DIC_model, col = "red")	
  abline(summary(lm_DIC_model)$coefficients[,1] + summary(lm_DIC_model)$coefficients[,2]*1.96, col = "red", lty = 3)	
  abline(summary(lm_DIC_model)$coefficients[,1] - summary(lm_DIC_model)$coefficients[,2]*1.96, col = "red", lty = 3)	
  text(x = 55, y = 750, "Autotrophy", cex = 1.1)	
  text(x = -55, y = -650, "Heterotrophy", cex = 1.1)	
  legend(	
    "topleft",	
    legend = "B",	
    bty = 'n', 	
    cex = 1.4)	
}	

# Plot the fluxes against the NEP	


# plot the NEP for the ECS and the modeled flux	



{	
  par(mfrow = c(1,2))	
  par(mar = c(5, 5, 2, 2))	
  plot(NEP_Flux_ECS_week$NEP_day[-c(15,4)], NEP_Flux_ECS_week$mean_flux_ECS[-c(15,4)], ylim = c(-300, 100), xlim = c(-100, 300),	
       ylab = expression("Weekly averaged CO"[2]*" flux [mmol m"^"-2"*"d"^"-1"*"]"), xlab = expression("NEP [mmol m"^"-2"*"d"^"-1"*"]"), cex.lab = 1.3, cex.axis = 1.3, pch = 19, col = "blue")	
  rect(xleft = 0, xright = 350, ytop = 0, ybottom = -400, col = rgb(0.1,0.7,0.0, alpha = 0.1) , border = NA)	
  rect(xleft = -200, xright = 0, ytop = 350, ybottom = 0, col = rgb(0.7,0.1,0.0, alpha = 0.1), border = NA)	
  points(NEP_Flux_model_week$NEP_day[-16], NEP_Flux_model_week$mean_flux_MacIntyre[-16])	
  abline(h=0, v=0, lty = 3)	
  abline(lm_NEP_ECS, col = "red")	
  abline(summary(lm_NEP_ECS)$coefficients[,1] + summary(lm_NEP_ECS)$coefficients[,2]*1.96, col = "red", lty = 3)	
  abline(summary(lm_NEP_ECS)$coefficients[,1] - summary(lm_NEP_ECS)$coefficients[,2]*1.96, col = "red", lty = 3)	
  text(x = 270, y = -290, "Autotrophy", cex = 1.1)	
  text(x = -55, y = 90, "Heterotrophy", cex = 1.1)	
  legend(	
    "topright",	
    legend = "A",	
    bty = 'n', 	
    cex = 1.4)	
  legend(	
    "bottomleft",	
    legend = c(expression("Modeled CO"[2]~"flux"), expression("Measured CO"[2]~"flux")),	
    cex = 1.1,	
    pch = c(16,1),	
    col = c("blue", "black"),	
    bty = 'n'	
  )	
  
  par(mar = c(5, 5, 2, 2))	
  plot(NEP_Flux_model_DO_week$NEP_day[-c(16,10,5)], NEP_Flux_model_DO_week$mean_flux_DO[-c(16,10,5)], ylim = c(-100, 300), xlim = c(-100, 300),	
       ylab = expression("Weekly averaged O"[2]*" flux [mmol m"^"-2"*"d"^"-1"*"]"), xlab = expression("NEP [mmol m"^"-2"*"d"^"-1"*"]"), cex.lab = 1.3, cex.axis = 1.3, pch = 1, col = "black")	
  rect(xleft = 0, xright = 350, ytop = 350, ybottom = 0, col = rgb(0.1,0.7,0.0, alpha = 0.1) , border = NA)	
  rect(xleft = -150, xright = 0, ytop = 0, ybottom = -150, col = rgb(0.7,0.1,0.0, alpha = 0.1), border = NA)	
  abline(h=0, v=0, lty = 3)	
  abline(lm_NEP_DO, col = "red")	
  abline(summary(lm_NEP_DO)$coefficients[,1] + summary(lm_NEP_DO)$coefficients[,2]*1.96, col = "red", lty = 3)	
  abline(summary(lm_NEP_DO)$coefficients[,1] - summary(lm_NEP_DO)$coefficients[,2]*1.96, col = "red", lty = 3)	
  text(x = 240, y = 290, "Autotrophy", cex = 1.1)	
  text(x = -55, y = -90, "Heterotrophy", cex = 1.1)	
  legend(	
    "topright",	
    legend = "B",	
    bty = 'n', 	
    cex = 1.4)	
}	



# The DIC inflow and the flux	

{	
  par(mfrow = c(1,3))	
  par(mar = c(5, 5, 2, 2))	
  plot(DIC_InOut_Flux$DIC_in_out, DIC_InOut_Flux$Flux_CO2_ECS,ylim = c(-300,100),	
       ylab = expression("Measured CO"[2]~"flux [mmol m"^"-2"*"d"^"-1"*"]"), xlab = expression("DIC flow [mmol m"^"-2"*"d"^"-1"*"]"), cex.lab = 1.3, cex.axis = 1.3, pch = 19, col = "blue")	
  abline(h=0, v=0, lty = 3)	
  abline(lm_ECS_DICInOut, col = "red")	
  abline(summary(lm_ECS_DICInOut)$coefficients[,1] + summary(lm_ECS_DICInOut)$coefficients[,2]*1.96, col = "red", lty = 3)	
  abline(summary(lm_ECS_DICInOut)$coefficients[,1] - summary(lm_ECS_DICInOut)$coefficients[,2]*1.96, col = "red", lty = 3)	
  
  legend(	
    inset = c(-0.08,-0.02),	
    "topleft",	
    legend = "A",	
    bty = 'n', 	
    cex = 1.4	
  )	
  
  par(mar = c(5, 5, 2, 2))	
  plot(DIC_InOut_Flux$DIC_in_out, DIC_InOut_Flux$Flux_CO2_model, ylim = c(-300,100),	
       ylab = expression("Modeled CO"[2]~"flux [mmol m"^"-2"*"d"^"-1"*"]"), xlab = expression("DIC flow [mmol m"^"-2"*"d"^"-1"*"]"), cex.lab = 1.3, cex.axis = 1.3, pch = 1, col = "black")	
  abline(h=0, v=0, lty = 3)	
  abline(lm_model_DICInOut, col = "red")	
  abline(summary(lm_model_DICInOut)$coefficients[,1] + summary(lm_DIC_model)$coefficients[,2]*1.96, col = "red", lty = 3)	
  abline(summary(lm_model_DICInOut)$coefficients[,1] - summary(lm_DIC_model)$coefficients[,2]*1.96, col = "red", lty = 3)	
  legend(	
    inset = c(-0.08,-0.02),	
    "topleft",	
    legend = "B",	
    bty = 'n', 	
    cex = 1.4)	
  
  par(mar = c(5, 5, 2, 2))	
  plot(DIC_InOut_Flux$DIC_in_out, DIC_InOut_Flux$Flux_O2_model, ylim = c(-100,300),	
       ylab = expression("Modeled O"[2]~"flux [mmol m"^"-2"*"d"^"-1"*"]"), xlab = expression("DIC flow [mmol m"^"-2"*"d"^"-1"*"]"), cex.lab = 1.3, cex.axis = 1.3, pch = 1, col = "black")	
  abline(h=0, v=0, lty = 3)	
  abline(lm_O2_DICInOut, col = "red")	
  abline(summary(lm_O2_DICInOut)$coefficients[,1] + summary(lm_DIC_model)$coefficients[,2]*1.96, col = "red", lty = 3)	
  abline(summary(lm_O2_DICInOut)$coefficients[,1] - summary(lm_DIC_model)$coefficients[,2]*1.96, col = "red", lty = 3)	
  legend(	
    inset = c(-0.08,-0.02),	
    "topleft",	
    legend = "C",	
    bty = 'n', 	
    cex = 1.4)	
}	







# Check the average EP during before and after the peak 	
t.test(Metabolic_ECS$GPP_day[Metabolic_ECS$GPP_day < 0 & Metabolic_ECS$Date < as.Date("2016-08-10")]/1000, Metabolic_ECS$GPP_day[Metabolic_ECS$GPP_day < 0 & Metabolic_ECS$Date > as.Date("2016-08-10")]/1000)	
t.test(Metabolic_ECS$R_day[Metabolic_ECS$R_day > 0 & Metabolic_ECS$Date < as.Date("2016-08-10")]/1000, Metabolic_ECS$R_day[Metabolic_ECS$R_day > 0 & Metabolic_ECS$Date > as.Date("2016-08-10")]/1000)	
t.test(Metabolic_ECS$NEP_day[ Metabolic_ECS$Date < as.Date("2016-08-10")]/1000, Metabolic_ECS$NEP_day[Metabolic_ECS$Date > as.Date("2016-08-10")]/1000)	


t.test(Metabolic_ECS$GPP_day[Metabolic_ECS$GPP_day < 0 & Metabolic_ECS$Date < as.Date("2016-08-10")]/Metabolic_ECS$Chl_lab[Metabolic_ECS$GPP_day < 0 & Metabolic_ECS$Date < as.Date("2016-08-10")], 	
       Metabolic_ECS$GPP_day[Metabolic_ECS$GPP_day < 0 & Metabolic_ECS$Date > as.Date("2016-08-10")]/Metabolic_ECS$Chl_lab[Metabolic_ECS$GPP_day < 0 & Metabolic_ECS$Date > as.Date("2016-08-10")])	
t.test(Metabolic_ECS$R_day[Metabolic_ECS$R_day > 0 & Metabolic_ECS$Date < as.Date("2016-08-10")]/Metabolic_ECS$Chl_lab[Metabolic_ECS$R_day > 0 & Metabolic_ECS$Date < as.Date("2016-08-10")], 	
       Metabolic_ECS$R_day[Metabolic_ECS$R_day > 0 & Metabolic_ECS$Date > as.Date("2016-08-10")]/Metabolic_ECS$Chl_lab[Metabolic_ECS$R_day > 0 & Metabolic_ECS$Date > as.Date("2016-08-10")])	
t.test(Metabolic_ECS$NEP_day[Metabolic_ECS$GPP_day < 0 & Metabolic_ECS$Date < as.Date("2016-08-10")]/Metabolic_ECS$Chl_lab[Metabolic_ECS$GPP_day < 0 & Metabolic_ECS$Date < as.Date("2016-08-10")], 	
       Metabolic_ECS$NEP_day[Metabolic_ECS$GPP_day < 0 & Metabolic_ECS$Date > as.Date("2016-08-10")]/Metabolic_ECS$Chl_lab[Metabolic_ECS$GPP_day < 0 & Metabolic_ECS$Date > as.Date("2016-08-10")])	




t.test(Metabolic_model$GPP_day[Metabolic_model$GPP_day < 0 & Metabolic_model$Date < as.Date("2016-08-10")]/1000, Metabolic_model$GPP_day[Metabolic_model$GPP_day < 0 & Metabolic_model$Date > as.Date("2016-08-10")]/1000)	
t.test(Metabolic_model$R_day[Metabolic_model$R_day > 0 & Metabolic_model$Date < as.Date("2016-08-10")]/1000, Metabolic_model$R_day[Metabolic_model$R_day > 0 & Metabolic_model$Date > as.Date("2016-08-10")]/1000)	
t.test(Metabolic_model$NEP_day[ Metabolic_model$Date < as.Date("2016-08-10")]/1000, Metabolic_model$NEP_day[Metabolic_model$Date > as.Date("2016-08-10")]/1000)	


t.test(Metabolic_model$GPP_day[Metabolic_model$GPP_day < 0 & Metabolic_model$Date < as.Date("2016-08-10")]/Metabolic_model$Chl_lab[Metabolic_model$GPP_day < 0 & Metabolic_model$Date < as.Date("2016-08-10")], 	
       Metabolic_model$GPP_day[Metabolic_model$GPP_day < 0 & Metabolic_model$Date > as.Date("2016-08-10")]/Metabolic_model$Chl_lab[Metabolic_model$GPP_day < 0 & Metabolic_model$Date > as.Date("2016-08-10")])	
t.test(Metabolic_model$R_day[Metabolic_model$R_day > 0 & Metabolic_model$Date < as.Date("2016-08-10")]/Metabolic_model$Chl_lab[Metabolic_model$R_day > 0 & Metabolic_model$Date < as.Date("2016-08-10")], 	
       Metabolic_model$R_day[Metabolic_model$R_day > 0 & Metabolic_model$Date > as.Date("2016-08-10")]/Metabolic_model$Chl_lab[Metabolic_model$R_day > 0 & Metabolic_model$Date > as.Date("2016-08-10")])	
t.test(Metabolic_model$NEP_day[Metabolic_model$GPP_day < 0 & Metabolic_model$Date < as.Date("2016-08-10")]/Metabolic_model$Chl_lab[Metabolic_model$GPP_day < 0 & Metabolic_model$Date < as.Date("2016-08-10")], 	
       Metabolic_model$NEP_day[Metabolic_model$GPP_day < 0 & Metabolic_model$Date > as.Date("2016-08-10")]/Metabolic_model$Chl_lab[Metabolic_model$GPP_day < 0 & Metabolic_model$Date > as.Date("2016-08-10")])	



t.test(Metabolic_model_DO$GPP_day[Metabolic_model_DO$GPP_day < 0 & Metabolic_model_DO$Date < as.Date("2016-08-10")]/1000, Metabolic_model_DO$GPP_day[Metabolic_model_DO$GPP_day < 0 & Metabolic_model_DO$Date > as.Date("2016-08-10")]/1000)	
t.test(Metabolic_model_DO$R_day[Metabolic_model_DO$R_day > 0 & Metabolic_model_DO$Date < as.Date("2016-08-10")]/1000, Metabolic_model_DO$R_day[Metabolic_model_DO$R_day > 0 & Metabolic_model_DO$Date > as.Date("2016-08-10")]/1000)	
t.test(Metabolic_model_DO$NEP_day[ Metabolic_model_DO$Date < as.Date("2016-08-10")]/1000, Metabolic_model_DO$NEP_day[Metabolic_model_DO$Date > as.Date("2016-08-10")]/1000)	


t.test(Metabolic_model_DO$GPP_day[Metabolic_model_DO$GPP_day < 0 & Metabolic_model_DO$Date < as.Date("2016-08-10")]/Metabolic_model_DO$Chl_lab[Metabolic_model_DO$GPP_day < 0 & Metabolic_model_DO$Date < as.Date("2016-08-10")], 	
       Metabolic_model_DO$GPP_day[Metabolic_model_DO$GPP_day < 0 & Metabolic_model_DO$Date > as.Date("2016-08-10")]/Metabolic_model_DO$Chl_lab[Metabolic_model_DO$GPP_day < 0 & Metabolic_model_DO$Date > as.Date("2016-08-10")])	
t.test(Metabolic_model_DO$R_day[Metabolic_model_DO$R_day > 0 & Metabolic_model_DO$Date < as.Date("2016-08-10")]/Metabolic_model_DO$Chl_lab[Metabolic_model_DO$R_day > 0 & Metabolic_model_DO$Date < as.Date("2016-08-10")], 	
       Metabolic_model_DO$R_day[Metabolic_model_DO$R_day > 0 & Metabolic_model_DO$Date > as.Date("2016-08-10")]/Metabolic_model_DO$Chl_lab[Metabolic_model_DO$R_day > 0 & Metabolic_model_DO$Date > as.Date("2016-08-10")])	
t.test(Metabolic_model_DO$NEP_day[Metabolic_model_DO$GPP_day < 0 & Metabolic_model_DO$Date < as.Date("2016-08-10")]/Metabolic_model_DO$Chl_lab[Metabolic_model_DO$GPP_day < 0 & Metabolic_model_DO$Date < as.Date("2016-08-10")], 	
       Metabolic_model_DO$NEP_day[Metabolic_model_DO$GPP_day < 0 & Metabolic_model_DO$Date > as.Date("2016-08-10")]/Metabolic_model_DO$Chl_lab[Metabolic_model_DO$GPP_day < 0 & Metabolic_model_DO$Date > as.Date("2016-08-10")])	




# check the difference between methods	
t.test(Metabolic_model_DO$NEP_day, Metabolic_ECS$NEP_day)	
t.test(Metabolic_model_DO$NEP_day, Metabolic_model$NEP_day)	
t.test(Metabolic_model$NEP_day, Metabolic_ECS$NEP_day)	



#### Check the regression with Chl and the difference on Chl over time	

{	
  Metabolic_ECS$diff_chl[2:length(Metabolic_ECS$Chl_lab)] <- diff(Metabolic_ECS$Chl_lab)	
  Metabolic_ECS$change_chl[2:length(Metabolic_ECS$Chl_lab)] <- diff(log(Metabolic_ECS$Chl_lab))	
  Metabolic_model$diff_chl[2:length(Metabolic_model$Chl_lab)] <- diff(Metabolic_model$Chl_lab)	
  Metabolic_model$change_chl[2:length(Metabolic_model$Chl_lab)] <- diff(log(Metabolic_model$Chl_lab))	
  Metabolic_model_DO$diff_chl[2:length(Metabolic_model_DO$Chl_lab)] <- diff(Metabolic_model_DO$Chl_lab)	
  Metabolic_model_DO$change_chl[2:length(Metabolic_model_DO$Chl_lab)] <- diff(log(Metabolic_model_DO$Chl_lab))	
  
  library(vegan)	
  
  rda <- Metabolic_ECS[complete.cases(Metabolic_ECS[,c("change_chl", "Chl_lab", "NEP_day", "GPP_day", "R_day", "diff_chl", "T_C", "I_60")]),]	
  plot(rda(rda[,c(8:10)]~ rda$change_chl * rda$Chl_lab * rda$diff_chl))	
  plot(rda(rda[,c(8,9)]~ rda$change_chl + rda$Chl_lab + rda$diff_chl + rda$T_C + rda$I))	
  
  
  
  # NEP 	
  plot(Metabolic_ECS$Chl_lab, Metabolic_ECS$NEP_day)	
  par(mfrow = c(2,2))	
  plot(lm(Metabolic_ECS$NEP_day[-c(49,126,9)] ~ Metabolic_ECS$change_chl[-c(49,126,9)]))	
  summary(lm(Metabolic_ECS$NEP_day[-c(49,126,9)] ~ Metabolic_ECS$change_chl[-c(49,126,9)]))	
  
  plot(lm(Metabolic_model$NEP_day ~ Metabolic_model$Chl_lab))	
  summary(lm(Metabolic_model$NEP_day ~ Metabolic_model$Chl_lab))	
  
  Metabolic_ECS$NEP_dayTrans <- sign(Metabolic_ECS$NEP_day) * log(abs(Metabolic_ECS$NEP_day))	
  
  
  # Change forms a better model 	
  # lm_NEP_Chl <- lm(Metabolic$NEP_dayTrans ~ Metabolic$diff_chl * Metabolic$Chl_lab)	
  
  # There is no interaction The AIC isn't better but the Adjusted R2 gets lower	
  lm_NEP_Chl <- lm(Metabolic_ECS$NEP_day[-c(49,126,9)] ~ Metabolic_ECS$change_chl[-c(49,126,9)] + Metabolic_ECS$Chl_lab[-c(49,126,9)])	
  
  
  # The model is very much skewed. Can't get it better as it is	
  plot(lm_NEP_Chl)	
  # check for heteroscedasticity	
  lmtest::bptest(lm_NEP_Chl)	
  # No heteroscedasticity	
  
  shapiro.test(residuals(lm_NEP_Chl))	
  
  
  vif(lm(log(Metabolic_ECS$NEP_day*-1) ~ Metabolic_ECS$diff_chl + Metabolic_ECS$Chl_lab))	
  summary(lm_NEP_Chl)	
  AIC(lm_NEP_Chl)	
  
  
  
  
  # GPP	
  
  
  
  # measured	
  lm_GPP_Chl_ECS <- lm(log(GPP_day*-1) ~ change_chl + Chl_lab, data = Metabolic_ECS)	
  
  par(mfrow = c(2,2))	
  plot(lm_GPP_Chl_ECS)	
  # check for heteroscedasticity	
  lmtest::bptest(lm_GPP_Chl_ECS)	
  # No heteroscedasticity	
  shapiro.test(residuals(lm_GPP_Chl_ECS))	
  hist(residuals(lm_GPP_Chl_ECS))	
  vif(lm(log(Metabolic_ECS$GPP_day*-1) ~ Metabolic_ECS$diff_chl + Metabolic_ECS$Chl_lab))	
  AIC(lm_GPP_Chl_ECS)	
  
  
  summary(lm_GPP_Chl_ECS)	
  
  
  boot(data = Metabolic_ECS, statistic = bootReg, formula = log(GPP_day*-1) ~ change_chl + Chl_lab, R = 5000)	
  
  # modeled	
  
  
  
  lm_GPP_Chl_model <- lm(log(Metabolic_model$GPP_day*-1) ~ Metabolic_model$change_chl + Metabolic_model$Chl_lab)	
  lm_GPP_Chl_model <- lm(log(Metabolic_model$GPP_day[-c(137)]*-1) ~ Metabolic_model$change_chl[-c(137)] + Metabolic_model$Chl_lab[-c(137)])	
  
  
  
  
  par(mfrow = c(2,2))	
  plot(lm_GPP_Chl_model)	
  # check for heteroscedasticity	
  lmtest::bptest(lm_GPP_Chl_model)	
  # No heteroscedasticity	
  shapiro.test(residuals(lm_GPP_Chl_model))	
  
  vif(lm(log(Metabolic_model$GPP_day*-1) ~ Metabolic_model$diff_chl + Metabolic_model$Chl_lab))	
  AIC(lm_GPP_Chl_model)	
  
  
  summary(lm_GPP_Chl_model)	
  
  
  # modeled O2	
  
  lm_GPP_Chl_model_DO <- lm(log(Metabolic_model_DO$GPP_day[-c(137)]*-1) ~ Metabolic_model_DO$change_chl[-c(137)] + Metabolic_model_DO$Chl_lab[-c(137)])	
  
  par(mfrow = c(2,2))	
  plot(lm_GPP_Chl_model_DO)	
  # check for heteroscedasticity	
  lmtest::bptest(lm_GPP_Chl_model_DO)	
  # No heteroscedasticity	
  shapiro.test(residuals(lm_GPP_Chl_model_DO))	
  
  vif(lm(log(Metabolic_model$GPP_day*-1) ~ Metabolic_model$diff_chl + Metabolic_model$Chl_lab))	
  AIC(lm_GPP_Chl_model_DO)	
  
  
  summary(lm_GPP_Chl_model_DO)	
  
  
  
  # R	
  
  #measured flux	
  summary(lm(log(Metabolic_ECS$R_day) ~ Metabolic_ECS$change_chl))	
  summary(lm(log(Metabolic_ECS$R_day) ~ Metabolic_ECS$Chl_lab))	
  summary(lm(log(Metabolic_ECS$R_day) ~ Metabolic_ECS$change_chl * Metabolic_ECS$Chl_lab))	
  
  
  lm_R_Chl_ECS <- lm(log(Metabolic_ECS$R_day) ~ Metabolic_ECS$change_chl + Metabolic_ECS$Chl_lab)	
  
  par(mfrow = c(2,2))	
  plot(lm_R_Chl_ECS)	
  # check for heteroscedasticity	
  lmtest::bptest(lm_R_Chl_ECS)	
  # No heteroscedasticity	
  
  shapiro.test(residuals(lm_R_Chl_ECS))	
  hist(residuals(lm_R_Chl_ECS))	
  
  AIC(lm_R_Chl_ECS)	
  
  vif(lm(log(Metabolic_ECS$GPP_day*-1) ~ Metabolic_ECS$diff_chl + Metabolic_ECS$Chl_lab))	
  
  
  
  summary(lm_R_Chl_ECS)	
  
  
  
  # modeled flux	
  
  summary(lm(log(Metabolic_model$R_day) ~ Metabolic_model$change_chl))	
  summary(lm(log(Metabolic_model$R_day) ~ Metabolic_model$Chl_lab))	
  summary(lm(log(Metabolic_model$R_day) ~ Metabolic_model$change_chl * Metabolic_model$Chl_lab))	
  
  
  lm_R_Chl_model <- lm(log(Metabolic_model$R_day[-c(9,142)]) ~ Metabolic_model$change_chl[-c(9,142)] + Metabolic_model$Chl_lab[-c(9,142)])	
  
  
  par(mfrow = c(2,2))	
  plot(lm_R_Chl_model)	
  # check for heteroscedasticity	
  lmtest::bptest(lm_R_Chl_model)	
  # No heteroscedasticity	
  
  shapiro.test(residuals(lm_R_Chl_model))	
  hist(residuals(lm_R_Chl_model))	
  
  AIC(lm_R_Chl_model)	
  
  vif(lm(log(Metabolic_model$GPP_day*-1) ~ Metabolic_model$diff_chl + Metabolic_model$Chl_lab))	
  
  
  
  summary(lm_R_Chl_model)	
  
  
  # chech for DO model	
  
  lm_R_Chl_DO_model <- lm(log(Metabolic_model_DO$R_day[-26]) ~ Metabolic_model_DO$change_chl[-26] + Metabolic_model_DO$Chl_lab[-26])	
  
  
  par(mfrow = c(2,2))	
  plot(lm_R_Chl_DO_model)	
  # check for heteroscedasticity	
  lmtest::bptest(lm_R_Chl_DO_model)	
  # No heteroscedasticity	
  
  shapiro.test(residuals(lm_R_Chl_DO_model))	
  hist(residuals(lm_R_Chl_DO_model))	
  
  AIC(lm_R_Chl_DO_model)	
  
  vif(lm(log(Metabolic_model$GPP_day*-1) ~ Metabolic_model$diff_chl + Metabolic_model$Chl_lab))	
  
  
  
  summary(lm_R_Chl_DO_model)	
  
  
  
  
  # Check when R exceeds GPP	
  
  
  # For the modeled CO2 flux	
  
  # assuming  average growth	
  P <- function(T)  (-exp(4.0373757)* exp(5.6280641*0.04) * exp(0.0059871*T) + exp(4.1713080) *exp(3.4333259*0.04) * exp(0.0048901*T)) 	
  TP <-function(T,P0) P(T)-P0 	
  plot(-400:400,P(-400:400)) 	
  abline(h=0)	
  uniroot(TP,c(-400,400),P=0) 	
  
  
  #s.e.	
  
  coef(summary(lm_R_Chl_ECS))[,2]	
  coef(summary(lm_GPP_Chl_ECS))[,2]	
  
  
  P <- function(T)  (-exp(4.0373757-0.1282699) * exp((5.6280641-0.8644973)*0.04) * exp((0.0059871-0.0005672)*T) + exp(4.1713080-0.1327674) * exp((3.4333259 - 0.7883963)*0.04) * exp((0.0048901-0.0005901)*T)) 	
  TP <-function(T,P0) P(T)-P0 	
  plot(1:400,P(1:400)) 	
  uniroot(TP,c(0,400),P=0)	
  
  # s.e. = 	
  42.06269 - 39.90473	
  
  
  
  # Check when the NEP becomes positve at pop decline	
  
  P <- function(T)  (-exp(4.0373757)* exp(5.6280641*-0.03) * exp(0.0059871*T) + exp(4.1713080) *exp(3.4333259*-0.03) * exp(0.0048901*T)) 	
  
  TP <-function(T,P0) P(T)-P0 	
  plot(1:400,P(1:400)) 	
  uniroot(TP,c(0,400),P=0) 	
  abline(h=0, col = "red")	
  
  #s.e.	
  
  coef(summary(lm_R_Chl))[,2]	
  coef(summary(lm_GPP_Chl))[,2]	
  
  P <- function(T)  (-exp(4.0373757-0.1282699) * exp((5.6280641-0.8644973)*-0.03) * exp((0.0059871-0.0005672)*T) + exp(4.1713080-0.1327674) * exp((3.4333259 - 0.7883963)*-0.03) * exp((0.0048901-0.0005901)*T)) 	
  TP <-function(T,P0) P(T)-P0 	
  plot(1:400,P(1:400)) 	
  uniroot(TP,c(0,400),P=0)	
  
  
  # s.e. 	
  182.1098 - 172.3314	
  
  
  # based on the regression models NEP would be 	
  
  
  #for the modeled O2 flux	
  
  # assuming average growth	
  P <- function(T)  (exp(4.3506687)* exp(4.4629046*0.04) * exp(0.0041343*T) - exp(4.4344539) *exp(1.6260009 * 0.04) * exp(0.0025672*T)) 	
  TP <-function(T,P0) P(T)-P0 	
  plot(1:400,P(1:400)) 	
  uniroot(TP,c(0,400),P=0) 	
  
  
  #s.e.	
  
  coef(summary(lm_R_Chl_model))[,2]	
  coef(summary(lm_GPP_Chl_model))[,2]	
  
  P <- function(T)  (-exp(3.88+0.1522318887)* exp((11.31+1.3898583157)*0.04) * exp((0.004995 + 0.0006534179) *T) + exp(4.0009524 + 0.149710852) *exp((8.8796114 + 1.354280386)* 0.04) * exp((0.0036338 + 0.000640467)*T)) 	
  TP <-function(T,P0) P(T)-P0 	
  plot(1:400,P(1:400)) 	
  uniroot(TP,c(0,400),P=0)	
  
  # s.e. = 	
  17.43818 - 14.40358	
  
  
  
  # Check when the NEP becomes negative at decline rate	
  
  P <- function(T)  (-exp(3.88)* exp(11.31 * -0.03) * exp(0.004995*T) + exp(4.0009524) *exp(8.8796114 * -0.03) * exp(0.0036338*T)) 	
  
  TP <-function(T,P0) P(T)-P0 	
  plot(1:400,P(1:400)) 	
  uniroot(TP,c(0,400),P=0) 	
  abline(h=0, col = "red")	
  
  #s.e.	
  
  coef(summary(lm_R_Chl))[,2]	
  coef(summary(lm_GPP_Chl))[,2]	
  
  P <- function(T)  (-exp(3.88+0.1522318887)* exp((11.31+1.3898583157)*-0.03) * exp((0.004995 + 0.0006534179) *T) + exp(4.0009524 + 0.149710852) *exp((8.8796114 + 1.354280386)* -0.03) * exp((0.0036338 + 0.000640467)*T)) 	
  TP <-function(T,P0) P(T)-P0 	
  plot(1:400,P(1:400)) 	
  uniroot(TP,c(0,400),P=0)	
  
  
  # s.e. 	
  142.4214 - 140.0213	
  
  
  # based on the regression models NEP would be 	
  
  
  
  
  
}	






{	
  
  
  # Check the model in combination with light and temperature	
  
  # GPP	
  
  
  
  plot(rda(rda[c(10:12)] ~ rda$Sun_hours + rda$Chl_lab + rda$T_C + rda$I + rda$change_chl))	
  plot(rda(rda[c(10:12)] ~ rda$Sun_hours * rda$Chl_lab * rda$T_C * rda$I * rda$change_chl))	
  
  PAR_dailyDose$Date <- as.Date(PAR_dailyDose$Date)	
  Metabolic_ECS <- join(Metabolic_ECS[c(1:5,7,8,10:14)], PAR_dailyDose, by = "Date")	
  Metabolic_model <- join(Metabolic_model[c(1:5,7,8,10:14)], PAR_dailyDose, by = "Date")	
  Metabolic_model_DO <- join(Metabolic_model_DO[c(1:5,7,8,10:14)], PAR_dailyDose, by = "Date")	
  
  
  
  
  lm_GPP_Chl_ECS <- lm(log(GPP_day*-1) ~ change_chl + Chl_lab * Sun_hours + I60  * Sun_hours + Chl_lab * I60 , data = Metabolic_ECS[-c(5),])	
  
  vif(lm(Metabolic$GPP_day~ Metabolic$Sun_hours + Metabolic$T_C + Metabolic$I))	
  
  par(mfrow = c(2,2))	
  plot(lm_GPP_Chl_ECS)	
  # check for heteroscedasticity	
  lmtest::bptest(lm_GPP_Chl_ECS)	
  # No heteroscedasticity	
  shapiro.test(residuals(lm_GPP_Chl_ECS))	
  
  vif(lm(log(Metabolic$GPP_day*-1) ~ Metabolic$diff_chl + Metabolic$Chl_lab))	
  AIC(lm_GPP_Chl_ECS)	
  
  
  summary(lm_GPP_Chl_ECS)	
  
  
  
  
  
  
  # for the modeled flux	
  
  lm_GPP_Chl_model <- lm(log(GPP_day*-1) ~ change_chl + Chl_lab * I60 + Sun_hours, data = Metabolic_model[-c(5),])	
  
  
  
  par(mfrow = c(2,2))	
  plot(lm_GPP_Chl_model)	
  # check for heteroscedasticity	
  lmtest::bptest(lm_GPP_Chl_model)	
  # No heteroscedasticity	
  shapiro.test(residuals(lm_GPP_Chl_model))	
  
  vif(lm(log(lm_GPP_Chl_model$GPP_day*-1) ~ lm_GPP_Chl_model$diff_chl + Metabolic$Chl_lab))	
  AIC(lm_GPP_Chl_model)	
  
  
  summary(lm_GPP_Chl_model)	
  
  
  # for O2	
  
  lm_GPP_Chl_model_DO <- lm(log(GPP_day*-1) ~ change_chl + Chl_lab * Sun_hours + I60, data = Metabolic_model_DO)	
  
  
  
  par(mfrow = c(2,2))	
  plot(lm_GPP_Chl_model_DO)	
  # check for heteroscedasticity	
  lmtest::bptest(lm_GPP_Chl_model_DO)	
  # No heteroscedasticity	
  shapiro.test(residuals(lm_GPP_Chl_model_DO))	
  
  vif(lm(log(lm_GPP_Chl_model$GPP_day*-1) ~ lm_GPP_Chl_model$diff_chl + Metabolic$Chl_lab))	
  AIC(lm_GPP_Chl_model_DO)	
  
  
  summary(lm_GPP_Chl_model_DO)	
  
  
  
  # R	
  
  summary(lm(log(Metabolic_ECS$R_day) ~ Metabolic_ECS$change_chl))	
  summary(lm(log(Metabolic_ECS$R_day) ~ Metabolic_ECS$Chl_lab))	
  summary(lm(log(Metabolic_ECS$R_day) ~ Metabolic_ECS$T_C))	
  summary(lm(log(Metabolic_ECS$R_day) ~ Metabolic_ECS$I_60))	
  summary(lm(log(Metabolic_ECS$R_day) ~ Metabolic_ECS$Sun_hours))	
  
  lm_R_Chl_ECS <- lm(log(GPP_day*-1) ~ change_chl *  Sun_hours +  Chl_lab + T_C, data = Metabolic_ECS[-c(129,130,126,131,5,9),])	
  
  
  
  
  par(mfrow = c(2,2))	
  plot(lm_R_Chl_ECS)	
  par(mfrow = c(1,1))	
  # check for heteroscedasticity	
  lmtest::bptest(lm_R_Chl_ECS)	
  # No heteroscedasticity	
  
  shapiro.test(residuals(lm_R_Chl_ECS))	
  hist(residuals(lm_R_Chl_ECS))	
  
  AIC(lm_R_Chl_ECS)	
  
  vif(lm(log(Metabolic_ECS$GPP_day*-1) ~ Metabolic_ECS$diff_chl + Metabolic_ECS$Chl_lab + Metabolic_ECS$T_C + Metabolic_ECS$Sun_hours))	
  
  
  
  summary(lm_R_Chl_ECS)	
  
  
  
  # for the modeled flux	
  
  summary(lm(log(Metabolic_model$R_day) ~ Metabolic_model$change_chl))	
  summary(lm(log(Metabolic_model$R_day) ~ Metabolic_model$Chl_lab))	
  summary(lm(log(Metabolic_model$R_day) ~ Metabolic_model$T_C))	
  summary(lm(log(Metabolic_model$R_day) ~ Metabolic_model$I_60))	
  summary(lm(log(Metabolic_model$R_day) ~ Metabolic_model$Sun_hours))	
  
  lm_R_Chl_model <- lm(log(GPP_day*-1) ~ change_chl +  Sun_hours +  Chl_lab + T_C, data = Metabolic_model[-c(5,8,140),])	
  
  
  
  
  par(mfrow = c(2,2))	
  plot(lm_R_Chl_model)	
  par(mfrow = c(1,1))	
  # check for heteroscedasticity	
  lmtest::bptest(lm_R_Chl)	
  # No heteroscedasticity	
  
  shapiro.test(residuals(lm_R_Chl_model))	
  hist(residuals(lm_R_Chl_model))	
  
  AIC(lm_R_Chl_model)	
  
  vif(lm(log(Metabolic_model$GPP_day*-1) ~ Metabolic_model$diff_chl + Metabolic_model$Chl_lab + Metabolic_model$I + Metabolic_model$T_C + Metabolic_model$Sun_hours))	
  
  
  
  summary(lm_R_Chl_model)	
  
  
  # For the modeled O2 flux	
  
  lm_R_Chl_modelDO <- lm(log(GPP_day*-1) ~ change_chl +  Sun_hours +  Chl_lab + T_C, data = Metabolic_model_DO[-c(8,111),])	
  
  
  
  
  par(mfrow = c(2,2))	
  plot(lm_R_Chl_modelDO)	
  par(mfrow = c(1,1))	
  # check for heteroscedasticity	
  lmtest::bptest(lm_R_Chl_modelDO)	
  # No heteroscedasticity	
  
  shapiro.test(residuals(lm_R_Chl_modelDO))	
  hist(residuals(lm_R_Chl_modelDO))	
  
  AIC(lm_R_Chl_modelDO)	
  
  vif(lm(log(Metabolic_model$GPP_day*-1) ~ Metabolic_model$diff_chl + Metabolic_model$Chl_lab + Metabolic_model$I + Metabolic_model$T_C + Metabolic_model$Sun_hours))	
  
  
  
  summary(lm_R_Chl_modelDO)	
  
  
  
}	





# Check the DIC metabolic by using the PI curve	
{	
  DIC_PI <- PI_curve(ChangeGas_Chl = DF_Metabolic$Tot_DIC_change_ECS_kal/DF_Metabolic$Chl_lab, Sample = DF_Metabolic$Days, PAR = DF_Metabolic$`I(60)`, Gas = "DIC")	
  DIC_PI <-  data.frame(Reduce(rbind,mapply(c, DIC_PI)))	
  DIC_PI$h_i[DIC_PI$r > 0] = NA	
  DIC_PI$P_max[DIC_PI$r > 0] = NA	
  DIC_PI$r[DIC_PI$r > 0] = NA	
  DIC_PI <- na.omit(DIC_PI)	
  DIC_PI <- DIC_PI*24	
  summary(DIC_PI)	
}	

#Based on the median values and the mean light irradiance at the surface 	
# In all cases the estimates based on the PI curve are higher than based on the diurnal changes in DIC	


#########################################################################################################################	
#	
# NEP DIC_ECS 	
#	
#########################################################################################################################	



DF_Metabolic_PI_ECS <- DF_Metabolic[c("Date","Tot_DIC_change_ECS_kal", "Days", "I(60)", "I (umol/m2/s)", "T_C_water", "ChlVoltNorm",  "Chl_lab")]	

DF_Metabolic_PI_ECS <- DF_Metabolic_PI_ECS[complete.cases(DF_Metabolic_PI_ECS["Tot_DIC_change_ECS_kal"]),]	



# Regression used in report	


{	
  
  
  
  
  PI_field_coef_DIC_L <- PI_curve(ChangeGas_Chl = DF_Metabolic_PI_ECS$Tot_DIC_change_ECS_kal/DF_Metabolic_PI_ECS$Chl_lab, Sample = DF_Metabolic_PI_ECS$Days, PAR = DF_Metabolic_PI_ECS$`I(60)`, Gas = "DIC")	
  
  count(DF_Metabolic_ECS$Date)	
  
  
  # Create a list for the minimum temperature	
  Temp_min <- list()	
  for (i in 1:max(DF_Metabolic_PI_ECS$Days)) {	
    Temp_min[[i]] <- min(DF_Metabolic_PI_ECS$T_C_water[DF_Metabolic_PI_ECS$Days == i])	
  }	
  
  
  Temp_max <- list()	
  for (i in 1:max(DF_Metabolic_PI_ECS$Days)) {	
    Temp_max[[i]] <- max(DF_Metabolic_PI_ECS$T_C_water[DF_Metabolic_PI_ECS$Days == i])	
  }	
  
  
  Temp_av <- list()	
  for (i in 1:max(DF_Metabolic_PI_ECS$Days)) {	
    Temp_av[[i]] <- mean(DF_Metabolic_PI_ECS$T_C_water[DF_Metabolic_PI_ECS$Days == i])	
  }	
  
  Chl_volt <- list()	
  for (i in 1:max(DF_Metabolic_PI_ECS$Days)) {	
    Chl_volt[[i]] <- mean(DF_Metabolic_PI_ECS$ChlVoltNorm[DF_Metabolic_PI_ECS$Days == i])	
  }	
  
  I60_max <- list()	
  for (i in 1:max(DF_Metabolic_PI_ECS$Days)) {	
    I60_max[[i]] <- mean(DF_Metabolic_PI_ECS$`I(60)`[DF_Metabolic_PI_ECS$Days == i])	
  }	
  
  I_max <- list()	
  for (i in 1:max(DF_Metabolic_PI_ECS$Days)) {	
    I_max[[i]] <- max(DF_Metabolic_PI_ECS$`I (umol/m2/s)`[DF_Metabolic_PI_ECS$Days == i])	
  }	
  
  Chl_lab <- list()	
  for (i in 1:max(DF_Metabolic_PI_ECS$Days)) {	
    Chl_lab[[i]] <- mean(DF_Metabolic_PI_ECS$Chl_lab[DF_Metabolic_PI_ECS$Days == i])	
  }	
  
  
  
  PI_field_coef_DIC <-  data.frame(Reduce(rbind,mapply(c, PI_field_coef_DIC_L, Temp_min, Temp_max, Temp_av, I60_max, I_max, Chl_lab, Chl_volt)))	
  
  
  
  names(PI_field_coef_DIC)[1] <- paste("r")	
  names(PI_field_coef_DIC)[2] <- paste("P_max")	
  names(PI_field_coef_DIC)[3] <- paste("h_i")	
  names(PI_field_coef_DIC)[4] <- paste("V4")	
  names(PI_field_coef_DIC)[5] <- paste("V5")	
  names(PI_field_coef_DIC)[6] <- paste("V6")	
  
  PI_field_coef_DIC$h_i[PI_field_coef_DIC$r > 0] = NA	
  PI_field_coef_DIC$P_max[PI_field_coef_DIC$r > 0] = NA	
  PI_field_coef_DIC$r[PI_field_coef_DIC$r > 0] = NA	
  PI_field_coef_DIC_ECS <- na.omit(PI_field_coef_DIC)	
  
  
  
  
  
  
  
  summary(PI_field_coef_DIC_ECS)	
  
  summary(PI_field_coef_DIC_ECS$r*24)	
  sd(PI_field_coef_DIC_ECS$r*24)/sqrt(length(PI_field_coef_DIC_ECS$r*24))	
  
  
  plot(PI_field_coef_DIC_ECS$P_max)	
  plot(PI_field_coef_DIC_ECS$h_i)	
  plot(PI_field_coef_DIC_ECS$r)	
  
  par(mfrow = c(2,2))	
  plot(PI_field_coef_DIC_ECS$h_i[PI_field_coef_DIC_ECS$h_i < 600], xlab = "Time steps (undefined)", ylab = "h_i")	
  #with(PI_field_coef_DIC, text(PI_field_coef_DIC$h_i, labels = row.names(PI_field_coef_DIC) ))	
  plot(PI_field_coef_DIC_ECS$V6[PI_field_coef_DIC_ECS$h_i < 600], PI_field_coef_DIC_ECS$h_i[PI_field_coef_DIC_ECS$h_i < 600], xlab = "Temperature", ylab = "h_i")	
  #with(PI_field_coef_DIC, text(PI_field_coef_DIC$V6, PI_field_coef_DIC$h_i, labels = row.names(PI_field_coef_DIC) ))	
  plot(PI_field_coef_DIC_ECS$X7[PI_field_coef_DIC_ECS$h_i < 2000], PI_field_coef_DIC_ECS$h_i[PI_field_coef_DIC_ECS$h_i < 2000], xlab = "Mean Chlorrophyll fluorescence", ylab = "h_i")	
  #with(PI_field_coef_DIC, text(PI_field_coef_DIC$V7, PI_field_coef_DIC$h_i, labels = row.names(PI_field_coef_DIC) ))	
  plot(PI_field_coef_DIC_ECS$X8[PI_field_coef_DIC_ECS$h_i < 600], PI_field_coef_DIC_ECS$h_i[PI_field_coef_DIC_ECS$h_i < 600], xlab = "max irradiance at 60 cm", ylab = "h_i")	
  #with(PI_field_coef_DIC, text(PI_field_coef_DIC$V8, PI_field_coef_DIC$h_i, labels = row.names(PI_field_coef_DIC) ))	
  par(mfrow = c(1,1))	
  
  
  
  #######################################	
  # P_max	
  #######################################	
  
  # There does not seem to be a clear trend in P_max and temp	
  
  Temp_range <- seq(steps = nrow(PI_field_coef_DIC_ECS),from = min(PI_field_coef_DIC_ECS$V4)-1, to = max(PI_field_coef_DIC_ECS$V5)+1 )	
  
  summary(PI_field_coef_DIC_ECS$P_max)*24	
  sd(PI_field_coef_DIC_ECS$P_max*24)/sqrt(length(PI_field_coef_DIC_ECS$P_max))	
  
  
  
  
  par(mfrow = c(2,2))	
  plot(lm(log(PI_field_coef_DIC_ECS$P_max[-64]) ~ PI_field_coef_DIC_ECS$V6[-64]))	
  par(mfrow = c(1,1))	
  summary(lm(log(PI_field_coef_DIC_ECS$P_max[-64]) ~ PI_field_coef_DIC_ECS$V6[-64]))	
  
  lm_Pmax_ECS <- lm(log(PI_field_coef_DIC_ECS$P_max) ~ PI_field_coef_DIC_ECS$V6)	
  
  Pmax_ECS_boot <- boot(statistic = bootReg, formula = log(P_max) ~ V6, data = PI_field_coef_DIC_ECS, R = 2000)	
  boot.ci(Pmax_ECS_boot, type = "bca", index = 1)	
  boot.ci(Pmax_ECS_boot, type = "bca", index = 2)	
  
  
  
  plot(PI_field_coef_DIC_ECS$V6, PI_field_coef_DIC_ECS$P_max, col = "orange", ylab = "Maximum photosynthesis [umol/L/h]", xlab = "Temperature [C]")	
  lines(Temp_range, exp(-1.370035) * exp(0.005546 * Temp_range), col = "blue")                                                    # Based on I60	
  lines(Temp_range, exp(-2.415) * exp(-0.0908 * Temp_range), col = "blue", lty = 3)  	
  lines(Temp_range, exp(0.728) * exp(0.0622 * Temp_range), col = "blue", lty = 3)	
  
  
  
  
  ####################################################33	
  # REspiration	
  ######################################################	
  
  
  plot(PI_field_coef_DIC_ECS[,4], PI_field_coef_DIC_ECS$r * -1, xlab = "Temperature", ylab = "Respiration", col = "orange")	
  points(PI_field_coef_DIC_ECS[,4], exp(47.9224)* exp(-1.2697/ (0.000086173303 * (PI_field_coef_DIC_ECS[,4] + 273.15))), col = "blue")	
  with(PI_field_coef_DIC_ECS, text(PI_field_coef_DIC_ECS[,4],PI_field_coef_DIC_ECS$r* -1, labels = row.names(PI_field_coef_DIC_ECS) ))	
  
  
  # Calculate the natural logarithm of r and the inverse of temp	
  PI_field_coef_DIC_ECS$r_ln <- log(PI_field_coef_DIC_ECS$r*-1)	
  PI_field_coef_DIC_ECS$invK <- 1/((PI_field_coef_DIC_ECS$V4+273.15)*0.000086173303)	
  
  par(mfrow = c(2,2))	
  plot(lm(PI_field_coef_DIC_ECS$r_ln ~ PI_field_coef_DIC_ECS$invK))	
  par(mfrow = c(1,1))	
  summary(lm(PI_field_coef_DIC_ECS$r_ln ~ PI_field_coef_DIC_ECS$invK))	
  
  
  summary(PI_field_coef_DIC_ECS$r)*24	
  sd(PI_field_coef_DIC_ECS$r*24)/sqrt(length(PI_field_coef_DIC_ECS$r))	
  
  lm_R_ECS <- summary(lm(PI_field_coef_DIC_ECS$r_ln ~ PI_field_coef_DIC_ECS$invK))	
  
  
  r_ECS_boot <- boot(statistic = bootReg, formula = r_ln ~ invK, data = PI_field_coef_DIC_ECS, R = 2000)	
  boot.ci(r_ECS_boot, type = "bca", index = 1)	
  boot.ci(r_ECS_boot, type = "bca", index = 2)	
  
  
  
  plot(PI_field_coef_DIC_ECS$V4, PI_field_coef_DIC_ECS$r*-1, col = "orange", ylab = "Respiration [umol/L/h]", xlab = "Temperature [C]", ylim = c(0,0.5))	
  lines(Temp_range, exp(49.1520)* exp(-1.3030/ (0.000086173303 * (Temp_range + 273.15))), col = "blue")                                                    # Based on I60	
  lines(Temp_range, exp(31.67)* exp(-1.786/ (0.000086173303 * (Temp_range + 273.15))), col = "blue", lty = 3)  	
  lines(Temp_range, exp(68.27) * exp(-0.865/ (0.000086173303 * (Temp_range + 273.15))), col = "blue", lty = 3)	
  
  
  # The intercept is off, compared to 1.96 * the s.e.	
  
  ##################################################	
  # h_i against chl votl	
  ##################################################	
  PI_field_coef_DIC_ECS_Chl <-  data.frame(Reduce(rbind,mapply(c, PI_field_coef_DIC_L, Chl_volt, I60_max, Chl_lab, Chl_volt, I60_max, Chl_lab)))	
  
  names(PI_field_coef_DIC_ECS_Chl)[1] <- paste("r")	
  names(PI_field_coef_DIC_ECS_Chl)[2] <- paste("P_max")	
  names(PI_field_coef_DIC_ECS_Chl)[3] <- paste("h_i")	
  
  
  PI_field_coef_DIC_ECS_Chl$h_i[PI_field_coef_DIC_ECS_Chl$r > 0] = NA	
  PI_field_coef_DIC_ECS_Chl$P_max[PI_field_coef_DIC_ECS_Chl$r > 0] = NA	
  PI_field_coef_DIC_ECS_Chl$r[PI_field_coef_DIC_ECS_Chl$r > 0] = NA	
  PI_field_coef_DIC_ECS_Chl$hi_ln <- log(PI_field_coef_DIC_ECS_Chl$h_i)	
  PI_field_coef_DIC_ECS_Chl <- PI_field_coef_DIC_ECS_Chl[complete.cases(PI_field_coef_DIC_ECS_Chl[,c("r","h_i", "V8", "V9", "hi_ln")]),]	
  
  
  
  
  summary(PI_field_coef_DIC_ECS_Chl$h_i)	
  sd(PI_field_coef_DIC_ECS$h_i)/sqrt(length(PI_field_coef_DIC_ECS$h_i))	
  
  PI_field_coef_DIC_ECS_Chl$hi_ln <- log(PI_field_coef_DIC_ECS_Chl$h_i)	
  
  
  
  PI_field_coef_DIC_ECS_Chl$I <-PI_field_coef_DIC_ECS_Chl$V8                           # The lm model fits better if the inverse Chlorophyll is used. However, at low chlorophyll the h_i becomes unrealistic high	
  
  
  par(mfrow = c(2,2))	
  plot(lm(PI_field_coef_DIC_ECS_Chl$hi_ln ~ PI_field_coef_DIC_ECS_Chl$V9))                	
  par(mfrow = c(1,1))	
  
  
  summary(lm(PI_field_coef_DIC_ECS_Chl$hi_ln ~ PI_field_coef_DIC_ECS_Chl$V9))	
  
  lm_hi_ECS <- lm(PI_field_coef_DIC_ECS_Chl$hi_ln ~ PI_field_coef_DIC_ECS_Chl$V9)	
  
  
  
  
  h_i_ECS_boot <- boot(statistic = bootReg, formula = hi_ln ~ V9, data = PI_field_coef_DIC_ECS_Chl, R = 2000)	
  boot.ci(h_i_ECS_boot, type = "bca", index = 1)	
  boot.ci(h_i_ECS_boot, type = "bca", index = 2)	
  
  ChlV_range <- seq(steps = nrow(PI_field_coef_DIC_ECS_Chl), from = min(PI_field_coef_DIC_ECS_Chl$V9, na.rm = TRUE)-10, to = max(PI_field_coef_DIC_ECS_Chl$V9, na.rm = TRUE))	
  
  plot(PI_field_coef_DIC_ECS_Chl$V9, PI_field_coef_DIC_ECS_Chl$h_i, col = "orange", pch = 16, ylim = c(0,1000),	
       ylab = "Half saturation constant [uE/m2/s]", xlab = "Chlorophyll"~alpha~" ["~mu~"g L"^"-1"*"]", cex.lab = 1.2, cex.axis = 1.2)	
  lines(ChlV_range,  exp(4.958625) * exp(-0.014497 * ChlV_range), col = "red", lwd = 2)                                                    # Based on I60	
  lines(ChlV_range, exp(4.958625 + 0.457345*1.96) * exp((-0.014497 + 0.002979*1.96) * ChlV_range), col = "red", lty = 3)  	
  lines(ChlV_range, exp(4.958625 - 0.457345*1.96) * exp((-0.014497 - 0.002979*1.96) * ChlV_range), col = "red", lty = 3)	
  
  
  # The bootstrap is compareble to the estimated sd by the model, therefore it may be assumed that is a genarel model applicable for other similar data set or models. 	
  
}	






################################################################################################################################################	
#	
#  NEP     DO_MacIntyre	
#	
################################################################################################################################################	

{	
  # correct for Chlorophyll	
  DF_Metabolic$Tot_DO_change_Chl <- (DF_Metabolic$delta.do + DF_Metabolic$`DOflux_umol/L/h` - DF_Metabolic$`CO2flux_benthic_umol/L/h`)/DF_Metabolic$Chl_lab	
  plot(DF_Metabolic$Tot_DO_change_Chl)	
  
  
  DF_Metabolic_PI <- DF_Metabolic[complete.cases(DF_Metabolic["Tot_DO_change_Chl"]),]	
  
  
  PI_field_coef_b <- PI_curve(ChangeGas_Chl = DF_Metabolic_PI$Tot_DO_change_Chl, Sample = DF_Metabolic_PI$Days, PAR = DF_Metabolic_PI$`I(60)`, Gas = "DO")	
  
  
  
  
  PI_field_coef <-  data.frame(Reduce(rbind,mapply(c, PI_field_coef_b, Temp_min, Temp_max, Temp_av, Chl_volt, I60_max, Chl_lab, SIMPLIFY=FALSE)))	
  PI_field_coef$P_max[PI_field_coef$r > 0] = NA	
  PI_field_coef$h_i[PI_field_coef$r > 0] = NA	
  PI_field_coef$r[PI_field_coef$r > 0] = NA	
  
  summary(PI_field_coef$P_max)*24	
  summary(PI_field_coef$h_i)	
  summary(PI_field_coef$r)*24	
  
  sd(PI_field_coef$P_max*24, na.rm = TRUE)/sqrt(length(PI_field_coef$P_max))	
  sd(PI_field_coef$h_i, na.rm = TRUE)/sqrt(length(PI_field_coef$h_i))	
  sd(PI_field_coef$r*24, na.rm = TRUE)/sqrt(length(PI_field_coef$r))	
  
  
  # remove missing values, these cause errors in the inverse tamperature	
  PI_field_coef <- PI_field_coef[complete.cases(PI_field_coef[c("r", "V6")]),]	
  
  
  plot(PI_field_coef[,4], PI_field_coef$r * -1, xlab = "Temperature", ylab = "Respiration", col = "orange")	
  points(PI_field_coef[,4], exp(-2.82402)* exp(-0.00564/ (0.000086173303 * (PI_field_coef[,4] + 273.15))), col = "blue")	
  with(PI_field_coef, text(PI_field_coef[,4],PI_field_coef$r* -1, labels = row.names(PI_field_coef) ))	
  
  
  # Calculate the natural logarithm of r and the inverse of temp	
  PI_field_coef$r_ln <- log(PI_field_coef$r*-1)	
  PI_field_coef$invK <- 1/((PI_field_coef[,4]+273.15)*0.000086173303)	
  
  
  
  
  
  # Extract the parameters with a lm	
  par(mfrow = c(2,2))	
  plot(lm(PI_field_coef$r_ln ~ PI_field_coef$invK))	
  par(mfrow = c(1,1))	
  
  summary(lm(PI_field_coef$r_ln ~ PI_field_coef$invK))                            # Non-significant model	
  
  r_DO_boot <- boot(statistic = bootReg, formula = r_ln ~ invK, data = PI_field_coef, R = 2000)	
  boot.ci(r_DO_boot, type = "bca", index = 1)	
  boot.ci(r_DO_boot, type = "bca", index = 2)	
  
  
  
  plot(PI_field_coef_DIC_ECS$V4, PI_field_coef_DIC_ECS$r*-1, col = "orange", ylab = "Respiration [umol/L/h]", xlab = "Temperature [C]", ylim = c(0,0.5))	
  lines(Temp_range, exp(-2.893979)* exp(-0.003888/ (0.000086173303 * (Temp_range + 273.15))), col = "blue")                                                    # Based on I60	
  lines(Temp_range, exp(-9.280)* exp(-0.1972/ (0.000086173303 * (Temp_range + 273.15))), col = "blue", lty = 3)  	
  lines(Temp_range, exp(4.848) * exp(0.1560/ (0.000086173303 * (Temp_range + 273.15))), col = "blue", lty = 3)	
  
  
  ###################################3	
  # P_max	
  ######################################	
  
  # There does not seem to be a clear trend in P_max and temp	
  
  summary(PI_field_coef$P_max)	
  plot(PI_field_coef$V6, PI_field_coef$P_max, col = "orange")	
  points(PI_field_coef$V6, exp(-0.41857) * exp(-0.03473 * PI_field_coef$V6), col= "blue")     # Based on I 60	
  points(PI_field_coef$V6, exp(-1.430003) * exp(0.001399 * PI_field_coef$V6), col= "blue")     # Based on I surface	
  
  
  par(mfrow = c(2,2))	
  plot(lm(log(PI_field_coef$P_max) ~ PI_field_coef$V6))	
  par(mfrow = c(1,1))	
  summary(lm(log(PI_field_coef$P_max) ~ PI_field_coef$V6))	
  
  lm_Pmax_DO <- lm(log(PI_field_coef$P_max) ~ PI_field_coef$V6)	
  
  
  
  Pmax_DO_boot <- boot(statistic = bootReg, formula = log(P_max) ~ V6, data = PI_field_coef, R = 2000)	
  boot.ci(Pmax_DO_boot, type = "bca", index = 1)	
  boot.ci(Pmax_DO_boot, type = "bca", index = 2)	
  
  plot(PI_field_coef$V6, PI_field_coef$P_max, col = "orange", ylab = "Maximum photosynthesis [umol/L/h]", xlab = "Temperature [C]")	
  lines(Temp_range, exp(-0.42634) * exp(-0.03424 * Temp_range), col = "blue")                                                    # Based on I60	
  lines(Temp_range, exp(-1.1199) * exp(-0.0666 * Temp_range), col = "blue", lty = 3)  	
  lines(Temp_range, exp(0.1960) * exp(0.0006 * Temp_range), col = "blue", lty = 3)	
  
  
  
  
  
  ##################################################	
  # h_i against chl votl	
  ##################################################	
  
  summary(PI_field_coef$h_i)	
  
  PI_field_coef$hi_ln <- log(PI_field_coef$h_i)	
  PI_field_coef$I <- PI_field_coef$V8                           # The lm model fits better if the inverse Chlorophyll is used. However, at low chlorophyll the h_i becomes unrealistic high	
  
  
  par(mfrow = c(2,2))	
  plot(lm(PI_field_coef$hi_ln ~ PI_field_coef$V9))                # This models shows non-linearity, however it fits in reality better than the inverse Chl model	
  par(mfrow = c(1,1))	
  
  
  summary(lm(PI_field_coef$hi_ln[PI_field_coef$hi_ln < 7] ~ PI_field_coef$V9[PI_field_coef$hi_ln < 7]))	
  
  lm_hi_DO <- lm(PI_field_coef$hi_ln[PI_field_coef$hi_ln < 7] ~ PI_field_coef$V9[PI_field_coef$hi_ln < 7])	
  
  
  
  
  
  h_i_DO_boot <- boot(statistic = bootReg, formula = hi_ln ~ V9, data = PI_field_coef, R = 2000)	
  boot.ci(h_i_DO_boot, type = "bca", index = 1)	
  boot.ci(h_i_DO_boot, type = "bca", index = 2)	
  
  
  
  par(mfrow =c(2,1))	
  par(mar = c(3,5,2,2))	
  plot(PI_field_coef$V9[PI_field_coef$hi_ln < 7], PI_field_coef$h_i[PI_field_coef$hi_ln < 7],  pch = 16, col = "orange", ylab = expression("Half saturation constant [ "~mu*"E m"^"-2"*" s"^-1*"]"), xlab = "", xaxt = 'n', cex.lab = 1.2, cex.axis = 1.2)	
  lines(light_range, exp(6.031456) * exp(-0.014622 * ChlV_range), col = "red", lwd = 2)                                                    # Based on I60	
  lines(light_range, exp(5.646) * exp(-0.0207 * ChlV_range), col = "red", lty = 3)  	
  lines(light_range, exp(6.851) * exp(-0.0101 * ChlV_range), col = "red", lty = 3)	
  legend(	
    "topright",	
    legend = "A",	
    bty = 'n',	
    cex = 1.4	
  )	
  
  ChlV_range <- seq(steps = nrow(PI_field_coef_DIC_ECS_Chl), from = min(PI_field_coef_DIC_ECS_Chl$V9, na.rm = TRUE)-10, to = max(PI_field_coef_DIC_ECS_Chl$V9, na.rm = TRUE))	
  par(mar = c(5,5,0,2))	
  plot(PI_field_coef_DIC_ECS_Chl$V9, PI_field_coef_DIC_ECS_Chl$h_i, col = "orange", pch = 16, ylim = c(0,1000),	
       ylab = expression("Half saturation constant [ "~mu*"E m"^"-2"*" s"^-1*"]"), xlab = expression("Chlorophyll"~paste(italic(a))*" ["~mu*"g L"^"-1"*"]"), cex.lab = 1.2, cex.axis = 1.2)	
  lines(ChlV_range,  exp(4.958625) * exp(-0.014497 * ChlV_range), col = "red", lwd = 2)                                                    # Based on I60	
  lines(ChlV_range, exp(4.158) * exp(-0.0199 * ChlV_range), col = "red", lty = 3)  	
  lines(ChlV_range, exp(5.827) * exp(-0.0097 * ChlV_range), col = "red", lty = 3)	
  legend(	
    "topright",	
    legend = "B",	
    bty = 'n',	
    cex = 1.4	
  )	
  
}	




###########################################################################################################################3	
#	
# NEP    DIC Mac_Intyre	
#	
############################################################################################################################	
# Not used	
{	
  DF_Metabolic$Tot_DIC_change_Chl <- (DF_Metabolic$delta.dic + DF_Metabolic$`CO2flux_umol/L/h` + DF_Metabolic$`CO2flux_benthic_umol/L/h` + DF_Metabolic$Ca_prec_umolL - DF_Metabolic$`DIC_IN-OUT_umol/h`)/DF_Metabolic$Chl_lab	
  plot(DF_Metabolic$Tot_DIC_change_Chl)	
  
  PI_field_coef_DIC_a <- PI_curve(ChangeGas_Chl = DF_Metabolic$Tot_DIC_change_Chl, Sample = DF_Metabolic$Days, PAR = DF_Metabolic$`I(60)`, Gas = "DIC")	
  
  
  
  
  PI_field_coef_DIC <-  data.frame(Reduce(rbind,mapply(c, PI_field_coef_DIC_a, Temp_min, Temp_max, Temp_av, Chl_volt, I60_max, Chl_lab, SIMPLIFY=FALSE)))	
  PI_field_coef_DIC$P_max[PI_field_coef_DIC$r > 0] = NA	
  PI_field_coef_DIC$h_i[PI_field_coef_DIC$r > 0] = NA	
  PI_field_coef_DIC$r[PI_field_coef_DIC$r > 0] = NA	
  PI_field_coef_DIC <- PI_field_coef_DIC[complete.cases(PI_field_coef_DIC[,c("r", "V6")]),]	
  
  summary(PI_field_coef_DIC$P_max)*24	
  summary(PI_field_coef_DIC$h_i)	
  summary(PI_field_coef_DIC$r)*24	
  
  
  sd(PI_field_coef_DIC$P_max*24)/sqrt(length(PI_field_coef_DIC$P_max))	
  sd(PI_field_coef_DIC$h_i)/sqrt(length(PI_field_coef_DIC$h_i))	
  sd(PI_field_coef_DIC$r*24)/sqrt(length(PI_field_coef_DIC$r))	
  
  
  par(mfrow = c(2,2))	
  plot(PI_field_coef_DIC$h_i[PI_field_coef_DIC$h_i < 600], xlab = "Time steps (undefined)", ylab = "h_i")	
  #with(PI_field_coef_DIC, text(PI_field_coef_DIC$h_i, labels = row.names(PI_field_coef_DIC) ))	
  plot(PI_field_coef_DIC$V6[PI_field_coef_DIC$h_i < 600], PI_field_coef_DIC$h_i[PI_field_coef_DIC$h_i < 600], xlab = "Temperature", ylab = "h_i")	
  #with(PI_field_coef_DIC, text(PI_field_coef_DIC$V6, PI_field_coef_DIC$h_i, labels = row.names(PI_field_coef_DIC) ))	
  plot(PI_field_coef_DIC$V7[PI_field_coef_DIC$h_i < 2000], PI_field_coef_DIC$h_i[PI_field_coef_DIC$h_i < 2000], xlab = "Mean Chlorrophyll fluorescence", ylab = "h_i")	
  #with(PI_field_coef_DIC, text(PI_field_coef_DIC$V7, PI_field_coef_DIC$h_i, labels = row.names(PI_field_coef_DIC) ))	
  plot(PI_field_coef_DIC$V8[PI_field_coef_DIC$h_i < 600], PI_field_coef_DIC$h_i[PI_field_coef_DIC$h_i < 600], xlab = "max irradiance at 60 cm", ylab = "h_i")	
  #with(PI_field_coef_DIC, text(PI_field_coef_DIC$V8, PI_field_coef_DIC$h_i, labels = row.names(PI_field_coef_DIC) ))	
  par(mfrow = c(1,1))	
  
  
  
  
  #######################################	
  #	
  # Respiration	
  #	
  ######################################	
  
  
  
  # Calculate the natural logarithm of r and the inverse of temp	
  PI_field_coef_DIC$r_ln <- log(PI_field_coef_DIC$r*-1)	
  PI_field_coef_DIC$invK <- 1/((PI_field_coef_DIC[,4]+273.15) * 0.000086173303)	
  
  
  # Extract the parameters with a lm	
  par(mfrow = c(2,2))	
  plot(lm(PI_field_coef_DIC$r_ln ~ PI_field_coef_DIC$invK))	
  par(mfrow = c(1,1))	
  
  summary(lm(PI_field_coef_DIC$r_ln ~ PI_field_coef_DIC$invK))	
  # Based on I surface	
  
  
  lm_r_DIC <- lm(PI_field_coef_DIC$r_ln ~ PI_field_coef_DIC$invK)	
  
  
  
  r_DIC_boot <- boot(statistic = bootReg, formula = r_ln ~ invK, data = PI_field_coef_DIC, R = 2000)	
  boot.ci(r_DIC_boot, type = "bca", index = 1)	
  boot.ci(r_DIC_boot, type = "bca", index = 2)	
  
  
  
  plot(PI_field_coef_DIC$V4, PI_field_coef_DIC$r*-1, col = "orange", ylab = "Respiration [umol/L/h]", xlab = "Temperature [C]", ylim = c(0,0.5))	
  lines(Temp_range, exp(1.9884500)* exp(-0.1164263/ (0.000086173303 * (Temp_range + 273.15))), col = "blue")                                                    # Based on I60	
  lines(Temp_range, exp(-6.761)* exp(-0.3495/ (0.000086173303 * (Temp_range + 273.15))), col = "blue", lty = 3)  	
  lines(Temp_range, exp(11.286) * exp(0.1004/ (0.000086173303 * (Temp_range + 273.15))), col = "blue", lty = 3)	
  
  
  
  
  ###################################3	
  #	
  # P_max	
  #	
  ######################################	
  
  # There does not seem to be a clear trend in P_max and temp	
  
  summary(PI_field_coef_DIC$P_max)	
  
  
  par(mfrow = c(2,2))	
  plot(lm(log(PI_field_coef_DIC$P_max[-c(89,73,13)]) ~ PI_field_coef_DIC$V6[-c(89,73,13)]))	
  par(mfrow = c(1,1))	
  summary(lm(log(PI_field_coef_DIC$P_max) ~ PI_field_coef_DIC$V6))	
  
  lm_Pmax_DIC <- lm(log(PI_field_coef_DIC$P_max) ~ PI_field_coef_DIC$V6)	
  
  
  
  Pmax_DIC_boot <- boot(statistic = bootReg, formula = log(P_max) ~ V6, data = PI_field_coef_DIC, R = 2000)	
  boot.ci(Pmax_DIC_boot, type = "bca", index = 1)	
  boot.ci(Pmax_DIC_boot, type = "bca", index = 2)	
  
  plot(PI_field_coef_DIC$V6, PI_field_coef_DIC$P_max, col = "orange", ylab = "Maximum photosynthesis [umol/L/h]", xlab = "Temperature [C]")	
  lines(Temp_range, exp(-1.168634059) * exp(-0.005904952 * Temp_range), col = "blue")                                                    # Based on I60	
  lines(Temp_range, exp(-2.065) * exp(-0.0726 * Temp_range), col = "blue", lty = 3)  	
  lines(Temp_range, exp(0.313) * exp(0.0414 * Temp_range), col = "blue", lty = 3)	
  
  
  
  ##################################################	
  #	
  # h_i against chl votl	
  #	
  ##################################################	
  
  summary(PI_field_coef_DIC$h_i)	
  
  PI_field_coef_DIC$hi_ln <- log(PI_field_coef_DIC$h_i)	
  PI_field_coef_DIC$I <- PI_field_coef_DIC$V8                           # The lm model fits better if the inverse Chlorophyll is used. However, at low chlorophyll the h_i becomes unrealistic high	
  
  
  par(mfrow = c(2,2))	
  plot(lm(PI_field_coef_DIC$hi_ln ~ PI_field_coef_DIC$V9))                # This models shows non-linearity, however it fits in reality better than the inverse Chl model	
  par(mfrow = c(1,1))	
  
  
  summary(lm(PI_field_coef_DIC$hi_ln ~ PI_field_coef_DIC$V9))	
  
  lm_hi_DIC <- lm(PI_field_coef_DIC$hi_ln ~ PI_field_coef_DIC$V9)	
  
  
  h_i_DIC_boot <- boot(statistic = bootReg, formula = hi_ln ~ I, data = PI_field_coef_DIC, R = 2000)	
  boot.ci(h_i_DIC_boot, type = "bca", index = 1)	
  boot.ci(h_i_DIC_boot, type = "bca", index = 2)	
  
  
  plot(PI_field_coef_DIC$V8, PI_field_coef_DIC$h_i, col = "orange", ylab = "Half saturation constant [umol/m2/s]", xlab = "Light intensity at 60 cm [umol/m2/s]", ylim = c(0,1000))	
  lines(light_range, exp(5.436087304) * exp(0.001602153 * light_range), col = "blue")                                                    # Based on I60	
  lines(light_range, exp(4.782) * exp(-0.0042 * light_range), col = "blue", lty = 3)  	
  lines(light_range, exp(6.117) * exp(0.0087 * light_range), col = "blue", lty = 3)	
  
  
  #####	
  # Compare means	
  
  t.test(PI_field_coef_DIC$P_max, PI_field_coef$P_max)	
  t.test(PI_field_coef_DIC$P_max,  PI_field_coef_DIC_ECS$P_max)	
  t.test(PI_field_coef$P_max, PI_field_coef_DIC_ECS$P_max)	
  
  # No difference between Pmax	
  
  
  t.test(PI_field_coef_DIC$h_i, PI_field_coef$h_i)	
  t.test(PI_field_coef_DIC$h_i,  PI_field_coef_DIC_ECS$h_i)	
  t.test(PI_field_coef$h_i, PI_field_coef_DIC_ECS$h_i)	
  
  # large differences between H_i	
  # The estimates for ECS are the lowest	
  
  
  t.test(PI_field_coef_DIC$r, PI_field_coef$r)	
  t.test(PI_field_coef_DIC$r,  PI_field_coef_DIC_ECS$r)	
  t.test(PI_field_coef$r, PI_field_coef_DIC_ECS$r)	
  
  
  mean(PI_field_coef_DIC$r[-62] +(PI_field_coef_DIC$P_max[-62] * PI_field_coef_DIC$V8[-62]) / (PI_field_coef_DIC$h_i[-62] + PI_field_coef_DIC$V8[-62]), na.rm = TRUE)*24	
  mean(PI_field_coef$r +(PI_field_coef$P_max* PI_field_coef$V8) / (PI_field_coef$h_i + PI_field_coef$V8), na.rm = TRUE)*24	
  mean(PI_field_coef_DIC_ECS$r +(PI_field_coef_DIC_ECS$P_max* PI_field_coef_DIC_ECS$V7) / (PI_field_coef_DIC_ECS$h_i + PI_field_coef_DIC_ECS$V7), na.rm = TRUE)*24	
  
  
  t.test(PI_field_coef_DIC$r + (PI_field_coef_DIC$P_max*  PI_field_coef_DIC$V8) / (PI_field_coef_DIC$h_i +  PI_field_coef_DIC$V8), PI_field_coef$r + (PI_field_coef$P_max*  PI_field_coef_DIC$V8) / (PI_field_coef$h_i +  PI_field_coef_DIC$V8))	
  t.test(PI_field_coef_DIC$r + (PI_field_coef_DIC$P_max*  PI_field_coef_DIC$V8) / (PI_field_coef_DIC$h_i +  PI_field_coef_DIC$V8),  PI_field_coef_DIC_ECS$r + (PI_field_coef_DIC_ECS$P_max*  PI_field_coef_DIC$V8) / (PI_field_coef_DIC_ECS$h_i +  PI_field_coef_DIC$V8))	
  t.test(PI_field_coef_DIC_ECS$r + (PI_field_coef_DIC_ECS$P_max*  PI_field_coef_DIC$V8) / (PI_field_coef_DIC_ECS$h_i +  PI_field_coef_DIC$V8), PI_field_coef$r + (PI_field_coef$P_max*  PI_field_coef_DIC$V8) / (PI_field_coef$h_i +  PI_field_coef_DIC$V8))	
  
  
  # The average values of DO in respiration is lower than that of DIC (model and ECS)	
  
  # The estimation for h_i ECS is lower, meaning more efficient, but R is higher. 	
  
  
}	










#########################################################################################################################################################	
#########################################################################################################################################################	
######################################                   Carbon budget             ######################################################################	
#########################################################################################################################################################	
#########################################################################################################################################################	
DF_data_macIntyre$`CO2flux_benthic_umol/L/h` <- ((10 * exp(0.0362 * DF_data_macIntyre$T_C_water + 1.635)/44.0095) * 1000 * -1)/1800/24        # conversion from mg/m2/d to umol/L/h multiplied by -1 because it relfects a negative flux, C into the the lake	



DF_CarbonBudget <- DF_Metabolic[complete.cases(DF_Metabolic[,c("delta.dic", "CO2flux_macIntyreR", "CO2flux_benthic_umol/L/h", "Ca_prec_umolL", "DIC_IN-OUT_umol/h", "CO2FluxECS_filter_umol_L_h")]),]	
DF_CarbonBudget$Tot_DIC_change <- (DF_CarbonBudget$delta.dic + DF_CarbonBudget$CO2flux_macIntyreR/1800 + DF_CarbonBudget$`CO2flux_benthic_umol/L/h` + DF_CarbonBudget$Ca_prec_umolL - DF_CarbonBudget$`DIC_IN-OUT_umol/h`)	

Moving_average(DF_CarbonBudget$Tot_DIC_change, 3)	

# CO2 flux	
# Tot_DIC change -> productivity 	
# benthic flux 	
# Ca precipitation 	
# Incoming water in the lake.	


# Create a datafile with the average values per day and average vlaues per hl bin	









# Check the daily change in carbon budget	




# Extract the ecosystem production 	
library(dplyr)	
DF_CarbonBudget %>%	
  group_by(Date_dusk,irr) %>%	
  summarise( mean = mean(Tot_DIC_change, na.rm=TRUE),	
             N_DIC = length(Tot_DIC_change),	
             sd_DIC   = sd(Tot_DIC_change, na.rm=TRUE)	
  )-> Daily_mean_NEP	
detach("package:dplyr", unload = TRUE)	


# create a file with the mean fluxes 	
library(dplyr)	
DF_CarbonBudget %>%	
  group_by(Date_dusk) %>%	
  summarise(mean_DIC_change = mean(Tot_DIC_change, na.rm=TRUE),	
            N_DIC_change = length(Tot_DIC_change),	
            sd_DIC_change   = sd(Tot_DIC_change, na.rm=TRUE),	
            se_DIC_change  = sd(Tot_DIC_change, na.rm=TRUE),	
            N_fluxAtm   = sum(!is.na(CO2flux_macIntyreR)),	
            mean_fluxAtmt  = mean(CO2flux_macIntyreR/1800, na.rm=TRUE),	
            sd_fluxAtm    = sd(CO2flux_macIntyreR/1800, na.rm=TRUE),	
            se_fluxAtm    = sd_fluxAtm  / sqrt(N_fluxAtm),	
            mean_fluxAtmtECS  = mean(CO2FluxECS_filter_umol_L_h, na.rm=TRUE),	
            N_fluxBent   = sum(!is.na(`CO2flux_benthic_umol/L/h`)),	
            mean_fluxBent  = mean(`CO2flux_benthic_umol/L/h`, na.rm=TRUE),	
            sd_fluxBent    = sd(`CO2flux_benthic_umol/L/h`, na.rm=TRUE),	
            se_fluxBent    = sd_fluxBent / sqrt(N_fluxBent),	
            N_CaPrec    = sum(!is.na(Ca_prec_umolL)),	
            mean_CaPrec  = mean(Ca_prec_umolL, na.rm=TRUE),	
            sd_CaPrec    = sd(Ca_prec_umolL, na.rm=TRUE),	
            se_CaPrec    = sd_CaPrec / sqrt(N_CaPrec ),	
            N_CInOut    = sum(!is.na(`DIC_IN-OUT_umol/h`)),	
            mean_CInOut  = mean(`DIC_IN-OUT_umol/h`, na.rm=TRUE),	
            sd_CInOut    = sd(`DIC_IN-OUT_umol/h`, na.rm=TRUE),	
            se_CInOut    = sd_CInOut,	
            mean_Chl_lab  = mean(Chl_lab, na.rm=TRUE),	
            N_Chl_lab = sum(!is.na(Chl_lab)),	
            sd_Chl_lab = sd(Chl_lab, na.rm = TRUE)	
  )-> Daily_mean_Cchange	
detach("package:dplyr", unload = TRUE)	




plot(Daily_mean_NEP$Date_dusk[Daily_mean_NEP$irr == 1], Daily_mean_NEP$mean[Daily_mean_NEP$irr == 1], col = "blue")	
points(Daily_mean_NEP$Date_dusk[Daily_mean_NEP$irr == 0], Daily_mean_NEP$mean[Daily_mean_NEP$irr == 0], col = "red")	
abline(h = 0)	



# Do the same per month	

library(dplyr)	
DF_CarbonBudget %>%	
  group_by(strftime(Date_dusk,"%m"),irr) %>%	
  summarise( mean = mean(Tot_DIC_change_ECS, na.rm=TRUE),	
             N_DIC = length(Tot_DIC_change_ECS),	
             sd_DIC   = sd(Tot_DIC_change_ECS, na.rm=TRUE)	
  )-> Monthly_NEP_mean	
detach("package:dplyr", unload = TRUE)	


library(dplyr)	
DF_CarbonBudget %>%	
  group_by(strftime(Date_dusk,"%m")) %>%	
  summarise(mean_DIC = mean(Tot_DIC_change_ECS, na.rm=TRUE),	
            N_DIC = length(Tot_DIC_change_ECS),	
            sd_DIC   = sd(Tot_DIC_change_ECS, na.rm=TRUE),	
            N_fluxAtm   = sum(!is.na(CO2flux_macIntyreR)),	
            mean_fluxAtmt  = mean(CO2flux_macIntyreR/1800, na.rm=TRUE),	
            sd_fluxAtm    = sd(CO2flux_macIntyreR/1800, na.rm=TRUE),	
            se_fluxAtm    = sd_fluxAtm  / sqrt(N_fluxAtm),	
            N_fluxBent   = sum(!is.na(`CO2flux_benthic_umol/L/h`)),	
            mean_fluxBent  = mean(`CO2flux_benthic_umol/L/h`, na.rm=TRUE),	
            sd_fluxBent    = sd(`CO2flux_benthic_umol/L/h`, na.rm=TRUE),	
            se_fluxBent    = sd_fluxBent / sqrt(N_fluxBent),	
            N_CaPrec    = sum(!is.na(Ca_prec_umolL)),	
            mean_CaPrec  = mean(Ca_prec_umolL, na.rm=TRUE),	
            sd_CaPrec    = sd(Ca_prec_umolL, na.rm=TRUE),	
            se_CaPrec    = sd_CaPrec / sqrt(N_CaPrec ), 	
            mean_Chl_fluo  = mean(ChlVoltNorm, na.rm=TRUE),	
            N_Chl_fluo = sum(!is.na(ChlVoltNorm)),	
            sd_Chl_fluo = sd(ChlVoltNorm, na.rm = TRUE)	
  )-> Monthly_mean_Cchange	
detach("package:dplyr", unload = TRUE)	


plot(Monthly_NEP_mean$`strftime(Date_dusk, "%m")`[Monthly_NEP_mean$irr == 1], Monthly_NEP_mean$mean[Monthly_NEP_mean$irr == 1], col = "blue")	
points(Monthly_NEP_mean$`strftime(Date_dusk, "%m")`[Monthly_NEP_mean$irr == 0], Monthly_NEP_mean$mean[Monthly_NEP_mean$irr == 0], col = "red")	
abline(h = 0)	


strftime(DF_Metabolic_ECS$Date,"%m")	


write.table(Daily_mean_NEP,"C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/AVP_balans/Carbon_balans/Daily_mean_NEP_umolL", sep = ";")	
write.table(Daily_mean_Cchange,"C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/AVP_balans/Carbon_balans/Daily_mean_Cchange_umolL", sep = ";")	

write.table(Monthly_NEP_mean,"C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/AVP_balans/Carbon_balans/Monthly_NEP_mean", sep = ";")	
write.table(Monthly_mean_Cchange,"C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/AVP_balans/Carbon_balans/Monthly_mean_Cchange", sep = ";")	




#########################	
# Plot the C balance	



# Note that all rates per hour 	
Daily_mean_Cchange <- read_csv("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/AVP_balans/Carbon_balans/Daily_mean_Cchange.csv", 	
                               col_types = cols(Date_dusk = col_date(format = "%m/%d/%Y")))	





Daily_Cchange_long <- Daily_mean_Cchange[c("Date_dusk", "mean_fluxAtmt_neg","mean_fluxAtmt_pos", "mean_fluxBent", "mean_CaPrec_neg", "mean_CaPrec_pos", 	
                                           "mean_CInOut_neg", "mean_CInOut_pos", "mean_Chl_lab", "GPP", "R")]	




Daily_Cchange_long <- reshape(Daily_Cchange_long, 	
                              varying = c("mean_fluxAtmt_neg","mean_fluxAtmt_pos", "mean_fluxBent", "mean_CaPrec_neg", "mean_CaPrec_pos", 	
                                          "mean_CInOut_neg", "mean_CInOut_pos", "GPP", "R"), 	
                              v.names = "Flow_rate_umolh", 	
                              timevar = "post",	
                              times = c("mean_fluxAtmt_neg","mean_fluxAtmt_pos", "mean_fluxBent", "mean_CaPrec_neg", "mean_CaPrec_pos", 	
                                        "mean_CInOut_neg", "mean_CInOut_pos",  "GPP", "R"),	
                              new.row.names = 1:10000,	
                              direction = "long")	



Daily_Cchange_long$Date_dusk <- as.POSIXct(Daily_Cchange_long$Date_dusk)	



library(ggplot2)	

{	
  CarbonBalance <- ggplot() + 	
    geom_area(data = Daily_Cchange_long, aes(as.POSIXct(Date_dusk), Flow_rate_umolh ,group=post, fill = post), position = 'stack')	
  
  
  CarbonBalance <- CarbonBalance + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 	
                                         panel.background = element_rect(fill="White", colour = "grey"), 	
                                         legend.position =c(0.5,0.06), legend.direction = "horizontal",	
                                         legend.background = element_rect(fill = NA),	
                                         legend.key = element_rect(size = 0),	
                                         axis.text=element_text(size=12),	
                                         axis.title=element_text(size=12),	
                                         legend.text = element_text(size=11)) +	
    labs(y = expression("Flow rate [m"^"3"*"d"^"-1"*"]"), x = "Months [2016]") 	
  
  
  CarbonBalance <- CarbonBalance + geom_line(aes(Date_dusk, mean_Chl_lab/10, col = "Chl"), Daily_Cchange_long)	
  
  
  CarbonBalance <- CarbonBalance + scale_fill_manual(values=c( "#009E73", "#56B4E9", "#0072B2",  "#E69F00",  "#D55E00",   "#F0E442","#000000","#56B4E9",   "#F0E442","#000000"),	
                                                     name="",	
                                                     breaks=c( "mean_fluxAtmt_neg","mean_fluxAtmt_pos", "mean_fluxBent", "mean_CaPrec_neg","mean_CaPrec_pos", "mean_CInOut_neg", "mean_CInOut_pos",  "GPP", "R"),	
                                                     labels=c("Atm-flux", "Atm-flux",  "Benth-flux",  "Ca-Prec", "Ca-Prec", "C inflow", "C outflow", "GPP", "R"),	
                                                     guide = guide_legend(ncol = 8, label.hjust = 0,	
                                                                          keywidth=0.7, keyheight=0.5, default.unit="cm"))	
  
  
  CarbonBalance <- CarbonBalance + scale_colour_manual(name="",	
                                                       values=c(Chl="black"),	
                                                       labels = c(expression("Chlorophyll")),	
                                                       guide = guide_legend(ncol = 1))	
  
  
  
  
  CarbonBalance <- CarbonBalance + guides(colour = guide_legend(override.aes = list(shape=32, linetype=1), ncol = 1))	
  
  number_ticks <- function(n) {function(limits) pretty(limits, n)}	
  CarbonBalance <- CarbonBalance + scale_y_continuous(expression("Flow rate [m"^"3"*"d"^"-1"*"]"),	
                                                      sec.axis = sec_axis(~ . *10, name = expression("Chloride [mg L"^"-1"*"]")),	
                                                      breaks=number_ticks(6))	
  
  CarbonBalance <- CarbonBalance + expand_limits(y=c(-50,50))	
  
  CarbonBalance	
}	







CarbonPie = ggplot(data=Daily_Cchange_long, 	
                   aes(x=factor(1),	
                       y=Flow_rate_umolh,	
                       fill = factor(post)	
                   )	
) 	












#####################################	
#	
# create a file with the averages over binned chl values	
#	
#####################################	



ind_a <- seq(length=17, from=0, to=400)	
ind_b <- seq(length=17, from=25, to=425)	
ind_c <- seq(1:17)	

DF_Metabolic_ECS$Chl_lab_bin <- DF_Metabolic_ECS$Chl_lab	
for (i in seq_along(ind_c)){	
  DF_Metabolic_ECS$Chl_lab_bin[DF_Metabolic_ECS$Chl_lab >= ind_a[i] & DF_Metabolic_ECS$Chl_lab < ind_b[i]] = ind_b[i]	
}	


library(dplyr)	
DF_Metabolic_ECS %>%	
  group_by(Chl_lab_bin,irr) %>%	
  summarise( mean = mean(Tot_DIC_change_ECS, na.rm=TRUE),	
             N_DIC = length(Tot_DIC_change_ECS),	
             sd_DIC   = sd(Tot_DIC_change_ECS, na.rm=TRUE)	
  )-> NEP_mean_Chl	
detach("package:dplyr", unload = TRUE)	


library(dplyr)	
DF_Metabolic_ECS %>%	
  group_by(Chl_lab_bin) %>%	
  summarise(mean_DIC = mean(Tot_DIC_change_ECS, na.rm=TRUE),	
            N_DIC = length(Tot_DIC_change_ECS),	
            sd_DIC   = sd(Tot_DIC_change_ECS, na.rm=TRUE),	
            mean_delta_DIC = mean(delta.dic, na.rm=TRUE),	
            N_delta_DIC = length(delta.dic),	
            sd_delta_DIC   = sd(delta.dic, na.rm=TRUE),	
            N_fluxAtm   = sum(!is.na(CO2FluxECS_filter_umol_m2_h)),	
            mean_fluxAtmt  = mean(CO2FluxECS_filter_umol_m2_h/1800, na.rm=TRUE),	
            sd_fluxAtm    = sd(CO2FluxECS_filter_umol_m2_h/1800, na.rm=TRUE),	
            se_fluxAtm    = sd_fluxAtm  / sqrt(N_fluxAtm),	
            N_fluxBent   = sum(!is.na(`CO2flux_benthic_umol/L/h`)),	
            mean_fluxBent  = mean(`CO2flux_benthic_umol/L/h`, na.rm=TRUE),	
            sd_fluxBent    = sd(`CO2flux_benthic_umol/L/h`, na.rm=TRUE),	
            se_fluxBent    = sd_fluxBent / sqrt(N_fluxBent),	
            N_CaPrec    = sum(!is.na(Ca_prec_umolL)),	
            mean_CaPrec  = mean(Ca_prec_umolL, na.rm=TRUE),	
            sd_CaPrec    = sd(Ca_prec_umolL, na.rm=TRUE),	
            se_CaPrec    = sd_CaPrec / sqrt(N_CaPrec ),	
            mean_Chl_fluo  = mean(ChlVoltNorm, na.rm=TRUE),	
            N_Chl_fluo = sum(!is.na(ChlVoltNorm)),	
            sd_Chl_fluo = sd(ChlVoltNorm, na.rm = TRUE)	
  )-> Mean_Cchange_chl	
detach("package:dplyr", unload = TRUE)	





write.table(Mean_Cchange_chl,"C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/AVP_balans/Carbon_balans/Mean_Cchange_chl", sep = ";")	
write.table(NEP_mean_Chl,"C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/AVP_balans/Carbon_balans/NEP_mean_Chl", sep = ";")	




# Checht the relation between incoming and outgoing DIC	

CarbonBudget_InOut <- read_delim("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/AVP_balans/Carbon_balans/CarbonBudget_InOut.csv", 	
                                 ";", escape_double = FALSE, trim_ws = TRUE)	
CarbonBudget_InOut <- CarbonBudget_InOut[,c(1:6)]	

CarbonBudget_InOut <- CarbonBudget_InOut[complete.cases(CarbonBudget_InOut),]	

par(mfrow = c(2,2))	
plot(lm(CarbonBudget_InOut$Incoming_MI~CarbonBudget_InOut$Outgoing_MI))	
plot(lm(CarbonBudget_InOut$Incoming_EC~CarbonBudget_InOut$Outgoing_EC))	
par(mfrow = c(1,1))	

summary(lm(CarbonBudget_InOut$Incoming_MI~CarbonBudget_InOut$Outgoing_MI))	
summary(lm(CarbonBudget_InOut$Incoming_EC~CarbonBudget_InOut$Outgoing_EC))	

