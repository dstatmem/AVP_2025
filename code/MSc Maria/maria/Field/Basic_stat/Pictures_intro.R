
par(mar = mar.2ysides)

inc <- seq(0.9,1.1,0.001666667)
DO <- sin(seq(0.1,12, 0.1))*50+250*inc
pH <- sin(seq(0.1,12,0.1))+8
DIC <- sin(seq(0.1,12, 0.1))*-100+1100/inc
hour <- seq(0.4,48, 0.4)



{
plot(hour,CO2, type = 'l', ylim = c(0,400), ylab = "Gas ["~mu~"mol L"^-1*"]", col = "orange", lwd = 3,
     cex.axis = 1.3, cex.lab = 1.3)
lines(DO, type = 'l', col = "darkblue", lwd = 3)
par(new = TRUE)
plot(pH, type = 'l', yaxt = 'n', ylab = "", ylim = c(6,11), col = "darkred", lwd = 3, lty = 3)
axis(side = 4, ylab = "pH", cex.axis = 1.3)
mtext(side = 4, line = 3, "pH", cex = 1.3)
legend(
  "topleft",
  horiz = TRUE,
  lty = c(1,1,3),
  lwd = 3,
  col = c("orange", "darkblue", "darkred"),
  legend = c(expression("CO"[2]), "DO", "pH"),
  cex = 1.3,
  bty = 'n'
)
}


{
  plot(hour, DO, type = 'l', ylim = c(0,500), xlim = c(-1,48), xlab = "Time [Hours]",ylab = "DO ["~mu~"mol L"^-1*"]", col = "darkblue", lwd = 3,
       cex.axis = 1.3, cex.lab = 1.3)
  rect(xleft = hour[16], 
       xright = hour[47],
       ybottom = -20, ytop = 600, col = rgb(0,0,0, alpha = 0.2), border = NA, mar = -5)
  rect(xleft = hour[78], 
       xright = hour[109],
       ybottom = -20, ytop = 600, col = rgb(0,0,0, alpha = 0.2), border = NA, mar = -5)
  par(new = TRUE)
  plot(hour,DIC, type = 'l', yaxt = 'n', xaxt = 'n', ylab = "", xlab = "", ylim = c(0,1500), col = "darkorange", lwd = 3, lty = 1)
  axis(side = 4, cex.axis = 1.3)
  mtext(side = 4, line = 3, "DIC ["~mu~"mol L"^-1*"]", cex = 1.3)
  legend(
    "topleft",
    horiz = FALSE,
    lty = c(1,1),
    lwd = 3,
    col = c("orange", "darkblue"),
    legend = c("DIC", "DO"),
    cex = 1.3,
    bty = 'n'
  )
  text(x=12.8, y=10, "night", cex = 1.3)
  text(x=37.5, y=10, "night", cex = 1.3)
  #text(x=12.8, y=1450, "R", cex = 1.3)
  #text(x=25.5, y=1450, "GPP", cex = 1.3)
  
}


# Graph to compare atmospheric flux and chlorphyll with own data


{
# Jones et al., 2016

# Jones used the summer flux which lasted 18 weeks
# The flux input is g C per summer 
# used the flux of COle at all might underestimate the flux rate


Jones <- data.frame(Chl = c(3.69047619	,
                            8.452380952	,
                            16.07142857	,
                            23.0952381	,
                            28.69047619	,
                            38.45238095	,
                            59.52380952	,
                            26.78571429	,
                            20.5952381	,
                            17.14285714	,
                            8.214285714	,
                            6.071428571	,
                            3.452380952	,
                            1.666666667	,
                            3.928571429), 
                    FLux = c(2.897727273	,
                             0.909090909	,
                             1.079545455	,
                             -1.079545455	,
                             -4.375	,
                             -6.647727273	,
                             -17.55681818	,
                             -8.579545455	,
                             -7.045454545	,
                             -5.738636364	,
                             -3.977272727	,
                             -3.75	,
                             -1.306818182	,
                             -4.090909091	,
                             -7.102272727))

Jones$FLuxmmolC_d <- Jones$FLux/18/7/44.01*1000

plot(Jones$Chl, Jones$FLuxmmolC_d)


# Finlay et al., 2010

# lakes with a high pH
# The inpt is given in g C per m2 for 120d

Finlay <- data.frame( Chl = c(4.8,
                              29.3,
                              12.6,
                              35.1,
                              22.9,
                              27.6),
                      Flux = c(1.838811505264772,
                               7.330963996379751,
                               -4.432573232378985,
                               -44.77163859517452,
                               -12.681839252996625,
                               -6.5119975276483615))


Finlay$FLuxmmolC_d <- Finlay$Flux/120/44.01*1000

plot(Finlay$Chl, Finlay$FLuxmmolC_d)


# Wilkinson et al., 2016

Wilkinson <- data.frame( Chl = c(2.1,
                              5.8,
                              11),
                      Flux = c(794.20,
                               196.57,
                               390.50))


Wilkinson$FLuxmmolC_d <- Wilkinson$Flux/100

plot(Wilkinson$Chl, Wilkinson$FLuxmmolC_d)



# Jonsson et al., 2008

Jonsson <- data.frame( Chl = 1.9,
                         Flux = c(221))


Jonsson$FLuxmmolC_d <- Jonsson$Flux/44.01

#PAcheco 2013

Pacheco <- data.frame( Chl = c(71,
                               44,
                               41,
                               138,
                               49),
                         Flux = c(10,
                                  -6,
                                  -12,
                                  -18,
                                  -6))


Pacheco$FLuxmmolC_d <- Pacheco$Flux/365/44.01*1000

plot(Pacheco$Chl, Pacheco$FLuxmmolC_d)




# Ouyang et al., 2017
# chl is in ln and flux in umol/m2/s

Ouyang <- data.frame( Chl = c(0.62	,
                              1.04	,
                              1.16	,
                              1.23	,
                              1.65	,
                              2.32	,
                              3.02	,
                              3.1	,
                              2.89	,
                              2.8	,
                              2.27	,
                              2.18	,
                              2.04	,
                              1.97	,
                              1.97	,
                              2.1	,
                              1.5	,
                              1.27	,
                              1.59	,
                              1.11	,
                              1.2	,
                              0.97),
                      Flux = c(0.81	,
                               1.15	,
                               1.11	,
                               0.59	,
                               0.1	,
                               0.26	,
                               0.15	,
                               -0.17	,
                               -0.67	,
                               -0.52	,
                               -0.42	,
                               -0.34	,
                               -0.3	,
                               -0.47	,
                               -0.99	,
                               -1.12	,
                               -0.35	,
                               0.03	,
                               -0.22	,
                               -0.15	,
                               -0.32	,
                               -0.51))

Ouyang$Chl <- exp(Ouyang$Chl)
Ouyang$FLuxmmolC_d <- Ouyang$Flux*60*60/1000

plot(Ouyang$Chl, Ouyang$FLuxmmolC_d)








# own Data
# the raw data is in umol/m2/h
# The data is the averaged data over a month


This_study <- data.frame( Chl = c(391.57744,    
                                  238.04304,    
                                  127.15327,    
                                  18.50942,     
                                  25.84215,     
                                  115.19038,    
                                  222.47849),
                       Flux = c( -2410.4369,
                                 -3723.4806,
                                 -3608.7804,
                                 -819.4526,
                                 450.6659,
                                 -1283.2122,
                                 -1392.5821))


This_study$FLuxmmolC_d <- This_study$Flux*24/1000

plot(This_study$Chl, This_study$FLuxmmolC_d)

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


par(mar = mar.default +1)
plot(This_study$Chl, This_study$FLuxmmolC_d, col = "blue", pch = 19, ylim = c(-90,18), xlim = c(-10,400),
     xlab = expression("Chlorophyll"~alpha~"["~mu~"g L"^"-1"*"]"), ylab = expression("CO"[2]*"flux [mmol m"^"-2"*"d"^"-1"*"]"),
     cex.lab = 1.4, cex.axis = 1.4, cex = 1.2)
abline(h=0, lty = 3, lwd = 2, col = "red")
points(Ouyang$Chl, Ouyang$FLuxmmolC_d, col = "#CC79A7", pch = 18, cex = 1.2)
points(Jones$Chl, Jones$FLuxmmolC_d, col = "#000000", pch = 2, cex = 1.2)
points(Finlay$Chl, Finlay$FLuxmmolC_d, col = "#E69F00", pch = 4, cex = 1.2)
points(Wilkinson$Chl, Wilkinson$FLuxmmolC_d, col = "#56B4E9", pch = 5, cex = 1.2)
points(Jonsson$Chl, Jonsson$FLuxmmolC_d, col = "#009E73", pch = 6, cex = 1.2)
points(Pacheco$Chl, Pacheco$FLuxmmolC_d, col = "#0072B2", pch = 7, cex = 1.2)
legend(
  "bottomleft",
  inset = c(0,0),
  bty = 'n',
  col = c("blue","#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73","#0072B2"),
  pch = c(19, 18, 2,4,5,6,7),
  cex = 1.1,
  legend = c("this study", "Ouyang (2017)","Jones (2016)", "Finlay (2010)", "Wilkinson (2016)", "Jonsson (2010)", "Pacheco (2013)")
)


par(mar = c(5,5,1,2))
plot(Ouyang$Chl, Ouyang$FLuxmmolC_d, col = "#CC79A7", pch = 18,  ylim = c(-10,10), xlim = c(-10,300),
     xlab = expression("Chlorophyll"~paste(italic(a))~"["*mu*"g L"^"-1"*"]"), ylab = expression("CO"[2]~"flux [mmol m"^"-2"*"d"^"-1"*"]"),
     cex.lab = 1.4, cex.axis = 1.4, cex = 1.2)
abline(h=0, lty = 3, lwd = 2, col = "red")
points(Ouyang$Chl, Ouyang$FLuxmmolC_d, col = "#CC79A7", pch = 18, cex = 1.2)
points(Jones$Chl, Jones$FLuxmmolC_d, col = "#000000", pch = 2, cex = 1.2)
points(Finlay$Chl, Finlay$FLuxmmolC_d, col = "#E69F00", pch = 4, cex = 1.2)
points(Wilkinson$Chl, Wilkinson$FLuxmmolC_d, col = "#56B4E9", pch = 5, cex = 1.2)
points(Jonsson$Chl, Jonsson$FLuxmmolC_d, col = "#009E73", pch = 6, cex = 1.2)
points(Pacheco$Chl, Pacheco$FLuxmmolC_d, col = "#0072B2", pch = 7, cex = 1.2)
legend(
  "bottomright",
  inset = c(0,0),
  bty = 'n',
  col = c("#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73","#0072B2"),
  pch = c(18, 2,4,5,6,7),
  cex = 1.1,
  legend = c("Ouyang (2017)","Jones (2016)", "Finlay (2010)", "Wilkinson (2016)", "Jonsson (2010)", "Pacheco (2013)")
)



}








# Plot the NEP 


#Finlay given in mg C m3 d1

NEP_Finlay <- data.frame( Chl = c(4.8,
                                   29.3,
                                   12.6,
                                   35.1,
                                   22.9,
                                   27.6),
                            GPP = c(17.6,
                                    124.9,
                                    83.7,
                                    276.5,
                                    54.2,
                                    74.7),
                             R = c(85.6,
                                   112.1,
                                   144.7,
                                   229.9,
                                   205,
                                   86.5),
                            NPP = c(-70.3,
                                    13.5,
                                    -59.1,
                                    46.5,
                                    -164,
                                    -10.3))


#convert to mmol/L/d
NEP_Finlay[,c(2,3,4)] <- NEP_Finlay[,c(2,3,4)]/44.01/1000



plot(NEP_Finlay$Chl, NEP_Finlay$NPP)
plot(NEP_Finlay$Chl, NEP_Finlay$NPP/NEP_Finlay$Chl)


# Hanson 2003

#Used the R DO and GPP DO to calculate the NEP

NEP_Hanson <- data.frame( Chl = c(4.4	,
                                  7.1	,
                                  9.4	,
                                  13.2	,
                                  13.8	,
                                  15.3	,
                                  17.9	,
                                  18.6	,
                                  19.9	,
                                  20.2	,
                                  21.3	,
                                  23.2	,
                                  25.2	,
                                  28	,
                                  28.1	,
                                  29.3	,
                                  34	,
                                  34.9	,
                                  39.8	,
                                  41.5	,
                                  42.5	,
                                  55.9	,
                                  78.3	,
                                  105.3	),
                          GPP = c(4.6	,
                                  7.6	,
                                  4.4	,
                                  10.6	,
                                  30.2	,
                                  3.1	,
                                  44.1	,
                                  31.6	,
                                  1.01	,
                                  3.6	,
                                  1	,
                                  25.4	,
                                  1.3	,
                                  8.6	,
                                  22.7	,
                                  39.4	,
                                  22.3	,
                                  75.3	,
                                  10.4	,
                                  48.5	,
                                  13.4	,
                                  4.1	,
                                  17.3	,
                                  2.3	),
                          R = c(3.5	,
                                9.5	,
                                7.2	,
                                12.4	,
                                33	,
                                83.5	,
                                59.3	,
                                50.7	,
                                24.7	,
                                14.1	,
                                20.3	,
                                38.8	,
                                1.3	,
                                5.6	,
                                28.3	,
                                33.1	,
                                93.8	,
                                119.8	,
                                9.8	,
                                48.8	,
                                14.3	,
                                5.8	,
                                14.6	,
                                0.5	),
                          NPP = c(1.1	,
                                  -1.9	,
                                  -2.8	,
                                  -1.8	,
                                  -2.8	,
                                  -80.4	,
                                  -15.2	,
                                  -19.1	,
                                  -23.69	,
                                  -10.5	,
                                  -19.3	,
                                  -13.4	,
                                  0	,
                                  3	,
                                  -5.6	,
                                  6.3	,
                                  -71.5	,
                                  -44.5	,
                                  0.6	,
                                  -0.3	,
                                  -0.9	,
                                  -1.7	,
                                  2.7	,
                                  1.8	 ))


#convert to mmol/L/d
NEP_Hanson[,c(2,3,4)] <- NEP_Hanson[,c(2,3,4)]/1000



plot(NEP_Hanson$Chl, NEP_Hanson$NPP)
plot(NEP_Hanson$Chl, NEP_Hanson$NPP/NEP_Finlay$Chl)


# This data is in mg O2 L d!!!!!

# Solomon 2013


NEP_Solomon <- data.frame( Chl = c(55.7	,
                                   2.3	,
                                   17.7	,
                                   2.6	,
                                   19.2	,
                                   1.8	,
                                   64.5	,
                                   5.3	,
                                   16.8	,
                                   3.8	,
                                   17.6	,
                                   33.5	,
                                   17	,
                                   4.2	,
                                   14.3	,
                                   1.3	,
                                   30.3	,
                                   1.9	,
                                   46	,
                                   1.7	,
                                   15	,
                                   41.2	,
                                   51	,
                                   11.6	),
                          NPP = c(4.45,
                                  -0.05,
                                  -0.97,
                                  -0.03,
                                  -0.41,
                                  -0.24,
                                  3.8,
                                  0.17,
                                  -0.88,
                                  1.95,
                                  0.18,
                                  1.75,
                                  0.14,
                                  -0.68,
                                  -0.43,
                                  0.08,
                                  0.41,
                                  0.68,
                                  0.64,
                                  -0.07,
                                  -0.14,
                                  1.03,
                                  0.01,
                                  -0.15))


#convert to mmol/L/d
NEP_Solomon[,c(2)] <- NEP_Solomon[,c(2)]/32



plot(NEP_Solomon$Chl, NEP_Solomon$NPP)
plot(NEP_Solomon$Chl, NEP_Solomon$NPP/NEP_Finlay$Chl)


NEP_Laas_O2 <- data.frame( Chl = c(45.4,
                                   21.8,
                                   47.8,
                                   44.2,
                                   63.6,
                                   44.4),
                           NEP = c(14.76/1000,
                                   6.05/1000,
                                   9.86/1000,
                                   2.62/1000,
                                   -13.31/1000,
                                   2.55/1000))

plot(NEP_Laas_O2$Chl, NEP_Laas_O2$NEP)
plot(NEP_Laas_O2$Chl, NEP_Laas_O2$NEP/NEP_Finlay$Chl)





# The file below was addapted in excell to include GPP 
Metabolic_ECS_P <- read_delim("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Metabolic_ECS.csv", 
                        ";", escape_double = FALSE, col_types = cols(Date = col_date(format = "%m/%d/%Y")), 
                        trim_ws = TRUE)

Metabolic_modelO_P <- read_delim("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/Metabolic_model_DO.csv", 
                              ";", escape_double = FALSE, col_types = cols(Date = col_date(format = "%d-%m-%Y")), 
                              trim_ws = TRUE)

# Metabolism in umol/L/h based on DIC 

Metabolic_ECS_P <- Metabolic_ECS_P[complete.cases(Metabolic_ECS_P[c("R", "NEP")]),]
Metabolic_modelO_P <- Metabolic_modelO_P[complete.cases(Metabolic_modelO_P[c("R", "NEP")]),]

library(dplyr)
Metabolic_ECS_P %>%
  group_by(months(Date)) %>%
  summarise(
    NEP_day = mean(NEP_day, na.rm = TRUE),
  GPP_day = mean(GPP_day, na.rm = TRUE),
  R_day = mean(R_day, na.rm = TRUE),
  Chl_lab = mean(Chl_lab, na.rm = TRUE)) -> NEP_month


Metabolic_modelO_P %>%
  group_by(months(Date)) %>%
  summarise(
    NEP_day = mean(NEP_day, na.rm = TRUE),
    GPP_day = mean(GPP_day, na.rm = TRUE),
    R_day = mean(R_day, na.rm = TRUE),
    Chl_lab = mean(Chl_lab, na.rm = TRUE)) -> NEP_month_O2

detach("package:dplyr", unload = TRUE)

{
par(mfrow = c(1,1))
par(mar = mar.default+1)
plot(NEP_month$Chl_lab, NEP_month$NEP_day*-1/1000, ylab = expression("Monthly averaged NEP [mmol L"^"-1"*"d"^"-1"*"]"), xlab = expression("Chlorophyll "~paste(italic(a))~"["*mu*"g L"^"-1"*"]"),
     pch = 19, col = "blue", cex = 1.1, ylim = c(-0.15,0.15), xlim = c(0,400), cex.axis = 1.2, cex.lab = 1.2)
abline(h=0, lty = 3, lwd = 2, col = "red")
points(NEP_Hanson$Chl, NEP_Hanson$NPP, col = "darkgreen", pch = 15, cex = 1)
points(NEP_Finlay$Chl, NEP_Finlay$NPP, col = "#E69F00", pch = 17, cex = 1.5)
points(NEP_month$Chl_lab, NEP_month$NEP_day*-1/1000, col = "blue", pch = 19, cex = 1.1)
points(NEP_month_O2$Chl_lab, NEP_month_O2$NEP_day*-1/1000, pch = 1, cex = 1.2, col = "blue")
points(NEP_Solomon$Chl, NEP_Solomon$NPP, col = "darkred", pch = 0, cex = 1.1)
points(NEP_Laas_O2$Chl, NEP_Laas_O2$NEP, col = "black", pch = 6, cex = 1.1)
legend(
  "bottomright",
  inset = c(0,0),
  bty = 'n',
  col = c(NA,"blue", "darkgreen", "#E69F00", NA, "blue", "darkred", "black"),
  pch = c(NA,19, 15,17,NA,1,0,6),
  cex = 1,
  legend = c("[Carbon]","This study", expression("Hanson "~paste(italic("et al"))*". (2003)"), expression("Finlay "~paste(italic("et al"))* ". (2010)"),"\n[Oxygen]","This study", expression("Solomon "~paste(italic("et al"))* ". (2013)"), expression("Laas "~paste(italic("et al"))* ". (2012)"))
)
}
