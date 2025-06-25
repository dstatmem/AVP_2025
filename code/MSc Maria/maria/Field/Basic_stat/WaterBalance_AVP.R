library(readr)




DF_WB <- read_csv("C:/Users/Maria/Dropbox/MSc Maria/Maria/Field/Data_files/Waterbalns_R.csv", 
                         col_types = cols(Date = col_date(format = "%m/%d/%Y"), 
                                          Date_Chloride = col_date(format = "%m/%d/%Y")))
View(DF_WB)


plot(DF_WB$TotIn_m3d, DF_WB$TotOut_m3d)



# Check the blance closure of Cl 
plot(DF_WB$Tot_Cl_measured, DF_WB$Tot_Cl_estimated)


# From the date at which data is availeble from the Hoornsloot
lm_Cl <- lm(DF_WB$Tot_Cl_estimated[DF_WB$Date_Chloride > "2013/01/01"] ~ DF_WB$Tot_Cl_measured[DF_WB$Date_Chloride > "2013/01/01"])
par(mfrow = c(2,2))
plot(lm_Cl)
par(mfrow = c(1,1))
lmtest::bptest(lm_Cl)

library("car")
summary(lm_Cl)
# test for linearity
shapiro.test(residuals(lm_Cl))
leveneTest(Tot_Cl_estimated ~ Tot_Cl_measured, data = DF_WB)

par(mar = mar.default+0.3)
#plot the regression of the choride ballance
plot(DF_WB$Tot_Cl_measured[DF_WB$Date_Chloride > "2013/01/01"], DF_WB$Tot_Cl_estimated[DF_WB$Date_Chloride > "2013/01/01"],
     ylab = expression("Estimated outflow Chloride [mg L"^"-1"*"]"), xlab = expression("Measured inflow Chloride [mg L"^"-1"*"]"), col = "darkblue", pch = 16,
     xlim = c(80,180), ylim = c(80,180))
    
abline(lm_Cl, col = "red", lwd = 2)
abline(coef(summary(lm_Cl))[,1] + coef(summary(lm_Cl))[,2], col = "red", lwd = 1, lty = 3)
abline(coef(summary(lm_Cl))[,1] - coef(summary(lm_Cl))[,2], col = "red", lwd = 1, lty = 3)
legend( 
  "topleft",
  bty = 'n',
  lty = c(NA, 1,3),
  pch = c(16, NA, NA),
  col = c("darkblue", "red", "red"),
  legend = c(expression("Cl [mg L"*"-1"*"]"), "Regression line", "95% confidence interval")
)



# Create a new data frame rearanged in the length by in put output

# uitgaande posten 
DF_WB$Verdamping_m3d <- DF_WB$Verdamping_m3d*-1
DF_WB$Wegzijging_m3d <- DF_WB$Wegzijging_m3d*-1
DF_WB$intrek_m3d <- DF_WB$intrek_m3d*-1
DF_WB$uitPolders_m3d <- (DF_WB$Legmeerdijksloot_m3d + DF_WB$Maalderij_m3d + DF_WB$Ouderkerkerlaan_m3d + DF_WB$Ouderkerkerlaan_m3d)*-1
DF_WB$uitlaat_m3d <- DF_WB$uitlaat_m3d*-1
DF_WB$Uitlaat_Tot <- DF_WB$uitlaat_m3d + DF_WB$uitPolders_m3d

#inkomende posten
DF_WB$Neerslag_m3d
DF_WB$uitspoeling_m3d
DF_WB$inlaat_m3d



DF_WB_long <- DF_WB[c("Date", "Verdamping_m3d", "Wegzijging_m3d", "intrek_m3d", 
                    "Uitlaat_Tot", "Neerslag_m3d", "uitspoeling_m3d", "inlaat_m3d", "Chloride_berekend_mgL", "Chloride_gemeten_mgL")]

DF_WB_long <- reshape(DF_WB_long, 
             varying = c("Verdamping_m3d", "Wegzijging_m3d", "intrek_m3d",
                         "Uitlaat_Tot", "Neerslag_m3d", "uitspoeling_m3d", "inlaat_m3d"), 
             v.names = "Amount_m3d", 
             timevar = "post",
             times = c("Verdamping_m3d", "Wegzijging_m3d", "intrek_m3d", 
                       "Uitlaat_Tot", "Neerslag_m3d", "uitspoeling_m3d", "inlaat_m3d"),
             new.row.names = 1:10000,
             direction = "long")


DF_WB_long <- DF_WB_long[order(DF_WB_long$Date),]
DF_WB_long$Date <- as.POSIXct(DF_WB_long$Date)





#plot the water balans



# Color blind friendly palet
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00" , "#CC79A7")

#Breaks for background rectangles
rects <- data.frame(xstart = as.POSIXct("2016-05-18 01:00:00"), 
                    xend = as.POSIXct("2016-11-17 01:00:00"))


library(ggplot2)

{
Waterbalans <- ggplot() + 
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf), alpha = 0.1) + 
            geom_area(data = DF_WB_long, aes(as.POSIXct(Date), Amount_m3d ,group=post, fill = post), position = 'stack') 



Waterbalans <- Waterbalans + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                   panel.background = element_rect(fill="White", colour = "black"), 
                                   legend.position =c(0.5,0.06), legend.direction = "horizontal",
                                   legend.background = element_rect(fill = NA),
                                   legend.key = element_rect(size = 0),
                                   axis.text=element_text(size=13),
                                   axis.title=element_text(size=14),
                                   legend.text = element_text(size=12)) +
                             labs(y = expression("Flow rate [m"^"3"*"d"^"-1"*"]"), x = "Time [Months]") 
                             

Waterbalans <- Waterbalans + geom_line(aes(Date, Chloride_berekend_mgL*100, col = "EstimatedCl"), DF_WB_long)
Waterbalans <- Waterbalans + geom_point(aes(Date, Chloride_gemeten_mgL*100, col = "MeasuredCl"), DF_WB_long)


Waterbalans <- Waterbalans + scale_fill_manual(values=c( "#009E73", "#D55E00", "#0072B2",  "#E69F00",  "#56B4E9",   "#F0E442","#000000", "#000000", "#000000"),
                                               name="",
                                               breaks=c("Neerslag_m3d", "uitspoeling_m3d", "inlaat_m3d", "Wegzijging_m3d", "Verdamping_m3d", "intrek_m3d", "Uitlaat_Tot", "Chloride_berekend_mgL", "Chloride_gemeten_mgL"),
                                               labels=c("Precipitation", "Run-off", "Inlet",  "Percolation", "Evaporation",  "Infiltration", "Outlet", "Cl", "Cl"),
                                               guide = guide_legend(ncol = 8, label.hjust = 0,
                                                                    keywidth=0.7, keyheight=0.5, default.unit="cm"))


Waterbalans <- Waterbalans + scale_colour_manual(name="",
                                                 values=c(EstimatedCl="black", MeasuredCl="blue"),
                                                 labels = c(expression("Estimated chloride"), expression("Measured chloride")),
                                                 guide = guide_legend(ncol = 2))


Waterbalans <- Waterbalans + scale_shape_manual(name="",
                                                values=c(EstimatedCl=2, MeasuredCl=1),
                                                labels = c(expression("Estimated chloride"), expression("Measured chloride")),
                                                guide = guide_legend(ncol = 2)) 




Waterbalans <- Waterbalans + guides(colour = guide_legend(override.aes = list(shape=c(32,16), linetype=c(1,0)), ncol = 2))

number_ticks <- function(n) {function(limits) pretty(limits, n)}
Waterbalans <- Waterbalans + scale_y_continuous(expression("Flow rate [m"^"3"*"d"^"-1"*"]"),
                                                sec.axis = sec_axis(~ . /100, name = expression("Chloride [mg L"^"-1"*"]")),
                                                breaks=number_ticks(6))
                                  
Waterbalans <- Waterbalans + expand_limits(y=c(-30000,30000))

Waterbalans
}





