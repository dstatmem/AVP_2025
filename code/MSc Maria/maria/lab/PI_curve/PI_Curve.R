library(car)
setwd("D:/UvA/Master/IBED stage/Jolanda/Stat/PI_curve")
PI_curves_R <- read.csv("D:/UvA/Master/IBED stage/Jolanda/Data/Lab_data/PI_curves/PI_curves_R.csv")
Iout_optode <- read.csv("D:/UvA/Master/IBED stage/Jolanda/Data/Lab_data/Iout_optode.csv")
PI_model_R <- read.csv("D:/UvA/Master/IBED stage/Jolanda/Stat/PI_curve/coef_PI_model_R.csv")
View(PI_curves_R)
names(PI_curves_R)


## 
##split up in runs

run <- list()

for (i in 1:max(PI_curves_R$run)) {
  run[[i]] <- PI_curves_R[which(PI_curves_R$run == i), ]
}

## if you want to output different files for the different runs 
#for (i in 1:max(PI_curves_R$run)) {
#  run[[i]] <- PI_curves_R[which(PI_curves_R$run == i), ]
#  fl <- paste0("run", i, ".txt")
#  write.table(run[[i]], file = fl, sep = "\t", row.name=F)
#}

PI_curves_R <- read.csv("D:/UvA/Master/IBED stage/Jolanda/Data/Lab_data/PI_curves/PI_curves_R.csv")

PI_curves_R <- read_csv("C:/Users/admin/Dropbox/MSc Maria/Maria/Lab/PI_curve/PI-Curves_R.csv")

##################
#Seperate values 

I1 <- Iout_optode$Ch1
I2 <- Iout_optode$Ch2
I3 <- Iout_optode$Ch3
I4 <- Iout_optode$Ch4


I <- PI_curves_R[c(1:7),c(4)]
PI1 <- data.frame(PI_curves_R[c(8:nrow(PI_curves_R)),c(1,4,5)]) #removed the first set of Ch1 gave an error
PI2 <- data.frame(PI_curves_R[,c(1,4,6)])
PI3 <- data.frame(PI_curves_R[,c(1,4,7)])
PI4 <- data.frame(PI_curves_R[,c(1,4,8)])

#create a list devided by run

PI1_list <- list()
for (i in 1:max(PI1$run)) {
  PI1_list[[i]] <- PI1[which(PI1$run == i), ]
}


PI2_list <- list()
for (i in 1:max(PI2$run)) {
  PI2_list[[i]] <- PI2[which(PI2$run == i), ]
}


PI3_list <- list()
for (i in 1:max(PI3$run)) {
  PI3_list[[i]] <- PI3[which(PI3$run == i), ]
}


PI4_list <- list()
for (i in 1:max(PI4$run)) {
  PI4_list[[i]] <- PI4[which(PI4$run == i), ]
}

##extract the third column the Ch(n) form the list of matrices.
## by using "[" you can extract colums or rows from a list. lapply(Mylist, "[", 1,) => row, lapply(Mylist, "[", , 2) => column, lapply(MyList,"[", 1, 2) => first row second collumn
## http://scl.io/fCPufvMf#gs.NPVmqrc
Ch1 <- lapply(PI1_list, "[", ,3)
Ch1[[65]] <- NULL
Ch1[[64]] <- NULL
Ch1[[63]] <- NULL
Ch1[[39]] <- NULL
Ch1[[36]] <- NULL
Ch1[[25]] <- NULL
Ch1[[24]] <- NULL
Ch2 <- lapply(PI2_list, "[", ,3)
Ch2[[65]] <- NULL
Ch2[[63]] <- NULL
Ch2[[61]] <- NULL
Ch3 <- lapply(PI3_list, "[", ,3)
Ch3[[65]] <- NULL
Ch3[[63]] <- NULL
Ch3[[62]] <- NULL
Ch3[[59]] <- NULL
Ch3[[19]] <- NULL
Ch4 <- lapply(PI4_list, "[", ,3)
Ch4[[65 ]] <- NULL
Ch4[[64]] <- NULL
Ch4[[63]] <- NULL
Ch4[[62]] <- NULL
Ch4[[11]] <- NULL
Ch4[[8]] <- NULL

## select only for delivered=True
## Y<-Mydf[Mydf$Delivered==TRUE,] 





#######################################
#             PI curves
#######################################


## assign the function 
Michaelis_Menten = function(par, I) {
  par[1] + (par[2]* I) / (par[3] + I)
}



########################
# fit the model

  
  
  
  #Ch1
  # ! be aware ! Run one is exually chamber two, etc
  coef_model_Ch1 <- list()
  coef_model_Ch1 <-  lapply(Ch1, function(x) (coef(nls(x ~ r + (P_max * I1) /(h_i + I1),  start=list(r=-0.1, P_max=0.6,h_i=400)))))  
  #nls(Ch1[[67]] ~ r + (P_max * I1) /(h_i + I1),  start=list(r=-0.1, P_max=0.6,h_i=400))
  # 24, 25, 36, 39, 63, 64, 65
  write.table(coef_model_Ch1,"coef_model_Ch1.txt", sep = "\t", row.names = FALSE)
  
  
  #Ch2
  coef_model_Ch2 <- list()
  coef_model_Ch2 <-  lapply(Ch2, function(x) (coef(nls(x ~ r + (P_max * I2) /(h_i + I2),  start=list(r=-0.1, P_max=0.6,h_i=400)))))
  #nls(Ch2[[67]] ~ r + (P_max * I2) /(h_i + I2),  start=list(r=-0.1, P_max=0.6,h_i=400))
  # 61, 63, 65
  write.table(coef_model_Ch2,"coef_model_Ch2.txt", sep = "\t", row.names = FALSE)
  
  #Ch3
  # ! be aware ! Run 19 are removed but numbering continues as mentioned
  coef_model_Ch3 <- list()
  coef_model_Ch3 <-  lapply(Ch3, function(x) (coef(nls(x ~ r + (P_max * I3) /(h_i + I3),  start=list(r=-0.1, P_max=0.6,h_i=400)))))
  #nls(Ch3[[67]] ~ r + (P_max * I3) /(h_i + I3),  start=list(r=-0.1, P_max=0.6,h_i=400))
  # 19, 59, 62, 63, 65, 
  write.table(coef_model_Ch3,"coef_model_Ch3.txt", sep = "\t", row.names = FALSE)
  
  #Ch4
  # ! be aware ! Run eight and ten are removed but numbering continues as mentioned
  coef_model_Ch4 <- list()
  coef_model_Ch4 <-  lapply(Ch4, function(x) (coef(nls(x ~ (r + (P_max * I4) /(h_i + I4)),  start=list(r=-0.1, P_max=0.6,h_i=400)))))
  #nls(Ch4[[67]] ~ r + (P_max * I4) /(h_i + I4),  start=list(r=-0.1, P_max=0.6,h_i=400))
  # 8, 11, 62, 63, 64, 65  
  write.table(coef_model_Ch4,"coef_model_Ch4.txt", sep = "\t", row.names = FALSE)
  

  
################################
# Plot the model
# and extract the values of the fitted line
  
  
  
# visualization code
# get the coefficients using the coef function

fit_line_Ch1 <- lapply(coef_model_Ch1, function(x) Michaelis_Menten(x,I1))
lapply(fit_line_Ch1, function(x) plot(x ~ I1,type="l"))
plot(fit_line_Ch1[[3]]~I1)
write.table(fit_line_Ch1,"
            ", sep = "\t", row.names = FALSE)
#lapply(Ch1, function(x) plot(x ~ I))

fit_line_Ch2 <- lapply(coef_model_Ch2, function(x) Michaelis_Menten(x,I2))
lapply(fit_line_Ch2, function(x) plot(x ~ I2,type="l"))
write.table(fit_line_Ch2,"fit_line_Ch2.txt", sep = "\t", row.names = FALSE)
#lapply(Ch2, function(x) plot(x ~ I))

fit_line_Ch3 <- lapply(coef_model_Ch3, function(x) Michaelis_Menten(x,I3))
lapply(fit_line_Ch3, function(x) plot(x ~ I3,type="l"))
write.table(fit_line_Ch3,"fit_line_Ch3.txt", sep = "\t", row.names = FALSE)
#lapply(Ch1, function(x) plot(x ~ I))

fit_line_Ch4 <- lapply(coef_model_Ch4, function(x) Michaelis_Menten(x,I4))
lapply(fit_line_Ch4, function(x) plot(x ~ I4,type="l"))
write.table(fit_line_Ch4,"fit_line_Ch4.txt", sep = "\t", row.names = FALSE)
#lapply(Ch4, function(x) plot(x ~ I))



##########################################
#
# Check for outliers 

old <- par(mfrow = c(2,2))
plot(glm(PI_model_R$comp_punt ~ PI_model_R$r....1))
plot(glm(PI_model_R$comp_punt ~ PI_model_R$h_i))
plot(glm(PI_model_R$comp_punt ~ PI_model_R$P_max))
plot(glm(PI_model_R$comp_punt ~ PI_model_R$Temp))
par(old)


summary(glm(PI_model_R$comp_punt ~ PI_model_R$r....1))
summary(glm(PI_model_R$comp_punt ~ PI_model_R$h_i))
summary(glm(PI_model_R$comp_punt ~ PI_model_R$P_max))
summary(glm(PI_model_R$comp_punt ~ PI_model_R$Temp))
