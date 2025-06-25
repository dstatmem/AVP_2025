

#' calculates the maximum photosynthesis, respiration rate and half saturation constant from the hourly change in a gas and the light intensity
#'
#' @param ChangeGas_Chl A list containing the rate of change of oxygen or carbon in umol/L, Change_perHour can be used
#' @param SampleNR Numeric sample coding
#' @param PAR The light intensity given in umol/m2/s
#' @param Gas Measured gas: DO = oxygen and DIC = carbon
#' @return A list per sampling NR containing the maximum photosynthesis umol/L, respiration umol/L, and the half saturation rate for light.
#' @examples
#' PI(DO_change, SampleNR, PAR)

PI_curve <- function(ChangeGas_Chl, SampleNR, PAR, Gas) {

if(gas == "DIC"){
  ChangeGas_Chl*-1
}

  DF_Metabolic_PI <-  ChangeGas_Chl
  DF_Metabolic_PI$change_Chl <-  ChangeGas_Chl
  DF_Metabolic_PI$Sample <-  Sample
  DF_Metabolic_PI$PAR <-  PAR

DF_Metabolic_PI <- ChangeGas_Chl[complete.cases(ChangeGas_Chl["ChangeGas_Chl"]),]




# create new lists
NEP <- list()
for (i in 1:max(DF_Metabolic_PI$Sample)) {
  NEP[[i]] <- DF_Metabolic_PI$change_Chl[DF_Metabolic_PI$Sample == i]
}

I <- list()
for (i in 1:max(DF_Metabolic_PI$Sample)) {
  I[[i]] <- DF_Metabolic_PI$PAR[DF_Metabolic_PI$Sample == i]
}



# Fit the non-linear model
PI_field_coef_DIC <- list()


for (i in 1:max(DF_Metabolic$Sample)){
  tryCatch({
    I_v <- I[[i]]
    PI_field_coef_DIC[[i]] <- coef(nls(NEP[[i]] ~ r + (P_max * I_v) / (h_i + I_v),  start=list(r=-0.04, P_max=0.16, h_i=168)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

return(PI_field_coef_DIC)

}

