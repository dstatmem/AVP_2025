#'
#' @title Calculates the moving average
#' @description Calculates the moving average over a number of data points given by steps
#' @param Data A list with numeric values
#' @param Steps The range of averaging given by natural number,
#' @return The moving average of \code{Data} over the number of \code{Steps}.
#' @examples
#' Moving_average(DIC, 4)

moving_average <- function(Data, Steps) {

  Steps1 <- rep(1/Steps, Steps)

  Out <- stats::filter(Data, Steps1, sides=2)

  return(Out)

}
