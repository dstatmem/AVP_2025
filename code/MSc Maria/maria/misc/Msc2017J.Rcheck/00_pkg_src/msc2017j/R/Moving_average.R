#' Calculates the moving average
#'
#' @param Data A list with numeric values
#' @param Steps A real number, range of averaging
#' @return The moving average of \code{Data} for \code{Steps} steps.
#' @examples
#' Moving_average(DIC, 4)

Moving_average <- function(Data, Steps) {

  Steps1 <- rep(1/Steps, Steps)

  Out <- filter(Data, Steps1, sides=2)

  return(Out)

}

