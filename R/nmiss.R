#' Amount of missing data per variable
#'
#' @description Amount of missing data per variable
#' @param dat Data frame that you want to inspect
#'
#' @return
#' @export
#'
#' @examples
#' data(mtcars)
#' mtcars[c(2:3), 4] <- NA #create NAs
#' nmiss(mtcars)
#'
nmiss <- function(dat){
  pmiss <- apply(dat, 2, function(x) sum(is.na(x))/length(x))
  pcomp <- sum(complete.cases(dat)) / nrow(dat)
  cat("Percent missing per variable:\n")
  print(pmiss)
  cat("\nPercent complete cases:", pcomp, "\n")
  cat("Number to impute:", round((1 - pcomp) * 100,0), '\n')
}
