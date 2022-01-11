#' Amount of missing data per variable
#'
#' @description Amount of missing data per variable
#' @param dat Data frame that you want to inspect.
#'
#' @export
#'
#' @examples
#' data(mtcars)
#' mtcars[c(2:3), 4] <- NA #create NAs
#' nmiss(mtcars)
#'
#'@return
#'By default, this function will print the following items to the console
#'\itemize{
#' \item The percent of missing data per variable.
#' \item The percent of complete cases (range: 0 to 1).
#' \item Suggested number of datasets to impute when using multiple imputation.
#'}
#'
nmiss <- function(dat){
  pmiss <- apply(dat, 2, function(x) sum(is.na(x))/length(x))
  pcomp <- sum(complete.cases(dat)) / nrow(dat)
  cat("Percent missing per variable:\n")
  print(pmiss)
  cat("\nPercent complete cases:", pcomp, "\n")
  cat("Number to impute:", round((1 - pcomp) * 100,0), '\n')
}
