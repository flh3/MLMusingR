#' Computes the group mean of a variable
#'
#'
#'
#' @param x Variable to compute the mean for (e.g., dataframe$varname)
#' @param grp Cluster/grouping variable (e.g., dataframe$cluster)
#'
#' @export
#' @import stats
#' @examples
#' data(mtcars)
#' #create a group centered variable
#' mtcars$mpg.barj <- group_mean(mtcars$mpg, mtcars$cyl)
group_mean <- function(x, grp) {
  #grp <- as.numeric(as.factor(grp))
  if(is.factor(x) == TRUE){
      ml <- length(levels(x))
      if(ml == 2){
      cat("Variable represents percent:", levels(x)[ml], '\n')
      x <- as.numeric(x) - 1
    } else {
      print("More than two levels in the factor.")
    }
  }
  return(ave(x, grp, FUN = function(x) mean(x, na.rm = TRUE)))
}
