#' Group-mean center a variable
#'
#' Also referred to as centering within cluster.
#'
#' @param x Variable to center (e.g., dataframe$varname)
#' @param grp Cluster/grouping variable (e.g., dataframe$cluster)
#'
#' @export
#'
#' @examples
#' data(mtcars)
#' #create a group centered variable
#' mtcars$mpg.gpc <- group_center(mtcars$mpg, mtcars$cyl)
group_center <- function(x, grp) {
  grp <- as.numeric(as.factor(grp))
  return(x - tapply(x, grp, mean, na.rm = TRUE)[grp])
}
