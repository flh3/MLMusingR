#' Group-mean center a variable
#'
#' Also referred to as centering within cluster.
#'
#' @param var Variable to center (e.g., dataframe$varname)
#' @param grp Cluster/grouping variable (e.g., dataframe$cluster)
#'
#' @return
#' @export
#'
#' @examples
#' data(mtcars)
#' #create a group centered variable
#' mtcars$mpg.gpc <- groupcenter(mtcars$mpg, mtcars$cyl)
groupcenter <- function(var, grp) {
  grp <- as.numeric(as.factor(grp))
  return(var - tapply(var, grp, mean, na.rm = TRUE)[grp])
}
