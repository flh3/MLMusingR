#' Group-mean center a variable
#'
#' Also referred to as centering within cluster (or within context) or demeaning the variable.
#' By default, uses \code{na.rm = TRUE} when computing group means.
#'
#' @param x Variable to center (e.g., \code{dataframe$varname}).
#' @param grp Cluster/grouping variable (e.g., \code{dataframe$cluster}).
#'
#' @export
#' @return
#' A vector of group-mean centered variables.
#'
#' @examples
#' data(mtcars)
#' #create a group centered variable
#' mtcars$mpg.gpc <- group_center(mtcars$mpg, mtcars$cyl)
group_center <- function(x, grp) {
  grp <- as.numeric(as.factor(as.character(grp)))
  return(as.numeric(x - tapply(x, grp, mean, na.rm = TRUE)[grp]))
}
#' @rdname group_center
#' @export
cwc <- group_center
