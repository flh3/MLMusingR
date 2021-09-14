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
group_mean <- function(x, grp) {
  #grp <- as.numeric(as.factor(grp))
  return(ave(x, grp, FUN = function(x) mean(x, na.rm = TRUE)))
}
