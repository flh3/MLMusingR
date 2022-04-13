#' Computes the group mean of a variable
#'
#' Computes the group means of a variable by a specified cluster/group. Can also be used with factors that have two levels.
#'
#' @param x Variable to compute the mean for (e.g., \code{dataframe$varname}).
#' @param grp Cluster/grouping variable (e.g., \code{dataframe$cluster}).
#' @param lm Compute reliability (lambda) adjusted means.
#' @return
#' Outputs a vector of group means.
#' @export
#' @import stats
#' @examples
#' data(mtcars)
#' #create a group mean aggregated variable
#' mtcars$mpg.barj <- group_mean(mtcars$mpg, mtcars$cyl)
group_mean <- function(x, grp, lm = FALSE) {
  #grp <- as.numeric(as.factor(grp))

  if(is.factor(x) == TRUE){
    ml <- length(levels(x))
    if(ml == 2){
      cat("Variable represents percent:", levels(x)[ml], '\n')
      x <- as.numeric(x) - 1
    } else {
      stop("More than two levels in the factor. Only works with two levels.")
    }
  }

  if (lm == FALSE){
    # these are the group means
    return(ave(x, grp, FUN = function(x) mean(x, na.rm = TRUE)))

  } else { #latent means :: test procedure

    g <- grp
    G <- length(table(g))
    freq <- data.frame(table(g))
    n <- length(g)
    scaling <- (n^2 - sum(freq$Freq^2)) / (n * (G - 1))

    ms <- ave(x, g, FUN = function(x) mean(x, na.rm = TRUE))
    cs <- x - ms #centered scores

    b.cov <- (var(ms, na.rm = T) * (n - 1)) / (G - 1) #group level cov matrix
    w.cov <- (var(cs, na.rm = T) * (n - 1)) / (n - G) #individual level cov matrix
    pb.cov <- (b.cov - w.cov) / scaling #estimate of pure/adjusted between cov matrix
    #tot.v <- w.cov + pb.cov
    #icc <- round((pb.cov / tot.v), 3)
    #print(icc)
    ns <- ave(x, g, FUN = length)

    lambda <- pb.cov / (pb.cov + (w.cov / ns)) #reliability
    ovm <- mean(x, na.rm = TRUE) #overall grand mean
    ms <- ave(x, grp, FUN = function(x) mean(x, na.rm = TRUE))
    #print(lambda)
    (1 - lambda) * ovm + (lambda * ms)

  }

}
