#' Scale of Sampling Weights
#'
#' Uses the \code{cluster} and \code{ecluster} (cluster size and effective cluster size) options specified in Mplus. See \href{https://statmodel.com/download/Scaling3.pdf}{note} from the Mplus website.
#' If there is no variation in weights within a cluster, the weights will scale to 1.
#' @param cluster The cluster variable.
#' @param data The original dataset.
#' @param wt The weight variable to scale.
#' @param type Either \code{cluster} or \code{ecluster}. See pdf from Mplus website.
#' @examples
#' data(pisa2012, package = 'MLMusingR')
#' pisa2012$clustwt <- wscale('schoolid', pisa2012, 'w_fschwt')
#' @export
wscale <- function(cluster, data, wt, type = 'cluster'){

  if(type != 'cluster' & type != 'ecluster') {warning("Invalid scaling type.")}
  if(sum(is.na((data[,c(cluster, wt)]))) > 0) warning('Missing value/s in cluster or weight variable. Inspect your data.')

  if(type == 'cluster'){
    ns <- as.numeric(ave(data[, cluster], data[, cluster], FUN = length)) #how many in cluster (numerator)
    swt <- ave(data[, wt], data[, cluster], FUN = sum) #sum of wij (denominator)
    swgt <- data[, wt] * (ns / swt) #wij x adjustment
  }

  if(type == 'ecluster'){
    num <- ave(data[, wt], data[, cluster], FUN = sum)
    num <- num^2
    den <- ave(data[, wt]^2, data[, cluster], FUN = sum)
    ess <- num / den
    totwgt <- ave(data[, wt], data[, cluster], FUN = sum)
    swgt <- data[, wt] * (ess / totwgt)
  }

  return(swgt) #scaled weight

}
