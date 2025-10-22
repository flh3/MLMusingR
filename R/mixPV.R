#' Fit Weighted Multilevel Models Using Plausible Values
#'
#' Helper function to fit multilevel models with plausible values
#' using weights at different levels using the mix function from the WeMix package (Bailey et al., 2023): see https://cran.r-project.org/web/packages/WeMix/WeMix.pdf.
#'
#' @importFrom stats formula residuals var coef pt model.matrix family weights fitted.values
#' @importFrom methods is
#' @importFrom WeMix mix
#' @importFrom parallel detectCores makeCluster clusterEvalQ clusterExport parLapply
#' @param fml The model formula. Multiple plausible values are specified using the form: \code{pv1 + pv2 + pv3 ~ x1} (depending how many PVs are present).
#' @param mc Option to use multiple cores to speed up processing (set to FALSE by default).
#' @param data Merged dataset to analyze (containing variables at different levels).
#' @param silent Option to show which plausible value is being analyzed (set to FALSE by default).
#' @param ... Options that are used by the mix function in the WeMix package.
#' @return A list object of \code{mix} results. Results are pooled using the summary function.
#'
#' @references
#' \cite{Huang, F. (2024). Using plausible values when fitting multilevel models with large-scale assessment data using R. Large-scale Assessments in Education, 12(7).
#' (\href{https://largescaleassessmentsineducation.springeropen.com/articles/10.1186/s40536-024-00192-0}{link})}
#'
#'
#' @author Francis Huang, \email{huangf@missouri.edu}
#'
#'
#' @examples
#' \dontrun{
#' data(pisa2012, package = 'MLMusingR')
#' m1 <- mixPV(pv1math + pv2math + pv3math + pv4math + pv5math ~ escs + (1|schoolid),
#' weights = c('w_fstuwt', 'w_fschwt'), data = pisa2012)
#' summary(m1)
#' }
#' @export
mixPV <- function(fml, data = NULL, mc = FALSE, silent = FALSE, ...){
  # require(WeMix)
  res <- list() #empty list
  xx <- deparse1(fml[[2]])
  xx <- gsub(" ", "", xx)
  outc <- unlist(strsplit(xx, split = '\\+'))
  #nout <- length(outc) #number of outcomes
  pred <- fml[[3]]

  if(mc == FALSE){
    res <- lapply(outc, function(x){
      if (silent == FALSE) cat("Analyzing plausible value:", x, "\n")
      newf <- reformulate(deparse(pred), response = x)
      WeMix::mix(newf, data = data, ...)}
    )
  } else { ## added multicore processing
    ### using parallel processing, this works with Windows
    # require(parallel)
    cores <- detectCores() - 1
    if(cores < 2) (stop("Unable to use multiple cores on this computer."))
    cat("Attempting to use", cores, "cores. Progress will not be displayed. Please wait.\n")
    cl1 <- makeCluster(cores)
    clusterEvalQ(cl1,
                 library(WeMix)
    )

    xx <- data
    clusterExport(cl1, varlist = "xx", envir = environment())
    res <- parLapply(cl1, outc, function(x){
      # print(x)
      newf <- reformulate(deparse(pred), response = x)
      WeMix::mix(newf, data = xx, ...)

    })
    parallel::stopCluster(cl1)

  }

  class(res) <- c("mixPV", 'list')
  return(res)
}
