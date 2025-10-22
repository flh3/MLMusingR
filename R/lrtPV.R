#' Likelihood Ratio Test with Model Results Using Plausible Values
#'
#' Compares two nested models (a full and a reduced model). Results in an F statistic (not the traditional chi-square) with a p-value (see Huang, 2024).
#' The full model must come first. Statistically significant results indicate that the full model
#' fits better than the reduced model. Uses computations shown by Li et al. (1991).
#'
#' @param mf The full model object fit using mixPV.
#' @param mr The reduced model object fit using mixPV.
#'
#' @references
#' \cite{Huang, F. (2024). Using plausible values when fitting multilevel models with large-scale assessment data using R. Large-scale Assessments in Education, 12(7).
#' (\href{https://largescaleassessmentsineducation.springeropen.com/articles/10.1186/s40536-024-00192-0}{link})}
#'
#' \cite{Li, K. H., Meng, X.L., Raghunathan, T. E., & Rubin, D. B. (1991). Significance levels from repeated p-values with multiply imputed data. Statistica Sinica, 65–92.
#' }
#'
#'
#' @examples
#' \dontrun{
#' data(pisa2012, package = 'MLMusingR')
#' reduced <- mixPV(pv1math + pv2math + pv3math + pv4math + pv5math ~
#'  escs + (1|schoolid), data = pisa2012,
#'  weights = c('w_fstuwt', 'w_fschwt'))
#'full <- mixPV(pv1math + pv2math + pv3math + pv4math + pv5math ~
#'  escs + (escs|schoolid), data = pisa2012,
#'  weights = c('w_fstuwt', 'w_fschwt'))
#'lrtPV(full, reduced)
#' }
#'
#' @export
lrtPV <- function(mf, mr){ #for mixPV
  nll <- sapply(mr, FUN = function(x) x$lnl)
  fll <- sapply(mf, FUN = function(x) x$lnl)
  a1 <- summary(mr[[1]])
  a2 <- summary(mf[[1]])
  k.r <- length(mr[[1]]$theta) + length(mr[[1]]$coef)
  k.f <- length(mf[[1]]$theta) + length(mf[[1]]$coef)
  k <- k.f - k.r
  lldif <- -2 * (nll - fll)
  m <- length(nll)
  dbar <- mean(lldif)
  r2 <- (1 + (1 / m)) * var(sqrt(lldif))
  Fs <- (dbar/k - ((m + 1)/(m - 1)) * r2) / (1 + r2) #from package
  v <- k^(-3 / m) * (m - 1) * (1 + r2^(-1))^2
  pv <- pf(Fs, k, v, lower.tail = FALSE)
  return(data.frame(F = Fs, df1 = k, df2 = v, r = r2, pv = round(pv, 4)))
}
