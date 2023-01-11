#' Tidy a CR2 object
#'
#' @param x A `CR2` object.
#' @param conf.int Logical indicating whether or not to include
#'   a confidence interval in the tidied output. Defaults to FALSE.
#' @param conf.level The confidence level to use for the confidence
#'   interval if conf.int = TRUE. Must be strictly greater than 0
#'   and less than 1. Defaults to 0.95, which corresponds to a
#'   95 percent confidence interval.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing component-level
#'   information about the model
#'
#' @importFrom generics tidy
#' @importFrom stats confint df qt
#' @export
tidy.CR2 <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {

  result <- x$results %>%
    tibble::as_tibble(rownames = "term") %>%
    dplyr::rename(estimate = Estimate,
                  std.error = `cr.se`,
                  statistic = `df`, #why the df
                  p.value = `p.val`)

  if (conf.int) {
    ci <- confint(x, level = conf.level)
    colnames(ci) <- c('conf.low', 'conf.high')
    #result <- dplyr::left_join(result, ci, by = "term")
    result <- dplyr::bind_cols(result, ci)

  }

  return(result)
}

#' @export
vcov.CR2 <- function(object, ...){
  object$vcov
}

#' Glance at goodness-of-fit statistics
#'
#' Helper function used to obtain supporting fit statistics for multilevel models. The R2s are computed using the `performance` package.
#'
#' @param x A `CR2` object.
#' @param ... Unused, included for generic consistency only.
#'
#' @return \code{glance} returns one row with the columns:
#'   \item{nobs}{the number of observations}
#'   \item{sigma}{the square root of the estimated residual variance}
#'   \item{logLik}{the data's log-likelihood under the model}
#'   \item{AIC}{Akaike Information Criterion}
#'   \item{BIC}{Bayesian Information Criterion}
#'   \item{r2.marginal}{marginal R2 based on fixed effects only using method of Nakagawa and Schielzeth (2013)}
#'   \item{r2.conditional}{conditional R2 based on fixed and random effects using method of Nakagawa and Schielzeth (2013)}
#'
#' @rawNamespace if(getRversion()>='3.3.0') importFrom(stats, sigma) else importFrom(lme4,sigma)
#'
#' @importFrom broom glance
#' @export
glance.CR2 <- function(x, ...) {

  tmp <- data.frame(
    sigma = sigma(x$orig),
    logLik = as.numeric(stats::logLik(x$orig)),
    AIC = stats::AIC(x$orig),
    BIC = stats::BIC(x$orig),
    nobs = stats::nobs(x$orig),
    r2.marginal = as.numeric(performance::r2_nakagawa(x$orig)[2]),
    r2.conditional = as.numeric(performance::r2_nakagawa(x$orig)[1])
    # )
  )
  return(tmp)
}


#' @export
confint.CR2 <- function(object, parm, level = 0.95, ...){
  z <- object
  k <- nrow(z$results)
  cf <- z$results$Estimate
  parm <- row.names(z$results)
  se <- z$results$cr.se
  a <- (1 - level) / 2

  crit <- qt(a, z$results$df)
  a <- c(a, 1 - a)
  pct <- sprintf("%0.1f%%", a * 100)
  ci <- array(NA, dim = c(k, 2),
              dimnames = list(parm, pct))

  for (i in 1:k){
    crit <- qt(a, z$results$df[i])
    ci[i, 1] <- cf[i] - abs(crit[1]) * se[i]
    ci[i, 2] <- cf[i] + abs(crit[2]) * se[i]
  }
  ci
}


