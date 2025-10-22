#' @export
#' @importFrom stats printCoefmat
print.CR2 <- function(x, ...){
  cat("\nStandard error type =", x$crtype, '\n')
  cat("Degrees of freedom =", x$df, '\n\n')
  printCoefmat(x$ttable, x$digits)
}

#' Compute the inverse square root of a matrix
#'
#' From Imbens and Kolesar (2016).
#' @param A The matrix object.
#' @export
MatSqrtInverse <- function(A) {

  ei <- eigen(A, symmetric = TRUE) #obtain eigenvalues and eigenvectors
  d <- pmax(ei$values, 10^-12) #set negatives values to zero
  #or near zero 10^-12
  d2 <- 1/sqrt(d) #get the inverse of the square root
  d2[d == 0] <- 0
  ei$vectors %*% (if (length(d2)==1) d2 else diag(d2)) %*% t(ei$vectors)
}

#' Helper function to allow use with modelsummary
#'
#' @param x The mixPV model object.
#' @param dfadj If set to TRUE (default), uses newer df computation. If FALSE, uses standard Rubin pooling formula.
#' @param ... Additional unspecified options.
#' @export
tidy.mixPV <- function(x, dfadj = TRUE, ...){
  m <- length(x)

  re <- sapply(x, FUN = function(x) x$vars)
  ns <- nrow(re)
  rese <- sapply(x, FUN = function(x) x$varDF$SEvcov)
  cfs <- rbind(sapply(x, coef), re)
  ses <- rbind(sapply(x, FUN = function(x) x$SE), rese[1:ns,]) #this is the SE

  cfs.res <- rowMeans(cfs)
  if(names(cfs.res)[1] == '') names(cfs.res)[1] <- "(Intercept)"
  B <- apply(cfs, 1, var) #Vb
  ubar <- apply(ses, 1, FUN = function(x) mean(x^2)) #Vw
  adj <- (1 + (1/m))
  combvar <- ubar + adj * (B)
  ses.res <- sqrt(combvar)
  tstat <- cfs.res / ses.res
  dof <- (m - 1) * (1 + ((m * ubar) / ((m + 1) * B)))^2
  #from Graham
  RIV <- (B + (B / m)) / ubar #same
  term <- names(cfs.res)

  if (dfadj){
    ###
    ns2 <- summary(x[[1]])$ngroups[1]
    k <- ns
    lambda <- RIV / (1 + RIV)
    adj2 <- ((ns2 - k) + 1) / ((ns2 - k) + 3)
    dfobs <- adj2 * ((ns2 - k) * (1 - lambda))
    dfold <- (m - 1) / lambda^2 # Rubin's way [no small sample adj]
    # dof <- (m - 1) * (1 + ((m * ubar) / ((m + 1) * B)))^2 #same
    ## adjusted / Barnard Rubin for small sample
    dof <- (dfold * dfobs) / (dfold + dfobs) #mice way; Barnard & Rubin (1999)
  }

  pv <- 2 * pt(-abs(tstat), dof)
  crit <- abs(qt(p = .025, df = dof)) #.05 / 2
  conf.low <- cfs.res - crit * ses.res
  conf.high <- cfs.res + crit * ses.res
  final <- data.frame(term = term,
                      estimate = cfs.res,
                      std.error = ses.res,
                      statistic = tstat,
                      dof = dof,
                      conf.low = conf.low,
                      conf.high = conf.high,
                      p.value = round(pv, 4),
                      RIV = RIV)
  row.names(final) <- NULL
  return(final)
}

#' @export
glance.mixPV <- function(x, ...){
  m <- length(x)
  ns <- summary(x[[1]])$ngroups[1]
  ll <- sapply(x, FUN = function(x) x$lnl)
  p <- (nrow(x[[1]]$varDF) + length(coef(x[[1]])))
  AICbar <- mean(-2 * ll) + (2 * p)
  BICbar <- mean(-2 * ll) + (log(ns) * p)
  gof <- data.frame(Nobs = ns, N.pv = m, AICbar = AICbar,
                    BICbar = BICbar)
  return(gof)
}

#' Use the summary function on a saved list of mixPV results
#' @param x The mixPV object.
#' @export
summary_all <- function(x){
  lapply(x, summary)
}

#' Pool plausible values using Rubin's rules
#'
#' @param Bs The regression coefficients.
#' @param SEs The standard errors.
#' @param ns2 The number of observations.
#' @param dfadj If set to TRUE (default), uses newer df computation. If FALSE, uses standard Rubin pooling formula.
#' @export
pool_pv <- function(Bs, SEs, ns2, dfadj = TRUE){
  cfs <- do.call(rbind, Bs)
  ses <- do.call(rbind, SEs)
  m <- nrow(cfs) #number of imputations
  cfs.res <- colMeans(cfs)
  ns <- ncol(cfs) #number of coefficients / estimates
  # does not inc number of covariances for rs models
  B <- apply(data.frame(cfs[,1:ns]), 2, var) #Vb
  ubar <- apply(data.frame(ses[,1:ns]), 2, FUN = function(x) mean(x^2)) #Vw
  adj <- (1 + (1 / m))
  combvar <- ubar + (adj * B)
  ses.res <- sqrt(combvar)
  tstat <- cfs.res / ses.res
  dof <- (m - 1) * (1 + ((m * ubar) / ((m + 1) * B)))^2

  #from Graham
  RIV <- (B + (B / m)) / ubar #same
  term <- names(cfs.res)

  ###
  if (dfadj){
    k <- ns #no of pred inc intercept
    lambda <- RIV / (1 + RIV)
    adj2 <- ((ns2 - k) + 1) / ((ns2 - k) + 3)
    dfobs <- adj2 * ((ns2 - k) * (1 - lambda))
    dfold <- (m - 1) / lambda^2 # Rubin's way same [no small sample adj]
    # dof <- (m - 1) * (1 + ((m * ubar) / ((m + 1) * B)))^2 #same
    ## adjusted / Barnard Rubin for small sample
    dof <- (dfold * dfobs) / (dfold + dfobs) #mice way; Barnard & Rubin (1999)
  }

  ###
  pv <- 2 * pt(-abs(tstat), dof)
  crit <- abs(qt(p = .025, df = dof)) #.05 / 2
  conf.low <- cfs.res - crit * ses.res
  conf.high <- cfs.res + crit * ses.res
  output <- cbind(
    estimate = cfs.res,
    std.error = ses.res,
    statistic = tstat,
    df = dof,
    conf.low = conf.low,
    conf.high = conf.high,
    p.value = round(pv, 4),
    RIV = RIV,
    "Pr(>t)" = round(pv, 4)
  )
  return(output)

}

#' Create summary output from the mixPV function
#'
#' @param object The mixPV object
#' @param dfadj If set to TRUE (default), uses newer df computation. If FALSE, uses standard Rubin pooling formula.
#' @param ... Additional unspecified options.
#' @export
summary.mixPV <- function(object, dfadj = TRUE, ...){
  # require(broom)
  m <- length(object)
  ns <- summary(object[[1]])$ngroups[1]
  #random effects
  re <- lapply(object, FUN = function(x) x$vars)
  rese <- lapply(object, FUN = function(x) x$varDF$SEvcov)
  re.final <- pool_pv(re, rese, ns, dfadj)

  #fixed effects
  cfs <- lapply(object, coef)
  ses <- lapply(object, FUN = function(x) x$SE) #this is the SE
  ttable <- pool_pv(cfs, ses, ns, dfadj) #added ns

  res <- list(ttable = ttable[drop = F],
              reff = re.final,
              gof = glance(object)
  )
  class(res) <- c("mPV", "list")
  return(res)
}

#' @export
print.mPV <- function(x, ...){
  cat("Results of multilevel analyses with", x$gof$N.pv, "plausible values.\n")
  cat("Number of observations:", x$gof$Nobs, '\n')
  cat("\nEstimates for random effects: \n")
  printCoefmat(x$reff[,-c(5:8), drop = FALSE], digits = 3)
  cat("\nEstimates for fixed effects: \n")
  printCoefmat(x$ttable[,-c(5:8), drop = FALSE], digits = 3)

}

