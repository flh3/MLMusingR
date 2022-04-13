#' Cluster robust standard errors with degrees of freedom adjustments for lmerMod/lme objects
#'
#' Function to compute the CR2/CR0 cluster
#' robust standard errors (SE) with Bell and McCaffrey (2002)
#' degrees of freedom (dof) adjustments. Suitable even with a low number of clusters.
#' The model based (mb) and cluster robust standard errors are shown for comparison purposes.
#'
#'
#' @importFrom stats nobs resid formula residuals var coef pt model.matrix family weights fitted.values
#' @param m1 The \code{lmerMod} or \code{lme} model object.
#' @param digits Number of decimal places to display.
#' @param type Type of cluster robust standard error to use ("CR2" or "CR0").
#' @param satt If Satterthwaite degrees of freedom are to be computed (if not, between-within df are used).
#' @param Gname Group/cluster name if more than two levels of clustering (does not work with lme).
#' @return A data frame (\code{results}) with the cluster robust adjustments with p-values.
#' \item{Estimate}{The regression coefficient.}
#' \item{mb.se}{The model-based (regular, unadjusted) SE.}
#' \item{cr.se}{The cluster robust standard error.}
#' \item{df}{degrees of freedom: Satterthwaite or between-within.}
#' \item{p.val}{p-value using CR0/CR2 standard error.}
#' \item{stars}{stars showing statistical significance.}
#'
#' @references
#' \cite{Bell, R., & McCaffrey, D. (2002). Bias reduction in standard errors for linear regression with multi-stage samples. Survey Methodology, 28, 169-182.
#' (\href{https://www150.statcan.gc.ca/n1/pub/12-001-x/2002002/article/9058-eng.pdf}{link})}
#'
#' Liang, K.Y., & Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. \emph{Biometrika, 73}(1), 13–22.
#' \doi{10.1093/biomet/73.1.13}
#'
#' @author Francis Huang, \email{huangf@missouri.edu}
#' @author Bixi Zhang, \email{bixizhang@missouri.edu}
#'
#' @examples
#' require(lme4)
#' data(sch29, package = 'MLMusingR')
#' robust_mixed(lmer(math ~ male + minority + mses + mhmwk + (1|schid), data = sch29))
#' @export
robust_mixed <- function(m1, digits = 3, type = 'CR2', satt = TRUE, Gname = NULL){

  ### for lmer
  if(class(m1) %in%  c('lmerMod', 'lmerModLmerTest')){ #if lmer
    dat <- m1@frame
    X <- model.matrix(m1) #X matrix
    B <- fixef(m1) #coefficients
    y <- m1@resp$y #outcome
    Z <- getME(m1, 'Z') #sparse Z matrix
    b <- getME(m1, 'b') #random effects

    if (is.null(Gname)){
      Gname <- names(getME(m1, 'l_i')) #name of clustering variable
      if (length(Gname) > 1) {
        stop("lmer: Can only be used with non cross-classified data. If more than two levels, specify highest level using Gname = 'clustername'")
      }
    }

    js <- table(dat[, Gname]) #how many observation in each cluster
    G <- bdiag(VarCorr(m1)) #G matrix

    #re <- as.numeric(y - (X %*% B + Z %*% b)) #not used, just checking
    #data.frame(re, resid(m1)) #the same
    #cor(re, resid(m1)) #1
    #qq <- getME(m1, 'q') #columns in RE matrix

    NG <- getME(m1, 'l_i') #number of groups :: ngrps(m1)
    NG <- NG[length(NG)]

    gpsv <- dat[, Gname] #data with groups

    # { #done a bit later than necessary but that is fine
    #   if(is.unsorted(gpsv)){
    #    # stop("Data are not sorted by cluster. Please sort your data first by cluster, run the analysis, and then use the function.\n")
    #   }
    # }

    getV <- function(x) {
      lam <- data.matrix(getME(x, "Lambdat"))
      var.d <- crossprod(lam)
      Zt <- data.matrix(getME(x, "Zt"))
      vr <- sigma(x)^2
      var.b <- vr * (t(Zt) %*% var.d %*% Zt)
      sI <- vr * diag(nobs(x))
      var.y <- var.b + sI
    }
    Vm <- getV(m1)
  }

  ## for nlme
  if(class(m1) == 'lme'){ #if nlme
    dat <- m1$data
    fml <- formula(m1)
    X <- model.matrix(fml, data = dat)
    B <- fixef(m1)
    NG <- m1$dims$ngrps[[1]]
    if (length(m1$dims$ngrps) > 3) {stop("Can only be used with two level data.")}
    Gname <- names(m1$groups)
    y <- dat[,as.character(m1$terms[[2]])]
    gpsv <- dat[,Gname]
    js <- table(gpsv)

    {#done a bit later than necessary but that is fine
      if(is.unsorted(gpsv)){
        stop("Data are not sorted by cluster. Please sort your data first by cluster, run the analysis, and then use the function.\n")
      }
    }

    ml <- list()
    for (j in 1:NG){
      test <- getVarCov(m1, individuals = j, type = 'marginal')
      ml[[j]] <- test[[1]]
    }

    Vm <- as.matrix(Matrix::bdiag(ml)) #to work with other funs
  }

  ### robust computation :: once all elements are extracted
  rr <- y - X %*% B #residuals with no random effects
  #n <- nobs(m1)
  cdata <- data.frame(cluster = gpsv, r = rr)
  k <- ncol(X) #number of predictors (inc intercept)
  gs <- names(table(cdata$cluster)) #name of the clusters
  u <- matrix(NA, nrow = NG, ncol = k) #LZ
  uu <- matrix(NA, nrow = NG, ncol = k) #CR2

  cnames <- names(table(gpsv))

  #cpx <- solve(crossprod(X))
  #cpx <- chol2inv(qr.R(qr(X))) #using QR decomposition, faster, more stable?
  #cdata <- data.frame(cluster = dat[,Gname])
  #NG <- length(cnames)

  ### quicker way, doing the bread by cluster
  tmp <- split(X, cdata$cluster)
  XX <- lapply(tmp, function(x) matrix(x, ncol = k)) #X per clust

  # to get Vc^-1 per cluster
  aa <- function(x){
    sel <- which(cdata$cluster == x)
    chol2inv(chol(Vm[sel, sel])) #this is V^-1
    #solve(Vm[sel, sel])
  }

  Vm2 <- lapply(cnames, aa) #Vc^-1 used

  a2 <- function(x){
    sel <- which(cdata$cluster == x)
    Vm[sel, sel] #this is V

  }

  Vm3 <- lapply(cnames, a2) #Vc used
  names(Vm2) <- names(Vm3) <- cnames #naming
  #Vm2 is the inverse, Vm3 is just the plain V matrix

  Vinv <- as.matrix(Matrix::bdiag(Vm2))
  # to get X V-1 X per cluster
  bb <- function(x){
    t(XX[[x]]) %*% Vm2[[x]] %*% XX[[x]]
  }

  dd <- lapply(cnames, bb)
  br <- solve(Reduce("+", dd)) #bread

  #rrr <- split(rr, getME(m1, 'flist'))
  rrr <- split(rr, cdata$cluster)

#### Meat matrix

  if (type == "CR2"){
  tXs <- function(s) {

    Ijj <- diag(nrow(XX[[s]]))
    Hjj <- XX[[s]] %*% br %*% t(XX[[s]]) %*% Vm2[[s]]
    IHjj <- Ijj - Hjj
    #MatSqrtInverse(Ijj - Hjj) #early adjustment / valid
    V3 <- chol(Vm3[[s]]) #based on MBB
    Bi <- V3 %*% IHjj %*% Vm3[[s]] %*% t(V3)
    t(V3) %*% MatSqrtInverse(Bi) %*% V3

  } # A x Xs / Need this first

  tX <- lapply(cnames, tXs)

  cc2 <- function(x){
    rrr[[x]] %*% tX[[x]] %*% Vm2[[x]] %*% XX[[x]]
  }

  u <- t(sapply(1:NG, cc2)) #using 1:NG instead

} else { ## meat matrix for CR0
  # residual x inverse of V matrix x X matrix
  cc0 <- function(x){
    rrr[[x]] %*% Vm2[[x]] %*% XX[[x]]
  }

  u <- t(sapply(cnames, cc0))
}

  ## e'(Vg)-1 Xg ## CR0
  ## putting the pieces together
  #br2 <- solve(t(X) %*% Vinv %*% X) #bread
  mt <- t(u) %*% u #meat
  clvc2 <- br %*% mt %*% br #variance covariance matrix
  rse <- sqrt(diag(clvc2)) #standard errors

  ### DEGREES OF FREEDOM :::::::::::::::::


  if (satt == TRUE){
    dfn <- satdf(m1, Gname = Gname, type = type)
  } else {

  ### HLM dof
  chk <- function(x){
    vrcheck <- sum(tapply(x, gpsv, var), na.rm = T) #L1,
    # na needed if only one observation with var = NA
    y <- 1 #assume lev1 by default
    if (vrcheck == 0) (y <- 2) #if variation, then L2
    return(y)
  }

  levs <- apply(X, 2, chk) #all except intercept
  # levs[1] <- 1 #intercept

  tt <- table(levs)
  l1v <- tt['1']
  l2v <- tt['2']

  l1v[is.na(l1v)] <- 0
  l2v[is.na(l2v)] <- 0

  ####
  n <- nobs(m1)

  df1 <- n - l1v - length(js)
  df2 <- NG - l2v

  dfn <- rep(df1, length(levs)) #naive
  dfn[levs == '2'] <- df2
}

  robse <- as.numeric(rse)
  FE_auto <- fixef(m1)
  cfsnames <- names(FE_auto)
  statistic <- FE_auto / rse
  p.values = round(2 * pt(abs(statistic), df = dfn, lower.tail = F), digits)
  stars <- as.character(cut(p.values, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   labels = c("***", "**", "*", ".", " "), include.lowest = TRUE)
  )

robs <- rse
pv <- p.values
vc <- clvc2


  #gams <- solve(t(X) %*% solve(Vm) %*% X) %*% (t(X) %*% solve(Vm) %*% y)
  #SEm <- as.numeric(sqrt(diag(solve(t(X) %*% solve(Vm) %*% X)))) #X' Vm-1 X
  #SE <- as.numeric(sqrt(diag(vcov(m1)))) #compare standard errors

  SE <- as.numeric(sqrt(diag(br)))
  ttable <- cbind(
    Estimate = round(FE_auto, digits),
    mb.se = round(SE, digits),
    robust.se = round(robs, digits),
    t.stat = round(FE_auto / robs, digits),
    df = round(dfn, 1),
    "Pr(>t)" = pv
  )
  results <- data.frame(
    Estimate = FE_auto,
    mb.se = SE,
    cr.se = robse,
    t.stat = FE_auto / robs,
    df = dfn,
    p.val = pv,
    stars
  )

  type <-  ifelse(type == 'CR2', 'CR2', 'CR0')
  dft <- ifelse(satt == TRUE, "Satterthwaite", "Between-within")

  res <- list(ttable = ttable,
              results = results,
              crtype = type,
              df = dft,
              digits = digits,
              vcov = vc,
              orig = m1)

  class(res) <- 'CR2'
  return(res)
}


