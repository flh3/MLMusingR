#' Robust standard errors for mixed models
#'
#' @param m1 lme4 or nlme model object
#' @param digits Number of digits for output
#' @param satt Satterthwaite degrees of freedom approximation; TRUE or FALSE
#' @param Gname Group/cluster name if more than two levels of clustering
#'
#' If there are more than two levels of clustering, the clustering variable should
#' set at the highest level
#'
#' @return
#' @export
robust_mixed <- function(m1, digits = 4, satt = TRUE, Gname = NULL){

  if(class(m1) %in%  c('lmerMod', 'lmerModLmerTest')){ #if lmer

    X <- model.matrix(m1) #X matrix
    B <- fixef(m1) #coefficients
    y <- m1@resp$y #outcome
    Z <- getME(m1, 'Z') #sparse Z matrix
    b <- getME(m1, 'b') #random effects

    if (is.null(Gname)){
    Gname <- names(getME(m1, 'l_i')) #name of clustering variable
    if (length(Gname) > 1) {
      stop("lmer: Can only be used with non cross-classified data. If more than two levels, specify Gname = 'clustername'")
    }
    }

    js <- table(m1@frame[, Gname]) #how many observation in each cluster
    G <- bdiag(VarCorr(m1)) #G matrix

    #re <- as.numeric(y - (X %*% B + Z %*% b)) #not used, just checking
    #data.frame(re, resid(m1)) #the same
    #cor(re, resid(m1)) #1

    qq <- getME(m1, 'q') #columns in RE matrix
    NG <- getME(m1, 'l_i') #number of groups :: ngrps(m1)
    NG <- NG[length(NG)]
    # nre <- getME(m1, 'p_i') #qq/NG --> number of random effects
    # inde <- cumsum(js) #number per group summed to create an index


    ### The following is used to create the V matrix
    ### Probably other (better and faster) ways to to do this but I think this is the most transparent
    ### and works with my basic knowledge of matrices in R

    gpsv <- m1@frame[, Gname] #data with groups

    # { #done a bit later than necessary but that is fine
    #   if(is.unsorted(gpsv)){
    #     stop("Data are not sorted by cluster. Please sort your data first by cluster, run the analysis, and then use the function.\n")
    #   }
    # }

  #   ml <- list() #empty list to store matrices
  #   cols <- seq(1, (NG * nre), by = nre) #dynamic columns for Z matrix
  #
  #   for (i in 1:length(inde)){
  #
  #     if (i == 1) {
  #       st = 1} else {
  #         st = inde[i - 1] + 1}
  #     end = st + js[i] - 1
  #
  #     nc <- cols[i]
  #     ncend <- cols[i] + (nre - 1)
  #
  #     Zi <- data.matrix(Z[st:end, nc:ncend]) #depends on how many obs in a cluster and how many rand effects
  #     ml[[i]] <- Zi %*% G %*% t(Zi) + diag(sigma(m1)^2, nrow = js[i]) #ZGZ' + r
  #   }
  #
  #   Vm <- bdiag(ml) #makes a block diagonal weighting matrix
  # }

    getV <- function(x){
      lam <- data.matrix(getME(x, 'Lambdat'))
      var.d <- crossprod(lam)
      # var.d <- t(lam) %*% lam
      Zt <- data.matrix(getME(x, "Zt"))
      vr <- sigma(x)^2
      var.b <- vr * (t(Zt) %*% var.d %*% Zt)
      sI <- vr * diag(nobs(x))
      var.y <- var.b + sI
    }

    Vm <- getV(m1)

  }

  if(class(m1) == 'lme'){ #if nlme
    require(nlme)
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

    { #done a bit later than necessary but that is fine
      if(is.unsorted(gpsv)){
       stop("Data are not sorted by cluster. Please sort your data first by cluster, run the analysis, and then use the function.\n")
      }
    }

    ml <- list()
    for (j in 1:NG){
      test <- getVarCov(m1, individuals = j, type = 'marginal')
      ml[[j]] <- test[[1]]
    }

    Vm <- bdiag(ml)
  }



  ### robust computation :: once all elements are extracted
  rr <- y - X %*% B #residuals with no random effects
  cdata <- data.frame(cluster = gpsv, r = rr)
  k <- ncol(X) #
  gs <- names(table(cdata$cluster)) #name of the clusters
  u <- matrix(NA, nrow = NG, ncol = k)

  # correct
  # for(i in 1:NG){
  #   tmp <- js[i] #how many in group
  #   u[i,] <- as.numeric(t(cdata$r[cdata$cluster == gs[i]]) %*% solve(ml[[i]]) %*% X[gpsv == gs[i], 1:k])
  # }

  for(i in 1:NG){ #2021.06.18
    sel <- cdata$cluster == gs[i] #selection per cluster
    u[i,] <- as.numeric(t(cdata$r[sel]) %*% solve(Vm[sel, sel]) %*% X[gpsv == gs[i], 1:k])
  }



  ## e' (Zg)-1 Xg
  ## putting the pieces together
  Vinv <- solve(Vm)
  br2 <- solve(t(X) %*% Vinv %*% X) #bread
  mt <- t(u) %*% u #meat :: t(u) %*% u
  clvc2 <- br2 %*% mt %*% br2
  rse <- sqrt(diag(clvc2))

  ### HLM dof
  chk <- function(x){
    vrcheck <- sum(tapply(x, gpsv, var), na.rm = T) #L1,
    # na needed if only one observation with var = NA
    y <- 1 #assume lev1 by default
    if (vrcheck == 0) (y <- 2) #if variation, then L2
    return(y)
  }

  levs <- apply(X, 2, chk) #all except intercept
  levs[1] <- 1 #intercept

  tt <- table(levs)
  l1v <- tt['1']
  l2v <- tt['2']

  l1v[is.na(l1v)] <- 0
  l2v[is.na(l2v)] <- 0


  ####
  n <- nobs(m1)
  #ns <- nobs(mod)
  df1 <- n - l1v - l2v
  df2 <- NG - l2v - 1

  dfn <- rep(df1, length(levs)) #naive
  dfn[levs == '2'] <- df2

  if (satt == T){
    dfn <- satdf(m1)
  }

  robse <- as.numeric(rse)
  FE_auto <- fixef(m1)
  statistic <- FE_auto / robse
  p.values = round(2 * pt(-abs(statistic), df = dfn), digits)

  stars <- cut(p.values, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
               labels = c("***", "**", "*", ".", " "), include.lowest = TRUE)

  ################# COMPARE RESULTS
  # gams <- solve(t(X) %*% solve(Vm) %*% X) %*% (t(X) %*% solve(Vm) %*% y)

  gams <- br2 %*% t(X) %*% Vinv %*% y
  # SEm <- as.numeric(sqrt(diag(solve(t(X) %*% solve(Vm) %*% X)))) #X' Vm-1 X
  SEm <- as.numeric(sqrt(diag(br2))) #X' Vm-1 X
  # SE <- as.numeric(sqrt(diag(data.matrix(vcov(m1))))) #compare standard errors
  return(data.frame(
    estimate = as.numeric(gams),
    #FE_auto,
    mb.se = SEm,
    #SE_auto = SE,
    robust.se = robse,
    df = round(dfn, 1),
    p.values,
    Sig = stars
  )
  )

}

#' Title
#'
#' @param A xxx
#'
#' @return
#' @export
#'
MatSqrtInverse <- function(A) {
  ##  Compute the inverse square root of a matrix
  ei <- eigen(A) #obtain eigenvalues and eigenvectors
  d <- pmax(ei$values, 0) #set negatives values to zero
  d2 <- 1/sqrt(d) #get the inverse of the square root
  d2[d == 0] <- 0
  ei$vectors %*% (if (length(d2)==1) d2 else diag(d2)) %*% t(ei$vectors)
}


## empirical DOF

#' Title
#'
#' @param mod xxx
#'
#' @return
#' @export
#'
satdf <- function(mod){
  if(class(mod) == 'lme') {
    dat <- mod$data
    fml <- formula(mod)
    X <- model.matrix(fml, data = dat)
    Gname <- names(mod$groups)
    gpsv <- dat[,Gname]

  } else if (class(mod) %in% c('lmerMod', 'lmerModLmerTest')){ #if lmer
    dat <- mod@frame
    X <- model.matrix(mod)
    xx <- length(getME(mod, 'l_i')) #number of groups :: ngrps(m1)

    Gname <- names(getME(mod, 'l_i'))[xx]
    gpsv <- mod@frame[, Gname]
  } else {
    stop("Type of object is not an lmer or lme object.")
  }

  cnames <- names(table(gpsv))
  cpx <- solve(crossprod(X))
  cdata <- data.frame(cluster = dat[,Gname])
  NG <- length(cnames)

  ## STEP 1
  tXs <- function(s) {
    Xs <- X[cdata$cluster == s, , drop = F]
    MatSqrtInverse(diag(NROW(Xs)) - Xs %*% cpx %*% t(Xs)) %*%
      Xs
  } # A x Xs / Need this first

  tX <- lapply(cnames, tXs)

  ## STEP 2
  tHs <- function(s) {
    Xs <- X[cdata$cluster == s, , drop = F]
    index <- which(cdata$cluster == s)
    ss <- matrix(0, nrow = n, ncol = length(index)) #all 0, n x G
    ss[cbind(index, 1:length(index))] <- 1 #indicator
    ss - X %*% cpx %*% t(Xs) #overall X x crossprod x Xs'
  }

  n <- nobs(mod)
  tH <- lapply(cnames, tHs) #per cluster

  ## STEP 3
  k <- ncol(X)
  id <- diag(k) #number of coefficients // for different df
  degf <- numeric(k) #vector for df // container

  for (j in 1:k){ #using a loop since it's easier to see

    Gt <- sapply(seq(NG), function(i) tH[[i]] %*%
                   tX[[i]] %*% cpx %*% id[,j])
    #already transposed because of sapply: this is G'
    #ev <- eigen(Gt %*% t(Gt))$values #eigen values: n x n
    ev <- eigen(t(Gt) %*% Gt)$values #much quicker this way, same result: p x p:
    degf[j] <- (sum(ev)^2) / sum(ev^2) #final step to compute df
  }

  return(degf)
}
