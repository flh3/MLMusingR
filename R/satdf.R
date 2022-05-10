#' Compute Satterthwaite degrees of freedom
#'
#' Function to compute empirical degrees of freedom
#' based on Bell and McCaffrey (2002).
#'
#'
#' @importFrom stats nobs resid formula residuals var coef pt model.matrix family weights fitted.values
#' @importFrom methods is
#' @param m1 The \code{lmerMod} or \code{lme} model object.
#' @param type The type of cluster robust correction used (i.e., CR2 or none).
#' @param Vinv2 Inverse of the variance matrix.
#' @param Vm2 The variance matrix.
#' @param br2 The bread component.
#' @param Gname The group (clustering variable) name'
#'
#' @author Francis Huang, \email{huangf@missouri.edu}
#' @author Bixi Zhang, \email{bixizhang@missouri.edu}
#'
#' @export
## empirical DOF
satdf <- function(m1, type = 'none', Vinv2, Vm2, br2, Gname = NULL){

  #require(Matrix)

  #if(class(m1) == 'lme'){ #if nlme
  if(is(m1, 'lme')){
    dat <- m1$data
    fml <- formula(m1)
    X <- model.matrix(fml, data = dat)
    B <- fixef(m1)
    NG <- m1$dims$ngrps[[1]]
    if (length(m1$dims$ngrps) > 3) {stop("Can only be used with two level data (for now).")}
    Gname <- names(m1$groups)
    y <- dat[,as.character(m1$terms[[2]])]
    gpsv <- dat[,Gname]
    js <- table(gpsv)
    K <- ncol(X)

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

  ### for lmer
  #if(class(m1) %in%  c('lmerMod', 'lmerModLmerTest')){ #if lmer
  if(is(m1, 'lmerMod')){
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


  Vm <- Matrix::drop0(Vm) #make a sparse matrix, if not already
  Vinv <- Matrix::solve(Vm)

  #Vinv <- chol2inv(chol(Vm))
  cpx <- solve(t(X) %*%  Vinv %*% X) #solve(Vm)
  ns <- nobs(m1)
  Im <- diag(ns) #identity matrix
  Hm <- X %*% cpx %*% t(X) %*% Vinv #Overall hat matrix
  IH <- as.matrix(Im - Hm) #difference

  nms <- names(table(dat[,Gname])) #names of clusters
  K <- ncol(X) #number of vars
  dd <- diag(K)
  NG <- length(nms)

  ### adjustments

  if (type == 'CR2') {

    tHs <- function(x) { #working CR2
      ind <- which(dat[,Gname] == x)
      Xs <- X[ind, ,drop = F]
      Vs <- Vm[ind, ind, drop = F]
      U <- chol(Vs) #with the cholesky matrix
      adj <- Xs %*% cpx %*% t(Xs) %*% chol2inv(U) #solve(Vs)

      Ws <- Vinv[ind, ind, drop = F] #Wj, Vinv in clusters
      ih <- IH[ind, , drop = F] #asymmetric, need rows(ind) here

      ng <- nrow(Xs)
      cr <- diag(ng) - adj

      t(ih) %*% t(U) %*% MatSqrtInverse(U %*% cr %*% Vs %*% t(U)) %*%
        U  %*% Ws %*% Xs %*% cpx ### this has the adjustment in the matsqrtinv & U
      # A(adjust matrix) is t(U) %*% MatSqrtInverse(U %*% cr %*% Vs %*% t(U)) %*% U
      #IHjj <- Ijj - Hjj
      #Bi <- chol(V3) %*% IHjj %*% V3 %*% t(chol(V3))
      #Ai <- t(chol(V3)) %*% MatSqrtInverse(Bi) %*% chol(V3)
    }

  } else {

    tHs <- function(x) { #CR0

      ind <- which(dat[,Gname] == x)
      Xs <- X[ind, ,drop = F]
      ih <- IH[ind, , drop = F]
      Ws <- Vinv[ind, ind, drop = F]
      t(ih) %*% Ws %*% Xs %*% cpx ### NO ADJUSTMENT but with Ws
    }
  }

  tmp <- lapply(nms, tHs)

  #Gm = do.call('cbind', tmp) #bind them together
  degf <- numeric() #container

  #Wm <- MatSqrtInverse(Vm) #as per Tipton 2015 -- this is new
  Wm <- Vm #W is just Vm or target variance in our case

  for (j in 1:K){ #using a loop since it's easier to see
    sel <- dd[j, ]
    Gt <- lapply(seq(NG), function(i) tmp[[i]] %*% sel)
    Gt <- as.matrix(do.call('cbind', Gt))
    #ev <- eigen(Wm %*% Gt %*% t(Gt) %*% Wm)$values
    #degf[j] <- (sum(ev)^2) / sum(ev^2) #final step to compute df
    #GG <- Wm %*% Gt %*% t(Gt) %*% Wm #avoids using eigen; from Kolesar
    #degf[j] <- sum(diag(GG))^2 / sum(GG * GG)

    GG <- t(Gt) %*% Wm %*% Gt  #from Pustejovsky and Tipton 2018 eq.11
    GGd <- GG[row(GG) == col(GG)] #just diag(GG)

    #degf[j] <- sum(diag(GG))^2 / sum(GG * GG) #lme issues?
    degf[j] <- sum(GGd)^2 / sum(GG * GG)
  }

  degf #manual computation for CR2 dof
}
