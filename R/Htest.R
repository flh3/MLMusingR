#' Test for homoskedasticity at level one
#'
#' @param newdata data to be used
#' @param fml level 1 formula
#' @param group grouping variable (in quotes)
#'
#' Based on Raudenbush and Bryk (2002).
#' A statistically significant Chisq indicates heteroskedasticity.
#' Output shows the H statistic, degrees of freedom, and p value.
#'
#' @importFrom stats formula complete.cases lm nobs pchisq
#' @export
Htest <- function(newdata, fml, group){

  gps <- names(table(newdata[,group]))
  gno <- length(gps)

  df <- matrix(NA, nrow = gno, ncol = 3)
  for (i in 1:gno){

    ss <- newdata[newdata[,group] == gps[i], ]

    Xm <- model.matrix(fml, data = ss)
    vr <- apply(Xm, 2, var)[-1] #detect variability  in predictors
    if (any(vr != 0)){
      tmp <- (lm(formula = fml, data = ss))
      df[i, 1] <- as.numeric(gps[i])
      df[i, 2] <- summary(tmp)$sigma
      df[i, 3] <- nobs(tmp)
    }
  }
  df2 <- data.frame(df)

  names(df2) <- c(group, 'rmse', 'n')
  df2 <- subset(df2, !is.na(rmse))

  deg <- tmp$rank  #includes intercept
  tst2 <- df2
  tst2$df <- tst2$n - deg
  tst2 <- subset(tst2, df > 9) #NOTE: RB 264 uses n :
  #Hoffman uses df
  tst2$lnv <- log(tst2$rmse^2) * tst2$df
  tot <- sum(tst2$lnv) / sum(tst2$df) #total

  tst2$d <- sqrt(tst2$df/2) * ((log(tst2$rmse ^2) - tot))
  tst2$d2 <- tst2$d^2
  ngps <- nrow(tst2) - 1 #used in the summation
  Hind <- sum(tst2$d2)
  data.frame(H = Hind, df = ngps, p = pchisq(Hind, ngps, lower.tail = F))
}
