#' Test for homoskedasticity at level one
#'
#' @param newdata data to be used
#' @param fml level 1 formula
#' @param group grouping variable (in quotes)
#'
#' Based on Raudenbush and Bryk (2002).
#' A statistically significant X2 indicates heteroskedasticity.
#'
#' @return
#' @export
Htest <- function(newdata, fml, group){

  gps <- newdata[,group]
  gps <- gps[!duplicated(gps)]
  gno <- length(gps)
  #err <- numeric()
  err <- 0
  df <- matrix(NA, nrow = gno, ncol = 3)
  for (i in 1:gno){
    cont <<- 0
    ss <- newdata[newdata[,group] == gps[i], ]
    # print(i)
    tmp <- tryCatch(lm(formula = fml, data = ss),
                    error = function(cond){
                      print("error")
                      err <<- err + 1
                      cont <<- 1
                    }
    )

    if (cont == 0){
      df[i, 1] <- gps[i]
      df[i, 2] <- summary(tmp)$sigma
      df[i, 3] <- nobs(tmp)
    }
    #cat(i,'::')
  }
  df2 <- data.frame(df)

  names(df2) <- c(group, 'rmse', 'n')
  df2 <- subset(df2, rmse > 0)

  deg <- tmp$rank  #includes intercept
  # print(deg)
  #deg <- 2 ###
  tst2 <- df2
  tst2$df <- tst2$n - deg
  tst2 <- subset(tst2, df > 9) #NOTE: RB 264 uses n :
  #Hoffman uses df
  tst2$lnv <- log(tst2$rmse^2) * tst2$df
  #print(tst2)
  tot <- sum(tst2$lnv) / sum(tst2$df) #total

  tst2$d <- sqrt(tst2$df/2) * ((log(tst2$rmse ^2) - tot))
  tst2$d2 <- tst2$d^2
  ngps <- nrow(tst2) - 1 #used in the summation
  Hind <- sum(tst2$d2)
  res <- data.frame(H = Hind, df = ngps, p = pchisq(Hind, ngps - 1, lower.tail = F))
  # print(ngps)
  # print(range(tst2$d))
  # print(pchisq(Hind, ngps - 1, lower.tail = F))
  return(res)
}