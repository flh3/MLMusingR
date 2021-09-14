

Htest <- function(newdata, fml, group){

  gps <- newdata[,group]
  gps <- gps[!duplicated(gps)]
  gno <- length(gps)
  #err <- numeric()
  err <- 0
  df <- matrix(NA, nrow = gno, ncol = 3)
  cont <- 0 #error count
  for (i in 1:gno){

    ss <- newdata[newdata[,group] == gps[i], ]
    # print(i)
    tmp <- tryCatch(lm(formula = fml, data = ss),
                    error = function(cond){
                      print("error")
                      #err <- err + 1
                      cont <<- cont + 1 #will not work if just <-
                    }
    )

    if (cont == 0){
      df[i, 1] <- gps[i]
      df[i, 2] <- summary(tmp)$sigma
      df[i, 3] <- nobs(tmp)
    }
    #cat(i,'::')
  }

  cat('cont', cont, '\n')
  #cat('error', err, '\n')
  df2 <- data.frame(df)

}

mtcars2 <- mtcars
mtcars2$mpg[mtcars2$cyl == 4] <- NA #remove outcome
mtcars2$mpg[mtcars2$cyl == 8] <- NA #remove outcome
res <- Htest(mtcars2, 'mpg ~ am', 'cyl')
res

