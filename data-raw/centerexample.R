
rm(list=ls())
library(lattice)
library(lmerTest)
library(anscombiser)
library(ggplot2)
library(dplyr)
set.seed(123)
ns <- 20
x1 <- rnorm(ns, m = 1)
y <- 5.5 + .5 * x1 + rnorm(ns)
cor(y, x1)
x1 = x1 + 5

#x2 <- rnorm(ns) #noise
dat <- data.frame(y, x = x1, cluster = 1)
cor(dat[,-3])

mat <- matrix(c(1, .445, .445, 1), 2, 2)
dat2 <- mimic(dat[,-3], correlation = mat)
dat2 <- data.frame(dat2)

dat3 <- mimic(dat2, correlation = mat)
dat3 <- data.frame(dat3)

names(dat2) <- names(dat3) <- c('y', 'x')

dat2$y <- dat2$y + 2
dat2$x <- dat2$x + 2
dat3$y <- dat3$y + 4
dat3$x <- dat3$x + 4


dat2$cluster <- 2
dat3$cluster <- 3

comb <- rbind(dat, dat2, dat3)
comb$x[comb$x < 0] <- .25
comb$cluster <- factor(comb$cluster)

cdata.ex <- comb
