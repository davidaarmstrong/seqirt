library(tibble)
library(dplyr)
library(pscl)
library(emIRT)
set.seed(43901)
mp <- function(x, inv=TRUE){
  if(inv)x <- 1/x
  x/sum(x)
}
sig <- diag(3)
sig[1,2] <- sig[2,1] <- .25
sig[1,3] <- sig[3,1] <- .25
sig[3,2] <- sig[2,3] <- .25
D <- diag(3)*.25
sig2 <-  D%*%sig%*%D

st <- sample(1:7, 100, prob=mp(1:7), replace=TRUE)
en <- sapply(st, function(x)sample(x:7, 1, prob = mp(x:7, inv=FALSE)))
en <- ifelse(st == 7, 7, en)
df <- NULL
for(i in 1:100){
  y <- matrix(NA, nrow=7, ncol=3)
  y[st[i],] <- MASS::mvrnorm(1,c(0,0,0),sig)
  if(st[i] != en[i]){
  for(t in (st[i]+1):en[i]){
    y[t, ] <- MASS::mvrnorm(1, y[(t-1),], sig2)
  }
  }
  tmp <- tibble(y1 = y[,1], y2=y[,2], y3=y[,3])
  tmp$time <- 1:7
  tmp$obs <- i
  df <- rbind(df, tmp)
}

nitems <- matrix(sample(25:75, 21, replace=TRUE), ncol=3, nrow=7)

mats <- vector(mode="list", length=7)

for(i in 1:7){
  a <- rbind(
    c(runif(nitems[i,1], 0, 1.5), rep(0, sum(nitems[i,2:3]))),
    c(rep(0, nitems[i,1]), runif(nitems[i,2], 0, 1.5), rep(0, nitems[i,3])),
    c(rep(0, sum(nitems[i,1:2])), runif(nitems[i,3], 0, 1.5)))
  a <- rbind(a, -colSums(a) * runif(sum(nitems[i,]), -1.5, 1))
  mats[[i]] <- a
}

df <- df %>% arrange(time, obs)
df$ins <- as.numeric(!is.na(df$y1))
X <- NULL
y1 <- df %>% filter(time == 1) %>% select(-time, -obs, -ins) %>% as.matrix
y1 <- cbind(y1, 1)
p1 <- pnorm(y1 %*% mats[[1]])
X <- cbind(X, array(rbinom(length(c(p1)), 1, p1), dim=dim(p1)))

y2 <- df %>% filter(time == 2) %>% select(-time, -obs, -ins) %>% as.matrix
y2 <- cbind(y2, 1)
p2 <- pnorm(y2 %*% mats[[2]])
X <- cbind(X,  array(rbinom(length(c(p2)), 1, p2), dim=dim(p2)))

y3 <- df %>% filter(time == 3) %>% select(-time, -obs, -ins) %>% as.matrix
y3 <- cbind(y3, 1)
p3 <- pnorm(y3 %*% mats[[3]])
X <- cbind(X, array(rbinom(length(c(p3)), 1, p3), dim=dim(p3)))

y4 <- df %>% filter(time == 4) %>% select(-time, -obs, -ins) %>% as.matrix
y4 <- cbind(y4, 1)
p4 <- pnorm(y4 %*% mats[[4]])
X <- cbind(X, array(rbinom(length(c(p4)), 1, p4), dim=dim(p4)))

y5 <- df %>% filter(time == 5) %>% select(-time, -obs, -ins) %>% as.matrix
y5 <- cbind(y5, 1)
p5 <- pnorm(y5 %*% mats[[5]])
X <- cbind(X, array(rbinom(length(c(p5)), 1, p5), dim=dim(p5)))

y6 <- df %>% filter(time == 6) %>% select(-time, -obs, -ins) %>% as.matrix
y6 <- cbind(y6, 1)
p6 <- pnorm(y6 %*% mats[[6]])
X <- cbind(X, array(rbinom(length(c(p6)), 1, p6), dim=dim(p6)))

y7 <- df %>% filter(time == 7) %>% select(-time, -obs, -ins) %>% as.matrix
y7 <- cbind(y7, 1)
p7 <- pnorm(y7 %*% mats[[7]])
X <- cbind(X, array(rbinom(length(c(p7)), 1, p7), dim=dim(p7)))

