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
    y[t, ] <- MASS::mvrnorm(1, y[(t-1),], matrix(0, nrow=3, ncol=3))
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
X1 <- array(rbinom(length(c(p1)), 1, p1), dim=dim(p1))
X <- cbind(X, X1)

y2 <- df %>% filter(time == 2) %>% select(-time, -obs, -ins) %>% as.matrix
y2 <- cbind(y2, 1)
p2 <- pnorm(y2 %*% mats[[2]])
X2 <- array(rbinom(length(c(p2)), 1, p2), dim=dim(p2))
X <- cbind(X, X2)

y3 <- df %>% filter(time == 3) %>% select(-time, -obs, -ins) %>% as.matrix
y3 <- cbind(y3, 1)
p3 <- pnorm(y3 %*% mats[[3]])
X3 <- array(rbinom(length(c(p3)), 1, p3), dim=dim(p3))
X <- cbind(X, X3)

y4 <- df %>% filter(time == 4) %>% select(-time, -obs, -ins) %>% as.matrix
y4 <- cbind(y4, 1)
p4 <- pnorm(y4 %*% mats[[4]])
X4 <- array(rbinom(length(c(p4)), 1, p4), dim=dim(p4))
X <- cbind(X, X4)

y5 <- df %>% filter(time == 5) %>% select(-time, -obs, -ins) %>% as.matrix
y5 <- cbind(y5, 1)
p5 <- pnorm(y5 %*% mats[[5]])
X5 <- array(rbinom(length(c(p5)), 1, p5), dim=dim(p5))
X <- cbind(X, X5)

y6 <- df %>% filter(time == 6) %>% select(-time, -obs, -ins) %>% as.matrix
y6 <- cbind(y6, 1)
p6 <- pnorm(y6 %*% mats[[6]])
X6 <- array(rbinom(length(c(p6)), 1, p6), dim=dim(p6))
X <- cbind(X, X6)

y7 <- df %>% filter(time == 7) %>% select(-time, -obs, -ins) %>% as.matrix
y7 <- cbind(y7, 1)
p7 <- pnorm(y7 %*% mats[[7]])
X7 <- array(rbinom(length(c(p7)), 1, p7), dim=dim(p7))
X <- cbind(X, X7)

nbills <- rowSums(nitems)
term <- rep(1:7, nbills)

xd1 <- X1[,1:nitems[1,1]]
xd2 <- X2[,1:nitems[2,1]]
xd3 <- X3[,1:nitems[3,1]]
xd4 <- X4[,1:nitems[4,1]]
xd5 <- X5[,1:nitems[5,1]]
xd6 <- X6[,1:nitems[6,1]]
xd7 <- X7[,1:nitems[7,1]]

d1x <- cbind(xd1, xd2, xd3, xd4, xd5, xd6, xd7)
d1term <- rep(1:7, nitems[,1])

Xrc <- pscl::rollcall(d1x, yea=1, nay=0, missing=NA)
cX <- convertRC(Xrc)
rc <- cX$votes

rg <- t(apply(rc, 1, function(x)c(min(which(x != 0)), max(which(x != 0)))))
sl <- matrix(d1term[rg[,1]], ncol=1)
el <- matrix(d1term[rg[,2]], ncol=1)

sim.dat <- list(
  rc = cX$votes,
  startlegis = matrix(sl-1, ncol=1),
  endlegis=matrix(el-1, ncol=1),
  bill.session = matrix(d1term-1, ncol=1),
  T = 7
)

myd <- rc %*% t(rc)
which(myd == min(c(myd)), arr.ind=TRUE)


xm0 <- rep(0, nrow(rc))
xm0[88] <- 2
xm0[64] <- -2
xs0 <- ifelse(xm0 == 0, 1, .1)

priors <- list(
  beta.mu = matrix(c(0,0), ncol=1),
  beta.sigma = diag(2),
  x.mu0 = matrix(xm0, ncol=1),
  x.sigma0 = matrix(xs0, ncol=1),
  omega2 = matrix(rep(.1, nrow(rc)), ncol=1))

startleg <- matrix(0, ncol=7, nrow=nrow(rc))
for(i in 1:100){
  startleg[i, sl[i]:el[i]] <- rnorm((el[i]-sl[i] + 1), 0, 1)
}
starts <- list(
  alpha = matrix(0, nrow=ncol(rc), ncol=1),
  beta = matrix(1, nrow=ncol(rc), ncol=1),
  x = startleg
)

lout <- dynIRT(.data = sim.dat,
               .starts = starts,
               .priors = priors)


               .control = list(
                 threads = 1,
                 verbose = TRUE,
                 thresh = 1e-6,
                 maxit=500
               ))


ymat <- matrix(df$y1, ncol=7)







