Xrc <- pscl::rollcall(X, yea=1, nay=0, missing=NA)
cX <- convertRC(Xrc)
rc <- cX$votes

out_oc <- oc(Xrc, dims = 3, polarity=c(60,30, 10))
res <- out_oc$rollcalls
mycat <- apply(res[,grep("norm", colnames(res))], 1, function(x)which.min(abs(x)))



rg <- t(apply(rc, 1, function(x)c(min(which(x != 0)), max(which(x != 0)))))
sl <- matrix(term[rg[,1]], ncol=1)
el <- matrix(term[rg[,2]], ncol=1)

out <- bin_screen(X, Dim=3, thresh=.4, type="absolute")
item_dim <- rep(rep(1:3, 7), c(t(nitems)))

sim.dat <- list(
  rc = cX$votes,
  startlegis = matrix(sl-1, ncol=1),
  endlegis=matrix(el-1, ncol=1),
  bill.session = matrix(term-1, ncol=1),
  T = 7
)

myd <- rc %*% t(rc)
w <- which(myd == min(c(myd)), arr.ind=TRUE)



#screen <- bin_screen(X, 3)


xm0 <- rep(0, nrow(rc))
xm0[w[1,1]] <- 2
xm0[w[1,2]] <- -2
xs0 <- ifelse(xm0 == 0, 1, .1)
priors <- list(
  beta.mu = matrix(c(0,0), ncol=1),
  beta.sigma = diag(2),
  x.mu0 = matrix(xm0, ncol=1),
  x.sigma0 = matrix(xs0, ncol=1),
  omega2 = matrix(rep(.1, nrow(rc)), ncol=1))

startleg <- matrix(0, ncol=7, nrow=nrow(rc))
for(i in 1:nrow(rc)){
  startleg[i, sl[i]:el[i]] <- rnorm((el[i]-sl[i] + 1), 0, 1)
}
starts <- list(
  alpha = matrix(0, nrow=ncol(rc), ncol=1),
  beta = matrix(1, nrow=ncol(rc), ncol=1),
  x = startleg
)

lout <- dynIRT(.data = sim.dat,
               .starts = starts,
               .priors = priors,
               .control = list(
                  threads = 1,
                  verbose = TRUE,
                  thresh = 1e-6,
                  maxit=500
                ))
obj <- lout
data <- rc
data[which(data == 0, arr.ind=TRUE)] <- NA
data[which(data == -1, arr.ind=TRUE)] <- 0
term_preds <- list()
for(i in 1:(ncol(obj$means$x))){
  tmp <- pnorm(cbind(1, obj$means$x[,i]) %*% t(cbind(obj$means$alpha, obj$means$beta)[which(term == i),]))
  term_preds[[i]] <- apply(tmp, 2, function(x)as.numeric(x > .5))
}
term_preds <- do.call(cbind, term_preds)

sames <- term_preds == data
pcp <- colMeans(sames, na.rm=TRUE)
pmc <- colMeans(data, na.rm=TRUE)
pmc <- ifelse(pmc > .5, pmc, 1-pmc)
pre <- (pcp-pmc)/(1-pmc)



lout <- dynIRT(.data = toronto.dat
               .starts = mq_data$cur.mq,
               .priors = mq_data$priors.mq,
               .control = {list(
                 threads = 1,
                 verbose = TRUE,
                 thresh = 1e-6,
                 maxit=500
               )})
