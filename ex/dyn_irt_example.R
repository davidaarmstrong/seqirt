Xrc <- pscl::rollcall(X, yea=1, nay=0, missing=NA)
cX <- convertRC(Xrc)
rc <- cX$votes

rg <- t(apply(rc, 1, function(x)c(min(which(x != 0)), max(which(x != 0)))))
sl <- matrix(term[rg[,1]], ncol=1)
el <- matrix(term[rg[,2]], ncol=1)

sim.dat <- list(
  rc = cX$votes,
  startlegis = sl-1,
  endlegis=el-1,
  bill.session = matrix(term-1, ncol=1),
  T = 7
)

screen <- bin_screen(X, 3)


xm0 <- rep(0, nrow(rc))
xm0[47] <- 2
xm0[10] <- -2
xs0 <- ifelse(xm0 == 0, 1, .1)
toronto.priors <- list(
  beta.mu = c(0,0),
  beta.sigma = diag(2),
  x.mu0 = matrix(xm0, ncol=1),
  x.sigma0 = matrix(xs0, ncol=1),
  omega2 = matrix(rep(.1, nrow(rc)), ncol=1))

#TODO
1. make starting values (using MDS?)

make_coef <- function(x, y){
  tmp.y <- y
  tmp.y <- case_when(y == 1 ~ 1L, y==-1 ~ 0L, y==0 ~ NA_integer_)
  tdf <- data.frame(y=tmp.y, x=x)
  coef(glm(y ~ x, data=tdf, family=binomial(link="probit")))
}

coefs <- apply(rc, 2, function(x)make_coef(out$legislators[,7], x))
coefs <- t(coefs)
coefs[,1] <- ifelse(abs(coefs[,1]) > 2, sign(coefs[,1])*2, coefs[,1])
coefs[,2] <- ifelse(abs(coefs[,2]) > 2, sign(coefs[,2])*2, coefs[,2])

lout <- dynIRT(.data = toronto.dat
               .starts = mq_data$cur.mq,
               .priors = mq_data$priors.mq,
               .control = {list(
                 threads = 1,
                 verbose = TRUE,
                 thresh = 1e-6,
                 maxit=500
               )})
