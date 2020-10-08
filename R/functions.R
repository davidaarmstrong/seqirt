## TODO
## dynamic IRT screener for dynamic IRT models
##  - binary screener doesn't work as well.



##' K-means Clustering With Some Known Cluster Memberships
##'
##' Estimates a k-means cluster analysis  on a single variable
##' using Lloyd's Method for some known cluster memberships
##' and some unknown-estimated cluster memberships.
##'
##' @param x variable to cluster
##' @param clusters a vector where known cluster memberships are coded
##' as either `1` or `2` and unknown are coded as `NA`.
##' @param centers optional argument providing the starting values
##' for the centers.
##' @param tol tolerance for sum of squared center differences to stop
##' iterating.
##' @param maxit maximum number of iterations.
##' @param ... not used.
##'
##' @return a list containing the centers and cluster memberships at the
##' final iterations.
##' @export
##'
km2 <- function(x, clusters, centers=NULL, tol=1e-05, maxit=15,...){
  init.clust <- clusters
  init.clust <- ifelse(is.na(clusters),
                       sample(1:2, sum(is.na(clusters)), replace=TRUE),
                       clusters)
  if(is.null(centers)){
    cdum <- cbind(as.numeric(init.clust == 1), as.numeric(init.clust == 2))
    cdum <- apply(cdum, 2, function(x)x/sum(x))
    centers <- c(x %*% cdum)
  }
  diff <- 1
  k <- 1
  while(diff > tol & k <= maxit){
    old.centers <- centers
    d <- abs(outer(x, centers, "-"))
    tmp.c <- apply(d, 1, which.min)
    new.clust <- ifelse(is.na(clusters), tmp.c, clusters)
    cdum <- cbind(as.numeric(new.clust == 1), as.numeric(new.clust == 2))
    cdum <- apply(cdum, 2, function(x)x/sum(x))
    centers <- c(x %*% cdum)
    diff <- sum((centers - old.centers)^2)
  }
  return(list(centers=centers, cluster=new.clust))
}

##' Proportional Reduction in Error for Binary IRT Models
##'
##' Calculates the proportional reduction in error for
##' binary IRT models estimated with `binIRT()` from the
##' `emIRT` package. `
##' @param obj Object of class `binIRT`.
##' @param data The data used to estimate `obj`.
##' @param threshold Threshold for PRE values - if all
##' PRE values are greater than the threshold, all
##' variables will be identified as in the same cluster.
##' @param type Should clustering be done on the PREs or
##' should the values greater than a pre-defined threshold
##' be returned?
##'
##' @return A list with the group for each variable in
##'  `cluster` and the cluster number with the high PRE
##'  values in `high.clust`.
##'
##' @export
##' @importFrom stats pnorm kmeans
##' @importFrom magrittr %>%
##' @importFrom dplyr group_by summarise arrange slice select pull tibble
##' @importFrom rlang .data


cluster.pre <- function(obj, data, threshold, type=c("cluster", "absolute")){
  type <- match.arg(type)
  data[which(data == 0, arr.ind=TRUE)] <- NA
  data[which(data == -1, arr.ind=TRUE)] <- 0
  preds <- pnorm(cbind(1, obj$means$x) %*% t(obj$means$beta))
  predx <- apply(preds, 2, function(x)as.numeric(x > .5))
  sames <- predx == data
  pcp <- colMeans(sames, na.rm=TRUE)
  pmc <- colMeans(data, na.rm=TRUE)
  pmc <- ifelse(pmc > .5, pmc, 1-pmc)
  pre <- (pcp-pmc)/(1-pmc)
  if(!all(pre >= threshold) & type == "cluster"){
    ins <- km2(pre, ifelse(pre < 0, 1, NA))
    tmp <- tibble(pre = pre, cluster=ins$cluster)
    high.clust <- tmp %>%
      group_by(.data$cluster) %>%
      summarise(mpre = mean(pre)) %>%
      arrange(.data$mpre) %>%
      slice(2) %>%
      select(.data$cluster) %>%
      pull %>%
      unname(.data)
    list(clusters = ins$cluster, high.clust=high.clust, pre = pre)
  }else{
    list(clusters = as.numeric(pre > threshold), high.clust=1, pre = pre)
  }
}

##' Binary IRT Screener
##'
##' Implements a screening pass of a cross-sectional binary IRT to identify
##' variables to be included for each dimension.
##'
##' @param X Roll call matrix
##' @param Dim Number of dimensions to estimate
##' @param thresh Threshold for the K-means clustering in the final
##' stage.  If all PRE values are above `thresh`, then no clustering will
##' be done and all and all remaining variables will be retained.
##' ##' @param type Should clustering be done on the PREs or
##' should the values greater than a pre-defined threshold
##' be returned?
##' @param control List of control arguments for the `binIRT()` function from the `emIRT` package.
##'
##' @export

bin_screen <- function(X, Dim = 2, thresh=.35, type=c("cluster", "absolute"),  control=NULL){
type = match.arg(type)
reuse <- 1:ncol(X)
keeps <- vector(mode="list", length=Dim)
if(is.null(control)){
  ctrl <- list(threads = 1,
               verbose = FALSE,
               thresh = 1e-6)
}else{
 ctrl <- control
}
for(i in 1:(Dim-1)){
  rcX <- rollcall(X[, reuse])
  cX <- convertRC(rcX)
  starts <- getStarts(cX$n, cX$m, 1)
  priors <- makePriors(cX$n, cX$m, 1)

  tmp <- binIRT(.rc=cX,
                .starts = starts,
                .priors = priors,
                .D = 1,
                .control = ctrl)
  cl <- cluster.pre(tmp, cX$votes, threshold=thresh, type=type)
  keeps[[i]] <- reuse[which(cl$clusters == cl$high.clust)]
  reuse <- reuse[-which(cl$clusters == cl$high.clust)]
}
rcX <- rollcall(X[, reuse])
cX <- convertRC(rcX)
starts <- getStarts(cX$n, cX$m, 1)
priors <- makePriors(cX$n, cX$m, 1)
tmp <- binIRT(.rc=cX,
              .starts = starts,
              .priors = priors,
              .D = 1,
              .control = ctrl)

cl <- cluster.pre(tmp, cX$votes, threshold=thresh, type=type)
keeps[[Dim]] <- reuse[which(cl$clusters == cl$high.clust)]

return(keeps)
}


##' Sequential MIRT Estimation
##'
##' Estimates a multidimensional IRT model by sequentially
##' identifying variables that correspond with each dimension.
##' Once identified, the models  are re-estimated with only the
##' relevant variables.
##'
##' @param .data Data frame with only variables that
##' get passed to the IRT function.
##' @param Dim Number of dimensions to estimate
##' @param thresh Threshold for the K-means clustering in the final
##' stage.  If all PRE values are above `thresh`, then no clustering will
##' be done and all and all remaining variables will be retained.
##' @param ... Arguments to be passed down to the `binIRT`
##' function including `.data`, `.starts` and `.priors`.
##'
##' @return A list with the output from each final `binIRT` model.
##'
##' @import emIRT
##' @importFrom pscl rollcall
##' @export
##'

seqIRT <- function(.data, Dim, thresh=.35, ...){
  X <- as.matrix(.data)
  dots <- list(...)
  if(".control" %in% names(dots)){
    ctrl <- dots$.control
  }else{
    ctrl <- list(threads = 1,
              verbose = FALSE,
              thresh = 1e-6)
  }
  screen <- bin_screen(X, Dim, ctrl)
  mods <- vector(mode="list", length=Dim)
  for(i in 1:Dim){
    rcX <- rollcall(X[, screen$keeps[[i]]])
    cX <- convertRC(rcX)
    starts <- getStarts(cX$n, cX$m, 1)
    priors <- makePriors(cX$n, cX$m, 1)

    mods[[i]] <- binIRT(.rc=cX,
                       .starts = starts,
                       .priors = priors,
                       .D = 1,
                       .control = ctrl)

  }
  return(list(mods = mods, keeps = screen$keeps))
}

##' Sequential Dynamic MIRT Estimation
##'
##' Estimates a multidimensional dynamic IRT model by sequentially
##' identifying variables that correspond with each dimension.
##' Once identified, the models  are re-estimated with only the
##' relevant variables.
##'
##' @param .data Data frame with only variables that
##' get passed to the IRT function.
##' @param Dim Number of dimensions to estimate
##' @param thresh Threshold for the K-means clustering in the final
##' stage.  If all PRE values are above `thresh`, then no clustering will
##' be done and all and all remaining variables will be retained.
##' @param ... Arguments to be passed down to the `binIRT`
##' function including `.data`, `.starts` and `.priors`.
##'
##' @return A list with the output from each final `binIRT` model.
##'
##' @import emIRT
##' @importFrom pscl rollcall
##' @export
##'

seq_dynIRT <- function(.data, Dim, thresh=.35, ...){
  X <- as.matrix(.data)
  dots <- list(...)
  if(".control" %in% names(dots)){
    ctrl <- dots$.control
  }else{
    ctrl <- list(threads = 1,
                 verbose = FALSE,
                 thresh = 1e-6)
  }
  screen <- bin_screen(X, Dim, ctrl)
  mods <- vector(mode="list", length=Dim)
  for(i in 1:Dim){
    rcX <- rollcall(X[, keeps[[i]]])
    cX <- convertRC(rcX)
    starts <- getStarts(cX$n, cX$m, 1)
    priors <- makePriors(cX$n, cX$m, 1)

    mods[[i]] <- binIRT(.rc=cX,
                        .starts = starts,
                        .priors = priors,
                        .D = 1,
                        .control = ctrl)

  }
  return(list(mods = mods, keeps = keeps))
}
