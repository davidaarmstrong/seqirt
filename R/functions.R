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


cluster.pre <- function(obj, data, threshold){
  data <- ifelse(data == -1, 0, data)
  preds <- pnorm(cbind(1, obj$means$x) %*% t(obj$means$beta))
  predx <- apply(preds, 2, function(x)as.numeric(x > .5))
  sames <- predx == data
  pcp <- colMeans(sames)
  pmc <- colMeans(data)
  pmc <- ifelse(pmc > .5, pmc, 1-pmc)
  pre <- (pcp-pmc)/(1-pmc)
  if(!all(pre >= threshold)){
    ins <- kmeans(pre, 2)
    tmp <- tibble(pre = pre, cluster=ins$cluster)
    high.clust <- tmp %>%
      group_by(.data$cluster) %>%
      summarise(mpre = mean(pre)) %>%
      arrange(.data$mpre) %>%
      slice(2) %>%
      select(.data$cluster) %>%
      pull %>%
      unname(.data)
    list(clusters = ins$cluster, high.clust=high.clust)
  }else{
    list(clusters = rep(1, length(pre)), high.clust=1)
  }
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
  reuse <- 1:ncol(X)
  keeps <- vector(mode="list", length=Dim)
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
    cl <- cluster.pre(tmp, cX$votes, threshold=1)
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

  cl <- cluster.pre(tmp, cX$votes, threshold=thresh)
  keeps[[Dim]] <- reuse[which(cl$clusters == cl$high.clust)]

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
  reuse <- 1:ncol(X)
  keeps <- vector(mode="list", length=Dim)
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
    cl <- cluster.pre(tmp, cX$votes, threshold=1)
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

  cl <- cluster.pre(tmp, cX$votes, threshold=thresh)
  keeps[[Dim]] <- reuse[which(cl$clusters == cl$high.clust)]

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
