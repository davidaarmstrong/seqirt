##' Proportional Reduction in Error for Binary IRT Models
##'
##' Calculates the proportional reduction in error for
##' binary IRT models estimated with `binIRT()` from the
##' `emIRT` package. `
##' @param obj Object of class `binIRT`.
##' @param data The data used to estimate `obj`.
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


cluster.pre <- function(obj, data){
  data <- ifelse(data == -1, 0, data)
  preds <- pnorm(cbind(1, obj$means$x) %*% t(obj$means$beta))
  predx <- apply(preds, 2, function(x)as.numeric(x > .5))
  sames <- predx == data
  pcp <- colMeans(sames)
  pmc <- colMeans(data)
  pmc <- ifelse(pmc > .5, pmc, 1-pmc)
  pre <- (pcp-pmc)/(1-pmc)
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
}


##' Sequential MIRT Estimation
##'
##' Estimates a multidimensional IRT model by sequentially
##' identifying variables that correspond with each dimension.
##' Once identified, the models  are re-estimated with only the
##' relevant variables.
##'
##' @param
