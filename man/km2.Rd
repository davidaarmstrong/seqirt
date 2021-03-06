% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{km2}
\alias{km2}
\title{K-means Clustering With Some Known Cluster Memberships}
\usage{
km2(x, clusters, centers = NULL, tol = 1e-05, maxit = 15, ...)
}
\arguments{
\item{x}{variable to cluster}

\item{clusters}{a vector where known cluster memberships are coded
as either \code{1} or \code{2} and unknown are coded as \code{NA}.}

\item{centers}{optional argument providing the starting values
for the centers.}

\item{tol}{tolerance for sum of squared center differences to stop
iterating.}

\item{maxit}{maximum number of iterations.}

\item{...}{not used.}
}
\value{
a list containing the centers and cluster memberships at the
final iterations.
}
\description{
Estimates a k-means cluster analysis  on a single variable
using Lloyd's Method for some known cluster memberships
and some unknown-estimated cluster memberships.
}
