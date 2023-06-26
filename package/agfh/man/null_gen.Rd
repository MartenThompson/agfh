\name{null_gen}
\alias{null_gen}
\title{Generate Data with Normal Sampling Errors}
\description{
  The Fay-Herriot small area model has a Normal latent variable and Normal observed response. This generates data according to that specification.
}
\usage{
  null_gen (M, p, D, lambda)
}
\arguments{
  \item{M}{number of areal units}
  \item{p}{dimension of regressors i.e. \eqn{x \in R^p}}
  \item{D}{vector of precisions for response, length \code{M}}
  \item{lambda}{value of latent variance}
}
\value{
  A list containing
  \item{D}{copy of argument `D`}
  \item{beta}{vector of length `p` latent coefficients}
  \item{lambda}{copy of argument `lambda`}
  \item{X}{matrix of independent variables}
  \item{theta}{vector of latent effects}
  \item{Y}{vector of responses}
  \item{err}{vector of sampling errors}
  \item{name}{name of sampling error distribution}
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  M <- 50
  p <- 3
  D <- rep(0.1, M)
  lamb <- 1/2
  dat <- null_gen(M, p, D, lamb)
}
