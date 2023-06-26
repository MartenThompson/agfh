\name{gamma_err_gen}
\alias{gamma_err_gen}
\title{Generate Data with Gamma Sampling Errors}
\description{
  The traditional Fay-Herriot small area model has a Normal latent variable and Normal observed response errors. This method generates data with Normal latent variables and Gamma errors on the response. Note that the sampling errors are transformed so their mean and variance match the the first two moments of the traditional model.
}
\usage{
  gamma_err_gen (M, p, D, lambda, shape, rate)
}
\arguments{
  \item{M}{number of areal units}
  \item{p}{dimension of regressors i.e. \eqn{x \in R^p}}
  \item{D}{vector of precisions for response, length \code{M}}
  \item{lambda}{value of latent variance}
  \item{shape}{shape parameter of Gamma distribution}
  \item{rate}{rate parameter of Gamma distribution}
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
  \item{name}{name of sampling error distribution, including shape and rate parameters}
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  M <- 50
  p <- 3
  D <- rep(0.1, M)
  lamb <- 1/2
  dat <- gamma_err_gen(M, p, D, lamb, 1/2, 10)
}
