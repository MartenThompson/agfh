\name{RM_theta_eblup}
\alias{rm_theta_eblup}
\alias{RM_theta_eblup}

\title{Traditional EBLUP Estimator of Theta}
\description{
  Traditional EBLUP Estimator of Theta
}
\usage{
  RM_theta_eblup(X, Y, D, theta.var.est)
}
\arguments{
  \item{X}{observed independent data to be analyzed}
  \item{Y}{observed dependent data to be analyzed}
  \item{D}{known precisions of response \code{Y}}
  \item{theta.var.est}{estimate of variance term for latent model; if \code{NA}, will automatically use method-of-moments}
}
\details{
  Traditional EBLUP estimator of latent values theta.
}
\value{
  Returns a vector of estimates of theta.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  X <- matrix(1:10, ncol=1)
  Y <- 2*X + rnorm(10, sd=1.1)
  D <- rep(1, 10)
  th.var.est <- 0.1
  RM_theta_eblup(X, Y, D, th.var.est)

  RM_theta_eblup(X, Y, D)
}
