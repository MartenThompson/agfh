\name{RM_beta_eblue}
\alias{rm_beta_eblue}
\alias{RM_beta_eblue}

\title{Traditional EBLUE Estimator of Beta}
\description{
  Traditional EBLUE Estimator of Beta
}
\usage{
  RM_beta_eblue(X, Y, D, theta_var_est)
}
\arguments{
  \item{X}{observed independent data to be analyzed}
  \item{Y}{observed dependent data to be analyzed}
  \item{D}{known precisions of response \code{Y}}
  \item{theta_var_est}{estimate of variance term for latent model}
}
\details{
  Traditional EBLUE estimator of beta.
}
\value{
  Returns a vector estimate of beta.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  X <- matrix(1:10, ncol=1)
  Y <- 2*X + rnorm(10, sd=1.1)
  D <- rep(1, 10)
  th.var.est <- 0.1
  RM_beta_eblue(X, Y, D, th.var.est)
}
