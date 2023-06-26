\name{RM_theta_var_moment_est}
\alias{rm_theta_var_moment_est}
\alias{RM_theta_var_moment_est}

\title{Moment-Based Estimator of Latent Model Variance}
\description{
  Simple moment-based estimator of the variance of the latent model.
}
\usage{
  RM_theta_var_moment_est(X, Y, D)
}
\arguments{
  \item{X}{observed independent data to be analyzed}
  \item{Y}{observed dependent data to be analyzed}
  \item{D}{known precisions of response \code{Y}}
}
\details{
  Simple moment-based estimator of the variance of the latent model.
}
\value{
  Returns a scalar estimate of variance.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  X <- matrix(1:10, ncol=1)
  Y <- 2*X + rnorm(10, sd=1.1)
  D <- rep(1, 10)
  RM_theta_var_moment_est(X, Y, D)
}
