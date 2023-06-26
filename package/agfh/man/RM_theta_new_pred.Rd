\name{RM_theta_new_pred}
\alias{rm_theta_new_pred}
\alias{RM_theta_new_pred}

\title{Traditional EBLUP Estimator of Theta for new \code{X} values}
\description{
  Traditional EBLUP Estimator of Theta for new \code{X} values
}
\usage{
  RM_theta_new_pred(X.new, beta.est)
}
\arguments{
  \item{X.new}{new independent data to be analyzed}
  \item{beta.est}{estimate of regression term for latent model}
}
\details{
  Simply \code{X}'\code{beta.est}
}
\value{
  Returns a vector of estimates of theta.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  X <- matrix(1:10, ncol=1)
  b <- 1
  RM_theta_new_pred(X, b)
}
