\name{resid_likelihood_theta_var_maker}
\alias{resid_likelihood_theta_var_maker}

\title{Maker Function: Residual Likelihood of Latent Variance}
\description{
  A maker function that returns a function. The returned function is the (non-adjusted) residual likelihood of the data for a given (latent) variance, from Yoshimori & Lahiri (2014).
}
\usage{
  resid_likelihood_theta_var_maker(X, Y, D)
}
\arguments{
  \item{X}{observed independent data to be analyzed}
  \item{Y}{observed dependent data to be analyzed}
  \item{D}{known precisions of response \code{Y}}
}
\value{
  Returns the (non-adjusted) residual likelihood as a function of the variance term in the latent model.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  X <- matrix(1:10, ncol=1)
  Y <- 2*X + rnorm(10, sd=1.1)
  D <- rep(1, 10)
  resid.lik <- resid_likelihood_theta_var_maker(X, Y, D)
  resid.lik(0.5)
}
