\name{make_gibbs_sampler}
\alias{make_gibbs_sampler}

\title{Maker Function: Traditional Fay-Herriot Gibbs Sampler}
\description{
  A maker function that returns a function. The returned function is a Gibbs sampler for the traditional Fay-Herriot model.
}
\usage{
  make_gibbs_sampler(X, Y, D, var_gamma_a=1, var_gamma_b=1)
}
\arguments{
  \item{X}{observed independent data to be analyzed}
  \item{Y}{observed dependent data to be analyzed}
  \item{D}{known precisions of response \code{Y}}
  \item{var_gamma_a}{latent variance prior parameter, \code{rgamma} shape}
  \item{var_gamma_b}{latent variance prior parameter, \code{rgamma} rate}
}
\value{
  Returns a Gibbs sampler, itself a function of initial parameter values (a list with values for \eqn{\beta, \theta}, and latent variance of \eqn{\theta}), number of samples, and thinning rate.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  n <- 10
  X <- matrix(1:n, ncol=1)
  Y <- 2*X + rnorm(n, sd=1.1)
  D <- rep(1, n)
  gib <- make_gibbs_sampler(X, Y, D)

  params.init <- list(
    beta=1,
    theta=rep(0,n),
    theta.var=1
  )
  gib.out <- gib(params.init, 5, 1)
}
