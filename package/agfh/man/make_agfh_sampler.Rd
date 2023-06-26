\name{make_agfh_sampler}
\alias{make_agfh_sampler}
\title{Maker Function: Agnostic Fay-Herriot Sampler}
\description{
  A maker function that returns a function. The returned function is a sampler for the agnostic Fay-Herriot model.
}
\arguments{
  \item{X}{observed independent data to be analyzed}
  \item{Y}{observed dependent data to be analyzed}
  \item{D}{known precisions of response \code{Y}}
  \item{var_gamma_a}{latent variance prior parameter, \code{rgamma} shape}
  \item{var_gamma_b}{latent variance prior parameter, \code{rgamma} rate}
  \item{S}{vector of starting support values for \eqn{g(\cdot)}}
  \item{kern.a0}{scalar variance parameter of GP kernel}
  \item{kern.a1}{scalar lengthscale parameter of GP kernel}
  \item{kern.fuzz}{scalar noise variance of kernel}
}
\details{
  Creates a Metropolis-within-Gibbs sampler of the agnostic Fay-Herriot model (AGFH).
}
\value{
  Returns a sampler, itself a function of initial parameter values (a list with values for \eqn{\beta, \theta}, the latent variance of \eqn{\theta}, and starting values for \eqn{g(.)}, typically zeros), number of samples, thinning rate, and scale of Metropolis-Hastings jumps for \eqn{\theta} sampling.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  n <- 10
  X <- matrix(1:n, ncol=1)
  Y <- 2*X + rnorm(n, sd=1.1)
  D <- rep(1, n)
  ag <- make_agfh_sampler(X, Y, D)

  params.init <- list(
    beta=1,
    theta=rep(0,n),
    theta.var=1,
    gamma=rep(0,n)
  )
  ag.out <- ag(params.init, 5, 1, 0.1)
}
