\name{theta_var_est_grid}
\alias{theta_var_est_grid}

\title{Basic Grid Optimizer for Likelihood}
\description{
  A basic grid search optimizer. Here, used to estimate the variance in the latent model by maximum likelihood.
}
\usage{
  theta_var_est_grid(likelihood_theta_var)
}
\arguments{
  \item{likelihood_theta_var}{some flavor of likelihood function in terms of latent variance}
}
\details{
  \code{likelihood_theta_var} may be created using \code{adj_resid_likelihood_theta_var_maker} or similar.

  We recommended implementing a more robust optimizer.
}
\value{
  The scalar value that optimizes \code{likelihood_theta_var}, or an error if this value is on the search boundary \eqn{[10^{-6}, 10^2]}.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  X <- matrix(1:10, ncol=1)
  Y <- 2*X + rnorm(10, sd=1.1)
  D <- rep(1, 10)
  adj.lik <- adj_resid_likelihood_theta_var_maker(X, Y, D)
  theta_var_est_grid(adj.lik)
}
