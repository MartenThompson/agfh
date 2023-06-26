\name{agfh_theta_new_pred}
\alias{agfh_theta_new_pred}

\title{Agnostic Fay-Herriot Hierarchical Bayesian Model Predictions at Latent Level}
\description{
  Find predictions of \eqn{\theta} using posterior samples from the AGFH model
}
\usage{
  agfh_theta_new_pred(X_new, beta_samples, theta_var_samples)
}
\arguments{
  \item{X_new}{single new independent data to be analyzed}
  \item{beta_samples}{posterior samples of latent regression parameter}
  \item{theta_var_samples}{posterior samples of latent variance parameter}
}
\details{
  \code{X_new} should be \eqn{1} \code{x} \eqn{p} shaped.

  \code{beta_samples} and \code{theta_var_samples} should contain the same number of samples (columns for the former, length of the latter).
}
\value{
  Vector containing n samples-many estimates of \eqn{\theta} at \code{X_new}.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  p <- 3
  n.post.samp <- 10
  X.new <- matrix(rep(1,p), nrow=1)
  beta.samp <- matrix(rnorm(n.post.samp*p, mean=2, sd=0.1), ncol=n.post.samp)
  thvar.samp <- runif(n.post.samp, 0.1, 1)

  th.preds <- agfh_theta_new_pred(X.new, beta.samp, thvar.samp)
}
