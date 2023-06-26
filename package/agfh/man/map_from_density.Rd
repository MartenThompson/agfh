\name{map_from_density}
\alias{map_from_density}
\title{Calculate the MAP Estimate from Posterior Samples}
\description{
  Find maximum a posteriori estimate using posterior samples
}
\usage{
  map_from_density(param.ts, plot=FALSE)
}
\arguments{
  \item{param.ts}{vector of scalar samples}
  \item{plot}{boolean, plot or not}
  }
\details{
  Finds location of max of density from samples.
}
\value{
  Scalar MAP estimate.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  n.post.samp <- 10
  beta.samp <- rnorm(n.post.samp, 0, 1/2)

  map_from_density(beta.samp)
}
