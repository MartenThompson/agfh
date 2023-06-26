\name{test_u_normal}
\alias{test_u_normal}

\title{Normality Test}
\description{
  Test a sample against the null hypothesis that it comes from a standard Normal distribution with the specified test.
}
\usage{
  test_u_normal(samples, test)
}
\arguments{
  \item{samples}{vector of values to be tested}
  \item{test}{name of test, one of \code{'SW'}, \code{'KS'}, \code{'CM'}, \code{'AD'}}
}
\details{
  Convenience function for consistent syntax in calling \code{shapiro_wilk}, \code{kolmogorov_smirnov}, \code{cramer_vonmises}, and \code{anderson_darling} tests.
}
\value{
  A list containing
  \item{name}{authors of normality test applied}
  \item{statistic}{scalar value of test statistics}
  \item{p.value}{corresponding p-value from test}
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  sample <- rnorm(100)
  test_u_normal(sample, 'SW')
}
