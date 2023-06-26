\name{cramer_vonmises}
\alias{cramer_vonmises}

\title{Cramer-Von Mises Normality Test}
\description{
  Test a sample against the null hypothesis that it comes from a standard Normal distribution.
}
\usage{
  cramer_vonmises(samples)
}
\arguments{
  \item{samples}{vector of values to be tested}
}
\details{
  Wrapper function for corresponding functionality in \code{goftest}. Originally developed in Cramer (1928), Mises (1931), and Smirnov (1936).
}
\value{
  A list containing
  \item{name}{authors of normality test applied i.e. 'Cramer von Mises'}
  \item{statistic}{scalar value of test statistics}
  \item{p.value}{corresponding p-value of the test}
}
\source{
  Cramer (1928), Mises (1931), and Smirnov (1936) via \code{goftest}.
}
\examples{
  sample <- rnorm(100)
  cramer_vonmises(sample)
}
