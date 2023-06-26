\name{kolmogorov_smirnov}
\alias{kolmogorov_smirnov}

\title{Kolmogorov-Smirnov Normality Test}
\description{
  Test a sample against the null hypothesis that it comes from a standard Normal distribution.
}
\usage{
  kolmogorov_smirnov(samples)
}
\arguments{
  \item{samples}{vector of values to be tested}
}
\details{
  Wrapper function for corresponding functionality in \code{stats}. Originally, from Kolmogorov (1933).
}
\value{
  A list containing
  \item{name}{name of normality test applied i.e. 'Komogorov Smirnov'}
  \item{statistic}{scalar value of test statistics}
  \item{p.value}{corresponding p-value from test}
}
\source{
  Kolmogorov (1933) via \code{stats}.
}
\examples{
  sample <- rnorm(100)
  kolmogorov_smirnov(sample)
}
