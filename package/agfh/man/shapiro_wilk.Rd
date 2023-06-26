\name{shapiro_wilk}
\alias{shapiro_wilk}

\title{Shaprio-Wilk Normality Test}
\description{
  Test a sample against the null hypothesis that it comes from a standard Normal distribution.
}
\usage{
  shapiro_wilk(samples)
}
\arguments{
  \item{samples}{vector of values to be tested}
}
\details{
  Wrapper function for corresponding functionality in \code{stats}. Originally, from Shapiro and Wilk (1975).
}
\value{
  A list containing
  \item{name}{authors of normality test applied i.e. 'Shapiro Wilk'}
  \item{statistic}{scalar value of test statistics}
  \item{p.value}{corresponding p-value of the test}
}
\source{
  Shapiro and Wilk (1975) via \code{stats}.
}
\examples{
  sample <- rnorm(100)
  shapiro_wilk(sample)
}
