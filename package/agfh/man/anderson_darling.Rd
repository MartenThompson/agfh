\name{anderson_darling}
\alias{anderson_darling}

\title{Anderson-Darling Normality Test}
\description{
  Test a sample against the null hypothesis that it comes from a standard Normal distribution.
}
\usage{
  anderson_darling(samples)
}
\arguments{
  \item{samples}{vector of values to be tested}
}
\details{
  Wrapper function for corresponding functionality in \code{goftest}. Originally, from Anderson and Darling (1954).
}
\value{
  A list containing
  \item{name}{authors of normality test applied i.e. 'Anderson Darling'}
  \item{statistic}{scalar value of test statistics}
  \item{p.value}{corresponding p-value of the test}
}
\source{
  Anderson and Darling (1954) via \code{goftest}.
}
\examples{
  sample <- rnorm(100)
  anderson_darling(sample)
}
