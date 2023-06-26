\name{mse}
\alias{mse}
\title{Calculate the Mean Squared Error Between two Vectors}
\description{
  Merely wanted to use such a function by name; nothing fancy
}
\usage{
  mse(x,y)
}
\arguments{
  \item{x}{vector of values}
  \item{y}{vector of values}
}
\value{
  A scalar: the MSE between \code{x} and \code{y}.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  mse(seq(1:10), seq(10:1))
}
