% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zipfpssMoments.R
\name{zipfpssMoments}
\alias{zipfpssMoments}
\title{Distribution Moments.}
\usage{
zipfpssMoments(k, alpha, lambda, isTruncated = FALSE,
  tolerance = 10^(-4))
}
\arguments{
\item{k}{Order of the moment to compute.}

\item{alpha}{Value of the \eqn{\alpha} parameter (\eqn{\alpha > k + 1}).}

\item{lambda}{Value of the \eqn{\lambda} parameter (\eqn{\lambda > 0}).}

\item{isTruncated}{Logical; if TRUE, the truncated version of the distribution is returned.}

\item{tolerance}{Tolerance used in the calculations (default = \eqn{10^{-4}}).}
}
\value{
A positive real value corresponding to the k-th moment of the distribution.
}
\description{
General function to compute the k-th moment of the Zipf-PSS distribution for any integer value \eqn{k \geq 1},
when it exists. The k-th moment exists if and only if  \eqn{\alpha > k + 1}.
}
\details{
The k-th moment of the Zipf-PSS distribution is finite for \eqn{\alpha} values
strictly greater than \eqn{k + 1}.
It is computed by calculating the partial sums of the serie, and stopping when two
consecutive partial sums differ less than the \code{tolerance} value.
The value of the last partial sum is returned.
}
\examples{
zipfpssMoments(1, 2.5, 2.3)
zipfpssMoments(1, 2.5, 2.3, TRUE)
}
