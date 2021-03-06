% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moezipf.R
\name{moezipf}
\alias{moezipf}
\alias{dmoezipf}
\alias{pmoezipf}
\alias{qmoezipf}
\alias{rmoezipf}
\title{The Marshal-Olkin Extended Zipf Distribution (MOEZipf).}
\usage{
dmoezipf(x, alpha, beta, log = FALSE)

pmoezipf(q, alpha, beta, log.p = FALSE, lower.tail = TRUE)

qmoezipf(p, alpha, beta, log.p = FALSE, lower.tail = TRUE)

rmoezipf(n, alpha, beta)
}
\arguments{
\item{x, q}{Vector of positive integer values.}

\item{alpha}{Value of the \eqn{\alpha} parameter (\eqn{\alpha > 1} ).}

\item{beta}{Value of the \eqn{\beta} parameter (\eqn{\beta > 0} ).}

\item{log, log.p}{Logical; if TRUE, probabilities p are given as log(p).}

\item{lower.tail}{Logical; if TRUE (default), probabilities are \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.}

\item{p}{Vector of probabilities.}

\item{n}{Number of random values to return.}
}
\value{
{
\code{dmoezipf} gives the probability mass function,
\code{pmoezipf} gives the cumulative distribution function,
\code{qmoezipf} gives the quantile function, and
\code{rmoezipf} generates random values from a MOEZipf distribution.}
}
\description{
Probability mass function, cumulative distribution function, quantile function and random number
generation for the MOEZipf distribution with parameters \eqn{\alpha} and \eqn{\beta}. The support of the MOEZipf
distribution are the strictly positive integer numbers large or equal than one.
}
\details{
The \emph{probability mass function} at a positive integer value \eqn{x} of the MOEZipf distribution with
parameters \eqn{\alpha} and \eqn{\beta} is computed as follows:

\deqn{p(x | \alpha, \beta) = \frac{x^{-\alpha} \beta \zeta(\alpha) }{[\zeta(\alpha) - \bar{\beta} \zeta (\alpha, x)] [\zeta (\alpha) - \bar{\beta} \zeta (\alpha, x + 1)]},\, x = 1,2,...,\, \alpha > 1, \beta > 0, }

where \eqn{\zeta(\alpha)} is the Riemann-zeta function at \eqn{\alpha}, \eqn{\zeta(\alpha, x)}
is the Hurtwitz zeta function with arguments \eqn{\alpha} and x, and \eqn{\bar{\beta} = 1 - \beta}.

The \emph{cumulative distribution function}, at a given positive integer value \eqn{x},
 is computed as \eqn{F(x) = 1 - S(x)}, where the survival function \eqn{S(x)} is equal to:
\deqn{S(x) = \frac{\beta\, \zeta(\alpha, x + 1)}{\zeta(\alpha) - \bar{\beta}\,\zeta(\alpha, x + 1)},\, x = 1, 2, .. }

The quantile of the MOEZipf\eqn{(\alpha, \beta)} distribution of a given probability value p
is equal to the quantile of the Zipf\eqn{(\alpha)} distribution at the value:
\deqn{p\prime = \frac{p\,\beta}{1 + p\,(\beta - 1)}}

The quantiles of the Zipf\eqn{(\alpha)} distribution are computed by means of the \emph{tolerance}
package.

To generate random data from a MOEZipf one applies the \emph{quantile} function over \emph{n} values randomly generated
from an Uniform distribution in the interval (0, 1).
}
\examples{
dmoezipf(1:10, 2.5, 1.3)
pmoezipf(1:10, 2.5, 1.3)
qmoezipf(0.56, 2.5, 1.3)
rmoezipf(10, 2.5, 1.3)

}
\references{
{
Casellas, A. (2013) \emph{La distribució Zipf Estesa segons la transformació Marshall-Olkin}. Universitat Politécnica de Catalunya.

Devroye L. (1986) Non-Uniform Random Variate Generation. Springer, New York, NY.

Duarte-López, A., Prat-Pérez, A., & Pérez-Casany, M. (2015). \emph{Using the Marshall-Olkin Extended Zipf Distribution in Graph Generation}. European Conference on Parallel Processing, pp. 493-502, Springer International Publishing.

Pérez-Casany, M. and Casellas, A. (2013) \emph{Marshall-Olkin Extended Zipf Distribution}. arXiv preprint arXiv:1304.4540.

Young, D. S. (2010). \emph{Tolerance: an R package for estimating tolerance intervals}. Journal of Statistical Software, 36(5), 1-39.
}
}
