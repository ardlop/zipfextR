% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getInitialValues.R
\name{getInitialValues}
\alias{getInitialValues}
\title{Calculates initial values for the parameters of the models.}
\usage{
getInitialValues(data, model = "zipf")
}
\arguments{
\item{data}{Matrix of count data.}

\item{model}{Specify the model that requests the initial values (default='zipf').}
}
\value{
Returns the initial values of the parameters for a given distribution.
}
\description{
The selection of appropiate initial values to compute the maximum likelihood estimations
reduces the number of iterations which in turn, reduces the computation time.
The initial values proposed by this function are computed using the first two empirical
frequencies.
}
\details{
The argument \code{data} is a two column matrix with the first column containing the observations and
the second column containing their frequencies. The argument \code{model} refers to the selected model of those
implemented in the package. The possible values are: \emph{zipf}, \emph{moezipf}, \emph{zipfpe},
\emph{zipfpss} or its zero truncated version \emph{zt_zipfpss}. By default, the selected model is the Zipf one.

For the MOEZipf, the Zipf-PE and the zero truncated Zipf-PSS models that contain the Zipf model as
a particular case, the \eqn{\beta} value will correspond to the one of the Zipf model (i.e. \eqn{\beta = 1} for the MOEZipf,
\eqn{\beta = 0} for the Zipf-PE and \eqn{\lambda = 0} for the zero truncated Zipf-PSS model) and the initial value for \eqn{\alpha}
is set to be equal to:
\deqn{\alpha_0 = log_2 \big (\frac{f_r(1)}{f_r(2)} \big),}
where \eqn{f_r(1)} and \eqn{f_r(2)} are the empirical relative frequencies of one and two.
This value is obtained equating the two empirical probabilities to their theoritical ones.

For the case of the Zipf-PSS the proposed initial values are obtained equating the empirical probability of zero
to the theoretical one which gives:
\deqn{\lambda_0 = -log(f_r(0)),}
where \eqn{f_r(0)} is the empirical relative frequency of zero. The initial value of \eqn{\alpha} is obtained
equating the ratio of the theoretical probabilities at zero and one to the empirical ones. This gives place to:
\deqn{\alpha_0 = \zeta^{-1}(\lambda_0 * f_r(0)/f_r(1)),}
where \eqn{f_r(0)} and \eqn{f_r(1)} are the empirical relative frequencies associated to the values 0 and 1 respectively.
The inverse of the Riemman Zeta function is obtained using the \code{optim} routine.
}
\examples{
data <- rmoezipf(100, 2.5, 1.3)
data <- as.data.frame(table(data))
data[,1] <- as.numeric(levels(data[,1])[data[,1]])
initials <- getInitialValues(data, model='zipf')
}
\references{
{ Güney, Y., Tuaç, Y., & Arslan, O. (2016). Marshall–Olkin distribution: parameter estimation and
application to cancer data. Journal of Applied Statistics, 1-13.}
}
