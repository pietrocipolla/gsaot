% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_functions.R
\name{ishi_homma_fun}
\alias{ishi_homma_fun}
\title{Ishigami-Homma function evaluation}
\usage{
ishi_homma_fun(N, A = 2, B = 1)
}
\arguments{
\item{N}{Number of input samples to generate.}

\item{A}{(default: \code{2}) Numeric, amplitude of the second sine component .}

\item{B}{(default: \code{1}) Numeric, coefficient of the interaction term.}
}
\value{
A list with two elements:
\itemize{
\item \code{x}: a numeric matrix of size \code{N x 8} containing the input samples.
\item \code{y}: a numeric vector of length \code{N} with the corresponding function outputs.
}
}
\description{
Evaluates the Ishigami-Homma function.
Input samples are drawn from a uniform distribution over \eqn{[-\pi, \pi]^3}
}
\details{
The Ishigami-Homma function is defined as:
\deqn{Y = \sin(X_1) + A \cdot \sin^2(X_2) + B \cdot X_3^4 \cdot \sin(X_1)}
where \eqn{X_i \sim \mathcal{U}(-\pi, \pi)}.
}
\examples{
result <- ishi_homma_fun(1000)
head(result$x)
head(result$y)

}
\seealso{
\code{\link{sobol_fun}}, \code{\link{gaussian_fun}}
}
