% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute-cis.R
\name{compute_rate_ci}
\alias{compute_rate_ci}
\title{Compute two-sided confidence intervals for binomial parameters.}
\usage{
compute_rate_ci(o, n, conf.level = 0.9, midp = TRUE,
  value_if_n0 = NA_real_, warn_if_n0 = TRUE)
}
\arguments{
\item{o}{the observed number of events (a vector of non-negative integers)}

\item{n}{total number of cases (a vector of non-negative integers that
satisfies \code{o <= n})}

\item{conf.level}{the confidence level (a vector of numbers between 0 and 1,
with default value 0.90)}

\item{midp}{a vector of logicals.  If \code{TRUE} (the default), then compute
mid-p confidence intervals instead of exact confidence intervals.}

\item{value_if_n0}{This parameter determines the confidence interval for the
case \code{o = n = 0}.  If \code{value_if_n0} consists of a single
numerical value (default: \code{NA_real_}), both boundaries of the
confidence interval are set to this value.  If \code{value_if_n0} is a
vector of two values, these values are taken as lower and upper boundary
of the confidence interval.  If \code{NULL}, the interval is [0, 1].}

\item{warn_if_n0}{logical of length one.  If \code{TRUE} (the default), then
print a warning if the combination \code{o = n = 0} occurs in the input.}
}
\value{
A \code{data.frame} with two columns, the \code{lower} and the
\code{upper} bound of the confidence interval.
}
\description{
Compute two-sided confidence intervals for a binomial parameter \eqn{p},
under the assumption that the underlying distribution of \eqn{o} is binomial
with parameters \eqn{p} and \eqn{n}.
}
\details{
The confidence interval at confidence level \code{conf.level} is computed by
inverting two one-sided tests, each with significance level
(1 - \code{conf.level})/2 (see \code{\link{compute_rate_pvalue}} for
details).

The case \eqn{o=n=0} deserves particular attention.  Formally, the
confidence interval is [0, 1] in this case.  However, when \eqn{n=0}, then
there is no data available to estimate the binomial parameter, and so the
confidence interval is meaningless.  Therefore, the default behaviour is to
assign a value of \code{NA_real_} to both boundaries of the confidence
interval and to print a warning.  The arguments \code{value_if_n0} and
\code{warn_if_n0} can be used to change this behaviour.

The arguments \code{o}, \code{n}, \code{conf.level} and \code{midp} may have
arbitrary lengths.  If the lengths differ, arguments are recycled according
to the usual \R rules (warnings may appear if lenghts are not multiples of
each other).  If any of these arguments has length zero, the result has zero
rows.  If any of these argument is \code{NA}, the result will also be
\code{NA}.
}
\examples{
compute_rate_ci(o = 1, n = 2)

}
\seealso{
Other confidence interval functions: \code{\link{compute_oe_ci}}

Other functions for rate indicators: \code{\link{compute_rate_pvalue}}
}
