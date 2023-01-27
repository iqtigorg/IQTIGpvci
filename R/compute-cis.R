# confidence interval functions for package IQTIGpvci
# Copyright (C) 2018 IQTIG
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Compute two-sided confidence intervals for binomial parameters.
#'
#' Compute two-sided confidence intervals for a binomial parameter \eqn{p},
#' under the assumption that the underlying distribution of \eqn{o} is binomial
#' with parameters \eqn{p} and \eqn{n}.
#'
#' The confidence interval at confidence level \code{conf.level} is computed by
#' inverting two one-sided tests, each with significance level
#' (1 - \code{conf.level})/2 (see \code{\link{compute_rate_pvalue}} for
#' details).
#'
#' The case \eqn{o=n=0} deserves particular attention.  Formally, the
#' confidence interval is [0, 1] in this case.  However, when \eqn{n=0}, then
#' there is no data available to estimate the binomial parameter, and so the
#' confidence interval is meaningless.  Therefore, the default behaviour is to
#' assign a value of \code{NA_real_} to both boundaries of the confidence
#' interval and to print a warning.  The arguments \code{value_if_n0} and
#' \code{warn_if_n0} can be used to change this behaviour.
#'
#' The arguments \code{o}, \code{n}, \code{conf.level} and \code{midp} may have
#' arbitrary lengths.  If the lengths differ, arguments are recycled according
#' to the usual \R rules (warnings may appear if lenghts are not multiples of
#' each other).  If any of these arguments has length zero, the result has zero
#' rows.  If any of these argument is \code{NA}, the result will also be
#' \code{NA}.
#'
#' @param o the observed number of events (a vector of non-negative integers)
#' @param n total number of cases (a vector of non-negative integers that
#'  satisfies \code{o <= n})
#' @param conf.level the confidence level (a vector of numbers between 0 and 1,
#'   with default value 0.90)
#' @param midp a vector of logicals.  If \code{TRUE} (the default), then compute
#'   mid-p confidence intervals instead of exact confidence intervals.
#' @param value_if_n0 This parameter determines the confidence interval for the
#'   case \code{o = n = 0}.  If \code{value_if_n0} consists of a single
#'   numerical value (default: \code{NA_real_}), both boundaries of the
#'   confidence interval are set to this value.  If \code{value_if_n0} is a
#'   vector of two values, these values are taken as lower and upper boundary
#'   of the confidence interval.  If \code{NULL}, the interval is [0, 1].
#' @param warn_if_n0 logical of length one.  If \code{TRUE} (the default), then
#'   print a warning if the combination \code{o = n = 0} occurs in the input.
#'
#' @return A \code{data.frame} with two columns, the \code{lower} and the
#' \code{upper} bound of the confidence interval.
#'
#' @export
#' @importFrom stats uniroot setNames
#'
#' @examples
#' compute_rate_ci(o = 1, n = 2)
#'
#' @family confidence interval functions
#' @family functions for rate indicators
compute_rate_ci <- function(o, n, conf.level = 0.90, midp = TRUE,
                            value_if_n0 = NA_real_, warn_if_n0 = TRUE) {
  assert_numeric_vector(o)
  assert_numeric_vector(n)
  assert_numeric_vector(conf.level)
  assert_between(conf.level, 0, 1)
  assert_integral(o, 0L, n)
  assert_integral(n)
  assert_logical(midp)
  if (!is.null(value_if_n0)) {
    assert_numeric_vector(value_if_n0)
    assert_length(value_if_n0, lmin = 1L, lmax = 2L)
  }
  assert_logical(warn_if_n0)
  assert_scalar(warn_if_n0)

  # empty input gives empty output
  if (min(length(o), length(n), length(conf.level), length(midp)) == 0L) {
     return(data.frame(lower = numeric(0L), upper = numeric(0L)))
  }

  warn_zero <- FALSE
  if (length(value_if_n0) == 1L) value_if_n0 <- rep(value_if_n0, 2L)
  results <- mapply(function(o, n, conf.level, midp) {
    if (is.na(o) || is.na(n) || is.na(conf.level) || is.na(midp)) {
      lower <- upper <- NA_real_
    } else {
      if (o == 0L) {
        lower <- 0
      } else {
        f <- function(t) {
          compute_rate_pvalue(o = o, n = n, t = t, alternative = "greater",
                              midp = midp) - (1 - conf.level) / 2
        }
        lower <- uniroot(f, interval = c(0, 1),
                         tol = .Machine$double.eps^0.75)$root
      }
      if (o == n) {
        upper <- 1
      } else {
        f <- function(t) {
          compute_rate_pvalue(o = o, n = n, t, alternative = "less",
                              midp = midp) - (1 - conf.level) / 2
        }
        upper <- uniroot(f, interval = c(0, 1),
                         tol = .Machine$double.eps^0.75)$root
      }
      if (n == 0L) {
        if (warn_if_n0) warn_zero <<- TRUE
        if (!is.null(value_if_n0)) {
          lower <- value_if_n0[1L]; upper <- value_if_n0[2L]
        }
      }
    }
    setNames(c(lower, upper), c("lower", "upper"))
  }, o, n, conf.level, midp)
  if (warn_zero) warning("n = 0 encountered.")
  data.frame(t(results))
}

#' Compute two-sided confidence intervals for Poisson parameters.
#'
#' Compute two-sided confidence intervals for the parameter \eqn{\lambda}, under
#' the assumption that the underlying distribution of \eqn{o} is Poisson with
#' parameter \eqn{e\cdot\lambda}{e·\lambda}.
#'
#' The confidence interval at confidence level \code{conf.level} is computed by
#' inverting two one-sided tests, each with significance level
#' (1 - \code{conf.level})/2 (see \code{\link{compute_oe_pvalue}} for details).
#'
#' The case \eqn{e=0} deserves particular attention.  In this degenerated case,
#' the parameter \eqn{\lambda} and the confidence interval become meaningless.
#' Therefore, the default behaviour is to assign a value of \code{NA_real_} to
#' both boundaries of the confidence interval and to print a warning.  The
#' arguments \code{value_if_e0} and \code{warn_if_e0} can be used to change this
#' behaviour.
#'
#' The arguments \code{o}, \code{e}, \code{conf.level} and \code{midp} may have
#' arbitrary lengths.  If the lengths differ, arguments are recycled according
#' to the usual \R rules (warnings may appear if lenghts are not multiples of
#' each other).  If any of these arguments has length zero, the result has
#' zero rows.  If any of these arguments is \code{NA}, the result will also be
#' \code{NA}.
#'
#' @param o the observed number of events (a vector of non-negative integers)
#' @param e the expected number of events (a vector of non-negative numbers)
#' @param conf.level the confidence level (a vector of numbers between 0 and 1,
#'   with default value 0.90)
#' @param midp a vector of logicals.  If \code{TRUE} (the default), then compute
#'   mid-p confidence intervals instead of exact confidence intervals.
#' @param value_if_e0 This parameter determines the confidence interval for the
#'   case \code{e = 0}.  If \code{value_if_e0} consists of a single
#'   numerical value (default: \code{NA_real_}), both boundaries of the
#'   confidence interval are set to this value.  If \code{value_if_e0} is a
#'   vector of two values, these values are taken as lower and upper boundary
#'   of the confidence interval.  If \code{NULL}, the interval is
#'   [0, \eqn{\infty}].
#' @param warn_if_e0 logical of length one.  If \code{TRUE} (the default), then
#'   print a warning if \code{e = 0} occurs in the input.
#'
#' @return A \code{data.frame} with two columns, the \code{lower} and the
#' \code{upper} bound of the confidence interval.
#'
#' @export
#' @importFrom stats uniroot setNames
#'
#' @examples
#' compute_oe_ci(o = 2, e = 1.8)
#'
#' @family confidence interval functions
#' @family functions for o/e indicators
compute_oe_ci <- function(o, e, conf.level = 0.90, midp = TRUE,
                          value_if_e0 = NA_real_, warn_if_e0 = TRUE) {
  assert_numeric_vector(o)
  assert_numeric_vector(e, xmin = 0)
  assert_numeric_vector(conf.level)
  assert_between(conf.level, 0, 1)
  assert_integral(o, 0L)
  assert_logical(midp)
  if (!is.null(value_if_e0)) {
    assert_numeric_vector(value_if_e0)
    assert_length(value_if_e0, lmin = 1L, lmax = 2L)
  }
  assert_logical(warn_if_e0)
  assert_scalar(warn_if_e0)

  # empty input gives empty output
  if (min(length(o), length(e), length(conf.level), length(midp)) == 0L) {
    return(data.frame(lower = numeric(0L), upper = numeric(0L)))
  }

  warn_zero <- FALSE
  if (length(value_if_e0) == 1L) value_if_e0 <- rep(value_if_e0, 2L)
  results <- mapply(function(o, e, conf.level, midp) {
    if (is.na(o) || is.na(e) || is.na(conf.level) || is.na(midp)) {
      lower <- upper <- NA_real_
    } else {
      if (e == 0) {
        if (warn_if_e0) warn_zero <<- TRUE
        if (is.null(value_if_e0)) {
          lower <- 0; upper <- Inf
        } else {
          lower <- value_if_e0[1L]; upper <- value_if_e0[2L]
        }
      } else {
        if (o == 0L) {
          lower <- 0
        } else {
          f <- function(t) {
            compute_oe_pvalue(o = o, e = e, t_smr = t,
                              alternative = "greater",
                              midp = midp) - (1 - conf.level) / 2
          }
          lower <- uniroot(f, interval = c(0, e),
                           tol = .Machine$double.eps^0.75,
                           extendInt = "upX")$root
        }
        f <- function(t) {
          compute_oe_pvalue(o = o, e = e, t_smr = t, alternative = "less",
                            midp = midp) - (1 - conf.level) / 2
        }
        upper <- uniroot(f, interval = c(0, 2 * e),
                         tol = .Machine$double.eps^0.75,
                         extendInt = "downX")$root
      }
    }
    setNames(c(lower, upper), c("lower", "upper"))
  }, o, e, conf.level, midp)
  if (warn_zero) warning("e = 0 encountered.")
  data.frame(t(results))
}
