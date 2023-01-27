# p-value functions for package IQTIGpvci
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

#' Compute (mid)p-values of exact one-sided binomial tests.
#'
#' Compute (mid)p-values of exact one-sided binomial tests.
#' \ifelse{latex}{If \code{alternative = "greater"}, then null hypothesis and
#'  alternative hypothesis are given as
#'  \deqn{H_0: p \le t \quad\mbox{and}\quad H_1: p > t.}
#'  Conversely, if \code{alternative = "less"}, then null hypothesis and the
#'  alternative hypothesis are given as
#'  \deqn{H_0: p \ge t \quad\mbox{and}\quad H_1: p < t.}
#' }{
#'  If \code{alternative = "greater"}, then the null hypothesis is
#'  \deqn{H_0: p \le t}
#'  with alternative hypothesis
#'  \eqn{H_1: p > t.}
#'  Conversely, if \code{alternative = "less"}, then the null hypothesis is
#'  \deqn{H_0: p \ge t}
#'  with alternative hypothesis
#'  \eqn{H_1: p < t.}
#' }
#'
#' When \code{alternative = "greater"}, the p-value is defined as
#' the probability \eqn{\mbox{Pr}(X\ge o)}{Pr(X \ge o)}, where \eqn{X} is a
#' binomial random variable with \eqn{n} trials and success probability \eqn{t}.
#' The corresponding mid-p-value is defined as
#' \deqn{\mbox{Pr}(X > o) + \frac12\mbox{Pr}(X = o).}{Pr(X > o) + ½ Pr(X = o).}
#'
#' The case \eqn{o=n=0} deserves particular attention.  Formally, the
#' p-value is 1, and the mid-p-value is 0.5 in this case.  However, when
#' \eqn{n=0}, then there is no data available to distinguish between the two
#' hypotheses, and so the (mid)p-value is meaningless.  Therefore, the default
#' behaviour is to assign a value of \code{NA_real_} to this case and to print a
#' warning.  The arguments \code{value_if_n0} and \code{warn_if_n0} can be
#' used to change this behaviour.
#'
#' The arguments \code{o}, \code{n}, \code{t} and \code{alternative} may have
#' arbitrary lengths.  If the lengths differ, arguments are recycled according
#' to the usual \R rules (warnings may appear if lenghts are not multiples of
#' each other).  If any of these argument has length zero, the result also has
#' length zero.  If any of these argument is \code{NA}, the result will also be
#' \code{NA}.
#'
#' @param o the observed number of events (a vector of non-negative integers)
#' @param n total number of cases (a vector of non-negative integers that
#'  satisfies \code{o <= n})
#' @param t the threshold value of the null hypothesis (a vector of numbers
#'  between 0 and 1)
#' @param alternative direction of the alternative (a character vector with
#'  elements \code{"greater"} for \eqn{H_1: p>t} or
#'  \code{"less"} for \eqn{H_1: p<t}).  The default value is \code{"greater"}.
#' @param midp a vector of logicals.  If \code{TRUE} (the default), then compute
#'   mid-p values instead of exact p-values.
#' @param value_if_n0 The (mid-)p-value for the case \code{o = n = 0} (a numeric
#'  of length one or \code{NULL}).  The default  value is \code{NA_real_}.
#'  If equal to \code{NULL}, then the p-value is equal to 1 and the mid-p-value
#'  is equal to 0.5.
#' @param warn_if_n0 logical of length one.  If \code{TRUE} (the default), then
#'  print a warning if the combination \code{o = n = 0} occurs in the input.
#'
#' @return A vector of (mid)p-values.
#'
#' @export
#' @importFrom stats pbinom dbinom
#'
#' @examples
#' compute_rate_pvalue(o = 1, n = 2, t = 0.95, alternative = "less")
#'
#' @family p-value functions
#' @family functions for rate indicators
compute_rate_pvalue <- function(o, n, t, alternative = "greater", midp = TRUE,
                                value_if_n0 = NA_real_, warn_if_n0 = TRUE) {
  assert_numeric_vector(o)
  assert_numeric_vector(n)
  assert_numeric_vector(t, 0, 1)
  assert_integral(o, 0L, n)
  assert_integral(n)
  assert_choices(alternative, c("greater", "less"))
  assert_logical(midp)
  if (!is.null(value_if_n0)) {
    assert_numeric_vector(value_if_n0)
    assert_scalar(value_if_n0)
  }
  assert_logical(warn_if_n0)
  assert_scalar(warn_if_n0)

  # empty input gives empty output
  if (min(length(o), length(n), length(t),
          length(alternative), length(midp)) == 0L) {
    return(numeric())
  }

  # the length of the output:
  len <- max(length(o), length(n), length(t), length(alternative), length(midp))
  # recycle all arguments that appear in ifelse statements to output length:
  alternative <- rep(alternative, length.out = len)
  midp <- rep(midp, length.out = len)
  # recycle n to output length, for two reasons:
  # 1. as argument to pbinom/dbinom, to ensure correct recycling there
  # 2. to correctly check whether o <= n
  n <- rep(n, length.out = len)

  # precompute probabilities, in order not to compute them twice in ifelse
  pbg <- pbinom(q = o, size = n, prob = t, lower.tail = FALSE)
  db <- dbinom(x = o, size = n, prob = t)

  result <-
    ifelse(alternative == "greater", {
      pbg + # the probability of x > o
        db * ifelse(midp, 0.5, 1)   # mid-p-correction
    }, ifelse(alternative == "less", {
      1 - pbg + # the probability of x <= o
        ifelse(midp, - 0.5 * db, 0) # mid-p-correction
    },
    NA_real_))
  if (!is.null(value_if_n0)) {
    zeros <- !is.na(n) & (n == 0L) & !is.na(o) & (o == 0L)
    if (any(zeros, na.rm = TRUE)) {
      if (warn_if_n0) warning("n = 0 encountered.")
      result[zeros] <- value_if_n0
    }
  }
  result
}

#' Compute the p-values of exact one-sided Poisson tests.
#'
#' Compute the p-values of exact one-sided Poisson tests.  Suppose that the
#' underlying distribution of \eqn{o} is Poisson with parameter
#' \eqn{e\cdot\lambda}{e·\lambda}.
#' \ifelse{latex}{If \code{alternative = "greater"}, then null hypothesis and
#'  alternative hypothesis are given as
#'  \deqn{H_0: \lambda \le t_{\mbox{SMR}}
#'    \quad\mbox{and}\quad H_1: \lambda > t_{\mbox{SMR}}.}
#'  Conversely, if \code{alternative = "less"}, then null hypothesis and
#'  alternative hypothesis are given as
#'  \deqn{H_0: \lambda \ge t_{\mbox{SMR}}
#'   \quad\mbox{and}\quad H_1: \lambda < t_{\mbox{SMR}}.}
#' }{If \code{alternative = "greater"}, then the null hypothesis is
#'  \deqn{H_0: \lambda \le t_smr}
#'  with alternative hypothesis
#'  \eqn{H_1: \lambda > t_smr.}
#'  Conversely, if \code{alternative = "less"}, then the null hypothesis is
#'  \deqn{H_0: \lambda \ge t_smr}
#'  with alternative hypothesis
#'  \eqn{H_1: \lambda < t_smr.}
#' }
#'
#' When \code{alternative = "greater"}, the p-value is defined as
#' the probability \eqn{\mbox{Pr}(X\ge o)}{Pr(X \ge o)}, where \eqn{X} is a
#' Poisson random variable with expected value \eqn{e\cdot\lambda}{e·\lambda}.
#' The corresponding mid-p-value is defined as
#' \deqn{\mbox{Pr}(X > o) + \frac12\mbox{Pr}(X = o).}{Pr(X > o) + ½ Pr(X = o).}
#'
#' The case \eqn{e=0} deserves particular attention.  In this case, \eqn{o}
#' should also be zero; when \eqn{e=0} and \eqn{o>0}, then the null hypothesis
#' can be immediately rejected, which corresponds to a (mid)p-value of zero.
#' When \eqn{o=e=0}, then, formally, the p-value is 1, and the mid-p-value is
#' 0.5.  However, when \eqn{e=0}, this usually means that there is no data
#' available to distinguish between the two hypotheses, and so the (mid)p-value
#' is meaningless.  Therefore, the default behaviour is to assign a value of
#' \code{NA_real_} to this case and to print a warning.  The arguments
#' \code{value_if_e0} and \code{warn_if_e0} can be used to change this
#' behaviour.
#'
#' The arguments \code{o}, \code{e}, \code{t_smr} and \code{alternative} may have
#' arbitrary lengths.  If the lengths differ, arguments are recycled according
#' to the usual \R rules (warnings may appear if lenghts are not multiples of
#' each other).  If any of these arguments has length zero, the result also has
#' length zero.  If any of these arguments is \code{NA}, the result will also be
#' \code{NA}.
#'
#' @param o the observed number of events (a vector of non-negative integers)
#' @param e the expected number of events (a vector of non-negative numbers)
#' @param t_smr the threshold value of the null hypothesis (a vector of
#'  non-negative numbers)
#' @param alternative direction of the alternative (a character vector with
#'  elements \code{"greater"} for
#'  \eqn{H_1: \lambda > t_{\mbox{SMR}}}{H_1: \lambda > t_smr} or
#'  \code{"less"} for
#'  \eqn{H_1: \lambda < t_{\mbox{SMR}}}{H_1: \lambda < t_smr}).
#'  The default value is \code{"greater"}.
#' @param midp a vector of logicals.  If \code{TRUE} (the default), then compute
#'  mid-p values instead of exact p-values
#' @param value_if_e0 The (mid)pvalue for the case \code{e = 0} (a numeric of
#'  length one or \code{NULL}).  The default  value is \code{NA_real_}.
#'  If equal to \code{NULL}, then the p-value is equal to 1 and the mid-p-value
#'  is equal to 0.5 when \eqn{o = e = 0}; otherwise, when \eqn{o > e = 0}, the
#'  (mid)p-value is 0.
#' @param warn_if_e0 logical of length one.  If \code{TRUE} (the default), then
#'  print a warning if \code{e = 0} occurs in the input.
#'
#' @return A vector of (mid)p-values.
#'
#' @export
#' @importFrom stats ppois dpois
#'
#' @examples
#' compute_oe_pvalue(o = 2, e = 1.8, t_smr = 2.18)
#'
#' @family p-value functions
#' @family functions for o/e indicators
compute_oe_pvalue <- function(o, e, t_smr, alternative = "greater", midp = TRUE,
                              value_if_e0 = NA_real_, warn_if_e0 = TRUE)
{
  assert_numeric_vector(o)
  assert_numeric_vector(e, xmin = 0)
  assert_numeric_vector(t_smr, xmin = 0)
  assert_integral(o, 0L)
  assert_choices(alternative, c("greater", "less"))
  assert_logical(midp)
  if (!is.null(value_if_e0)) {
    assert_numeric_vector(value_if_e0)
    assert_scalar(value_if_e0)
  }
  assert_logical(warn_if_e0)
  assert_scalar(warn_if_e0)

  # empty input gives empty output
  if (min(length(o), length(e), length(t_smr),
          length(alternative), length(midp)) == 0L) {
    return(numeric())
  }

  # the length of the output:
  len <- max(length(o), length(e), length(t_smr),
             length(alternative), length(midp))
  # recycle all arguments that appear in ifelse statements to output length:
  alternative <- rep(alternative, length.out = len)
  midp <- rep(midp, length.out = len)
  # recycle e to output length to ensure correct recycling in ppois/dpois
  e <- rep(e, length.out = len)

  # precompute probabilities, in order not to compute them twice in ifelse
  ppg <- ppois(q = o, lambda = e * t_smr, lower.tail = FALSE)
  dp <- dpois(x = o, lambda = e * t_smr)

  result <-
    ifelse(alternative == "greater", {
      ppg + # the probability of x > o
        dp * ifelse(midp, 0.5, 1)   # mid-p-correction
    }, ifelse(alternative == "less", {
      1 - ppg + # the probability of x <= o
        ifelse(midp, - 0.5 * dp, 0) # mid-p-correction
    },
    NA_real_))
  if (!is.null(value_if_e0)) {
    zeros <- !is.na(e) & (e == 0L)
    if (any(zeros, na.rm = TRUE)) {
      if (warn_if_e0) warning("e = 0 encountered.")
      result[zeros] <- value_if_e0
    }
  }
  result
}
