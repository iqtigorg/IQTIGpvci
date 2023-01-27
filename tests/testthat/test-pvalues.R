# Tests of p-value functions
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

context("p-value functions")

#### compute_rate_pvalue
test_that("compute_rate_pvalue returns the same values as exactci::binom.exact and respects NAs", {
  for (it in 1:sample_size) {
    # Add information to note in case the test fails
    seed <- sample.int(n = .Machine$integer.max, size = 1)
    set.seed(seed)
    add_seed <- function(note) {
      paste0(note, "\nRandom Seed: ", seed)
    }

    o <- rintvec(rintvec(1))
    # n needs to be larger than o and at least 1
    # for simplicity, we only test the case that min(n) >= max(o)
    # and that length(n) >= length(o)
    n <- pmax(1,
              (if (length(o) > 0) max(o) else 0) + rintvec(rintvec(1)))
    t <- runif(rintvec(1))

    alternative <- sample(x = c("greater", "less", NA), size = rintvec(1),
                          prob = c(0.49, 0.50, 0.01), replace = TRUE)
    midp <- sample(x = c(FALSE, TRUE, NA), size = rintvec(1),
                   prob = c(0.49, 0.50, 0.01), replace = TRUE)

    NAo <- sample(x = c(FALSE, TRUE), size = length(o),
                  prob = c(0.95, 0.05), replace = TRUE)
    NAn <- sample(x = c(FALSE, TRUE), size = length(n),
                  prob = c(0.98, 0.02), replace = TRUE)
    NAt <- sample(x = c(FALSE, TRUE), size = length(t),
                  prob = c(0.99, 0.01), replace = TRUE)

    o[NAo] <- NA
    n[NAn] <- NA
    t[NAt] <- NA

    # apply function
    # suppress warnings due to argument recycling
    lengths <- c(
      length(o), length(n), length(midp), length(t), length(alternative))
    len <- max(lengths)
    if (all(len %% lengths == 0, na.rm = TRUE))
      p <- compute_rate_pvalue(o = o, n = n, t = t,
                            alternative = alternative, midp = midp)
    else
      suppressWarnings(
        p <- compute_rate_pvalue(o = o, n = n, t = t,
                                 alternative = alternative, midp = midp))
        # regexp = "is not a multiple of",
        # info = add_seed("No warning for incongruent lengths of vectors."),
        # all = TRUE)

    # check output
    if (min(lengths) == 0) {
      expect_true(length(p) == 0,
                  add_seed("Empty input, but non-empty output."))
    } else {
      # recycle by hand:
      o <- rep(o, length.out = len)
      n <- rep(n, length.out = len)
      t <- rep(t, length.out = len)
      midp <- rep(midp, length.out = len)
      alternative <- rep(alternative, length.out = len)

      noNA <- !is.na(o) & !is.na(n) & !is.na(midp) &
        !is.na(t) & !is.na(alternative)
      expect_true(all(is.na(p[!noNA])),
                  add_seed("NA input, but output not NA."))
      if (any(noNA)) {
        alt <- compute_rate_pvalue_alt_vect(o = o[noNA], n = n[noNA],
                                            t = t[noNA],
                                            alternative = alternative[noNA],
                                            midp = midp[noNA])
        expect_true(max(abs(alt - p[noNA])) <= tol,
                    add_seed("Output differs from exactci::binom.exact."))
      }
    }
  }
})

test_that("Testing arguments value/warn_if_n0 to compute_rate_pvalue", {
  expect_warning(res <- compute_rate_pvalue(0, 0, 0.5),
                 regexp = "n = 0")
  expect_true(is.na(res))
  expect_silent(res <- compute_rate_pvalue(0, 0, 0.5, warn_if_n0 = FALSE))
  expect_true(is.na(res))
  expect_silent(res <- compute_rate_pvalue(0, 0:10, 0.5, value_if_n0 = 7,
                                           warn_if_n0 = FALSE))
  expect_equal(res, c(7, compute_rate_pvalue(0, 1:10, 0.5)))
  expect_equal(compute_rate_pvalue(0:2, 0:2, 0.5, value_if_n0 = NULL,
                                   warn_if_n0 = FALSE),
               c(1/2, 1/4, 1/8))
  expect_equal(compute_rate_pvalue(0:2, seq(0, 4, 2), 0.5, midp = FALSE,
                                   value_if_n0 = NULL, warn_if_n0 = FALSE),
               c(1, 3/4, 11/16))
})

test_that("Input handling of compute_rate_pvalue", {
  expect_error(compute_rate_pvalue(c("1", "2"), 40, 0.5),
               regexp = "o.*not.*numeric")
  expect_error(compute_rate_pvalue(list(1, 2), 40, 0.5),
               regexp = "o.*not.*numeric")
  expect_error(compute_rate_pvalue(c(0.5, 1:30), 40, 0.5),
               regexp = "o.*not.*integ")
  expect_error(compute_rate_pvalue(10:-1, 40, 0.5),
               regexp = "o.*small.*0")
  expect_error(
    suppressWarnings(compute_rate_pvalue(1:10, 10:4, 0.5)),
               regexp = "o.*larg.*n")
  expect_error(compute_rate_pvalue(0, c("1", "2"), 0.5),
               regexp = "n.*not.*numeric")
  expect_error(compute_rate_pvalue(0, list(1, 2), 0.5),
               regexp = "n.*not.*numeric")
  expect_error(compute_rate_pvalue(0, c(0.5, 1:30), 0.5),
               regexp = "n.*not.*integ")
  expect_error(compute_rate_pvalue(2, 5, c("0.1", "0.2")),
               regexp = "t.*not.*numeric")
  expect_error(compute_rate_pvalue(2, 5, list(0.1, 0.2)),
               regexp = "t.*not.*numeric")
  expect_error(compute_rate_pvalue(2, 5, c(0.2, -0.1, 0.3)),
               regexp = "t.*small.*0")
  expect_error(compute_rate_pvalue(2, 5, c(0.2, 1.1, 0.5)),
               regexp = "t.*large.*1")
  expect_error(compute_rate_pvalue(2, 5, 0.3, alternative = "asdf"),
               regexp = "alternative.*not in")
  expect_error(compute_rate_pvalue(2, 5, 0.3, midp = 0),
               regexp = "midp.*not.*logic")
  expect_error(compute_rate_pvalue(2, 5, 0.3, value_if_n0 = "yes"),
               regexp = "value_if_n0.*not.*numeric")
  expect_error(compute_rate_pvalue(2, 5, 0.3, value_if_n0 = c(1, 2)),
               regexp = "value_if_n0.*not.*length")
  expect_error(compute_rate_pvalue(2, 5, 0.3, warn_if_n0 = 0),
               regexp = "warn_if_n0.*not.*logic")
  expect_error(compute_rate_pvalue(2, 5, 0.3, warn_if_n0 = c(TRUE, TRUE)),
               regexp = "warn_if_n0.*not.*length")
})

#### compute_oe_pvalue

test_that("compute_oe_pvalue returns the same values as exactci::binom.exact and respects NAs", {
  for (it in 1:sample_size) {
    # Add information to note in case the test fails
    seed <- sample.int(n = .Machine$integer.max, size = 1)
    set.seed(seed)
    add_seed <- function(note) {
      paste0(note, "\nRandom Seed: ", seed)
    }

    o <- rintvec(rintvec(1))
    e <- rexp(rintvec(1))
    t <- runif(rintvec(1))
    alternative <- sample(x = c("greater", "less", NA), size = rintvec(1),
                          prob = c(0.49, 0.50, 0.01), replace = TRUE)
    midp <- sample(x = c(FALSE, TRUE, NA), size = rintvec(1),
                   prob = c(0.49, 0.50, 0.01), replace = TRUE)

    NAo <- sample(x = c(FALSE, TRUE), size = length(o),
                  prob = c(0.95, 0.05), replace = TRUE)
    NAe <- sample(x = c(FALSE, TRUE), size = length(e),
                  prob = c(0.98, 0.02), replace = TRUE)
    NAt <- sample(x = c(FALSE, TRUE), size = length(t),
                  prob = c(0.99, 0.01), replace = TRUE)

    o[NAo] <- NA
    e[NAe] <- NA
    t[NAt] <- NA

    # apply function
    # suppress warnings due to argument recycling
    lengths <- c(
      length(o), length(e), length(midp), length(t), length(alternative))
    len <- max(lengths)
    if (all(len %% lengths == 0, na.rm = TRUE))
      p <- compute_oe_pvalue(o = o, e = abs(e), t_smr = t,
                             alternative = alternative, midp = midp)
    else
      suppressWarnings(
        p <- compute_oe_pvalue(o = o, e = abs(e), t_smr = t,
                               alternative = alternative, midp = midp))

    # check output
    if (min(lengths) == 0) {
      expect_true(length(p) == 0,
                  add_seed("Empty input, but non-empty output."))
    } else {
      # recycle by hand:
      o <- rep(o, length.out = len)
      e <- rep(e, length.out = len)
      t <- rep(t, length.out = len)
      midp <- rep(midp, length.out = len)
      alternative <- rep(alternative, length.out = len)

      noNA <- !is.na(o) & !is.na(e) & !is.na(midp) &
        !is.na(t) & !is.na(alternative)
      expect_true(all(is.na(p[!noNA])),
                  add_seed("NA input, but output not NA."))
      if (any(noNA)) {
        alt <- compute_oe_pvalue_alt_vect(o = o[noNA], e = abs(e[noNA]),
                                          t_smr = t[noNA],
                                          alternative = alternative[noNA],
                                          midp = midp[noNA])
        expect_true(max(abs(alt - p[noNA])) <= tol,
                    add_seed("Output differs from exactci::binom.exact."))
      }
    }
  }
})

test_that("Testing arguments value/warn_if_e0 to compute_oe_pvalue", {
  expect_warning(res <- compute_oe_pvalue(c(0, 1), 0, 0.5),
                 regexp = "e = 0")
  expect_true(all(is.na(res)))
  expect_silent(res <- compute_oe_pvalue(c(0, 3), 0, 0.5, warn_if_e0 = FALSE))
  expect_true(all(is.na(res)))
  expect_silent(res <- compute_oe_pvalue(0, 0:10, 0.5, value_if_e0 = 7,
                                           warn_if_e0 = FALSE))
  expect_equal(res, c(7, compute_oe_pvalue(0, 1:10, 0.5)))
  expect_equal(compute_oe_pvalue(0:3, c(0:2, 0), 0.5, value_if_e0 = NULL,
                                   warn_if_e0 = FALSE),
               c(1/2, compute_oe_pvalue(1:2, 1:2, 0.5), 0))
  expect_equal(compute_oe_pvalue(0:3, c(0, 4, 2, 0), 0.5, midp = FALSE,
                                   value_if_e0 = NULL, warn_if_e0 = FALSE),
               c(1, compute_oe_pvalue(1:2, c(4, 2), 0.5, midp = FALSE), 0))
})

test_that("Input handling of compute_oe_pvalue", {
  expect_error(compute_oe_pvalue(c("1", "2"), 40, 0.5),
               regexp = "o.*not.*numeric")
  expect_error(compute_oe_pvalue(list(1, 2), 40, 0.5),
               regexp = "o.*not.*numeric")
  expect_error(compute_oe_pvalue(c(0.5, 1:30), 40, 0.5),
               regexp = "o.*not.*integ")
  expect_error(compute_oe_pvalue(10:-1, 40, 0.5),
               regexp = "o.*small.*0")
  expect_error(compute_oe_pvalue(0, c("1", "2"), 0.5),
               regexp = "e.*not.*numeric")
  expect_error(compute_oe_pvalue(0, list(1, 2), 0.5),
               regexp = "e.*not.*numeric")
  expect_error(compute_oe_pvalue(0, c(5:-5), 0.5),
               regexp = "e.*smaller than 0")
  expect_error(compute_oe_pvalue(2, 5, c("0.1", "0.2")),
               regexp = "t_smr.*not.*numeric")
  expect_error(compute_oe_pvalue(2, 5, list(0.1, 0.2)),
               regexp = "t_smr.*not.*numeric")
  expect_error(compute_oe_pvalue(2, 5, c(0.2, -0.1, 0.3)),
               regexp = "t_smr.*small.*0")
  expect_error(compute_oe_pvalue(2, 5, 0.3, alternative = "asdf"),
               regexp = "alternative.*not in")
  expect_error(compute_oe_pvalue(2, 5, 0.3, midp = 0),
               regexp = "midp.*not.*logic")
  expect_error(compute_oe_pvalue(2, 5, 0.3, value_if_e0 = "yes"),
               regexp = "value_if_e0.*not.*numeric")
  expect_error(compute_oe_pvalue(2, 5, 0.3, value_if_e0 = c(1, 2)),
               regexp = "value_if_e0.*not.*length")
  expect_error(compute_oe_pvalue(2, 5, 0.3, warn_if_e0 = 0),
               regexp = "warn_if_e0.*not.*logic")
  expect_error(compute_oe_pvalue(2, 5, 0.3, warn_if_e0 = c(TRUE, TRUE)),
               regexp = "warn_if_e0.*not.*length")
})
