# Tests of confidence interval functions
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

context("confidence interval functions")

#### compute_rate_ci
test_that("compute_rate_ci returns the same values as exactci::binom.exact, fulfills definition and respects NAs", {
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
    conf.level <- runif(rintvec(1))
    midp <- sample(x = c(FALSE, TRUE, NA), size = rintvec(1),
                   prob = c(0.49, 0.50, 0.01), replace = TRUE)

    NAo <- sample(x = c(FALSE, TRUE), size = length(o),
                  prob = c(0.95, 0.05), replace = TRUE)
    NAn <- sample(x = c(FALSE, TRUE), size = length(n),
                  prob = c(0.98, 0.02), replace = TRUE)
    NAc <- sample(x = c(FALSE, TRUE), size = length(conf.level),
                  prob = c(0.99, 0.01), replace = TRUE)

    o[NAo] <- NA
    n[NAn] <- NA
    conf.level[NAc] <- NA

    # apply function
    # suppress warnings due to argument recycling
    lengths <- c(
      length(o), length(n), length(conf.level), length(midp))
    len <- max(lengths)
    if (all(len %% lengths == 0, na.rm = TRUE))
      ci <- compute_rate_ci(o = o, n = n, conf.level = conf.level, midp = midp)
    else
      suppressWarnings(
        ci <- compute_rate_ci(o = o, n = n, conf.level = conf.level,
                              midp = midp))

    # check output
    if (min(lengths) == 0) {
      expect_true(nrow(ci) == 0, add_seed("Empty input, but non-empty output."))
    } else {
      # recycle by hand:
      o <- rep(o, length.out = len)
      n <- rep(n, length.out = len)
      midp <- rep(midp, length.out = len)
      conf.level <- rep(conf.level, length.out = len)

      noNA <- !is.na(o) & !is.na(n) & !is.na(midp) & !is.na(conf.level)
      expect_true(all(is.na(ci[!noNA, ])),
                  add_seed("NA input, but output not NA."))
      if (any(noNA)) {
        alt <- data.frame(t(
          compute_rate_ci_alt_vect(o = o[noNA], n = n[noNA],
                                   conf.level = conf.level[noNA],
                                   midp = midp[noNA])))

        # compare with exactci (precision of exactci is quite bad):
        expect_true(max(abs(alt - ci[noNA, ])) <= 1e-04,
                    add_seed("Output differs from exactci::binom.exact."))

        # compare with definition:
        ois0 <- (noNA & o == 0)
        oisn <- (noNA & o == n)
        if (any(ois0))
          expect_true(max(abs(ci$lower[ois0])) <= tol,
                      info = add_seed("Lower bound not zero for o == 0."))
        if (any(oisn))
          expect_true(max(abs(ci$upper[oisn] - 1)) <= tol,
                      info = add_seed("Upper bound not one for o == n."))
        regi <- noNA & !ois0  # check lower bound where o is not zero and not NA
        if (any(regi)) {
          pl <- mapply(compute_rate_pvalue,
                       o = o[regi], n = n[regi], t = ci$lower[regi],
                       alternative = "greater", midp = midp[regi])
          expect_true(max(abs(pl - (1 - conf.level[regi]) / 2),
                          na.rm = TRUE) <= tol,
                      info = add_seed("p-value of lower bound not as expected."))
        }
        regi <- noNA & !oisn  # check upper bound where o is not n and not NA
        if (any(regi)) {
          pu <- mapply(compute_rate_pvalue,
                       o = o[regi], n = n[regi], t = ci$upper[regi],
                       alternative = "less", midp = midp[regi])
          expect_true(max(abs(pu - (1 - conf.level[regi]) / 2),
                          na.rm = TRUE) <= tol,
                      info = add_seed("p-value of upper bound not as expected."))
        }
      }
    }
  }
})

test_that("Testing arguments value/warn_if_n0 to compute_rate_ci", {
  expect_warning(res <- compute_rate_ci(0, 0), regexp = "n = 0")
  expect_true(is.na(res$lower), is.na(res$upper))
  expect_silent(res <- compute_rate_ci(0, 0, warn_if_n0 = FALSE))
  expect_true(is.na(res$lower), is.na(res$upper))
  expect_silent(res <- compute_rate_ci(0, 0:10, value_if_n0 = 7,
                                       warn_if_n0 = FALSE))
  expect_equal(res, rbind(data.frame(lower = 7, upper = 7),
                          compute_rate_ci(0, 1:10)))
  ci0 <- data.frame(lower = 0, upper = 1)
  expect_equal(compute_rate_ci(c(0, 1, 1, 0), c(0:2, 0), value_if_n0 = NULL,
                               warn_if_n0 = FALSE),
               rbind(ci0, compute_rate_ci(1, 1:2), ci0))
  expect_equal(compute_rate_ci(c(0, 1, 2, 0), c(0, 4, 2, 0), midp = FALSE,
                               value_if_n0 = NULL, warn_if_n0 = FALSE),
               rbind(ci0, compute_rate_ci(1:2, c(4, 2), midp = FALSE), ci0))
})

test_that("Input handling of compute_rate_ci", {
  expect_error(compute_rate_ci(c("1", "2"), 40, 0.5),
               regexp = "o.*not.*numeric")
  expect_error(compute_rate_ci(list(1, 2), 40, 0.5),
               regexp = "o.*not.*numeric")
  expect_error(compute_rate_ci(c(0.5, 1:30), 40, 0.5),
               regexp = "o.*not.*integ")
  expect_error(compute_rate_ci(10:-1, 40, 0.5),
               regexp = "o.*small.*0")
  expect_error(
    suppressWarnings(compute_rate_ci(1:10, 10:4, 0.5)),
    regexp = "o.*larg.*n")
  expect_error(compute_rate_ci(0, c("1", "2"), 0.5),
               regexp = "n.*not.*numeric")
  expect_error(compute_rate_ci(0, list(1, 2), 0.5),
               regexp = "n.*not.*numeric")
  expect_error(compute_rate_ci(0, c(0.5, 1:30), 0.5),
               regexp = "n.*not.*integ")
  expect_error(compute_rate_ci(2, 5, c("0.1", "0.2")),
               regexp = "conf.level.*not.*numeric")
  expect_error(compute_rate_ci(2, 5, list(0.1, 0.2)),
               regexp = "conf.level.*not.*numeric")
  expect_error(compute_rate_ci(2, 5, c(0.2, -0.1, 0.3)),
               regexp = "conf.level.*not greater than 0")
  expect_error(compute_rate_ci(2, 5, c(0.2, 1.1, 0.5)),
               regexp = "conf.level.*not less than 1")
  expect_error(compute_rate_ci(2, 5, 0.3, midp = 0),
               regexp = "midp.*not.*logic")
  expect_error(compute_rate_ci(2, 5, 0.3, value_if_n0 = "yes"),
               regexp = "value_if_n0.*not.*numeric")
  expect_error(compute_rate_ci(2, 5, 0.3, value_if_n0 = c(1, 2, 3)),
               regexp = "value_if_n0.*length")
  expect_error(compute_rate_ci(2, 5, 0.3, warn_if_n0 = 0),
               regexp = "warn_if_n0.*not.*logic")
  expect_error(compute_rate_ci(2, 5, 0.3, warn_if_n0 = c(TRUE, TRUE)),
               regexp = "warn_if_n0.*not.*length")
})

#### compute_oe_ci
test_that("compute_oe_ci returns the same values as exactci::binom.exact, fulfills definition and respects NAs", {
  for (it in 1:sample_size) {
    # Add information to note in case the test fails
    seed <- sample.int(n = .Machine$integer.max, size = 1)
    set.seed(seed)
    add_seed <- function(note) {
      paste0(note, "\nRandom Seed: ", seed)
    }

    o <- rintvec(rintvec(1))
    e <- rexp(rintvec(1))
    conf.level <- runif(rintvec(1))
    midp <- sample(x = c(FALSE, TRUE, NA), size = rintvec(1),
                   prob = c(0.49, 0.50, 0.01), replace = TRUE)

    NAo <- sample(x = c(FALSE, TRUE), size = length(o),
                  prob = c(0.95, 0.05), replace = TRUE)
    NAe <- sample(x = c(FALSE, TRUE), size = length(e),
                  prob = c(0.98, 0.02), replace = TRUE)
    NAc <- sample(x = c(FALSE, TRUE), size = length(conf.level),
                  prob = c(0.99, 0.01), replace = TRUE)

    o[NAo] <- NA
    e[NAe] <- NA
    conf.level[NAc] <- NA

    # apply function
    # suppress warnings due to argument recycling
    lengths <- c(
      length(o), length(e), length(conf.level), length(midp))
    len <- max(lengths)
    if (all(len %% lengths == 0, na.rm = TRUE))
      ci <- compute_oe_ci(o = o, e = e, conf.level = conf.level, midp = midp)
    else
      suppressWarnings(
        ci <- compute_oe_ci(o = o, e = e, conf.level = conf.level, midp = midp))

    # check output
    if (min(lengths) == 0) {
      expect_true(nrow(ci) == 0, add_seed("Empty input, but non-empty output."))
    } else {
      # recycle by hand:
      o <- rep(o, length.out = len)
      e <- rep(e, length.out = len)
      midp <- rep(midp, length.out = len)
      conf.level <- rep(conf.level, length.out = len)

      noNA <- !is.na(o) & !is.na(e) & !is.na(midp) & !is.na(conf.level)
      expect_true(all(is.na(ci[!noNA, ])),
                  add_seed("NA input, but output not NA."))
      if (any(noNA)) {
        # remove NAs, test finite values.
        o <- o[noNA]
        e <- e[noNA]
        midp <- midp[noNA]
        conf.level <- conf.level[noNA]
        ci <- ci[noNA, ]

        alt <- data.frame(t(compute_oe_ci_alt_vect(o = o, e = e,
                                                   conf.level = conf.level,
                                                   midp = midp)))

        # compare with exactci:
        # The precision of the confidence intervals computed byexactci is quite
        # bad.  Thus, if our confidence interval differs notably from the one by
        # exactci, then we check the definition by computing the pvalue.  The
        # pvalue is computed using exactci, in order not to confirm artefacts
        # in the ci computations that arise from errors in pvalue computations.
        reldiff <- abs(alt$lower - ci$lower) /
          ifelse(alt$lower > 1, alt$lower, 1)
        expect(max(reldiff) < 1e-2,
               add_seed("Lower bound differs from exactci::binom.exact."))
        bad_indices <- reldiff > 1e-6
        if (any(bad_indices)) {
          expect(
            all(abs(mapply(compute_oe_pvalue_alt, o[bad_indices],
                           e[bad_indices], ci$lower[bad_indices],
                           alternative = "greater", midp = midp[bad_indices]) -
                      0.5 * (1 - conf.level[bad_indices])) <
                  0.1 * abs(mapply(compute_oe_pvalue_alt,
                                   o[bad_indices], e[bad_indices],
                                   alt$lower[bad_indices],
                                   alternative = "greater",
                                   midp = midp[bad_indices]) -
                              0.5 * (1 - conf.level[bad_indices]))),
            add_seed("Lower bound differs from exactci::binom.exact."))
        }
        reldiff <- abs(alt$upper - ci$upper) /
          ifelse(alt$upper > 1, alt$upper, 1)
        expect(max(reldiff) < 1e-2,
               add_seed("Upper bound differs from exactci::binom.exact."))
        bad_indices <- reldiff > 1e-6
        if (any(bad_indices)) {
          expect(
            all(abs(mapply(compute_oe_pvalue_alt, o[bad_indices],
                           e[bad_indices], ci$upper[bad_indices],
                           alternative = "less", midp = midp[bad_indices]) -
                      0.5 * (1 - conf.level[bad_indices])) <
                  0.1 * abs(mapply(compute_oe_pvalue_alt, o[bad_indices],
                                   e[bad_indices], alt$upper[bad_indices],
                                   alternative = "less",
                                   midp = midp[bad_indices]) -
                              0.5 * (1 - conf.level[bad_indices]))),
            add_seed("Upper bound differs from exactci::binom.exact."))
        }

        # compare with definition:
        ois0 <- (o == 0)
        if (any(ois0))
          expect_true(max(abs(ci$lower[ois0])) <= tol,
                      info = add_seed("Lower bound not zero for o == 0."))
        if (any(!ois0)) {  # check lower bound where o is not zero (and no NA)
          pl <- mapply(compute_oe_pvalue,
                       o = o[!ois0], e = e[!ois0], t = ci$lower[!ois0],
                       alternative = "greater", midp = midp[!ois0])
          expect_true(max(abs(pl - (1 - conf.level[!ois0]) / 2),
                          na.rm = TRUE) <= tol,
                      info = add_seed("p-value of lower bound not as expected."))
        }
        pu <- mapply(compute_oe_pvalue, o = o, e = e, t = ci$upper,
                     alternative = "less", midp = midp)
        expect_true(max(abs(pu - (1 - conf.level) / 2), na.rm = TRUE) <= tol,
                    info = add_seed("p-value of upper bound not as expected."))
      }
    }
  }
})

test_that("Testing arguments value/warn_if_e0 to compute_oe_ci", {
  expect_warning(res <- compute_oe_ci(c(1, 0), 0), regexp = "e = 0")
  expect_true(all(is.na(res$lower)), all(is.na(res$upper)))
  expect_silent(res <- compute_oe_ci(c(0, 5), 0, warn_if_e0 = FALSE))
  expect_true(all(is.na(res$lower)), all(is.na(res$upper)))
  expect_silent(res <- compute_oe_ci(0, 0:10, value_if_e0 = 7,
                                         warn_if_e0 = FALSE))
  expect_equal(res, rbind(data.frame(lower = 7, upper = 7),
                          compute_oe_ci(0, 1:10)))
  ci0 <- data.frame(lower = 0, upper = Inf)
  expect_equal(compute_oe_ci(0:3, c(0:2, 0), value_if_e0 = NULL,
                                 warn_if_e0 = FALSE),
               rbind(ci0, compute_oe_ci(1:2, 1:2), ci0))
  expect_equal(compute_oe_ci(0:3, c(0, 4, 2, 0), midp = FALSE,
                             value_if_e0 = NULL, warn_if_e0 = FALSE),
               rbind(ci0, compute_oe_ci(1:2, c(4, 2), midp = FALSE), ci0))
})

test_that("Input handling of compute_oe_ci", {
  expect_error(compute_oe_ci(c("1", "2"), 40, 0.5),
               regexp = "o.*not.*numeric")
  expect_error(compute_oe_ci(list(1, 2), 40, 0.5),
               regexp = "o.*not.*numeric")
  expect_error(compute_oe_ci(c(0.5, 1:30), 40, 0.5),
               regexp = "o.*not.*integ")
  expect_error(compute_oe_ci(10:-1, 40, 0.5),
               regexp = "o.*small.*0")
  expect_error(compute_oe_ci(0, c("1", "2"), 0.5),
               regexp = "e.*not.*numeric")
  expect_error(compute_oe_ci(0, list(1, 2), 0.5),
               regexp = "e.*not.*numeric")
  expect_error(compute_oe_ci(0, 10:-1, 0.5),
               regexp = "e.*smaller.*0")
  expect_error(compute_oe_ci(2, 5, c("0.1", "0.2")),
               regexp = "conf.level.*not.*numeric")
  expect_error(compute_oe_ci(2, 5, list(0.1, 0.2)),
               regexp = "conf.level.*not.*numeric")
  expect_error(compute_oe_ci(2, 5, c(0.2, -0.1, 0.3)),
               regexp = "conf.level.*not greater than 0")
  expect_error(compute_oe_ci(2, 5, c(0.2, 1.1, 0.5)),
               regexp = "conf.level.*not less than 1")
  expect_error(compute_oe_ci(2, 5, 0.3, midp = 0),
               regexp = "midp.*not.*logic")
  expect_error(compute_oe_ci(2, 5, 0.3, value_if_e0 = "yes"),
               regexp = "value_if_e0.*not.*numeric")
  expect_error(compute_oe_ci(2, 5, 0.3, value_if_e0 = c(1, 2, NA_real_)),
               regexp = "value_if_e0.*length")
  expect_error(compute_oe_ci(2, 5, 0.3, warn_if_e0 = 0),
               regexp = "warn_if_e0.*not.*logic")
  expect_error(compute_oe_ci(2, 5, 0.3, warn_if_e0 = c(TRUE, TRUE)),
               regexp = "warn_if_e0.*not.*length")
})
