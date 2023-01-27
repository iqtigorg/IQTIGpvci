# Helper functions for checking inputs
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

# Stop with an error message.
#
# Add information about the calling function.
#
# Must not be called interactively, since otherwise sys.call(-2) makes an error.
raise_error <- function(msg) {
  sc <- as.list(sys.call(-2))
  if (length(sc) > 0) {
    calling_function <- sc[[1]]
    fm2 <- paste0(" (in function ", calling_function, ")")
  } else {
    fm2 <- ""
  }
  stop(msg, fm2, ".", call. = FALSE)
}

# Generate a string from the value of a variable
vs <- function(x) {
  value <- deparse(x)
  if (length(value) > 1)
    value <- paste0(value[1], "...")
  value
}

# Check whether argument is a numeric vector
#
# Accept NA, NA_integer_, NA_real_, but not NA_character_, NA_complex_, ...
assert_numeric_vector <- function(x, xmin = NULL, xmax = NULL) {
  if (!(is.numeric(x) || (all(is.na(x)) && is.logical(x)))) {
    raise_error(paste0(deparse(substitute(x)), " = ", deparse(x),
                       " is not a numeric vector"))
  }
  if (! is.null(xmin) && any(x < xmin, na.rm = TRUE)) {
    raise_error(paste0(deparse(substitute(x)), " = ", deparse(x),
                       " is smaller than ", deparse(substitute(xmin))))
  }
  if (! is.null(xmax) && any(x > xmax, na.rm = TRUE)) {
    raise_error(paste0(deparse(substitute(x)), " = ", deparse(x),
                       " is larger than ", deparse(substitute(xmax))))
  }
}

# Check whether x is a vector of integers.
#
# If xmin and/or xmax are given, check whether xmin <= x <= xmax
#
# We do not use is.integer, since this throws an error when the argument is of
# type numeric (e.g. 10 instead of 10L). Instead, we check whether applying
# as.integer changes the value
assert_integral <- function(x, xmin = NULL, xmax = NULL) {
  if (any(!(is.na(x) | (x == as.integer(x))))) {
    raise_error(paste0(deparse(substitute(x)), " = ", vs(x),
                       " is not integral"))
  }
  if (! is.null(xmin) && any(x < xmin, na.rm = TRUE)) {
    raise_error(paste0(deparse(substitute(x)), " = ", vs(x),
                       " is smaller than ", deparse(substitute(xmin))))
  }
  if (! is.null(xmax) && any(x > xmax, na.rm = TRUE)) {
    raise_error(paste0(deparse(substitute(x)), " = ", vs(x),
                       " is larger than ", deparse(substitute(xmax))))
  }
}

# Check whether x is a scalar (i.e. has length one)
assert_scalar <- function(x) {
  if (length(x) != 1)
    raise_error(paste0(deparse(substitute(x)), " = ", vs(x),
                       " does not have length one (length is ", length(x), ")"))
}

# Check bounds on the length of x
assert_length <- function(x, lmin, lmax) {
  if (length(x) > lmax)
    raise_error(paste0(deparse(substitute(x)), " = ", vs(x),
                       " is of length greater than ", lmax,
                       " (length is ", length(x), ")"))
  if (length(x) < lmin)
    raise_error(paste0(deparse(substitute(x)), " = ", vs(x),
                       " is of length less than ", lmin,
                       " (length is ", length(x), ")"))
}

# Check whether argument lies strictly between a and b
assert_between <- function(x, xmin, xmax) {
  if (any(x <= xmin, na.rm = TRUE)) {
    raise_error(paste0(deparse(substitute(x)), " = ", vs(x),
                       " is not greater than ", deparse(substitute(xmin))))
  }
  if (any(x >= xmax, na.rm = TRUE)) {
    raise_error(paste0(deparse(substitute(x)), " = ", vs(x),
                       " is not less than ", deparse(substitute(xmax))))
  }
}

# Check whether argument belongs to a set of choices
assert_choices <- function(x, choices) {
  if (any(!(is.na(x) | x %in% choices))) {
    raise_error(paste0(deparse(substitute(x)), " = ", vs(x),
                       " is not in ", deparse(substitute(choices))))
  }
}

# Check whether argument is logical
assert_logical <- function(x) {
  if (!is.logical(x)) {
    raise_error(paste0(deparse(substitute(x)), " = ", vs(x),
                       " is not logical"))
  }
}
