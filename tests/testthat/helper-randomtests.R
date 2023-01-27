# Constants and helper functions for randomized testing
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

sample_size <- 100L             # number of random tests (per test)
imax <- 100L                    # maximal size of integers (100)
tol <- .Machine$double.eps^0.5  # precision when comparing different implementations

# generate random vectors
rintvec <- function(size, max = if (exists("imax")) imax else 100L, min = 0L) {
  sample.int(n = (max - min + 1L), size = size, replace = TRUE) + min - 1L
}
