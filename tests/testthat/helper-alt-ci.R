# Alternative implementations of p-value functions
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

########## alternative implementation, relying on package exactci

compute_rate_ci_alt <- function(o, n, conf.level = 0.90, midp = TRUE) {
  setNames(
    as.vector(
      exactci::binom.exact(x = o, n = n,
                           p = 0.5,  # need to specify dummy value for p
                           tsmethod = "central",
                           conf.level = conf.level,
                           alternative = "two.sided", midp = midp)$conf.int),
    c("lower", "upper"))
}

compute_oe_ci_alt <- function(o, e, conf.level = 0.90, midp = TRUE) {
  setNames(
    as.vector(
      exactci::poisson.exact(x = o, T = e,
                             tsmethod = "central",
                             conf.level = conf.level,
                             alternative = "two.sided", midp = midp)$conf.int),
    c("lower", "upper"))
}

compute_rate_ci_alt_vect <- Vectorize(compute_rate_ci_alt)
compute_oe_ci_alt_vect <- Vectorize(compute_oe_ci_alt)
