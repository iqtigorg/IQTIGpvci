# Alternative implementations of confidence interval functions
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

compute_rate_pvalue_alt <- function(o, n, t, alternative = "greater",
                                    midp = TRUE) {
  exactci::binom.exact(x = o, n = n, p = t, alternative = alternative,
                       midp = midp)$p.value
}

compute_oe_pvalue_alt <- function(o, e, t_smr, alternative = "greater",
                                  midp = TRUE) {
  exactci::poisson.exact(x = o, T = e * t_smr, alternative = alternative,
                         midp = midp)$p.value
}

compute_rate_pvalue_alt_vect <- Vectorize(compute_rate_pvalue_alt)
compute_oe_pvalue_alt_vect <- Vectorize(compute_oe_pvalue_alt)
