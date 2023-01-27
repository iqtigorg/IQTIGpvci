# IQTIGpvci: R Functions for Hospital Profiling

The R package `IQTIGpvci` provides statistical helper functions for hospital
profiling based on statistical hypothesis testing.  The package is part of the
scientific work in statistical methodology at the Federal Institute for Quality
Assurance and Transparency in Healthcare (IQTIG).  Moreover, the package
increases transparency about the methods used in the plan. QI procedure for
German hospital profiling.  In particular, it provides implementations of the
functions `compute_rate_pvalue` and `compute_oe_pvalue` mentioned in the
[plan. QI directive](https://www.g-ba.de/informationen/richtlinien/91/), which
perform computations of exact binomial and Poisson tests using mid-p-values.

For official details about the plan. QI procedure (in German), see the
[plan. QI directive](https://www.g-ba.de/informationen/richtlinien/91/) and the
[plan. QI final development report](https://iqtig.org/downloads/berichte/2016/IQTIG_Planungsrelevante-Qualitaetsindikatoren_Abschlussbericht.pdf).

The package provides four functions to compute (mid)p-values and confidence
intervals for rate indicators and for o/e (smr) indicators:
```R
compute_rate_pvalue
compute_rate_ci
compute_oe_pvalue
compute_oe_ci
```

See the documentation of these functions and the package vignette
(`vignette("IQTIGpvci")`) for more information.

### License and disclaimer

The directive that governs the plan. QI project may be amended and modified each year.  For current information on the project, we refer to the websites of the [Federal Joint Committee](https://www.g-ba.de) and the [IQTIG](https://www.iqtig.org).  The examples in this document are based on the directive that was published in 2016 and that applies to data collected in 2017.
Note that p-values computed using the package may deviate slightly from the
p-values used in official reports due to, for example, numerical reasons.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.
