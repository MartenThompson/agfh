
<!-- README.md is generated from README.Rmd.
    Re-render README.Rmd before release with devtools::build_readme()-->

# agfh

<!-- badges: start -->
<!--[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/agfh)](https://cran.r-project.org/package=agfh) -->
<!-- badges: end -->

## Overview

Small area statistics concerns estimation techniques for sub-populations
when direct estimation would be unreliable. The `agfh` package
implements the Agnostic Fay-Herriot model (AGFH), an extension of the
traditional small area model. In place of normal sampling errors, the
sampling error distribution is modified by a Gaussian process to
accommodate a broader class of distributions.

This flexibility is most useful in the presence of bounded, bi-modal, or
heavily skewed sampling errors. Practitioners should consider the AGFH
model when they have evidence of such departures from the traditional
methods

## Installation

Install the official version from CRAN:

``` r
install.packages('agfh')
```

Next, consult the accompanying paper for a thorough background (under
review), or the vignette within this package for an end-to-end
illustration of the package.

## Getting Started

The AGFH model is implemented as a Metropolis-within-Gibbs sampler; use
`make_agfh_sampler()` to instantiate a sampler. Doing so requires
supplying the observed response (as an
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")-vector
of univariate values), accompanying covariates (as an
![n \\times p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20p "n \times p")
matrix of values), and sampling error precision (again an
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")-vector
of univariate values). Additionally, prior hyperparameters can be
supplied.

`make_agfh_sampler()` creates a sampler function; calling it will
produce MCMC samples targeting the posterior. It requires starting
values for the Gibbs components as well as the desired number of steps
and thinning rate. Note, `n.mcmc=100` and `n.thin=10` would make 1000
MCMC steps and return every tenth.

The sampler returns a list of relevant samples and summary values.
Typically, the contents of `param.samples.list` are most interesting;
these are the posterior samples from the AGFH model. The convenience
method `map_from_density()` may be used to get a maximum a posteriori
point estimate.

Parallel analysis with the traditional Fay-Herriot model is also
possible with `agfh`, as detailed in the vignette. In particular,
`make_gibbs_sampler()` returns a Gibbs sampler of the traditional model
that can be used in the same manner as `make_agfh_sampler()`.
