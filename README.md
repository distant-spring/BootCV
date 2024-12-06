
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BootCV: Bootstrapping the Cross-Validation Estimate

<!-- badges: start -->
<!-- badges: end -->

We provide fast bootstrap methods to quantify the uncertainty of
cross-validation estimate, which quickly estimates the standard error of
the cross-validation estimate and produces valid confidence intervals
for a population parameter measuring average model performance. The
method overcomes the computational challenge inherent in bootstrapping
the cross-validation estimate by estimating the variance component
within a random effects model and is flexible for evaluating the
performance of general predictive models. Details can be found in Cai et
al. ([2023](https://doi.org/10.48550/arXiv.2307.00260)).

## Installation

You can install the development version of `BootCV` from
[GitHub](https://github.com/distant-spring/BootCV) with:

``` r
# install.packages("devtools")
devtools::install_github("distant-spring/BootCV")
```

To use `BootCV`, you also need to install the required package `lme4`
from [CRAN](https://cran.r-project.org/web/packages/lme4/index.html).

## Help Document

The [help
document](https://github.com/distant-spring/BootCV/blob/main/BootCV%20Document.pdf)
of `BootCV` provides detailed guidance and examples on how to use this
package.

## References

Bryan Cai, Fabio Pellegrini, Menglan Pang, Carl De Moor, Changyu Shen,
Vivek Charu and Lu Tian, 2023. “Bootstrapping the Cross-Validation
Estimate.” <https://doi.org/10.48550/arXiv.2307.00260>.
