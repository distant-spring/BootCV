
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Bootstrapping the Cross-Validation Estimate

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

You can install the released version of `BootCV` from
[CRAN](https://cran.r-project.org/web/packages/available_packages_by_name.html#available-packages-G)
with:

``` r
# install.packages("BootCV")
```

And the development version from [GitHub](https://github.com) with:

``` r
# install.packages("devtools")
# devtools::install_github("distant-spring/BootCV")
```

To use `BootCV`, we also need to install the required package `lme4`
from [CRAN](https://cran.r-project.org/web/packages/lme4/index.html).

## Quick Start

This section gives users a general sense of the package with a basic
example about using BootCV to generate cross-validation estimate and the
corresponding confidence interval. We will briefly go over the key
steps, basic operations, inputs and outputs.

First we load `BootCV` (and required package `lme4`)：

``` r
# library(lme4) # load required package
# library(BootCV)
```

We generate the data for illustration, which is a data matrix of
dimension nobs x nvars:

``` r
set.seed(1)
# data generation
n <- 90
p <- 10
x <- matrix(rnorm(n * p), ncol = p)
beta <- rnorm(p)
y <- x %*% beta + rnorm(n)
data <- cbind(y, x)
```

Users need to specify the the summary statistics L for evaluation of
model performance, which is calculated as a function of training data
and testing data. Here we use the mean square error in linear regression
as an example:

``` r
# summary statistics
L <- function(train.data, test.data) {
  y <- train.data[, 1]
  x <- train.data[, -1]
  yt <- test.data[, 1]
  xt <- test.data[, -1]

  fit <- lm(y ~ x)
  beta <- fit$coef

  return(mean((yt - cbind(1, xt) %*% beta)^2))
}
```

After specifying the sample size of training set in cross-validation, we
can use `Boot.CV` to estimate the standard error of cross-validation
estimate of summary statistics L:

``` r
m <- 50 # training set size
# boot = Boot.CV(data, L, m)
# result with adjustment
# result1 = CV.confint(boot,data,L,m,method='Boot.CV',adj=T,print=T)
# result without adjustment
# result2 = CV.confint(boot,data,L,m,method='Boot.CV',adj=F,print=T)
```

`boot` is a list that contains all the relevant information for
calculating the cross-validation estimate and confidence interval. Users
can adjust the number of bootstraps, the number of cross-validations and
the tuning parameter in determining the adjusted sample size of training
set flexibly. Based on `boot`, users can use `CV.confint` to generate
cross-validation estimation and confidence interval. Users can determine
whether to adjust for the reduced training sample size and whether to
print the results with the bool inputs `adj` and `print` and confidence
level with `alpha`.

Sometimes, training the prediction model can be very expensive in terms
of computation, and it may not be feasible to conduct many times
bootstraps. Users may consider use `Boot.Cal` in this case. One thing to
note is that we need to specify the confidence level when using
`Boot.Cal` and keep it same as that used in `CV.confint` (default is
0.95).

``` r
m <- 50 # training set size
# boot = Boot.Cal(data, L, m)
# result with adjustment
# result1 = CV.confint(boot,data,L,m,method='Boot.Cal',adj=T,print=T)
# result without adjustment
# result2 = CV.confint(boot,data,L,m,method='Boot.Cal',adj=F,print=T)
```

## More Examples

The algorithm quickly estimates the standard error of cross-validation
estimate
$$\widehat{Err}^{CV}_m=\frac{1}{B_{CV}}\sum_{b=1}^{B_{CV}} L\left\{D_{test}^b, \hat{\psi}(D_{train}^b)\right\},$$
where data $D=D_{train}^b\cup D_{test}^b$ represents the b-th split in
cross-validation and summary statistics L evaluates the performance of
estimation $\hat{\psi}$ fitted with training set in testing set. Many
cross-validation applications can fit into this very general framework.
This section provides several typical examples about how to construct
the summary statistics L in different settings:

### Example 1 (Precision Medicine)

When the prediction model is a precision medicine strategy which is a
binary classification rule to recommend a treatment to a patient based
on his or her baseline characteristics to maximize the treatment
benefit, we can construct an individualized treatment response (ITR)
score by minimizing a loss function based on a training dataset
$D_{train},$
$$\sum_{X_i\in D_{train}} \left\{Y_i-\gamma'\tilde{Z}_i-(G_i-\pi)\beta'\tilde{Z}_i \right\}^2 ,$$
where $X_i=(Z_i,G_i,Y_i),$ $Y_i$ is the response of interest with a
higher value being desirable, $G_i \in \{0, 1\}$ is a binary treatment
indicator and independent of the baseline covariate $Z_i$ (i.e., the
treatment is randomly assigned to patients in the training set),
$\tilde{Z}_i=(1, Z_i')'$ and $\pi=\Pr(G_i=1).$ Let the resulting
minimizers of $\gamma$ and $\beta$ be $\hat{\gamma}(D_{train})$ and
$\hat{\beta}(D_{train})$ (Tian et al.,
[2014](https://doi.org/10.1080/01621459.2014.951443)). The sample
average of cross-validated treatment effect estimators is our final
cross-validation estimator measuring the performance of the treatment
recommendation system and the summary statistics measuring the
prediction performance
is:$$L\left(D_{test}, \psi\right)=\frac{\sum_{X_i\in D_{test}}Y_iG_iI(\psi'Z_i>0)}{\sum_{X_i\in D_{test}} G_i I(\psi'Z_i>0)}-\frac{\sum_{X_i\in D_{test}}Y_i(1-G_i)I(\psi'Z_i>0)}{\sum_{X_i\in D_{test}}(1-G_i)I(\psi'Z_i>0)}.$$
In this example, $X=(Z,G,Y)$ with $Z$, $G$ and $Y$ being predictors, the
treatment assignment indicator, and a binary outcome, respectively,
$\hat{\psi}(D_{train})=\hat{\beta}(D_{train})$.

### Example 2 (Binary Outcomes)

When the prediction model is constructed via fitting a logistic
regression model, i.e., calculating a regression coefficient vector
$\hat{\beta}(D_{train})$ by maximizing the log-likelihood function

$$ \sum_{(Z_i,Y_i)\in D_{train} } \left[ \beta'\tilde{Z}_i Y_i-\log\left\{1+\exp(\beta'\tilde{Z}_i)  \right\}\right],$$
based on a training dataset $D_{train},$ where $Z_i$ is the predictor,
$\tilde{Z}_i=(1, Z_i')',$ and $Y_i\in \{0, 1\}$ is the binary outcome.
If the dimension of $Z_i$ is high, a lasso-regularization can be used in
estimating $\beta$. In any case, the c-index in a testing set $D_{test}$
can be calculated as
$$\hat{\theta}(D_{train}, D_{test})=\frac{1}{\tilde{n}_{test, 0}\tilde{n}_{test, 1}}\sum_{X_i \in D_{test}(0)}\sum_{X_j\in D_{test1}(1)}I\left(\hat{\beta}(D_{train})'\tilde{Z}_i<\hat{\beta}(D_{train})'\tilde{Z}_j \right),$$
where $\tilde{n}_{test,g}$ is the number of observations in the set
$D_{test}(g)=\{X_i=(Z_i,Y_i)\in D_{test}: Y_i=g\}, g\in \{0, 1\}.$ The
sample average of cross-validated c-index estimators is our final
cross-validation estimator measuring the classification performance of
the logistic regression and the summary statistics measuring the
prediction performance
is:$$L\left(D_{test}, \psi\right)=\frac{1}{\tilde{n}_{test, 0}\tilde{n}_{test, 1}}\sum_{X_i\in D_{test}(0)}\sum_{X_j\in D_{test}(1)} I\left(\psi'\tilde{Z}_i<\psi'\tilde{Z}_j \right).$$
In this example, $X=(Z, Y)$ with $Z$ and $Y$ being the predictor and a
binary outcome of interest, respectively,
$\hat{\psi}(D_{train})=\hat{\beta}(D_{train})$, which can also be fitted
with other methods such as random forest.

### Example 3 Mean Absolute Prediction Error (MAPE)

When the prediction model is constructed by fitting a standard linear
regression model, i.e., calculating a regression coefficient vector
$\hat{\beta}(D_{train})$ by minimizing a $L_2$ loss function
$$ \sum_{X_i\in D_{train} } \left(Y_i- \beta'\widetilde{Z}_i \right)^2,$$
based on a training dataset $D_{train}$, where $Z_i$ is the baseline
covariate for the $i$th patient and $\widetilde{Z}_i=(1, Z_i')'$
including an intercept. When the prediction error in testing set
$D_{test}$ is calculated by mean absolute prediction error (Tian et al.,
[2007](https://doi.org/10.1093/biomet/asm036))$$\hat{\theta}(D_{train}, D_{test})=\frac{1}{n_{test}}\sum_{X_i \in D_{test}}|Y_i-\hat{\beta}(D_{train})'\tilde{Z}_i|,$$
where $n_{test}$ is the number of observations in the testing set. The
sample average of cross-validated mean absolute prediction error
estimators is our final estimator measuring the predictive performance
of the linear model and the summary statistics measuring the prediction
performance
is:$$L\left(D_{test}, \psi\right)=\frac{1}{n_{test}}\sum_{X_i\in D_{test}}|Y_i-\psi'\widetilde{Z}_i|.$$
In this example, $X=(Z, Y)$ with $Z$ and $Y$ being the predictor and
outcome of interest, respectively,
$\hat{\psi}(D_{train})=\hat{\beta}(D_{train})$. If nonlinear prediction
models are considered, then one may construct the prediction model via a
more flexible machine learning algorithm such as random forest or neural
network and the summary statistic is similar.

## References

Bryan Cai, Fabio Pellegrini, Menglan Pang, Carl De Moor, Changyu Shen,
Vivek Charu and Lu Tian, 2023. “Bootstrapping the Cross-Validation
Estimate.” <https://doi.org/10.48550/arXiv.2307.00260>.

Lu Tian, Ash A Alizadeh, Andrew J Gentles and Robert Tibshirani. “A
simple method for estimating interactions between a treatment and a
large number of covariates.” Journal of the American Statistical
Association 109 1517–1532.
<https://doi.org/10.1080/01621459.2014.951443>.

Lu Tian, Tianxi Cai, Els Goetghebeur, L. J. Wei. “Model evaluation based
on the sampling distribution of estimated absolute prediction error.”
Biometrika, Volume 94, Issue 2, June 2007, Pages 297–311.
<https://doi.org/10.1093/biomet/asm036>.
