## Packages
library(lme4)

#' Bootstrap Calibration
#'
#' Boot.Cal implements the proposed algorithm (algorithm 3 in reference paper) with small B.bt and B.cv and returns the bootstrap standard error
#' of cross-validation estimator, cutoff for confidence interval and inflation factor.
#'
#' @details
#' The algorithm quickly estimates the standard error of cross-validation estimate
#' \deqn{\widehat{Err}^{CV}_m=\frac{1}{B_{CV}}\sum_{b=1}^{B_{CV}} L\left\{D_{test}^b, \hat{\psi}(D_{train}^b)\right\},}{}
#' where data \eqn{D=D_{train}^b\cup D_{test}^b} represents the b-th split in cross-validation and
#' summary statistics L evaluates the performance of estimation \eqn{\hat{\psi}} fitted with training set in testing set.
#' The adjusted sample size of training set \eqn{m_{adj}} is determined by minimizing
#' \eqn{(\frac{m_{adj}}{m/0.632}-1)^2+\lambda_0 (\frac{n-m}{n-m_{adj}}-1)^2.}
#' @param data data matrix of dimension nobs x nvars; each row is an observation vector.
#' @param L summary statistics with two inputs, train.data and test.data,
#' which are in the same form as data; it evaluates the performance of estimation fitted with train.data in test.data.
#' @param m training set size.
#' @param B.bt the number of bootstraps (default is 20).
#' @param B.cv the number of cross-validations (default is 50).
#' @param Brm.bt the number of bootstraps for variance estimator (default is 1000)
#' @param alpha 1-confidence level (default is 0.05).
#' @param lambda0 tuning parameter in determining the adjusted sample size of training set (default is 0.368).
#'
#' @return
#' \item{sigma}{the bootstrap standard error of cross-validation estimator.}
#' \item{cutoff}{the quantile in computation of confidence interval.}
#' \item{factor}{inflation factor to account for the reduced sample size in bootstrapped training set.}
#' @export
#'
#' @examples
#' # library(lme4) # load required package
#' library(BootCV)
#'
#' set.seed(1)
#' # data generation
#' n <- 90
#' p <- 10
#' x <- matrix(rnorm(n * p), ncol = p)
#' beta <- rnorm(p)
#' y <- x %*% beta + rnorm(n)
#' data <- cbind(y, x)
#'
#' # summary statistics
#' L <- function(train.data, test.data) {
#'   y <- train.data[, 1]
#'   x <- train.data[, -1]
#'   yt <- test.data[, 1]
#'   xt <- test.data[, -1]
#'
#'   fit <- lm(y ~ x)
#'   beta <- fit$coef
#'
#'   return(mean((yt - cbind(1, xt) %*% beta)^2))
#' }
#'
#' m <- 50 # training set size
#' boot <- Boot.Cal(data, L, m)
#' result1 <- CV.confint(boot, data, L, m, method = "Boot.Cal", adj = T, print = T)
#' result2 <- CV.confint(boot, data, L, m, method = "Boot.Cal", adj = F, print = T)
Boot.Cal <- function(data, L, m, B.bt = 20, B.cv = 50, Brm.bt = 1000, alpha = 0.05, lambda0 = 0.368) {
  # data: data points in n x p matrix
  # L: summary statistics with two parameters, training set and testing set
  # m: training set size
  # B.bt: the number of bootstraps (default is 20)
  # B.cv: the number of cross-validations (default is 50)
  # Brm.bt: the number of bootstraps for variance estimator (default is 1000)
  # alpha: 1-confidence level (default is 0.05)
  # lambda0: tuning parameter in determing the adjusted sample size of training set (default is 0.368)
  n <- nrow(data)
  weight.mat <- rmultinom(B.bt, size = n, prob = rep(1 / n, n)) # frequencies of bootstrapped dataset

  if (n - m < 100) {
    m.test <- m:(n - round(0.05 * n))
  } else {
    m.test <- round(seq(m, n - round(0.05 * n), length.out = 100))
  }
  m.adj <- m.test[which.min((m / (m.test * 0.632) - 1)^2 + lambda0 * ((n - m) / (n - m.test) - 1)^2)] # adjusted sample size of training set
  factor <- sqrt(n - m.adj + m.adj * 0.632) / sqrt(n)

  result.btcv <- matrix(NA, B.bt, B.cv)
  for (b.bt in 1:B.bt) {
    weightt <- weight.mat[, b.bt]
    for (b.cv in 1:B.cv) {
      set.seed(10^6 * b.bt + 10^3 * b.cv)
      id <- sample(n, m.adj, F)
      # split the data
      train.data <- data[id, ]
      test.data <- data[-id, ]
      train.weight <- weightt[id]
      test.weight <- weightt[-id]
      # construct the weighted dataset
      weighted.train.data <- train.data[rep(1:nrow(train.data), train.weight), ]
      weighted.test.data <- test.data[rep(1:nrow(test.data), test.weight), ]

      if (sum(test.weight) >= 1) {
        result.btcv[b.bt, b.cv] <- L(weighted.train.data, weighted.test.data)
      }
    }
  }
  result.vec.btcv <- as.vector(result.btcv)
  # random effects model
  id.bt <- rep(1:B.bt, B.cv)
  fit.result <- lmer(result.vec.btcv ~ (1 | id.bt))
  sigma.bt <- as.data.frame(VarCorr(fit.result))[1, 5]

  # bootstrap mixed effects model
  r.bt <- rep(NA, Brm.bt)
  for (brm.bt in 1:Brm.bt) {
    result.btcv.bt <- result.btcv[sample(B.bt, B.bt, T), ]
    result.vec.btcv.bt <- as.vector(result.btcv.bt)
    fit.result <- lmer(result.vec.btcv.bt ~ (1 | id.bt))
    r.bt[brm.bt] <- sigma.bt / as.data.frame(VarCorr(fit.result))[1, 5]
  }

  z.bt <- r.bt * rnorm(Brm.bt)
  return(list(sigma = sigma.bt, cutoff = quantile(abs(na.omit(z.bt)), 1 - alpha), factor = factor))
}
