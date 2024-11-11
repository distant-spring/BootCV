## Packages
library(lme4)

#' Cross-Validation estimation and confidence interval
#'
#' CV.confint returns the cross-validation estimation and confidence interval based on
#' the result of Boot.CV ot Boot.Cal
#'
#' @details
#' This function calculates (1) the cross-validation estimation
#' \deqn{\widehat{Err}^{CV}_m=\frac{1}{B_{CV}}\sum_{b=1}^{B_{CV}} L\left\{D_{test}^b, \hat{\psi}(D_{train}^b)\right\},}{}
#' where data \eqn{D=D_{train}^b\cup D_{test}^b} represents the b-th split in cross-validation and
#' summary statistics L evaluates the performance of estimation \eqn{\hat{\psi}} fitted with training set in testing set.
#' and (2) the corresponding confidence interval
#' \deqn{\left[\widehat{Err}^{CV}_m-z_{1-\alpha/2}\times \widehat{\sigma}_m^{CV},  \widehat{Err}^{CV}_m+z_{1-\alpha/2}\times \widehat{\sigma}_m^{CV} \right] (\text{for Boot.CV})}{}
#' \deqn{\left[\widehat{Err}^{CV}_m-c_{1-\alpha/2}\times \widehat{\sigma}_m^{CV},  \widehat{Err}^{CV}_m+c_{1-\alpha/2}\times \widehat{\sigma}_m^{CV} \right] (\text{for Boot.Cal})}{}
#' where \eqn{\widehat{\sigma}_m^{CV}} is the bootstrap standard error from Boot.CV or Boot.Cal, \eqn{z_{1-\alpha/2}} and \eqn{c_{1-\alpha/2}} are corresponding quantiles.
#' @param boot the return list from Boot.CV or Boot.Cal function.
#' @param data data matrix of dimension nobs x nvars; each row is an observation vector.
#' @param L summary statistics with two inputs, train.data and test.data,
#' which are in the same form as data; it evaluates the performance of estimation fitted with train.data in test.data.
#' @param m training set size.
#' @param method 'Boot.CV' or 'Boot.Cal'; use Boot.CV or Boot.Cal algorithms.
#' @param adj bool value; whether to adjust for the reduced training sample size (default is TRUE).
#' @param B.cv the number of cross-validations (default is 400).
#' @param alpha 1-confidence level (default is 0.05).
#' @param print bool value; whether to print the result (default is FALSE).
#'
#' @return
#' \item{est}{the cross-validation estimation.}
#' \item{CI}{the bootstrap cross-validation confidence interval.}
#' @export
#'
#' @examples
#' set.seed(1)
#' # data generation
#' n=90
#' p=10
#' x=matrix(rnorm(n*p), ncol=p)
#' beta=rnorm(p)
#' y=x%*%beta+rnorm(n)
#' data=cbind(y,x)
#'
#' m=50 # training set size
#' # summary statistics
#' L=function(train.data,test.data){
#'   y=train.data[,1]
#'   x=train.data[,-1]
#'   yt=test.data[,1]
#'   xt=test.data[,-1]
#'
#'   fit=lm(y~x)
#'   beta=fit$coef
#'
#'   return(mean((yt-cbind(1,xt)%*%beta)^2))
#' }
#'
#' boot=Boot.CV(data,L,m)
#' result1=CV.confint(boot,data,L,m,method='Boot.CV',adj=T,print=T)
#' result2=CV.confint(boot,data,L,m,method='Boot.CV',adj=F,print=T)
#'
#' boot=Boot.Cal(data,L,m)
#' result1=CV.confint(boot,data,L,m,method='Boot.Cal',adj=T,print=T)
#' result2=CV.confint(boot,data,L,m,method='Boot.Cal',adj=F,print=T)
CV.confint=function(boot,data,L,m,method=c('Boot.CV','Boot.Cal'),adj=TRUE,B.cv=400,alpha=0.05,print=FALSE){
  # boot: return from Boot.CV or Boot.Cal function
  # data: data points in n x p matrix
  # L: summary statistics with two parameters, training set and testing set
  # m: training set size
  # method: use Boot.CV or Boot.Cal algorithms
  # adj: whether to adjust for the reduced training sample size (default is TRUE)
  # B.cv: the number of cross-validations (default is 400)
  # alpha: 1-confidence level (default is 0.05)
  # print: whether to print the result (default is FALSE)
  n=nrow(data)
  result.cv=rep(NA, B.cv)

  for(b in 1:B.cv){
    set.seed(b)
    id=sample(n, m, F)
    train.data=data[id,]
    test.data=data[-id,]

    result.cv[b]=L(train.data,test.data)
  }

  est=mean(na.omit(result.cv))
  sigma=ifelse(adj,boot$sigma*boot$factor,boot$sigma)
  cutoff=ifelse(method=='Boot.Cal',boot$cutoff,qnorm(1-alpha/2))
  CI=c(est-cutoff*sigma,est+cutoff*sigma)

  if(print==TRUE){
    cat('The Cross-Validation estimation is ',est,'.\n',sep='')
    cat('The confidence interval of Bootstrap Cross-Validation is [',CI[1],',',CI[2],'].\n',sep='')
  }

  return(list(est=est,CI=CI))
}
