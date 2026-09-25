
# autotunegg

<!-- badges: start -->
<!-- badges: end -->

autotunegg provides a fast, cross validation-free tuning for the group lasso with least squares loss.

## Installation

You can install the development version of autotunegg from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("stevebroll/autotunegg")
```

## Example Usage:

``` r
# generate simple data
library(autotunegg)
set.seed(2026)
n <- 100
pg <- 10
p <- pg*n
groups <- rep(1:n, each = pg)
# 2 target groups (beta = 1)
s <- 2*pg
beta <- 0; beta[groups == sample(unique(groups), 2)] = 1
x <- matrix(rnorm(n*p,n,p))
snr <- 4
error.sd <- sqrt(sum(beta^2)/ snr)
err <- rnorm(n, 0, error.sd)
y <- x %*% beta + err
fit <- autotunegg(x, y, groups, active = T, trace_it = T)
# # plot beta's
plot(beta, pch = 1)
plot(fit$beta, pch = 16)
# print confusion matrix
table(beta == 1, fit$beta == 1)
# Estimated \eqn{\hat \sigma^2} path vs empirical
fit$CD.path.details$count_sig_beta
var(err)
```

