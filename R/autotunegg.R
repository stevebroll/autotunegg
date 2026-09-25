#' Fast and automatic tuning and fitting for group lasso penalized least squares
#'
#' @name autotunegg
#'
#' @description
#' Fits group lasso penalized least squares model over regulraization path
#' determined by intermediate estimates of noise variance \eqn{\hat \sigma^2}.
#'
#' @param x Matrix of predictors, of dimension \eqn{n \times p} with \eqn{n}.
#'  observations of \eqn{p} variables
#' @param y Vector of responses, length \eqn{n}.
#' @param group Vector of consecutive integer group labels, length \eqn{p}
#' @param alpha Significance level for sequential Wald tests, default is 0.01
#' @param standardize Whether to standardize prior to tuning and fitting, default
#' is \code{TRUE}.
#' Coefficients are transformed back to original scale after model convergence.
#' Strongly recommended as this algorithm does not implement standardized group
#' lasso.
#' @param standarize_response Whether Y is demeaned, default is \code{TRUE}.
#' @param intercept Whether to include an intercept, default is \code{TRUE}.
#' @param active Whether to implement active set selection. Default is \code{TRUE} when
#' all groups have 5 or more members, and \code{FALSE} otherwise.
#' @param trace_it Whether to print out iteration details, default is \code{FALSE}.
#' @param sigma_tolerance Convergence termination tolerance for the initial
#' tuning of \eqn{\sigma^2} and \eqn{\lambda}. Default is \eqn{1e-8}.
#' @param beta_tolerance Convergence termination tolerance after
#' \eqn{\sigma^2} and \eqn{\lambda} are obtained. Default is \eqn{1e-8}.
#' @param sigma_iter_max Maximum number of iterations allowed for the initial
#' tuning of \eqn{\sigma^2} and \eqn{\lambda}. Default is \code{100}.
#' @param beta_iter_max Maximum  number of iterations allowed after \eqn{\sigma^2}
#' and \eqn{\lambda} are obtained. Default is \code{100}.
#' @param active_iter_max If \code{active = TRUE}, limits the number of updates
#' to the active support set. Default is \code{100}.
#' @param tau Scaling for initialization of \eqn{\lambda}. Initial
#' \eqn{\lambda = \lambda_{\max} \times \tau}. Must be smaller than \code{1} or the
#' algorithm will get stuck. This scaling is removed once \eqn{\hat{\sigma}^2}
#' is smaller than \eqn{\text{Var}(Y)}.
#'
#' @returns A list with final model fit as well as regularization path details.
#'
#' @export
#'
#' @examples
#' # generate simple data
#' library(autotunegg)
#' set.seed(2026)
#' n <- 100
#' pg <- 10
#' p <- pg*n
#' groups <- rep(1:n, each = pg)
#' # 2 target groups (beta = 1)
#' s <- 2*pg
#' beta <- 0; beta[groups == sample(unique(groups), 2)] = 1
#' x <- matrix(rnorm(n*p,n,p))
#' snr <- 4
#' error.sd <- sqrt(sum(beta^2)/ snr)
#' err <- rnorm(n, 0, error.sd)
#' y <- x %*% beta + err
#' fit <- autotunegg(x, y, groups, active = T, trace_it = T)
#' # # plot beta's
#' plot(beta, pch = 1)
#' plot(fit$beta, pch = 16)
#' # print confusion matrix
#' table(beta == 1, fit$beta == 1)
#' # Estimated \eqn{\hat \sigma^2} path vs empirical
#' fit$CD.path.details$count_sig_beta
#' var(err)

autotunegg <- function(
    x,
    y,
    group,
    alpha = 0.01,
    standardize = T,
    standarize_response = T,
    intercept = T,
    active = NULL,
    trace_it = FALSE,
    sigma_tolerance = 1e-8,
    beta_tolerance = 1e-8,
    sigma_iter_max = 100,
    beta_iter_max = 100,
    active_iter_max = 100,
    tau = 0.5
) {

  if(is.null(active)){
    if(max(table(group >= 1))){
      active = TRUE
    } else{
      active = FALSE
    }
  }

  fit <- autotunegg_bcd_cpp(
    xin = x,
    yin = y,
    alpha = alpha,
    standardize = standardize,
    standarize_response = standarize_response,
    intercept = intercept,
    active = active,
    trace_it = trace_it,
    sigma_tolerance = sigma_tolerance,
    beta_tolerance = beta_tolerance,
    sigma_iter_max = sigma_iter_max,
    active_iter_max = active_iter_max,
    beta_iter_max = beta_iter_max,
    tau = tau
  )

  fit$x <- x
  fit$y <- y

  return(fit)
}
