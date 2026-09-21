// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List autotune_gg(arma::mat xin,
                 arma::vec yin,
                 arma::uvec group,
                 float alpha = 0.01,
                 bool standardize = true,
                 bool standardize_response = true,
                 bool intercept = true,
                 bool active = false,
                 float tau = 0.5,
                 bool trace_it = false,
                 double sigma_tolerance = 1e-8,
                 double beta_tolerance = 1e-8,
                 short int sigma_iter_max = 100,
                 short int active_iter_max = 100,
                 short int beta_iter_max = 100) {

  // CHECKS ----
  // Check if inputs xin, yin, and group are valid
  arma::vec y; arma::mat x;
  if(yin.is_vec()){
    y = yin;
  } else {
    stop("y must be a numeric vector");
  }

  if(!group.is_vec()){
    stop("group must be a numeric vector");
  }

  if(xin.is_colvec()){
    stop("x must be a numeric matrix with at least 2 variables");
  } else {
    x = xin;
  }


  // Check x, y dimensions
  int n = x.n_rows; int p = x.n_cols;
  if(y.n_elem != x.n_rows){
    stop("length of y does not match the number of rows in x");
  }

  // Check if group index is valid
  if(group.n_elem != x.n_cols){
    stop("Number of elements in group must equal the number of columns in x");
  } else if (!group.is_sorted("ascend") || group(0) != 1){
    stop("group should be a vector of non-decreasing consecutive integers starting with 1");
  }

  // Check valid significance level threshold
  if(alpha <= 0 || alpha >= 1){
    stop("alpha must be strictly between 0 and 1");
  }








  List cd_path_details = List::create(
    _["sorted_predictors"] = arma::conv_to<std::vector<double>>::from(active_indices + 1),
    _["sigma_sq_seq"] = sigma2_seq,
    _["no_of_iter_before_lambda_conv"] = sigma_iteration,
    _["no_of_iter_after_lambda_conv"] = beta_iteration,
    _["no_of_iterations"] = sigma_iteration + beta_iteration,
    _["support_set"] = support_set + 1,
    _["active_set"] = active_set + 1,
    _["count_sig_grps"] = vec_sig_beta_count,
    _["lambda0"] = lambda_value,
    _["null_support"] = null_support
  );

  return List::create(
    _["beta"] = beta,
    _["a0"] = intercept_estimate,
    _["lambda"] = lambda_effective / sqrt(n),
    _["sigma_sq"] = final_sigma,
    _["CD.path.details"] = cd_path_details
  );
}
