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
