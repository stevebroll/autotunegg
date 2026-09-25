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

  // STANDARDIZATION ----
  arma::field<arma::mat> Rfield(gmax), Qfield(gmax);
  if(standardize) {
    arma::mat Q, R;
    for(int g = 0; g < gmax; g++) {
      idx = arma::find(group == g);
      arma::mat Q, R;
      arma::qr_econ(Q, R, x.cols(idx));
      x.cols(idx) = Q;
      Qfield(g) = Q;
      Rfield(g) = R;
    }
  }

  double ymean = 0.0;
  if(standardize_response) {
    ymean = arma::mean(y);
    y = y - ymean;
  }

// OUTER LOOP ----

  while (error  > tolerance && iteration <= iter_max){
    lambda_effective = lambda_value * sigma2hat * tau;

    old_support_set = support_set;
    s = support_set.size();
    if(s > 0){
      support_set.erase(support_set.begin(), support_set.end());
    }
    old_beta = beta;

    // Group-wise iterative updates to beta and r
    for(int g = 0; g < gmax; g++){
      gactive = active_indices[g];
      idx = arma::find(group == gactive);

      arma::mat xg = x.cols(idx);
      r = r + (xg * beta(idx));
      beta_temp(idx) = xg.t() * r;
      double penscale = sqrt(pg[gactive]) * lambda_effective / (norm(xg * beta_temp(idx),2));
      beta(idx) = std::max((1 - penscale), 0.0) * beta_temp(idx);
      r = r - (xg * beta(idx));

    }

    // Partial residual ranking
      for(int g = 0; g < gmax; g++){
        idx = arma::find(group == g);
        arma::vec partial_res = r + (x.cols(idx) * beta(idx));
        partial_res_l2(g) = norm(partial_res, 2);
      }

    active_indices = seq_len(gmax) - 1;
    std::sort(active_indices.begin(), active_indices.end(),
              [&partial_res_l2](int a, int b) {
                return partial_res_l2[a] > partial_res_l2[b];
              });


    while((psum + pg[gactive]) < max_no_of_preds){
      k++;
      xg = x.cols(idx);
      arma::mat xtx = xg.t() * xg;
      betahat = xg.t() * ytemp;
      yhat = xg * betahat;
      arma::mat Sigma = arma::mat(pg[gactive], pg[gactive], arma::fill::eye) * sigma2hat;
      double f_stat = arma::as_scalar(betahat.t() * Sigma.i() * betahat);
      double cutoff = R::qf(1-alpha, pg[gactive], n - psum - pg[gactive], true, false) * pg[gactive];

      if(f_stat < cutoff){
        break;
      } else {
        ytemp = ytemp - yhat;
        support_set.push_back(gactive);
        psum = psum + pg[gactive];
        sigma2hat = pow(norm(ytemp,2),2)/(n - psum);
        gactive = active_indices[k];
        idx = arma::find(group == gactive);
      }
      iter++;

    }

    vec_sig_beta_count[iteration - 1] = iter;
    sigma2_seq[iteration - 1] = sigma2hat;


    beta_crit = abs(beta - old_beta)/ (1 + abs(beta));
    error = beta_crit.max();

    if (trace_it) {Rcout << "\rIteration: " << iteration << std::flush;}
    iteration++;


    if (setdiff(support_set, old_support_set).size() == 0) {
      flag--;
    } else {
      flag = 1;
    }

    if (flag == 0) {
      break;
    }

  }
  iteration--;




  // TRANSFORM BACK TO ORIGINAL BASIS ----

  if(standardize){
    x = xin;
    arma::vec beta_origscale = beta;
    for(int g = 0; g < gmax; g++) {
      idx = arma::find(group == g);
      if(norm(beta(idx), 2) > 0){
        arma::mat Q, R;
        R = Rfield(g);
        beta_origscale(idx) = R.i() * beta(idx);
      }
      beta = beta_origscale;
    }
  }

  if(standardize_response){
    y = y + ymean;
  }

  beta_iteration--;

  if (trace_it) {
    Rcout << "\nNo of predictor group significant for sigma estimation: " << vec_sig_beta_count[beta_iteration - 1] << std::endl;
  }

  arma::vec xmeans = mean(x,0).t();
  double intercept_estimate = 0.0;
  if(intercept) {
    intercept_estimate = ymean - arma::dot(beta, xmeans);
  }

  double final_sigma = rev(sigma2_seq)[0];




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
