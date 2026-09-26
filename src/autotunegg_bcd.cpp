// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List autotunegg_bcd_cpp(arma::mat xin,
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
                 short int beta_iter_max = 100,
                 short int active_iter_max = 100) {

  // Check xin, yin, group types
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
  int n = x.n_rows; int p = x.n_cols;

  // Check x, y observations
  if(y.n_elem != x.n_rows){
    stop("length of y does not match the number of rows in x");
  }

  // Check if group index is valid
  if(group.n_elem != x.n_cols){
    stop("Number of elements in group must equal the number of columns in x");
  } else if (!group.is_sorted("ascend")){
    stop("group should be a vector of non-decreasing consecutive integers");
  } else if(group(0) != 1){
    stop("group should start with 1");
  }
  int gmax = group.max();
  group = group - 1;
  arma::vec pg(gmax);
  int gactive;
  for(int g = 0; g < gmax; ++g){
    pg(g) = arma::accu(group == g);
  }
  // Check valid significance level threshold
  if(alpha <= 0 || alpha >= 1){
    stop("alpha must be strictly between 0 and 1");
  }

  // initialize parameters
  arma::vec beta, betahat, predmean, predsd; beta = predmean = predsd = arma::zeros(p);
  arma::mat xxt;
  arma::uvec idx;


  // Groupwise orthonormalization
  if(standardize) {
    for(int g = 0; g < gmax; g++) {
      idx = arma::find(group == g);
      arma::mat xg = x.cols(idx);
      arma::mat Q, R;
      arma::qr_econ(Q, R, xg);
      x.cols(idx) = Q;
    }
  }

  // standardize y
  double ymean = 0.0;
  if(standardize_response) {
    ymean = arma::mean(y);
    y = y - ymean;
  }


  // initalize variables for algorithm
  arma::vec r = y; double sigma2hat = arma::var(r);
  arma::ivec active_indices;
  IntegerVector support_set, old_support_set;
  active_indices = seq_len(gmax) - 1;
  NumericVector vec_sig_beta_count(sigma_iter_max), sigma2_seq(sigma_iter_max);
  // flagged for review in grp lasso case
  int max_no_of_preds = std::min((4*p) / 5, (4*n) / 5);
  double init_lambda, lambda_value;
  short int flag, s; flag = 1;
  bool null_support = FALSE;

  arma::vec temp(gmax);
  for (int g = 0; g < gmax; g++) {
    // temp(g) = norm(x.cols(arma::find(group == g)).t() * y, 2) / sqrt(pg[g]);
    temp(g) = norm(x.cols(arma::find(group == g)) * x.cols(arma::find(group == g)).t() * y, 2) / sqrt(pg[g]);

  }
  // init_lambda = max(temp) / sqrt(n);
  init_lambda = max(temp);

  lambda_value = init_lambda * (1.0  / (sigma2hat));


  int iteration = 1;
  double error = arma::datum::inf;
  arma::mat temp_stor(n,p);
  double lambda_effective, beta_l2norm;
  arma::vec old_beta(p), beta_temp(p), beta_crit(p), partial_res_l2(gmax);

  // Outer Loop

  while (error  > sigma_tolerance && iteration <= sigma_iter_max){
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

      // beta_temp(idx) = (x.cols(idx).t() * r) + beta(idx);
      // beta_l2norm = norm(beta_temp(idx),2) / sqrt(n);
      // double penscale = sqrt(pg[gactive]) * lambda_effective / beta_l2norm;
      // beta(idx) = std::max((1 - penscale), 0.0) * beta_temp(idx);
      // r = r - (x.cols(idx) * (beta(idx) - old_beta(idx)));

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

    arma::vec ytemp = y;
    int iter = 0;
    arma::mat xg;
    arma::vec yhat;

    arma::vec xk(n);

    int k = 0;
    gactive = active_indices[k];
    idx = arma::find(group == gactive);
    int psum = 1;

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

  if(support_set.size() == 0) {
    null_support = TRUE;
    if(iteration >= 2) {
      sigma2hat = sigma2_seq[iteration - 2];
    } else {
      sigma2hat = var(y) / 10;
    }
  }

  if (iteration < sigma_iter_max) {
    vec_sig_beta_count = vec_sig_beta_count[Range(0, iteration - 1)];
    sigma2_seq = sigma2_seq[Range(0, iteration - 1)];
  }

  lambda_effective = lambda_value * sigma2hat;
  int beta_iteration = 1, active_set_size = 0;
  error = R_PosInf;
  IntegerVector active_set = support_set;
  s = support_set.size();
  old_beta = beta;
  short int active_iterations = 1, sg_active = 0, iterations_finding_beta = 0;
  int psum = 0;
  NumericVector act_pred_count(active_iter_max);

  if(active){
    r = y, beta.zeros();
    for(int g = 0; g < s; g++){
      gactive = active_indices[g];
      idx = arma::find(group == gactive);
      beta(idx) = old_beta(idx);
    }
    r = r - x * (beta);

    while(active_iterations <= active_iter_max){
      beta_iteration = 1;
      error = R_PosInf;
      while(error > beta_tolerance && beta_iteration <= beta_iter_max) {
        old_beta = beta;
        psum = 0;
        for(int& g : active_set){
          gactive = active_indices[g];
          idx = arma::find(group == gactive);
          arma::mat xg = x.cols(idx);
          psum = psum + pg[gactive];

          r = r + (xg * beta(idx));
          beta_temp(idx) = xg.t() * r;
          double penscale = sqrt(pg[gactive]) * lambda_effective / (norm(xg * beta_temp(idx),2));
          beta(idx) = std::max((1 - penscale), 0.0) * beta_temp(idx);
          r = r - (xg * beta(idx));
        }
        beta_crit = abs(beta - old_beta)/ (1 + abs(beta));
        error = beta_crit.max();
        beta_iteration++;
      }
      iterations_finding_beta = iterations_finding_beta + --beta_iteration;
      act_pred_count[active_iterations - 1] = active_set.size();

      sg_active = 0;
      active_set_size = active_set.size();
      int j = active_set_size;
      gactive = active_indices[j];
      psum = psum + pg[gactive];

      while(psum < p){
        idx = arma::find(group == gactive);

        arma::mat xg = x.cols(idx);
        beta_l2norm = norm(xg * xg.t() * r,2);
        if(beta_l2norm >  (sqrt(n) * sqrt(pg[gactive]) * lambda_effective)){
          active_set.push_back(gactive);
          sg_active++;
        }

        j++;
        gactive = active_indices[j];
        psum = psum + pg[gactive];
      }
      if(sg_active == 0){
        break;
      }
      active_iterations++;
    }
  }else {
    while(error > beta_tolerance && beta_iteration <= beta_iter_max) {
      old_beta = beta;
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

      beta_crit = abs(beta - old_beta)/ (1 + abs(beta));
      error = beta_crit.max();
      if (trace_it) {Rcout << "\rLambda converged, Iteration: " << iteration + beta_iteration << std::flush;}
      beta_iteration++;
    }
  }


  if(standardize){
    arma::vec beta_origscale = beta;

    for(int g = 0; g < gmax; g++) {
      idx = arma::find(group == g);
      if(norm(beta(idx), 2) > 0){
        arma::mat xg = xin.cols(idx);
        arma::mat Q, R;
        arma::qr_econ(Q, R, xg);
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
    Rcout << "\nNo of predictor group significant for sigma estimation: " << vec_sig_beta_count[iteration - 1] << std::endl;
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
    _["no_of_iter_before_lambda_conv"] = iteration,
    _["no_of_iter_after_lambda_conv"] = beta_iteration,
    _["no_of_iterations"] = iteration + beta_iteration,
    _["support_set"] = support_set + 1,
    _["active_set"] = active_set + 1,
    _["count_sig_grps"] = vec_sig_beta_count,
    _["lambda0"] = lambda_value * tau,
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
