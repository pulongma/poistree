#ifndef _USE_Armadillo
#define _USE_Armadillo
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#endif

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <cmath>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
double base_mloglik(int n, double area, double a=1e-3, double b=1e-3){

	double loglik = lgamma(n+0.5) - (n+0.5)*log(area) - lgamma(n+1.0);

	return loglik;
}

double lbeta_ratio(int nL, int nR, double w0=1.0, double u0=.5){

	double lb1 = lgamma(w0*u0 + nL) + lgamma(w0*(1.0-u0)+nR) - lgamma(w0+nL+nR);
	double lb2 = lgamma(w0*u0) + lgamma(w0*(1.0-u0)) - lgamma(w0);
	return lb1 - lb2;
}

// [[Rcpp::export]]
double split_mloglik(int nL, int nR, double areaL, double areaR,
	double a=1e-3, double b=1e-3, double w0=1.0, double u0=.5){
	double lb_ratio = lbeta_ratio(nL, nR, w0, u0);
	double loglik_L, loglik_R;
	loglik_L = base_mloglik(nL, areaL, a, b);
	loglik_R = base_mloglik(nR, areaR, a, b);

	return lb_ratio + loglik_L + loglik_R;
}

double log_sum_exp(const arma::vec & x) {
  if(x.n_elem == 0) return -arma::datum::inf;
  if(x.n_elem == 1) return x(0);
  double c = x.max();
  return c + std::log(arma::sum(arma::exp(x - c)));
}

// [[Rcpp::export]]
double log_sum_exp_two(double x, double y){
	double result;

	if( (std::isinf(fabs(x)) == 1) && (std::isinf(fabs(y)) == 0) ){
		result = y;
	}else if((std::isinf(fabs(x)) == 0) && (std::isinf(fabs(y)) == 1)){
		result = x;
	}else if((std::isinf(fabs(x)) == 1) && (std::isinf(fabs(y)) == 1)){
		result = x;
	}else if(x-y>=100.0){
		result = x;
	}else if(x-y<=-100){
		result = y;
	}else{
		if(x>y){
			result = x + log(1.0 + exp(y-x));
		}else{
			result = y + log(1.0 + exp(x-y));
		}
	}

	return result;
}

arma::uvec in_region_nd(const arma::mat& x, const arma::mat& region) {
  int n = x.n_rows, d=x.n_cols;
  arma::uvec inside(n, arma::fill::ones);

  for(int j = 0; j < d; ++j) {
    double left = region(j, 0), right = region(j, 1);

    if(j < d-1) {
      inside = inside % (x.col(j) >= left) % (x.col(j) < right);
    } else {
      inside = inside % (x.col(j) >=left) % (x.col(j) <=right);
    }
  }
  return inside;
}

// Compute a type-1 quantile from sorted observations.
double quantile_type1(const arma::vec& sorted_x, double p) {
  int n = sorted_x.n_elem;
  if (n == 0) return NA_REAL;
  double h = (n - 1) * p + 1.0;
  int i = std::max(1, (int)std::floor(h));
  return sorted_x(i - 1);
}

// Evaluate log(1 + exp(x)) stably.
double log1pexp(double x)
{
    if (x <= 0.0)
        return std::log1p(std::exp(x));
    else
        return x + std::log1p(std::exp(-x));
}
