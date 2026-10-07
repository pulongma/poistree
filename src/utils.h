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

struct TreeNode {

    arma::mat region;
    arma::uvec idx;
    int depth;

    bool is_empty;
    bool is_leaf;

    int S;
    int J;
    double L;

    double prob_split;
    arma::vec prob_axis;
    arma::vec prob_cut;

    double lambda;

    TreeNode* clone() const {
        return new TreeNode(*this);
    }

    TreeNode() :
        depth(0), is_empty(true), is_leaf(true), S(-1), J(-1), L(NA_REAL),
        prob_split(NA_REAL)
    {}

    TreeNode(const arma::mat& region_,
             const arma::uvec& idx_,
             int depth_,
             bool is_empty_,
             bool is_leaf_,
             int S_ = -1,
             int J_ = -1,
             double L_ = NA_REAL,
             double prob_split_= NA_REAL,
             arma::vec prob_axis_=arma::vec(),
             arma::vec prob_cut_=arma::vec(),
             double lambda_ = 0
             )
      : region(region_), idx(idx_), depth(depth_),
        is_empty(is_empty_), is_leaf(is_leaf_),
        S(S_), J(J_), L(L_), prob_split(prob_split_), prob_axis(prob_axis_),
        prob_cut(prob_cut_), lambda(lambda_) {}

};

double base_mloglik(int n, double area, double a=1e-3, double b=1e-3);
double lbeta_ratio(int nL, int nR, double w0=1.0, double u0=.5);
double split_mloglik(int nL, int nR, double areaL, double areaR,
	double a=1e-3, double b=1e-3, double w0=1.0, double u0=.5);
double log_sum_exp(const arma::vec & x);
double log_sum_exp_two(double x, double y);
arma::uvec in_region_nd(const arma::mat& x, const arma::mat& region);
double quantile_type1(const arma::vec& sorted_x, double p);
double log1pexp(double x);
