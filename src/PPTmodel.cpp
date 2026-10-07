#ifndef _USE_Armadillo
#define _USE_Armadillo
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#endif

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <cmath>
#endif

using namespace Rcpp;

#include "PPT.h"
#include "SMCtree.h"

// [[Rcpp::export]]
Rcpp::List PPT_fit_PG(const arma::mat& pts, const arma::mat& grid, const arma::mat& region, double max_depth,
	int niter, int P, int min_leaf_n, double resample_thresh,
	double a, double b, bool verbose, double max_aspect_ratio, int cut_grid_n = 50)
{
  const int depth = ppt_checked_depth(max_depth);
  ppt_check_dense_storage(depth, P);
	int d = pts.n_cols;
	int n = pts.n_rows;

	double rho = 0.5, lam=1.0/d;

	bool force_mid_cut = false;

	PPT ref_tree;
	ref_tree.initialize(
        region, n, depth, min_leaf_n, a, b, rho, cut_grid_n,
        2.0, max_aspect_ratio
    );

    SMCtree smc(P, n, resample_thresh);
	smc.init_AS(ref_tree, region, n, depth,
		min_leaf_n, a, b, rho, cut_grid_n, max_aspect_ratio);
	smc.force_mid_cut = force_mid_cut;

	Rcpp::List PG = smc.PPT_PGAS(pts, grid, niter,
                          depth, min_leaf_n,
                          a, b,
                          rho, lam, cut_grid_n, max_aspect_ratio, verbose);

	return PG;

}

// Summarize weighted particle intensities and point-process diagnostics.

void summarize_particles(SMCtree&   Tsmc,
                             const arma::mat& grid,
                             const arma::vec& alphas,
                             arma::vec&      lam_mean,
                             arma::mat&      qmat,
                             double&         post_loglik,
                             double&         post_lppd,
                             arma::mat& lam_draws)
{
    int P  = Tsmc.P;
    int np = grid.n_rows;
    int K  = alphas.n_elem;

    lam_mean.zeros(np);
    post_loglik = 0.0;
    post_lppd   = 0.0;

    for (int p = 0; p < P; ++p)
    {
        Tsmc.particles[p].PPT_draw_lambda();
        Tsmc.particles[p].get_TreeLoglik();
        arma::vec lam = Tsmc.particles[p].predict_lambda(grid);

        lam_draws.col(p) = lam;

        lam_mean   += Tsmc.weights[p] * lam;
        post_loglik += Tsmc.weights[p] * Tsmc.particles[p].loglik;
        post_lppd   += Tsmc.particles[p].PPT_get_lppd(lam) * Tsmc.weights[p];
    }

    qmat.set_size(K, np);
	std::vector<std::pair<double,double>> vw(P);

    for (int j = 0; j < np; ++j)
    {

        for (int p = 0; p < P; ++p)
            vw[p] = { lam_draws(j, p), Tsmc.weights[p] };

        std::sort(vw.begin(), vw.end(),
                  [](auto& a, auto& b){ return a.first < b.first; });

        double csum = 0.0;
        int k = 0;
        for (const auto& pr : vw) {
            csum += pr.second;
            while (k < K && csum >= alphas(k)) {
                qmat(k, j) = pr.first;
                ++k;
            }
            if (k == K) break;
        }
        for (; k < K; ++k) qmat(k, j) = vw.back().first;
    }

	return;
}

// [[Rcpp::export]]
Rcpp::List PPT_fit_SMC(const arma::mat& pts, const arma::mat& grid, const arma::mat& region, double max_depth,
	int P, int min_leaf_n, double resample_thresh,
	double a, double b, double max_aspect_ratio, int cut_grid_n = 50)
{
  const int depth = ppt_checked_depth(max_depth);
  ppt_check_dense_storage(depth, P);
	int d = pts.n_cols;
	int n = pts.n_rows;
	int np = grid.n_rows;
	double rho = 0.5, lam=1.0/d;

	bool force_mid_cut = false;

	SMCtree pp_SMC1(P, n, resample_thresh);
	pp_SMC1.force_mid_cut = force_mid_cut;
	pp_SMC1.init(
        region, n, depth, min_leaf_n, a, b, rho, cut_grid_n,
        max_aspect_ratio
    );
	pp_SMC1.PPT_SMC(
        pts, depth, min_leaf_n, a, b, rho, lam, cut_grid_n,
        max_aspect_ratio
    );

	arma::vec alphas = {0.025, 0.5, 0.975};
	arma::vec lam_mean(np, arma::fill::zeros);
	arma::mat lam_draws(np, P);
	arma::mat qmat(3, np, arma::fill::zeros);
	double loglik_post, lppd_post;

	summarize_particles(pp_SMC1, grid, alphas,
                        lam_mean, qmat,
                        loglik_post, lppd_post, lam_draws);

	Rcpp::List PPparticle(P);
	for (int p = 0; p < P; ++p) {
	    PPparticle[p] = pp_SMC1.particles[p].to_R_list();
	}

	Rcpp::List lam_out = List::create(
		_["mean"] = lam_mean,
		_["median"] = qmat.row(1),
		_["lower95"] = qmat.row(0),
		_["upper95"] = qmat.row(2),
		_["draws"] = lam_draws
		);

	Rcpp::List SMCout = List::create(
  	_["loglik"] = loglik_post,
  	_["lppd"] = lppd_post,
  	_["lambda"] = lam_out,
  	_["particle"] = PPparticle,
  	_["weights"] = pp_SMC1.weights,
	_["resample_thresh"] = pp_SMC1.resample_thresh,
  	_["ESS"] = pp_SMC1.ESS_hist,
        _["logZ"] = pp_SMC1.logZ_hat,
    _["logZ_inc"] = pp_SMC1.logZ_inc,
    _["logZ_run"] = pp_SMC1.logZ_run
		);

	return SMCout;

}

// Inspect candidate cuts and their support constraints.
// [[Rcpp::export]]
Rcpp::List PPT_valid_cuts(const arma::mat& pts, const arma::mat& region,
                          int min_leaf_n, int cut_grid_n,
                          double max_aspect_ratio, bool force_mid_cut)
{
    PPT tree;
    tree.initialize(
        region, pts.n_rows, 1, min_leaf_n, 0.5, 0.0, 0.5,
        cut_grid_n, 2.0, max_aspect_ratio
    );
    std::vector<std::vector<double> > cuts = tree.find_valid_cuts(
        pts, region, force_mid_cut
    );
    Rcpp::List out(cuts.size());
    for (size_t j = 0; j < cuts.size(); ++j) out[j] = cuts[j];
    return out;
}

// Enumerate root transition probabilities for regression checks.

// [[Rcpp::export]]
Rcpp::List PPT_transition_probabilities(
    const arma::mat& pts, const arma::mat& region,
    int min_leaf_n, int cut_grid_n, double max_aspect_ratio)
{
    if (pts.n_rows == 0 || pts.n_cols == 0)
        Rcpp::stop("pts must be a non-empty matrix");
    if (region.n_rows != pts.n_cols || region.n_cols != 2)
        Rcpp::stop("region must be a d by 2 matrix");

    PPT tree;
    tree.initialize(
        region, pts.n_rows, 1, min_leaf_n, 0.5, 0.0, 0.5,
        cut_grid_n, 2.0, max_aspect_ratio
    );

    const TreeNode* parent = tree.nodes[0];
    std::vector<std::vector<double> > cuts = tree.find_valid_cuts(
        pts, region, false
    );

    PPT forward = tree.deep_copy();
    double log_increment = 0.0;
    forward.PPT_one_step_ahead(log_increment, 0, pts);
    const double proposal_split_probability = forward.nodes[0]->prob_split;

    std::vector<std::string> action;
    std::vector<int> axis;
    std::vector<double> cut;
    std::vector<double> log_probability;

    TreeNode stop = *parent;
    stop.S = 0;
    stop.J = -1;
    stop.L = NA_REAL;
    action.push_back("stop");
    axis.push_back(NA_INTEGER);
    cut.push_back(NA_REAL);
    log_probability.push_back(
        tree.PPT_log_transition_prob(parent, &stop, pts)
    );

    for (size_t j = 0; j < cuts.size(); ++j) {
        for (size_t k = 0; k < cuts[j].size(); ++k) {
            TreeNode split = *parent;
            split.S = 1;
            split.J = static_cast<int>(j);
            split.L = cuts[j][k];
            action.push_back("split");
            axis.push_back(static_cast<int>(j) + 1);
            cut.push_back(cuts[j][k]);
            log_probability.push_back(
                tree.PPT_log_transition_prob(parent, &split, pts)
            );
        }
    }

    Rcpp::NumericVector probability(log_probability.size());
    double replay_split_probability = 0.0;
    for (size_t k = 0; k < log_probability.size(); ++k) {
        probability[k] = std::exp(log_probability[k]);
        if (action[k] == "split") replay_split_probability += probability[k];
    }

    return Rcpp::List::create(
        Rcpp::_["action"] = action,
        Rcpp::_["axis"] = axis,
        Rcpp::_["cut"] = cut,
        Rcpp::_["log_probability"] = log_probability,
        Rcpp::_["probability"] = probability,
        Rcpp::_["proposal_split_probability"] = proposal_split_probability,
        Rcpp::_["replay_split_probability"] = replay_split_probability
    );
}

