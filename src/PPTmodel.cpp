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
Rcpp::List PPT_fit_PG(const arma::mat& pts, const arma::mat& grid, const arma::mat& region, int max_depth, 
	int niter=500, int P=1000, int min_leaf_n=2, double resample_thresh=0.5,
	double a=0.0, double b=0.0, bool verbose=true)
{
	int d = pts.n_cols;
	int n = pts.n_rows;
	// int np = grid.n_rows; 
	double rho = 0.5, lam=1.0/d; 
	int cut_grid_n = 30;

	bool force_mid_cut = false;
	// PPT 
	

	PPT ref_tree;
	ref_tree.initialize(region, n, max_depth, min_leaf_n, a, b, rho, cut_grid_n);
	

    SMCtree smc(P, n, resample_thresh);
	smc.init_AS(ref_tree, region, n, max_depth, 
		min_leaf_n, a, b, rho, cut_grid_n);	
	smc.force_mid_cut = force_mid_cut; 

	Rcpp::List PG = smc.PPT_PGAS(pts, grid, niter, 
                          max_depth, min_leaf_n, 
                          a, b,
                          rho, lam, cut_grid_n, verbose);



	return PG;

}


// Inputs
//------------------------------------------------------------------
//   pp_SMC1         : the completed SMC object
//   grid (np×d)     : prediction locations
//   alphas (K)      : e.g. {0.025, 0.50, 0.975}
//------------------------------------------------------------------
// Outputs
//------------------------------------------------------------------
//   lam_mean (np)          posterior mean of λ
//   qmat (K × np)          rows = requested α–quantiles
//   post_loglik            Σ w_p  ·  log p(y | tree_p )
//   post_lppd              Σ w_p  ·  lppd(tree_p)
//--------------------------------------------
void summarize_particles(SMCtree&   Tsmc,
                             const arma::mat& grid,          // np × d
                             const arma::vec& alphas,        // K
                             arma::vec&      lam_mean,       // ← np
                             arma::mat&      qmat,           // ← K × np
                             double&         post_loglik,
                             double&         post_lppd,
                             arma::mat& lam_draws)
{
    int P  = Tsmc.P;           // # particles
    int np = grid.n_rows;
    int K  = alphas.n_elem;

    // 1) containers ----------------------------------------------------
    // arma::mat lam_draws(np, P);                   // now location-major
    lam_mean.zeros(np);
    post_loglik = 0.0;
    post_lppd   = 0.0;


    // 2) collect λ(p,j) and accumulate means / log-sums ---------------
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

    // 3) weighted quantiles location-wise -----------------------------
    qmat.set_size(K, np);
	std::vector<std::pair<double,double>> vw(P);

    for (int j = 0; j < np; ++j)
    {
        // fill scratch with (λ,w) for fixed location j
        for (int p = 0; p < P; ++p)
            vw[p] = { lam_draws(j, p), Tsmc.weights[p] };

        std::sort(vw.begin(), vw.end(),
                  [](auto& a, auto& b){ return a.first < b.first; });

        double csum = 0.0;
        int k = 0;
        for (const auto& pr : vw) {
            csum += pr.second;
            while (k < K && csum >= alphas(k)) {
                qmat(k, j) = pr.first;          // row k, col j
                ++k;
            }
            if (k == K) break;
        }
        for (; k < K; ++k) qmat(k, j) = vw.back().first;
    }

	return;
}



// [[Rcpp::export]]
Rcpp::List PPT_fit_SMC(const arma::mat& pts, const arma::mat& grid, const arma::mat& region, int max_depth, 
	int P=1000, int min_leaf_n=2, double resample_thresh=0.5,
	double a=0.5, double b=0.0)
{
	int d = pts.n_cols;
	int n = pts.n_rows;
	int np = grid.n_rows; 
	double rho = 0.5, lam=1.0/d; 
	int cut_grid_n = 30;

	bool force_mid_cut = false;
	// PPT 
	SMCtree pp_SMC1(P, n, resample_thresh);
	pp_SMC1.force_mid_cut = force_mid_cut; 
	pp_SMC1.init(region, n, max_depth, min_leaf_n, a, b, rho, cut_grid_n);
	pp_SMC1.PPT_SMC(pts, max_depth, min_leaf_n, a, b, rho, lam, cut_grid_n);
	
	// arma::vec PPlambda_weighted(np, arma::fill::zeros);
	// double lppd=0.0;
	// double loglik =0.0; 

	// for(int p=0; p<P; p++){
	// 	pp_SMC1.particles[p].PPT_draw_lambda();
	// 	arma::vec lamp = pp_SMC1.particles[p].predict_lambda(grid);
	// 	PPlambda_weighted += pp_SMC1.weights[p] * lamp;
	// 	 pp_SMC1.particles[p].get_TreeLoglik(); 
	// 	loglik += pp_SMC1.particles[p].loglik * pp_SMC1.weights[p];
	// 	lppd += pp_SMC1.particles[p].PPT_get_lppd(lamp) * pp_SMC1.weights[p]; 
	// }

	arma::vec alphas = {0.025, 0.5, 0.975};
	arma::vec lam_mean(np, arma::fill::zeros);
	arma::mat lam_draws(np, P);
	arma::mat qmat(3, np, arma::fill::zeros);
	double loglik_post, lppd_post;

	summarize_particles(pp_SMC1, grid, alphas,
                        lam_mean, qmat,
                        loglik_post, lppd_post, lam_draws);

	// Serialise particles AFTER summarize_particles so each leaf carries its
	// drawn intensity (PPT_draw_lambda runs inside summarize_particles).
	// Previously this ran first, so returned leaf lambdas were all 0.
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
  	_["ESS"] = pp_SMC1.ESS_hist,
  	_["logZ"] = pp_SMC1.logZ_hat,      // scalar: SMC log marginal-evidence estimate, relative to root (final)
  	_["logZ_inc"] = pp_SMC1.logZ_inc,  // vector: per-step increment Delta_t (sum = logZ)
  	_["logZ_run"] = pp_SMC1.logZ_run   // vector: running cumulative estimate (last active entry = logZ)
		);

	return SMCout;

}


















