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

#ifndef _USE_Utils
#define _USE_Utils
#include "utils.h"
#endif

#ifndef _USE_ProgressBar
#define _USE_ProgressBar
#include <R_ext/Utils.h>   // interrupt the Gibbs sampler from R
#include <iostream>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#endif


using namespace Rcpp;

#include "PPT.h"
#include "SMCtree.h"


/************************************************************************/
/************************************************************************/
// SMC for treed PP

void SMCtree::PPT_SMC(const arma::mat& pts,int max_depth, 
  int min_leaf_n, 
  double a, double b,
  double rho, double lam, int cut_grid_n,
  double max_aspect_ratio,
  bool verbose)
{
  // if(verbose){
  //   Rcpp::Rcout<<"\n***********************************************\n";
  //   Rcpp::Rcout<<" cSMC starts ......\n";    
  // }

  int n = pts.n_rows;
  int d = pts.n_cols;
  // arma::vec lamvec = lam * arma::vec(d, arma::fill::ones);
  lam = 1.0 / d;
  int max_nodes = std::pow(2, max_depth) - 1; // number of nodes up to the finest level
  arma::mat region_root;
  if (!this->particles.empty() && !this->particles[0].nodes.empty() &&
      this->particles[0].nodes[0] != nullptr) {
      region_root = this->particles[0].nodes[0]->region;
  } else {
      region_root = arma::mat(d, 2, arma::fill::zeros);
      region_root.col(1).fill(1.0);
  }
  
  // --- Particle and ancestry storage ---
  this->particles.resize(P);
  this->weights = arma::vec(P, arma::fill::ones) / P; // initial weights = 1/P
                                                      // (was 1.0/ones == all 1s)

  arma::vec logw(P, arma::fill::zeros);
  arma::vec ESS_hist(max_nodes, arma::fill::zeros);
  arma::vec logZ_inc(max_nodes, arma::fill::zeros);   // per-step incremental log evidence Delta_t
  arma::vec logZ_run(max_nodes, arma::fill::zeros);   // running cumulative log-evidence estimate
  double logZ = 0.0;   // running log evidence relative to the root model (unbiased for Z/Q0(root))


  // Initialize particles
  for (int p = 0; p < P; ++p) {
      particles[p].clear();
      particles[p].initialize(
        region_root, n, max_depth, min_leaf_n, a, b, rho, cut_grid_n,
        2.0, max_aspect_ratio
      );
      // Initialize lambda (you may have a helper for this)
      const double root_area = arma::prod(
          region_root.col(1) - region_root.col(0)
      );
      particles[p].nodes[0]->lambda = R::rgamma(
          particles[p].a + n, 1.0 / (particles[p].b + root_area)
      );
      // this->particles[p] = tree->deep_copy();
  }

  arma::uvec idx = arma::linspace<arma::uvec>(0, P-1, P); // For resampling
  
  for (int id = 0; id < max_nodes; ++id){
    double lse_before = log_sum_exp(logw);   // for the marginal-likelihood estimate
    for (int p = 0; p < P; ++p) {
      auto& tree = this->particles[p];
      if (id >= (int)tree.size() || tree[id] == nullptr || tree[id]->is_empty) continue;
      if (tree[id]->depth > max_depth) continue;
      if (tree[id]->idx.n_elem < (unsigned int)min_leaf_n) continue;
      // Advance particle by one node
      double log_inc;
      tree.PPT_one_step_ahead(log_inc, id, pts);
      // Replace nodes
      logw[p] += log_inc;
    }
    // Unbiased SMC log marginal-likelihood increment (adaptive-resampling form):
    //   logZ += lse(logw_after) - lse(logw_before).  After a resample logw is reset
    //   to 0 so lse(logw_before) = log(P) on the next step, which is correct.
    double dlogZ = log_sum_exp(logw) - lse_before;   // per-step log-evidence increment Delta_t
    logZ += dlogZ;
    logZ_inc(id) = dlogZ;                    // per-step increment (sum -> logZ_hat)
    logZ_run(id) = logZ;                     // running cumulative log evidence (last active entry -> logZ_hat)

    // // draw lambda at leaf
    // for (int p = 0; p < P; ++p) {
    //     particles[p].PPT_draw_lambda();
    // }

    // --- Normalize weights (with NaN/negative guard) ---
    double maxw = logw.max();
    this->weights = arma::exp(logw - maxw);
    for (arma::uword i = 0; i < this->weights.n_elem; ++i)
        if (!std::isfinite(this->weights[i]) || this->weights[i] < 0)
            this->weights[i] = 0.0;
    double sumw = arma::sum(this->weights);
    if (sumw <= 0.0) this->weights.fill(1.0 / P);
    else             this->weights /= sumw;
    double ESS = 1.0 / arma::sum(arma::square(this->weights));
    ESS_hist(id) = ESS;
    idx = arma::linspace<arma::uvec>(0, P-1, P);
    // Rcpp::Rcout<<" ESS "<<ESS<<", ";

    if (ESS < resample_thresh * P) {
    // Rcpp::Rcout<<" Resampling"<<"\n";

      // Multinomial resampling
      IntegerVector sampled = Rcpp::sample(P, P, true, NumericVector(this->weights.begin(), this->weights.end()), false);
      // Propagate children via idx
      std::vector< PPT > new_particles(P);
      for (int j = 0; j < P; ++j) {
        new_particles[j] = this->particles[sampled[j]].deep_copy();
      }
      // free memory before reassigning 
      for(int p=0; p<P; ++p){
        this->particles[p].clear();
      }
      this->particles = new_particles;
      for(int p=0; p<P; ++p){
        new_particles[p].clear();
      }
      logw.fill(0.0); 
      this->weights.fill(1.0/P);
    }
    // Rcpp::Rcout<<"\n";

  }

  for(int p=0; p<P; p++){
    this->particles[p].get_TreeLoglik();
  }


  this->ESS_hist = ESS_hist;
  this->logZ_hat = logZ;                // SMC log-evidence estimate, relative to root (= last active logZ_run)
  this->logZ_inc = logZ_inc;      // per-step increments, aligned step-for-step with ESS_hist
  this->logZ_run = logZ_run;      // running cumulative log evidence (cumsum of logZ_inc)


  return;
}

// Conditional SMC for an exact Particle Gibbs update.
//
// The first particle is a prefix of ref_tree.  Its action at every scheduled
// heap node is forced, while the remaining particles use the ordinary
// one-step-ahead proposal.  Resampling is performed at fixed tree-level
// boundaries, which keeps the schedule independent of the random particle
// system.  The reference particle is retained during every conditional
// resampling step, and the returned trajectory is sampled from the final
// normalized weights.
//
// This is deliberately Particle Gibbs without ancestor sampling.  A valid
// ancestor sampler for this non-Markovian tree construction requires the full
// suffix target ratio, not merely the probability of the next reference
// action.  The previous local replay approximation did not define a PGAS
// kernel and has therefore been removed.
void SMCtree::PPT_cSMC(
  PPT* ref_tree, const arma::mat& pts,int max_depth, 
  int min_leaf_n, 
  double a, double b, 
  double rho, double lam, int cut_grid_n,
  double max_aspect_ratio,
  bool verbose)
{

  int n = pts.n_rows;
  int d = pts.n_cols;
  lam = 1.0 / d;
  int max_nodes = std::pow(2, max_depth) - 1; // number of nodes up to the finest level
  if (ref_tree == nullptr || ref_tree->nodes.empty() ||
      ref_tree->nodes[0] == nullptr) {
      Rcpp::stop("conditional SMC requires a nonempty reference tree");
  }
  arma::mat region_root = ref_tree->nodes[0]->region;
  if (region_root.n_rows != static_cast<arma::uword>(d) ||
      region_root.n_cols != 2) {
      Rcpp::stop("conditional-SMC reference region has the wrong dimension");
  }

  if (P < 1) Rcpp::stop("conditional SMC requires at least one particle");

  // Every particle starts from the root-only state.  In particular, do not
  // copy the complete reference tree into particle zero: doing so exposes
  // future decisions before their SMC time and allows the ordinary proposal
  // to overwrite the purportedly conditioned trajectory.
  for (auto& tree : this->particles) tree.clear();
  this->particles.clear();
  this->particles.resize(P);
  for (int p = 0; p < P; ++p) {
      particles[p].initialize(
        region_root, n, max_depth, min_leaf_n, a, b, rho, cut_grid_n,
        2.0, max_aspect_ratio
      );
  }

  arma::vec logw(P, arma::fill::zeros);
  this->weights = arma::vec(P, arma::fill::ones) / P;
  this->ESS_hist = arma::vec(max_nodes, arma::fill::zeros);

  for (int id = 0; id < max_nodes; ++id){
      for (int p = 0; p < P; ++p) {
        auto& tree = this->particles[p];
        if (id >= (int)tree.size() || tree[id] == nullptr || tree[id]->is_empty) continue;
        if (tree[id]->depth > max_depth) continue;
        if (tree[id]->idx.n_elem < (unsigned int)min_leaf_n) continue;

        double log_inc;
        if (p == 0)
          tree.PPT_force_reference_step(log_inc, id, *ref_tree, pts);
        else
          tree.PPT_one_step_ahead(log_inc, id, pts);
        logw[p] += log_inc;
      }

      // --- Normalize weights ---
      double maxw = logw.max();
      if (std::isfinite(maxw)) this->weights = arma::exp(logw - maxw);
      else                     this->weights.zeros(P);

      // Remove NaNs and negative values, and re-normalize
      for (arma::uword i = 0; i < this->weights.n_elem; ++i) {
          if (!std::isfinite(this->weights[i]) || this->weights[i] < 0)
              this->weights[i] = 0.0;
      }
      double sum_weights = arma::sum(this->weights);
      if (sum_weights <= 0.0) {
          Rcpp::stop("conditional SMC produced no finite-weight particle");
      } else {
          this->weights /= sum_weights;
      }

      const double ESS = 1.0 / arma::sum(arma::square(this->weights));
      this->ESS_hist(id) = ESS;

      // Fixed resampling at the end of every completed heap level.  End-of-
      // level indices are 0, 2, 6, 14, ... (id + 2 is a power of two).  Do not
      // resample after the terminal SMC step; its weights select the output.
      const unsigned int level_marker = static_cast<unsigned int>(id + 2);
      const bool end_of_level =
          (level_marker & (level_marker - 1U)) == 0U;
      if (end_of_level && id + 1 < max_nodes) {
          std::vector< PPT > new_particles(P);
          new_particles[0] = this->particles[0].deep_copy();

          if (P > 1) {
            IntegerVector sampled = Rcpp::sample(
                P, P - 1, true,
                NumericVector(this->weights.begin(), this->weights.end()),
                false
            );
            for (int j = 1; j < P; ++j)
              new_particles[j] =
                  this->particles[sampled[j - 1]].deep_copy();
          }

          for (int p = 0; p < P; ++p) this->particles[p].clear();
          this->particles = std::move(new_particles);
          logw.fill(0.0);
          this->weights.fill(1.0/P);
      }
  }

  // Select the terminal trajectory by its final importance weight.  Each
  // particle already contains its complete genealogy because resampling uses
  // deep copies; following an ancestor matrix a second time is incorrect.
  const int selected = Rcpp::sample(
      P, 1, false,
      Rcpp::NumericVector(this->weights.begin(), this->weights.end()),
      false
  )[0];
  PPT ref_traj = this->particles[selected].deep_copy();
  ref_tree->clear();
  *ref_tree = ref_traj.deep_copy();
}


// Particle Gibbs for treed PP, using the exact conditional-SMC kernel above.
// The historical C++ entry-point name is retained for API compatibility, but
// ancestor sampling is intentionally disabled until an exact full-suffix
// backward weight is available for this tree construction.
Rcpp::List SMCtree::PPT_PGAS(
  const arma::mat& pts, const arma::mat& grid, int niter, int max_depth, 
  int min_leaf_n, 
  double a, double b,
  double rho, double lam, int cut_grid_n,
  double max_aspect_ratio,
  bool verbose)
{
    int n = pts.n_rows, d = pts.n_cols;
    if (max_depth < 0) max_depth = std::floor(std::log2(n / double(min_leaf_n)));
    max_nodes = std::pow(2, max_depth) - 1;
    // if (Rcpp::NumericVector::is_na(lam)) lam = 1.0 / d;
    lam = 1.0 / d;
    // --- Initialization: run ordinary SMC and sample its weighted output ---
    
    this->particles.resize(P);
    this->ESS_hist = arma::vec(max_nodes, arma::fill::zeros);
    this->weights = arma::vec(P, arma::fill::ones); // dummy initial weights

    // Initial SMC initializes all particles internally.
    PPT_SMC(
      pts, max_depth, min_leaf_n, a, b, rho, lam, cut_grid_n,
      max_aspect_ratio, false
    );
    // Particles and weights now updated in-place

    int initial_index = Rcpp::sample(
        P, 1, false,
        Rcpp::NumericVector(this->weights.begin(), this->weights.end()),
        false
    )[0];
    PPT* ref_tree = new PPT(particles[initial_index].deep_copy());
    // Prepare storage for reference trajectories
    std::vector<Rcpp::List> ref_traj_list(niter);
    int np = grid.n_rows;
    arma::mat lambda_sample(np, niter, arma::fill::zeros); 
    arma::vec loglik(niter, arma::fill::zeros);
    arma::mat weight_samp(P, niter, arma::fill::zeros);
    arma::vec lppd(niter, arma::fill::zeros);

    // Rcpp::Rcout<<" Start PGAS: \n";
    Progress prog(niter, verbose);
    // --- Particle Gibbs main loop ---
    for (int iter = 0; iter < niter; ++iter) {
        if(Progress::check_abort()){
          return R_NilValue;
        }
        prog.increment();

        // cSMC sweep conditioned on current reference trajectory
        PPT_cSMC(
          ref_tree, pts, max_depth, min_leaf_n, a, b, rho, lam,
          cut_grid_n, max_aspect_ratio, false
        );
        // ref_tree is updated in place
        // sample the lambda 
        ref_tree->PPT_draw_lambda();
        ref_tree->get_TreeLoglik();
        loglik(iter) = ref_tree->loglik;

        // compute lambda based on reference tree 
        // draw lambda over leaf nodes 
        // ref_tree.draw_lambda_at_leaf(); already done in the above 
        // assign lambda over the prediction locations (grid / XX) and evaluate
        // the predictive log-density there (lppd over the prediction points).
        lambda_sample.col(iter) = ref_tree->predict_lambda(grid);
        lppd(iter) = ref_tree->PPT_get_lppd(lambda_sample.col(iter));

        // // free all non-reference particles 
        // for(int p=1; p<P; ++p){
        //   particles[p].clear();
        // }

        // Store reference tree as Rcpp::List
        ref_traj_list[iter] = ref_tree->to_R_list();
        weight_samp.col(iter) = this->weights;
    }
    
    
    // compute posterior summary of lambda 
    arma::vec lambda_mean = arma::mean(lambda_sample, 1); // row-wise mean
    arma::vec lambda_median = arma::median(lambda_sample, 1); // row-wise median
    arma::vec lambda_lower = arma::quantile(lambda_sample, arma::vec({0.025}), 1);
    arma::vec lambda_upper = arma::quantile(lambda_sample, arma::vec({0.975}), 1);

    Rcpp::List Rlam = Rcpp::List::create(
                        Rcpp::_["mean"]=lambda_mean,
                        Rcpp::_["median"]=lambda_median,
                        Rcpp::_["lower95"]=lambda_lower,
                        Rcpp::_["upper95"]=lambda_upper,
                        Rcpp::_["draws"]=lambda_sample
                        ); 

    Rcpp::List out = Rcpp::List::create(
        Rcpp::_["loglik"] = loglik,
        Rcpp::_["lppd"] = lppd,
        Rcpp::_["lambda"] = Rlam,
        // Rcpp::_["Plambda"] = Plam,
        // Rcpp::_["ESS"] = ESS_hist,
        Rcpp::_["weights"] = weight_samp,
        Rcpp::_["particles"] = ref_traj_list

    );

    delete ref_tree;

    return out;
}

/************************************************************************/
/************************************************************************/
