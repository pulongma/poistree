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


PPT SMCtree::get_trajectory(const arma::imat& ancestors, 
      const std::vector< PPT >& particles) {
  int n_nodes = ancestors.n_cols;
  int idx = 0;  // reference trajectory starts at first particle (0 in C++)
  for(int t = n_nodes - 1; t >= 0; --t) {
    idx = ancestors(idx, t);
  }
  // Return a copy of the reference trajectory
  PPT ref_tree;
  ref_tree = particles[idx].deep_copy();
  return ref_tree;
}

/************************************************************************/
/************************************************************************/
// SMC for treed PP

void SMCtree::PPT_SMC(const arma::mat& pts,int max_depth, 
  int min_leaf_n, 
  double a=.5, double b=0.0,
  double rho=0.5, double lam=1, int cut_grid_n=30,
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
  arma::mat region_root(d, 2, arma::fill::zeros);
  region_root.col(1).fill(1.0); // [0,1] hypercube
  
  // --- Particle and ancestry storage ---
  this->particles.resize(P);
  this->weights = arma::vec(P, arma::fill::ones) / P; // initial weights = 1/P
                                                      // (was 1.0/ones == all 1s)

  // std::vector< std::vector<TreeNode*> > particles(P);
  // arma::imat ancestors(P, max_nodes, arma::fill::ones);
  arma::vec logw(P, arma::fill::zeros);
  arma::vec ESS_hist(max_nodes, arma::fill::zeros);
  arma::vec logZ_inc(max_nodes, arma::fill::zeros);   // per-step incremental log evidence Delta_t
  arma::vec logZ_run(max_nodes, arma::fill::zeros);   // running cumulative log-evidence estimate
  double logZ = 0.0;   // running log evidence relative to the root model (unbiased for Z/Q0(root))


  // Initialize particles
  for (int p = 0; p < P; ++p) {
      particles[p].initialize(region_root, n, max_depth, min_leaf_n, a, b, rho, cut_grid_n);
      // Initialize lambda (you may have a helper for this)
      particles[p].nodes[0]->lambda = R::rgamma(particles[p].a + n, 1.0 / (particles[p].b + 1.0)); // Replace with correct area!
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

// SMC with AS
void SMCtree::PPT_cSMC(
  PPT* ref_tree, const arma::mat& pts,int max_depth, 
  int min_leaf_n, 
  double a=0.5, double b=0.0, 
  double rho=0.5, double lam=1.0, int cut_grid_n=30,
  bool verbose)
{

  int n = pts.n_rows;
  int d = pts.n_cols;
  // arma::vec lamvec = lam.isNotNull() ? as<arma::vec>(lam) : arma::vec(d, arma::fill::ones) / d;
  lam = 1.0 / d;
  arma::vec lamvec = lam * arma::vec(d, arma::fill::ones);
  int max_nodes = std::pow(2, max_depth) - 1; // number of nodes up to the finest level
  arma::mat region_root(d, 2, arma::fill::zeros);
  region_root.col(1).fill(1.0); // [0,1] hypercube
  
  // std::vector< std::vector<TreeNode*> > particles(P);
  

  // Initialize particles
  // deep copy ref_tree
  this->particles[0] = ref_tree->deep_copy();
  for (int p = 1; p < P; ++p) {
      particles[p].initialize(region_root, n, max_depth, min_leaf_n, a, b, rho, cut_grid_n);
      // Initialize lambda (you may have a helper for this)
      double alpha = particles[p].a * particles[p].get_area(0);
      double beta = particles[p].b + particles[p].get_area(0);
      particles[p].nodes[0]->lambda = R::rgamma(alpha + n, 1.0 / beta); // Replace with correct area!
  }
  

  // int count_resamp = 0;
  arma::uvec idx = arma::linspace<arma::uvec>(0, P-1, P); // For resampling
  // arma::vec weights = arma::vec(1/P, arma::fill::ones);
  arma::vec logw(P, arma::fill::zeros);

  for (int id = 0; id < max_nodes; ++id){

      for (int p = 0; p < P; ++p) {
        auto& tree = this->particles[p];
        if (id >= (int)tree.size() || tree[id] == nullptr || tree[id]->is_empty) continue;
        if (tree[id]->depth > max_depth) continue;
        if (tree[id]->idx.n_elem < (unsigned int)min_leaf_n) continue;
        // Advance particle by one node
        double log_inc;
        tree.PPT_one_step_ahead(log_inc, id, pts);
        // Replace nodes
        // this->particles[p] = tree.deep_copy(); // no need to do this 
        logw[p] += log_inc;
      }
      // --- Normalize weights ---
      double maxw = logw.max();
      weights = arma::exp(logw - maxw);
      weights /= arma::sum(weights);

      // Remove NaNs and negative values, and re-normalize
      for (int i = 0; i < weights.n_elem; ++i) {
          if (!std::isfinite(weights[i]) || weights[i] < 0) weights[i] = 0.0;
      }
      double sum_weights = arma::sum(weights);
      if (sum_weights <= 0.0) {
          // Defensive: if all zero, set uniform
          weights.fill(1.0 / weights.n_elem);
      } else {
          weights /= sum_weights;
      }

      double ESS = 1.0 / arma::sum(arma::square(weights));
      ESS_hist(id) = ESS;
      idx = arma::linspace<arma::uvec>(0, P-1, P);
      for (int p = 0; p < P; p++) this->ancestors(p, id) = p;
      // resample if necessary  
      if (ESS < resample_thresh * P) { 
          // count_resamp++;
          // Multinomial resampling
          IntegerVector sampled = Rcpp::sample(P, P, true, NumericVector(weights.begin(), weights.end()), false);
          
          // Propagate 
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
          // record ancestor indices
          for (int p = 0; p < P; p++) this->ancestors(p, id) = sampled[p];
      }

    // --- Ancestor sampling for the reference trajectory (PGAS) ---
    // The decision at node `id` is generated from node `id`'s own region/points
    // (see PPT_one_step_ahead), so the transition probability of reference node
    // `id` must be evaluated under particle p's node `id` -- NOT node `id-1`,
    // which in an array tree is a sibling/cousin, not the parent (the parent is
    // (id-1)/2). Backward weights are also formed in log-space with a
    // max-subtraction to avoid the over/underflow of exp(logw + logprob).
    // NOTE: this fixes a clear indexing/stability bug; the overall PGAS
    // ancestor-sampling scheme for this node-wise tree SMC should still be
    // unit-tested by the author (the primary, fully verified path is SMC).
      arma::vec log_bw(P, arma::fill::zeros);
      for (int p = 0; p < P; ++p) {
        const TreeNode* pnode =
            (id < (int)particles[p].nodes.size()) ? particles[p].nodes[id] : nullptr;
        double logprob = PPT_log_transition_prob(pnode, ref_tree->nodes[id], pts);
        log_bw[p] = logw[p] + logprob;
      }
      double max_lbw = log_bw.max();
      arma::vec backward_weights(P, arma::fill::zeros);
      if (std::isfinite(max_lbw))
          backward_weights = arma::exp(log_bw - max_lbw);

      for (arma::uword i = 0; i < backward_weights.n_elem; ++i)
          if (!std::isfinite(backward_weights[i]) || backward_weights[i] < 0)
              backward_weights[i] = 0.0;
      double sum_bw = arma::sum(backward_weights);
      if (sum_bw <= 0.0) backward_weights.fill(1.0 / backward_weights.n_elem);
      else               backward_weights /= sum_bw;

      int ancestor_idx = Rcpp::sample(P, 1, false,
                         Rcpp::NumericVector(backward_weights.begin(),
                         backward_weights.end()), false)[0];
      this->ancestors(0, id) = ancestor_idx; // only for p=0
  }

  // this->weights = weights;

  // --- Reconstruct new reference tree trajectory ---
  PPT ref_traj = this->get_trajectory(this->ancestors, particles);

  ref_tree->clear();
  *ref_tree = ref_traj.deep_copy();


  return;
}


// Particle Gibbs with ancestor sampling for treed PP
Rcpp::List SMCtree::PPT_PGAS(
  const arma::mat& pts, const arma::mat& grid, int niter, int max_depth, 
  int min_leaf_n, 
  double a=.5, double b=0.0,
  double rho=0.5, double lam=1, int cut_grid_n=30,
  bool verbose)
{
    int n = pts.n_rows, d = pts.n_cols;
    if (max_depth < 0) max_depth = std::floor(std::log2(n / double(min_leaf_n)));
    max_nodes = std::pow(2, max_depth) - 1;
    // if (Rcpp::NumericVector::is_na(lam)) lam = 1.0 / d;
    lam = 1.0 / d;
    arma::mat region_root(d, 2, arma::fill::zeros);
    region_root.col(1).fill(1.0);

    // --- Initialization: Run unconditional cSMC to get initial MAP reference tree ---
    
    this->particles.resize(P);
    this->ancestors = arma::imat(P, max_nodes, arma::fill::zeros); 
    this->ESS_hist = arma::vec(max_nodes, arma::fill::zeros);
    this->weights = arma::vec(P, arma::fill::ones); // dummy initial weights

    // Initialize all particles (single-tree initializers)
    for (int p = 0; p < P; ++p) {
        particles[p].initialize(region_root, n, max_depth, min_leaf_n, a, b, rho, cut_grid_n);
    }

    // initial cSMC
    PPT_SMC(pts, max_depth, min_leaf_n, a, b, rho, lam, cut_grid_n, false);
    // Particles and weights now updated in-place

    PPT* ref_tree = new PPT(particles[0].deep_copy());
    // Rcpp::Rcout<<" Initialization: MAP tree id="<<map_idx<<" is selected"<<"\n";
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

        for (int p = 0; p < P; ++p) particles[p].clear();
        particles.clear();
        particles.resize(P);

        // set reference trajectory and reinitialize 
        // particles[0] = ref_tree->deep_copy();
        // for(int p=1; p<P; ++p){
        //   particles[p].initialize(region_root, n, max_depth, min_leaf_n, a, b, rho, cut_grid_n);
        // } 
        // this->weights = arma::vec(P, arma::fill::ones) / P;  

        // cSMC sweep conditioned on current reference trajectory
        PPT_cSMC(ref_tree, pts, max_depth,min_leaf_n, a, b, rho, lam, cut_grid_n, false);
        // ref_tree is updated in place
        // sample the lambda 
        ref_tree->PPT_draw_lambda();
        ref_tree->get_TreeLoglik();
        loglik(iter) = ref_tree->loglik;

        // we do not need other particles except particle[0] 
        // since they are fully regenerated as in the next cSMC sweep 
        // free memory in SMC before overriding in the next iteration
        particles[0].clear();
        particles[0] = ref_tree->deep_copy();
        for(int p=1; p<P; ++p){
          particles[p].clear();
          particles[p].initialize(region_root, n, max_depth, min_leaf_n, a, b, rho, cut_grid_n);
        }

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
