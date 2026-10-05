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
  ppt_check_dense_storage(max_depth, P);
  int max_nodes = ppt_tree_steps(max_depth); // number of nodes up to the finest level
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
  arma::vec logZ_inc(max_nodes, arma::fill::zeros);   // per-step log relative-normalizer increment Delta_t
  arma::vec logZ_run(max_nodes, arma::fill::zeros);   // running log relative-normalizer estimate
  double logZ = 0.0;   // log relative normalizer estimate; exp(logZ), not logZ, is unbiased


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
    // Log of the relative SMC normalizer estimate (adaptive resampling):
    //   logZ += lse(logw_after) - lse(logw_before).  After a resample logw is reset
    //   to 0 so lse(logw_before) = log(P) on the next step, which is correct.
    double dlogZ = log_sum_exp(logw) - lse_before;   // per-step log relative-normalizer increment Delta_t
    logZ += dlogZ;
    logZ_inc(id) = dlogZ;                    // per-step increment (sum -> logZ_hat)
    logZ_run(id) = logZ;                     // running log relative normalizer (last active entry -> logZ_hat)

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

  this->ESS_hist = ESS_hist;
  this->logZ_hat = logZ;                // SMC log-normalizer estimate, relative to root (= last active logZ_run)
  this->logZ_inc = logZ_inc;      // per-step increments, aligned step-for-step with ESS_hist
  this->logZ_run = logZ_run;      // running log relative normalizer (cumsum of logZ_inc)


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
  ppt_check_dense_storage(max_depth, P);
  int max_nodes = ppt_tree_steps(max_depth); // number of nodes up to the finest level
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
    if (max_depth < 0) {
      if (n < 1 || min_leaf_n < 1) Rcpp::stop("Invalid automatic depth controls.");
      max_depth = ppt_checked_depth(std::floor(std::log2(n / double(min_leaf_n))));
    }
    ppt_check_dense_storage(max_depth, P);
    max_nodes = ppt_tree_steps(max_depth);
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
        prog.increment();
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


/************************************************************************/
/************************************************************************/
// Soft terminal-leaf PPT: shared-path SMC, conditional SMC with ancestor
// sampling, and Particle Gibbs.  See soft_smc.h for the model and notation and
// CLAUDE/shared_path_smc.tex for the algorithm.  Particles hold handles into
// `store`; the expansion of a coloured path node is computed once per distinct
// node (Phase 1), draws are made per particle in heap order (Phase 2), and
// resampling moves or copies record maps (Phase 3).

// Cache one-dimensional geometry, preserving the exact ordered path and the
// original integral routine. No quadrature approximation or RNG change is made.
double SoftSMCtree::cached_axis_log_integral(const arma::vec& gate,
    const std::vector<SoftGateStep>& path, int j,
    double extra_cut, int extra_side) {
  if ((int)axis_integrals.size() != M.d()) {
    axis_integrals.resize(M.d());
    integral_gate.set_size(M.d());
    integral_gate.fill(std::numeric_limits<double>::quiet_NaN());
  }
  auto& cache = axis_integrals[j];
  if (integral_gate[j] != gate[j]) {
    cache.clear();
    integral_gate[j] = gate[j];
  }
  AxisPathKey key;
  for (const SoftGateStep& g : path)
    if (g.axis == j) key.emplace_back(g.cut, g.side);
  if (extra_side != 0) key.emplace_back(extra_cut, extra_side);
  const auto found = cache.find(key);
  if (found != cache.end()) return found->second;
  const double value = soft_axis_log_integral(M, gate, path, j, extra_cut, extra_side);
  // Bound memory even for a fixed gate and a long chain. Eviction only affects speed.
  if (cache.size() >= 8192) cache.clear();
  cache.emplace(std::move(key), value);
  return value;
}

double SoftSMCtree::cached_log_exposure(const arma::vec& gate,
    const std::vector<SoftGateStep>& path) {
  double logH = 0.0;
  for (int j = 0; j < M.d(); ++j)
    logH += cached_axis_log_integral(gate, path, j);
  return logH;
}

void SoftSMCtree::build_gate_table(const arma::vec& gate) {
  const int d = M.d(), n = M.n();
  G.points = &M.X;
  G.gate = gate;
  G.width.set_size(d);
  for (int j = 0; j < d; ++j) G.width[j] = M.width(j);
  if (table_gate.n_elem != (arma::uword)d) {
    G.offset.assign(d + 1, 0);
    for (int j = 0; j < d; ++j) G.offset[j + 1] = G.offset[j] + (int)M.grid[j].size();
    G.n_cand = G.offset[d];
    G.cand_axis.resize(G.n_cand); G.cand_cut.resize(G.n_cand);
    for (int j = 0; j < d; ++j)
      for (size_t k = 0; k < M.grid[j].size(); ++k) {
        G.cand_axis[G.offset[j] + k] = j;
        G.cand_cut[G.offset[j] + k] = M.grid[j][k];
      }
    G.loglft.resize((size_t)n * G.n_cand);
    G.lft.resize((size_t)n * G.n_cand);
    table_gate.set_size(d);
    table_gate.fill(std::numeric_limits<double>::quiet_NaN());
  }
  for (int j = 0; j < d; ++j) {
    if (table_gate[j] == gate[j]) continue;
    for (int i = 0; i < n; ++i)
      for (int c = G.offset[j]; c < G.offset[j + 1]; ++c) {
        const double l = soft_log_left(M, gate, j, G.cand_cut[c], M.X(i, j));
        G.loglft[(size_t)i * G.n_cand + c] = l;
        G.lft[(size_t)i * G.n_cand + c] = std::exp(l);
      }
    table_gate[j] = gate[j];
  }
}

std::vector<SoftGateStep> SoftSMCtree::path_of(int v) const {
  std::vector<SoftGateStep> path;
  for (int u = v; store[u].parent >= 0; u = store[u].parent) path.push_back(store[u].gate);
  std::reverse(path.begin(), path.end());
  return path;
}

int SoftSMCtree::make_root() {
  SoftPathNode root;
  root.pts.resize(M.n());
  for (int i = 0; i < M.n(); ++i) root.pts[i] = i;
  root.axisH.resize(M.d());
  for (int j = 0; j < M.d(); ++j) root.axisH[j] = std::log(M.width(j));
  root.logH = std::accumulate(root.axisH.begin(), root.axisH.end(), 0.0);
  store.push_back(std::move(root));
  return (int)store.size() - 1;
}

// Child of store node `parent` for candidate `cand` on the given side, with the
// coloured points `pts`.  The exposure uses the parent's per-axis integrals.
int SoftSMCtree::make_child(int parent, int cand, int side, std::vector<int>&& pts) {
  const SoftPathNode& par = store[parent];
  SoftPathNode c;
  c.heap = 2 * par.heap + (side < 0 ? 0 : 1);
  c.depth = par.depth + 1;
  c.parent = parent;
  c.gate = {G.cand_axis[cand], G.cand_cut[cand], side};
  c.pts = std::move(pts);
  c.axisH = par.axisH;
  c.axisH[c.gate.axis] = side < 0 ? par.axisL[cand] : par.axisR[cand];
  c.logH = side < 0 ? par.logHL[cand] : par.logHR[cand];
  store.push_back(std::move(c));
  return (int)store.size() - 1;
}

// Return the shared reference for a cut, or a direct candidate reference when
// extreme scaled logits make the common shift unusable. The fallback is local.
std::shared_ptr<SoftExactAxis> SoftSMCtree::exact_cut_axis(SoftPathNode& A,
    int c, double& delta) {
  const int j = G.cand_axis[c], m = A.pts.size();
  auto axis = A.exact_axis.at(j);
  delta = cur_gate[j] * (G.cand_cut[c] - axis->refcut) / M.width(j);
  if (axis->finite && std::isfinite(delta) &&
      std::abs(delta) < std::numeric_limits<double>::max() / (m + 1.0)) return axis;
  std::vector<double> logits(m);
  for (int s = 0; s < m; ++s)
    logits[s] = cur_gate[j] * (G.cand_cut[c] - M.X(A.pts[s], j)) / M.width(j);
  delta = 0.0;
  return soft_exact_axis(logits, G.cand_cut[c], cur_gate[j], M.width(j), M);
}

// Phase 1: compute each coordinate once, retaining count coefficients for
// allocation draws. Root payloads survive sweeps; all other nodes are local.
void SoftSMCtree::expand_node(int v) {
  SoftPathNode& A = store[v];
  if (A.expanded) return;
  const std::vector<SoftGateStep> path = path_of(v);
  const int mA = A.pts.size(), d = M.d(), C = G.n_cand;
  const double HA = std::exp(A.logH), rho_d = M.rho_depth(A.depth);
  const bool exact = mA <= M.exact_max;
  const bool root = A.parent < 0;
  if (root && (int)root_axis_cache.size() != d) root_axis_cache.resize(d);
  A.logQA = M.logQ(mA, HA);
  A.axisL.resize(C); A.axisR.resize(C); A.logHL.resize(C); A.logHR.resize(C);
  A.count_norm.assign(C, std::numeric_limits<double>::quiet_NaN());
  A.score.assign(C + 1, 0.0); A.logprior.assign(C + 1, 0.0); A.logq.assign(C + 1, 0.0);
  A.score[0] = std::log(1.0 - rho_d) + A.logQA;
  A.logprior[0] = std::log(1.0 - rho_d);
  for (int j = 0; j < d; ++j) {
    const int begin = G.offset[j], end = G.offset[j + 1], Mj = end - begin;
    const double logprior = std::log(rho_d) - std::log((double)d) - std::log((double)Mj);
    const bool reuse = root && cache_root && root_axis_cache[j].ready &&
      root_axis_cache[j].gate == cur_gate[j] && root_axis_cache[j].exact == exact;
    if (reuse) {
      const RootAxisCache& cache = root_axis_cache[j];
      ++root_axis_hits;
      if (exact) A.exact_axis[j] = cache.exact_axis;
      for (int c = begin; c < end; ++c) {
        const int q = c - begin;
        A.axisL[c] = cache.axisL[q]; A.axisR[c] = cache.axisR[q];
        A.logHL[c] = cache.logHL[q]; A.logHR[c] = cache.logHR[q];
        A.count_norm[c] = cache.count_norm[q];
        A.logprior[c + 1] = logprior;
        A.score[c + 1] = logprior + cache.logPsi[q];
      }
      continue;
    }
    if (root) ++root_axis_builds;
    if (exact) {
      const auto bounds = std::minmax_element(M.grid[j].begin(), M.grid[j].end());
      const double refcut = *bounds.first + 0.5 * (*bounds.second - *bounds.first);
      std::vector<double> logits(mA);
      for (int s = 0; s < mA; ++s)
        logits[s] = cur_gate[j] * (refcut - M.X(A.pts[s], j)) / M.width(j);
      A.exact_axis[j] = soft_exact_axis(logits, refcut, cur_gate[j], M.width(j), M);
    }
    RootAxisCache cache;
    if (root && cache_root) cache.logPsi.reserve(Mj);
    for (int c = begin; c < end; ++c) {
      const double cut = G.cand_cut[c];
      A.axisL[c] = cached_axis_log_integral(cur_gate, path, j, cut, -1);
      A.axisR[c] = cached_axis_log_integral(cur_gate, path, j, cut, +1);
      const double otherH = A.logH - A.axisH[j];
      A.logHL[c] = otherH + A.axisL[c];
      A.logHR[c] = otherH + A.axisR[c];
      const double lbL = std::log(M.b + std::exp(A.logHL[c]));
      const double lbR = std::log(M.b + std::exp(A.logHR[c]));
      double logPsi;
      if (exact) {
        double delta;
        const auto axis = exact_cut_axis(A, c, delta);
        logPsi = soft_exact_log_psi(*axis, delta, M, lbL, lbR, &A.count_norm[c]);
      } else {
        logPsi = soft_log_psi_laplace(M, G, A.pts, c, lbL, lbR);
      }
      A.logprior[c + 1] = logprior;
      A.score[c + 1] = logprior + logPsi;
      if (root && cache_root) cache.logPsi.push_back(logPsi);
    }
    if (root && cache_root) {
      cache.ready = true; cache.exact = exact; cache.gate = cur_gate[j];
      cache.axisL.assign(A.axisL.begin() + begin, A.axisL.begin() + end);
      cache.axisR.assign(A.axisR.begin() + begin, A.axisR.begin() + end);
      cache.logHL.assign(A.logHL.begin() + begin, A.logHL.begin() + end);
      cache.logHR.assign(A.logHR.begin() + begin, A.logHR.begin() + end);
      cache.count_norm.assign(A.count_norm.begin() + begin, A.count_norm.begin() + end);
      if (exact) cache.exact_axis = A.exact_axis.at(j);
      root_axis_cache[j] = std::move(cache);
    }
  }
  A.logPhi = soft_lse(A.score);
  for (int k = 0; k <= C; ++k) {
    const double soft_part = std::log(1.0 - M.defensive) + A.score[k] - A.logPhi;
    A.logq[k] = M.defensive > 0.0 ? soft_lse({soft_part, std::log(M.defensive) + A.logprior[k]}) : soft_part;
  }
  A.log_inc = A.logPhi - A.logQA;
  A.expanded = true;
  ++n_expanded;
}

void SoftSMCtree::expand_level(int level) {
  const int lo = 1 << level, hi = 2 * lo;
  std::vector<int> frontier;
  for (const SoftParticleS& p : particles)
    for (auto it = p.rec.lower_bound(lo); it != p.rec.end() && it->first < hi; ++it)
      if (it->second.S == -1 && !store[it->second.node].expanded) frontier.push_back(it->second.node);
  std::sort(frontier.begin(), frontier.end());
  frontier.erase(std::unique(frontier.begin(), frontier.end()), frontier.end());
  for (int v : frontier) expand_node(v);
}

static inline unsigned long long soft_hash_bits(int cand, const std::vector<char>& B) {
  unsigned long long h = 1469598103934665603ULL ^ (unsigned long long)cand;
  for (char b : B) { h ^= (unsigned long long)(b + 1); h *= 1099511628211ULL; }
  return h;
}

// Phase 2 at heap position t: every particle with an undecided record at t
// draws (or, for particle 0 of a conditional SMC, is forced to) its action
// and colouring, creates or reuses the children, and updates its weight.
bool SoftSMCtree::sample_position(int t, SoftRef* ref) {
  bool advanced = false;
  for (int p = 0; p < P; ++p) {
    SoftParticleS& par = particles[p];
    auto it = par.rec.find(t);
    if (it == par.rec.end() || it->second.S != -1) continue;
    SoftRecord& R = it->second;
    const int v = R.node;
    SoftPathNode& A = store[v];
    if (A.depth >= M.Dmax) { R.S = 0; continue; }
    advanced = true;
    const bool forced = ref != nullptr && p == 0;
    const int mA = A.pts.size();
    const bool exact = mA <= M.exact_max;

    int act;
    if (!forced) act = soft_sample_log(A.logq);
    else {
      const SoftNodeP& D = ref_decision(*ref, t);
      act = 0;
      if (D.S == 1) {
        act = 1 + G.index(D.J, D.L);
        if (act <= 0) Rcpp::stop("conditional SMC: reference cut is not on the grid");
      }
    }
    if (act == 0) {
      R.S = 0;
      par.logw += exact && M.defensive == 0.0 ? A.log_inc : A.logprior[0] - A.logq[0];
      continue;
    }
    const int c = act - 1;
    const double HL = std::exp(A.logHL[c]), HR = std::exp(A.logHR[c]);
    std::vector<char> B(mA, 0);
    if (forced) for (int s = 0; s < mA; ++s) B[s] = ref_bit(*ref, A.pts[s], t);
    double logq_bits = 0.0;
    if (exact) {
      if (forced) {
        // q(B | cut) cancels from the exact importance ratio. Reference bits
        // only need routing; no count proposal, prefix table, or replay.
        ++exact_forced_routes;
      } else {
        double delta;
        const auto axis = exact_cut_axis(A, c, delta);
        const double lbL = std::log(M.b + HL), lbR = std::log(M.b + HR);
        const int k = soft_exact_sample_count(*axis, delta, lbL, lbR, A.count_norm[c]);
        if (k == 0 || k == mA) ++exact_boundary_draws;
        if (soft_exact_draw_bits(*axis, k, B)) ++exact_prefix_builds;
      }
    } else if (alloc_rates) {
      // auxiliary child rates (rate_explicit_smc_soft_ppt.tex, eq. for w_t):
      // lambda ~ q_lambda (Gamma posteriors at the expected allocation), then
      // independent Bernoulli allocations; the weight carries
      // q_lambda(lambda) prod Bernoulli / pi~(lambda | allocation).
      double Nbar = 0.0, Nref = 0.0;
      for (int s = 0; s < mA; ++s) { Nbar += G.r(A.pts[s], c); Nref += B[s]; }
      const double shL = M.a + (forced ? Nref : Nbar), shR = M.a + mA - (forced ? Nref : Nbar);
      const double lamL = R::rgamma(shL, 1.0 / (M.b + HL)), lamR = R::rgamma(shR, 1.0 / (M.b + HR));
      const double llL = std::log(lamL), llR = std::log(lamR);
      int N = 0;
      for (int s = 0; s < mA; ++s) {
        const int i = A.pts[s];
        const double wL = llL + G.left(i, c), wR = llR + G.right(i, c), D = soft_logaddexp(wL, wR);
        if (!forced) B[s] = std::log(R::unif_rand()) < wL - D ? 1 : 0;
        logq_bits += (B[s] ? wL : wR) - D;
        N += B[s];
      }
      logq_bits += R::dgamma(lamL, M.a + Nbar, 1.0 / (M.b + HL), 1) + R::dgamma(lamR, M.a + mA - Nbar, 1.0 / (M.b + HR), 1)
                 - R::dgamma(lamL, M.a + N, 1.0 / (M.b + HL), 1) - R::dgamma(lamR, M.a + mA - N, 1.0 / (M.b + HR), 1);
    } else {
      int mL = 0, mR = 0;
      for (int s = 0; s < mA; ++s) {
        const int i = A.pts[s];
        const double wL = G.left(i, c) + std::log(M.a + mL) - std::log(M.b + HL);
        const double wR = G.right(i, c) + std::log(M.a + mR) - std::log(M.b + HR);
        const double D = soft_logaddexp(wL, wR);
        if (!forced) B[s] = std::log(R::unif_rand()) < wL - D ? 1 : 0;
        logq_bits += (B[s] ? wL : wR) - D;
        if (B[s]) ++mL; else ++mR;
      }
    }

    std::vector<int> left, right;
    double lgate = 0.0;
    for (int s = 0; s < mA; ++s) {
      const int i = A.pts[s];
      if (B[s]) left.push_back(i); else right.push_back(i);
      if (!exact) lgate += B[s] ? G.left(i, c) : G.right(i, c);
    }
    if (exact) {
      // p0(c) Phi(c) / {Q(parent) q(c)}; with no defensive mixture all
      // actions have the same increment, avoiding cancellation roundoff.
      par.logw += M.defensive == 0.0 ? A.log_inc : A.score[act] - A.logQA - A.logq[act];
    } else {
      const int mL = left.size(), mR = right.size();
      const double log_target = A.logprior[act] + M.logQ(mL, HL) + M.logQ(mR, HR) + lgate - A.logQA;
      par.logw += log_target - A.logq[act] - logq_bits;
    }

    // children: reuse an identical coloured pair if some particle created it
    const unsigned long long key = soft_hash_bits(c, B);
    int lid = -1, rid = -1;
    auto range = A.children.equal_range(key);
    for (auto ch = range.first; ch != range.second; ++ch)
      if (store[ch->second.first].pts == left) { lid = ch->second.first; rid = ch->second.second; break; }
    if (lid < 0) {
      lid = make_child(v, c, -1, std::move(left));
      rid = make_child(v, c, +1, std::move(right));
      store[v].children.emplace(key, std::make_pair(lid, rid));
    }
    SoftRecord& R2 = par.rec[t];           // `A`/`R` may be invalidated by store growth
    R2.S = 1; R2.act = c;
    par.rec[2 * t] = SoftRecord{lid, -1, -1};
    par.rec[2 * t + 1] = SoftRecord{rid, -1, -1};
  }
  return advanced;
}

// Phase 3: conditional multinomial resampling of the record maps (first
// offspring by move, further ones by copy); with a reference and ancestor
// sampling, particle 0's ancestor is drawn from the ancestor weights.  The
// event takes place iff ESS <= ess_threshold * P: always with the default
// threshold 1, adaptively below it (e.g. 0.999 skips equal-weight events).
// The rule is a function of the weights, so the law of the sweep stays well
// defined and the Particle Gibbs argument is unchanged (shared_path_smc.tex).
void SoftSMCtree::resample(SoftRef* ref) {
  std::vector<double> lw(P);
  for (int p = 0; p < P; ++p) lw[p] = particles[p].logw;
  const double norm = soft_lse(lw);
  Rcpp::NumericVector W(P);
  double sumsq = 0.0;
  for (int p = 0; p < P; ++p) { W[p] = std::exp(lw[p] - norm); sumsq += W[p] * W[p]; }
  if (1.0 / sumsq > ess_threshold * P * (1.0 + 1e-9)) return;
  ++n_resampled;
  std::vector<int> anc(P);
  if (ref != nullptr) {
    int first = 0;
    if (use_as) {
      std::vector<double> alw(P);
      for (int p = 0; p < P; ++p) alw[p] = as_log_weight(particles[p], *ref);
      first = soft_sample_log(alw);
      if (first != 0) ++as_moved;
    }
    anc[0] = first;
    if (P > 1) {
      Rcpp::IntegerVector s = Rcpp::sample(P, P - 1, true, W, false);
      for (int p = 1; p < P; ++p) anc[p] = s[p - 1];
    }
  } else {
    Rcpp::IntegerVector s = Rcpp::sample(P, P, true, W, false);
    for (int p = 0; p < P; ++p) anc[p] = s[p];
  }
  std::vector<SoftParticleS> next(P);
  std::vector<int> first_use(P, -1);
  for (int p = 0; p < P; ++p) {
    const int a = anc[p];
    if (first_use[a] < 0) { first_use[a] = p; next[p].rec = std::move(particles[a].rec); }
    else next[p].rec = next[first_use[a]].rec;
    next[p].logw = 0.0;
  }
  particles.swap(next);
}

// Dense export of a particle (decisions, labels) as a reference trajectory.
void SoftSMCtree::export_particle(const SoftParticleS& p, SoftRef& out) {
  out.node.assign(M.n_nodes, SoftNodeP());
  out.z.assign(M.n(), 1);
  out.bits.clear();
  out.valid = true;
  for (const auto& kv : p.rec) {
    const SoftPathNode& v = store[kv.second.node];
    SoftNodeP& N = out.node[kv.first];
    N.active = true; N.depth = v.depth; N.m = v.pts.size(); N.logH = v.logH; N.path = path_of(kv.second.node);
    if (kv.second.S == 1) { N.S = 1; N.J = G.cand_axis[kv.second.act]; N.L = G.cand_cut[kv.second.act]; N.m = 0; }
    else { N.S = 0; for (int i : v.pts) out.z[i] = kv.first; }
  }
}

// Reference decision at heap node h; an inactive node is a prior draw.
const SoftNodeP& SoftSMCtree::ref_decision(SoftRef& ref, int h) {
  SoftNodeP& R = ref.node[h];
  if (R.S == -1) {
    R.depth = M.depth_of(h);
    if (R.depth >= M.Dmax) { R.S = 0; return R; }
    R.S = R::unif_rand() < M.rho_depth(R.depth) ? 1 : 0;
    if (R.S == 1) {
      R.J = std::min(M.d() - 1, (int)(R::unif_rand() * M.d()));
      R.L = M.grid[R.J][std::min((int)M.grid[R.J].size() - 1, (int)(R::unif_rand() * M.grid[R.J].size()))];
    }
  }
  return R;
}

// Reference bit of point i at node h: implied by the reference label when h is
// on the route of i, otherwise an inactive bit drawn from Bern(left gate).
int SoftSMCtree::ref_bit(SoftRef& ref, int i, int h) {
  int node = ref.z[i];
  while (node > h) {
    if (node / 2 == h) return node == 2 * h ? 1 : 0;
    node /= 2;
  }
  const long long key = (long long)i * (long long)M.n_nodes + h;
  auto it = ref.bits.find(key);
  if (it != ref.bits.end()) return it->second;
  const SoftNodeP& R = ref_decision(ref, h);
  const char bit = std::log(R::unif_rand()) < G.left(i, G.index(R.J, R.L)) ? 1 : 0;
  ref.bits[key] = bit;
  return bit;
}

// Sum of log Q over the glued leaves below node h when `pts` are routed
// through the reference suffix.
double SoftSMCtree::glued_logQ(SoftRef& ref, int h, const std::vector<int>& pts,
                               const std::vector<SoftGateStep>& path) {
  const SoftNodeP& R = ref_decision(ref, h);
  if (R.S == 0) return M.logQ(pts.size(), std::exp(cached_log_exposure(cur_gate, path)));
  std::vector<int> left, right;
  for (int i : pts) (ref_bit(ref, i, h) ? left : right).push_back(i);
  std::vector<SoftGateStep> pL = path, pR = path;
  pL.push_back({R.J, R.L, -1}); pR.push_back({R.J, R.L, +1});
  return glued_logQ(ref, 2 * h, left, pL) + glued_logQ(ref, 2 * h + 1, right, pR);
}

// Ancestor-sampling log weight; the glued quantity of a frontier node depends
// on the immutable coloured path node, gate, and reference only.  Cache it for
// the whole fixed-reference sweep.  Inactive reference decisions/bits are
// materialized lazily once and never changed during that sweep; evaluating a
// glued subtree materializes every reference value needed by that result.
double SoftSMCtree::as_log_weight(const SoftParticleS& p, SoftRef& ref) {
  double lw = p.logw;
  for (const auto& kv : p.rec) {
    if (kv.second.S != -1) continue;
    SoftPathNode& v = store[kv.second.node];
    if (v.depth >= M.Dmax) continue;
    if (v.glued_stamp != as_stamp) {
      v.glued = glued_logQ(ref, kv.first, v.pts, path_of(kv.second.node)) - M.logQ(v.pts.size(), std::exp(v.logH));
      v.glued_stamp = as_stamp;
    }
    lw += v.glued;
  }
  return lw;
}

// One SMC sweep (ref_in == nullptr) or one conditional SMC sweep with the
// reference in particle 0.
void SoftSMCtree::sweep(SoftRef* ref_in, const arma::vec& gate, SoftRef& ref_out) {
  cur_gate = gate;
  build_gate_table(gate);
  store.clear();
  // New path store and fixed reference: invalidate once per sweep, not at
  // every resampling event.  Avoid integer overflow in very long runs.
  as_stamp = (as_stamp == std::numeric_limits<int>::max()) ? 0 : as_stamp + 1;
  n_expanded = 0;
  as_moved = 0;
  n_resampled = 0;
  const int root = make_root();
  particles.assign(P, SoftParticleS());
  for (int p = 0; p < P; ++p) particles[p].rec[1] = SoftRecord{root, -1, -1};

  // candidate resampling events: after every heap position at which some
  // particle advanced (resample_node) or after every level; never after the
  // last position, whose weights select the output trajectory
  for (int level = 0; level < M.Dmax; ++level) {
    expand_level(level);
    const int lo = 1 << level, hi = 2 * lo;
    // A split creates only next-level children; resampling only copies current
    // particles.  Thus no new current-level position can appear after this
    // union is collected.  Retain heap order to preserve the random stream.
    std::vector<int> active_positions;
    for (const SoftParticleS& p : particles)
      for (auto it = p.rec.lower_bound(lo); it != p.rec.end() && it->first < hi; ++it)
        if (it->second.S == -1) active_positions.push_back(it->first);
    std::sort(active_positions.begin(), active_positions.end());
    active_positions.erase(std::unique(active_positions.begin(), active_positions.end()),
                           active_positions.end());
    for (int t : active_positions) {
      // A prior resampling may have removed this position from all particles.
      const bool advanced = sample_position(t, ref_in);
      const bool final_position = level == M.Dmax - 1 && t == hi - 1;
      if (resample_node && advanced && !final_position) resample(ref_in);
    }
    // Preserve the original level-boundary events even for an empty level.
    if (!resample_node && level < M.Dmax - 1) resample(ref_in);
  }

  std::vector<double> lw(P);
  for (int p = 0; p < P; ++p) lw[p] = particles[p].logw;
  const double norm = soft_lse(lw);
  weights.set_size(P);
  for (int p = 0; p < P; ++p) weights[p] = std::exp(lw[p] - norm);
  final_ess = 1.0 / arma::accu(arma::square(weights));
  export_particle(particles[soft_sample_log(lw)], ref_out);
}

// Gibbs sweep over the labels given the tree and the gate.
void SoftSMCtree::label_sweep(SoftRef& ref, const arma::vec& gate) {
  std::vector<int> leaves;
  for (int h = 1; h < M.n_nodes; ++h) if (ref.node[h].active && ref.node[h].S == 0) leaves.push_back(h);
  std::vector<double> lw(leaves.size());
  for (int i = 0; i < M.n(); ++i) {
    ref.node[ref.z[i]].m -= 1;
    for (size_t k = 0; k < leaves.size(); ++k) {
      const SoftNodeP& U = ref.node[leaves[k]];
      lw[k] = soft_log_phi(M, gate, U.path, i) + std::log(M.a + U.m) - std::log(M.b + std::exp(U.logH));
    }
    const int pick = leaves[soft_sample_log(lw)];
    ref.z[i] = pick;
    ref.node[pick].m += 1;
  }
}

void SoftSMCtree::refresh_exposures(SoftRef& ref, const arma::vec& gate) {
  for (int h = 1; h < M.n_nodes; ++h)
    if (ref.node[h].active) ref.node[h].logH = cached_log_exposure(gate, ref.node[h].path);
}

double SoftSMCtree::log_target_gate(const SoftRef& ref, const arma::vec& gate,
                                    const arma::vec& a_gate, const arma::vec& b_gate,
                                    const arma::vec& gate_min, bool shared) {
  if (!gate.is_finite() || arma::any(gate <= gate_min) || arma::any(gate <= 0.0))
    return -std::numeric_limits<double>::infinity();
  double out = 0.0;
  const arma::uword prior_n = shared ? 1 : gate.n_elem;
  for (arma::uword j = 0; j < prior_n; ++j) out += (a_gate[j] - 1.0) * std::log(gate[j]) - b_gate[j] * gate[j];
  for (int h = 1; h < M.n_nodes; ++h) {
    const SoftNodeP& U = ref.node[h];
    if (!U.active || U.S != 0) continue;
    out += M.logQ(U.m, std::exp(cached_log_exposure(gate, U.path)));
  }
  for (int i = 0; i < M.n(); ++i) out += soft_log_phi(M, gate, ref.node[ref.z[i]].path, i);
  return out;
}

// Log-random-walk Metropolis update of the gate (shared or per axis).
// One Metropolis proposal per coordinate (systematic scan), or a single
// proposal for a shared gate; returns the number of accepted proposals.
int SoftSMCtree::gate_update(SoftRef& ref, arma::vec& gate, const arma::vec& a_gate,
                             const arma::vec& b_gate, const arma::vec& sd_gate,
                             const arma::vec& gate_min, bool shared) {
  const int nup = shared ? 1 : (int)gate.n_elem;
  int accepted = 0;
  double log_cur = log_target_gate(ref, gate, a_gate, b_gate, gate_min, shared);
  for (int which = 0; which < nup; ++which) {
    arma::vec prop = gate;
    const double cur = gate[which], proposed = std::exp(std::log(cur) + R::rnorm(0.0, sd_gate[which]));
    if (shared) prop.fill(proposed); else prop[which] = proposed;
    const double log_prop = log_target_gate(ref, prop, a_gate, b_gate, gate_min, shared);
    if (std::log(R::unif_rand()) < log_prop - log_cur + std::log(proposed) - std::log(cur)) {
      gate = prop; log_cur = log_prop; ++accepted;
    }
  }
  if (accepted) refresh_exposures(ref, gate);
  return accepted;
}

// Particle Gibbs with ancestor sampling for S-PPT.  Gibbs cycle: conditional
// SMC over (tree, labels) | gate; label sweeps; gate update; rates drawn
// conjugately for output.  Rows of each state matrix: (heap id, axis, cut,
// lambda, NA, m) as consumed by ppt_eval_state().
Rcpp::List SoftSMCtree::PGAS(const arma::mat& grid, const arma::mat& xtest,
                             arma::vec gate, const arma::vec& a_gate, const arma::vec& b_gate,
                             const arma::vec& sd_gate, const arma::vec& gate_min, bool gate_shared,
                             int niter, int burn, int thin, int label_sweeps, bool update_gate,
                             bool verbose) {
  const int n = M.n();
  int ns = 0;
  for (int it = burn; it < niter; ++it) if ((it - burn) % thin == 0) ++ns;

  arma::mat draws(grid.n_rows, ns, arma::fill::zeros);
  arma::vec loglik(ns), loglik_test(xtest.n_rows > 0 ? ns : 0), integrated(ns), nleaf(ns), maxdepth(ns), ess(ns), as_rate(ns), expanded(ns), resampled(ns);
  arma::mat gates(ns, M.d());
  std::vector<arma::mat> states(ns);
  double gate_acc = 0.0, gate_tot = 0.0;

  SoftRef ref;
  sweep(nullptr, gate, ref);
  refresh_exposures(ref, gate);

  Progress prog(niter, verbose);
  int si = 0;
  for (int it = 0; it < niter; ++it) {
    if (Progress::check_abort()) return R_NilValue;

    for (int h = 1; h < M.n_nodes; ++h) if (!ref.node[h].active) ref.node[h] = SoftNodeP();
    ref.bits.clear();
    SoftRef next;
    sweep(&ref, gate, next);
    ref = next;
    const double as_frac = n_resampled > 0 ? (double)as_moved / n_resampled : NA_REAL;

    for (int k = 0; k < label_sweeps; ++k) label_sweep(ref, gate);
    if (update_gate) {
      gate_acc += gate_update(ref, gate, a_gate, b_gate, sd_gate, gate_min, gate_shared);
      gate_tot += gate_shared ? 1.0 : (double)gate.n_elem;
    }

    if (it >= burn && (it - burn) % thin == 0) {
      std::vector<int> active;
      for (int h = 1; h < M.n_nodes; ++h) if (ref.node[h].active) active.push_back(h);
      arma::mat st(active.size(), 6);
      std::unordered_map<int, double> lam;
      double comp = 0.0; int nl = 0, md = 0;
      for (size_t r = 0; r < active.size(); ++r) {
        const SoftNodeP& U = ref.node[active[r]];
        double lv = 0.0;
        if (U.S == 0) {
          const double H = std::exp(U.logH);
          lv = R::rgamma(M.a + U.m, 1.0 / (M.b + H));
          lam[active[r]] = lv; comp += lv * H; ++nl;
        }
        md = std::max(md, U.depth);
        st(r, 0) = active[r]; st(r, 1) = U.S == 1 ? U.J : -1;
        st(r, 2) = U.S == 1 ? U.L : NA_REAL; st(r, 3) = lv; st(r, 4) = NA_REAL; st(r, 5) = U.m;
      }
      auto intensity = [&](const arma::rowvec& x) {
        double out = 0.0;
        for (const auto& kv : lam) {
          double lp = 0.0;
          for (const SoftGateStep& g : ref.node[kv.first].path) lp += soft_log_gate(M, gate, g, x[g.axis]);
          out += kv.second * (lp < -745.0 ? 0.0 : std::exp(lp));
        }
        return std::max(out, 1e-300);
      };
      for (arma::uword g = 0; g < grid.n_rows; ++g) draws(g, si) = intensity(grid.row(g));
      double ll = -comp;
      for (int i = 0; i < n; ++i) ll += std::log(intensity(M.X.row(i)));
      loglik[si] = ll;
      if (xtest.n_rows > 0) {
        double lt = -comp;
        for (arma::uword i = 0; i < xtest.n_rows; ++i) lt += std::log(intensity(xtest.row(i)));
        loglik_test[si] = lt;
      }
      integrated[si] = comp; nleaf[si] = nl; maxdepth[si] = md; ess[si] = final_ess; as_rate[si] = as_frac; expanded[si] = n_expanded; resampled[si] = n_resampled;
      gates.row(si) = gate.t(); states[si] = st;
      ++si;
    }
    if ((it & 63) == 0) Rcpp::checkUserInterrupt();
    prog.increment();
  }

  Rcpp::List sn(ns);
  for (int s = 0; s < ns; ++s) sn[s] = states[s];
  return Rcpp::List::create(
    Rcpp::_["draws"] = draws, Rcpp::_["loglik"] = loglik, Rcpp::_["loglik_test"] = loglik_test,
    Rcpp::_["integrated_intensity"] = integrated, Rcpp::_["nleaf"] = nleaf,
    Rcpp::_["max_depth"] = maxdepth, Rcpp::_["ess"] = ess, Rcpp::_["as_rate"] = as_rate, Rcpp::_["expanded"] = expanded, Rcpp::_["resampled"] = resampled,
    Rcpp::_["gate"] = gates, Rcpp::_["gate_accept"] = gate_tot > 0 ? gate_acc / gate_tot : NA_REAL,
    Rcpp::_["state_nodes"] = sn, Rcpp::_["ndraws"] = ns);
}
