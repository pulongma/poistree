
#ifndef _USE_Armadillo
#define _USE_Armadillo
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#endif

#include <iomanip>       // std::setw / std::setprecision
#include <sstream>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <cmath>
#endif

using namespace Rcpp;


#include "PPT.h"


void PPT::split_node(int i, const arma::mat& pts, int axis, double cut, int min_leaf_n) {
    if (nodes[i] == nullptr || nodes[i]->is_empty) {
        // Rcpp::Rcout << "Node is null or empty, split skipped.\n";
        return;
    }

    arma::mat region = nodes[i]->region;
    arma::mat x = pts.rows(nodes[i]->idx);
    arma::mat left_region = region;
    left_region(axis, 1) = cut;
    arma::mat right_region = region;
    right_region(axis, 0) = cut;

    arma::uvec inside_left = in_region_nd(x, left_region);
    arma::uvec left_idx = nodes[i]->idx.elem(arma::find(inside_left));
    arma::uvec right_idx = nodes[i]->idx.elem(arma::find(inside_left == 0));
    int depth_child = nodes[i]->depth + 1;

    if (left_idx.n_elem < min_leaf_n || right_idx.n_elem < min_leaf_n) {
        return;
    }

    nodes[2 * i + 1] = new TreeNode(left_region, left_idx, depth_child, false, true);
    nodes[2 * i + 2] = new TreeNode(right_region, right_idx, depth_child, false, true);

    nodes[i]->is_leaf = false;
    nodes[i]->S = 1;
    nodes[i]->J = axis;
    nodes[i]->L = cut;
}



bool good_shape(const arma::mat& R, double lmin, double rmax)
{
    arma::vec side = R.col(1) - R.col(0);
    if (side.min() < lmin) return false;
    if (std::isinf(rmax) && rmax > 0) return true;
    if (!std::isfinite(rmax) || rmax < 1.0) return false;
    double asp = side.max() / side.min();
    return asp <= rmax;
}

std::vector<std::vector<double> > PPT::find_valid_cuts(
    const arma::mat& x,
    const arma::mat& region,
    bool force_mid_cut) 
{ 
  double buffer = 1e-3; 
  int n = x.n_rows, d = x.n_cols;
  std::vector<std::vector<double>> valid_cuts(d);

  arma::vec quantiles = arma::linspace(0.05, 0.95, cut_grid_n);

  for (int j = 0; j < d; ++j) {
    arma::vec xj = x.col(j);
    arma::rowvec axis_range = region.row(j);
    bool can_split = (n >= 2 * min_leaf_n);
    std::vector<double> out;

    if (!can_split) {
      valid_cuts[j] = out;
      continue;
    }

    if (force_mid_cut) {
      // Only allow mid-area cut
      double mid_cut = 0.5 * (axis_range(0) + axis_range(1));
      if ((mid_cut > axis_range(0) + buffer) && (mid_cut < axis_range(1) - buffer)) {
        arma::uvec left_idx  = arma::find(xj < mid_cut);
        arma::uvec right_idx = arma::find(xj >= mid_cut);

        arma::mat regionL = region;
        arma::mat regionR = region;
        regionL(j, 1) = mid_cut;
        regionR(j, 0) = mid_cut;

        double areaL = arma::prod(regionL.col(1) - regionL.col(0));
        double areaR = arma::prod(regionR.col(1) - regionR.col(0));

        if (good_shape(regionL, 1e-2, max_aspect_ratio) &&
            good_shape(regionR, 1e-2, max_aspect_ratio) &&
            left_idx.n_elem >= (unsigned)min_leaf_n &&
            right_idx.n_elem >= (unsigned)min_leaf_n &&
            areaL > 0 && areaR > 0) {
          out.push_back(mid_cut);
        }
      }
    } else {
      arma::vec uniq_xj = arma::unique(xj);
      uniq_xj = arma::sort(uniq_xj);
      arma::uvec keep = arma::find((uniq_xj > axis_range(0) + buffer) && (uniq_xj < axis_range(1) - buffer));
      arma::vec candidates = uniq_xj.elem(keep);

      if (candidates.n_elem >= 1) {
        arma::vec quantile_cuts(quantiles.n_elem);
        for (unsigned int q = 0; q < quantiles.n_elem; ++q) {
          quantile_cuts(q) = quantile_type1(candidates, quantiles(q));
        }
        quantile_cuts = arma::unique(quantile_cuts);

        for (unsigned int k = 0; k < quantile_cuts.n_elem; ++k) {
          double cut = quantile_cuts(k);
          arma::uvec left_idx  = arma::find(xj < cut);
          arma::uvec right_idx = arma::find(xj >= cut);

          arma::mat regionL = region;
          arma::mat regionR = region;
          regionL(j, 1) = cut;
          regionR(j, 0) = cut;

          if (!good_shape(regionL, 1e-2, max_aspect_ratio)) continue;
          if (!good_shape(regionR, 1e-2, max_aspect_ratio)) continue;

          double areaL = arma::prod(regionL.col(1) - regionL.col(0));
          double areaR = arma::prod(regionR.col(1) - regionR.col(0));
          if (left_idx.n_elem >= (unsigned)min_leaf_n &&
              right_idx.n_elem >= (unsigned)min_leaf_n &&
              areaL > 0 && areaR > 0) {
            out.push_back(cut);
          }
        }
      }

      // If no valid cut, try the median
      if (out.size() == 0) {
        double med_cut = arma::median(xj);
        if ((med_cut > axis_range(0) + buffer) && (med_cut < axis_range(1) - buffer)) {
          arma::uvec left_idx  = arma::find(xj < med_cut);
          arma::uvec right_idx = arma::find(xj >= med_cut);

          arma::mat regionL = region;
          arma::mat regionR = region;
          regionL(j, 1) = med_cut;
          regionR(j, 0) = med_cut;
          if (!good_shape(regionL, 1e-2, max_aspect_ratio)) continue;
          if (!good_shape(regionR, 1e-2, max_aspect_ratio)) continue;
          double areaL = arma::prod(regionL.col(1) - regionL.col(0));
          double areaR = arma::prod(regionR.col(1) - regionR.col(0));
          if (left_idx.n_elem >= (unsigned)min_leaf_n &&
              right_idx.n_elem >= (unsigned)min_leaf_n &&
              areaL > 0 && areaR > 0) {
            out.push_back(med_cut);
          }
        }
      }
    }
    valid_cuts[j] = out;
  }
  return valid_cuts;
}





// Convert a vector of pointers to TreeNode to an Rcpp::List
Rcpp::List PPT::to_R_list() {
    int n = nodes.size();
    Rcpp::List out(n);

    for (int i = 0; i < n; ++i) {
        if (nodes[i] == nullptr) {
            out[i] = R_NilValue; // missing node
        } else {
            Rcpp::List node = Rcpp::List::create(
                Rcpp::Named("region")   = nodes[i]->region,
                Rcpp::Named("idx")      = nodes[i]->idx + 1, // convert to 1-based R indices
                Rcpp::Named("depth")    = nodes[i]->depth,
                Rcpp::Named("is_empty") = nodes[i]->is_empty,
                Rcpp::Named("is_leaf")  = nodes[i]->is_leaf,
                Rcpp::Named("S")        = nodes[i]->S,
                Rcpp::Named("J")        = nodes[i]->J,
                Rcpp::Named("L")        = nodes[i]->L,
                Rcpp::Named("prob_split") = nodes[i]->prob_split,
                Rcpp::Named("prob_axis") = nodes[i]->prob_axis,
                Rcpp::Named("prob_cut") = nodes[i]->prob_cut,
                Rcpp::Named("lambda") = nodes[i]->lambda
            );
            out[i] = node;
        }
    }
    return out;
}



// Treed Poisson Process helper
/**********************************************************************************/


void PPT::PPT_one_step_ahead(double& log_inc, int i, 
    const arma::mat& pts)
{
    if (i < 0 || i >= (int)nodes.size() || nodes[i] == nullptr || nodes[i]->is_empty) {
        log_inc = 0.0;
        return;
    }
    log_inc = 0.0;
    TreeNode* node = nodes[i];
    arma::mat region = node->region;
    arma::uvec idx = node->idx;
    arma::mat x = pts.rows(idx);
    int d = x.n_cols;
    int n = idx.n_elem;

    arma::vec lamvec = lam * arma::ones(d);
    // arma::vec log_axis_prior = arma::log(lamvec);

    // depth-dependent split prior (Bayesian CART; Chipman et al. 1998, paper §2.2):
    //   P(split | depth d) = rho * (1 + d)^{-eta}
    // This regularises tree depth; with a constant split prob the tree over-grows
    // to max_depth and the piecewise-constant intensity becomes noisy.
    double rho_d = rho * std::pow(1.0 + node->depth, -eta);
    if (rho_d < 1e-12)       rho_d = 1e-12;
    if (rho_d > 1.0 - 1e-12) rho_d = 1.0 - 1e-12;

    // 1) No split marginal likelihood
    double area = arma::prod(region.col(1) - region.col(0));
    double loglik_nosplit; //
    loglik_nosplit = PPT_base_mloglik(n, area, a, b);
    // loglik_parent = PPT_base_mloglik(n, area, a, b); 
    // 2) Setup splits (use member function)
    std::vector<std::vector<double>> valid_cuts = find_valid_cuts(x, region, false);

    // 3) Count total cuts
    int total_cuts = 0;
    for (int j = 0; j < d; ++j) total_cuts += valid_cuts[j].size();

    if (total_cuts == 0) {
        node->S = 0;
        node->J = -1;
        node->L = NA_REAL;
        node->is_leaf = true;
        log_inc = 0.0;
        return;
    }

    // 4) Compute per-axis split marginal loglikelihood
    std::vector<arma::vec> loglik_splits(d);
    for (int axis = 0; axis < d; ++axis) {
        int K = valid_cuts[axis].size(); // # of cuts along a specific axis
        // if(K){lamvec(axis) = 0.0;}
        arma::vec loglik_per_cut(K);
        loglik_per_cut.fill(-arma::datum::inf);
        for (int j = 0; j < K; ++j) {
            double cut = valid_cuts[axis][j];
            arma::mat regionL = region; regionL(axis, 1) = cut;
            arma::mat regionR = region; regionR(axis, 0) = cut;
            arma::uvec inside_L = in_region_nd(x, regionL);
            int nL = arma::accu(inside_L);
            int nR = n - nL;
            if (nL < min_leaf_n || nR < min_leaf_n) continue;
            double areaL = arma::prod(regionL.col(1) - regionL.col(0));
            double areaR = arma::prod(regionR.col(1) - regionR.col(0));
            double mL = PPT_base_mloglik(nL, areaL, a, b);
            double mR = PPT_base_mloglik(nR, areaR, a, b);
            // loglik_per_cut(j) = loglik_nosplit - loglik_parent + mL + mR;
            loglik_per_cut(j) = mL + mR;
        }
        loglik_splits[axis] = loglik_per_cut;
    }

    lamvec /= arma::sum(lamvec);
    arma::vec log_axis_prior = arma::log(lamvec); 

    // 5) Collapse to per-axis marginal loglikelihood
    arma::vec loglik_per_axis(d); 
    loglik_per_axis.fill(-arma::datum::inf);
    for (int axis = 0; axis < d; ++axis) {
        arma::vec v = loglik_splits[axis];
        arma::uvec ok = arma::find_finite(v);
        if (ok.n_elem == 0) continue;
        double log_loc_prior = -std::log(ok.n_elem);
        arma::vec vals = v.elem(ok);
        loglik_per_axis(axis) = log_sum_exp(log_loc_prior + vals);
    }

    // 6) Posterior for S (depth-dependent split prior rho_d)
    arma::vec log_post_axis = log_axis_prior + loglik_per_axis;
    double log_post_split = std::log(rho_d) + log_sum_exp(log_post_axis);
    double log_post_nosplit = std::log(1.0 - rho_d) + loglik_nosplit;
    arma::vec log_post(2);
    log_post(0) = log_post_nosplit;
    log_post(1) = log_post_split;
    double Phi = log_sum_exp(log_post); // marginal likelihood
    arma::vec post_S = arma::exp(log_post - Phi);
    if (!post_S.is_finite() || arma::accu(post_S) == 0) { post_S = {1, 0}; }

    // 7) Sample S
    int S = R::rbinom(1, post_S(1));

    // Rcpp::Rcout<<"nonsplit="<<log_post_nosplit<<", split="<<log_post_split<<", S="<<S<<"\n";
    node->S = S;
    node->prob_split = post_S(1);

    // 8) If split
    arma::uvec isf = arma::find_finite(loglik_per_axis);
    if (S == 1 && isf.n_elem > 0) {
        // sample axis J
        arma::vec norm_log_post_axis = log_axis_prior + loglik_per_axis - log_sum_exp(log_post_axis);
        arma::vec post_axis = arma::exp(norm_log_post_axis);
        int J = Rcpp::sample(d, 1, false, Rcpp::NumericVector(Rcpp::wrap(post_axis)), false)[0];
        node->J = J;
        node->prob_axis = post_axis;

        // sample cut point L
        arma::vec prob_J = loglik_splits[J];
        arma::uvec valid_id = arma::find_finite(prob_J);
        if (valid_id.n_elem == 0) {
            node->S = 0;
            node->J = -1;
            node->L = NA_REAL;
            node->is_leaf = true;
            log_inc = 0.0;
            return;
        }
        arma::vec log_loc_prior(valid_id.n_elem); 
        log_loc_prior.fill(-std::log(valid_id.n_elem));
        arma::vec log_post_loc = log_loc_prior + prob_J.elem(valid_id);
        double norm_post_loc = log_sum_exp(log_post_loc);
        arma::vec post_loc = arma::exp(log_post_loc - norm_post_loc);
        int pick_id = Rcpp::sample(valid_id.n_elem, 1, false, Rcpp::NumericVector(Rcpp::wrap(post_loc)), false)[0];
        double L = valid_cuts[J][valid_id[pick_id]];
        node->L = L;
        node->prob_cut = post_loc;

        // --- Compute log prior and proposal for this action ---
        double log_prior = log(rho_d) + log_axis_prior(J) + log_loc_prior(pick_id);     // e.g., uniform = -log(#cuts)
        double log_prop = log(node->prob_split) + log(node->prob_axis(J)) + log(node->prob_cut(pick_id));

        log_inc += log_prior - log_prop;
        

        // cmopute log of Bayes factor 
        double logBF =  loglik_splits[J](valid_id[pick_id]) - loglik_nosplit;
        log_inc += logBF;
        // Use class member split_node
        this->split_node(i, pts, J, L, min_leaf_n);
    } else {
        node->is_leaf = true;
        node->S = 0;
        node->J = -1;
        node->L = NA_REAL;
        double log_prior = log(1.0-rho_d);
        double log_prop = log(post_S(0));
        log_inc = log_prior - log_prop;
    }
}


void PPT::PPT_force_reference_step(double& log_inc, int i,
                                   const PPT& ref_tree,
                                   const arma::mat& pts)
{
    log_inc = 0.0;

    const TreeNode* ref_node =
        (i >= 0 && i < static_cast<int>(ref_tree.nodes.size()))
            ? ref_tree.nodes[i] : nullptr;
    TreeNode* node =
        (i >= 0 && i < static_cast<int>(nodes.size())) ? nodes[i] : nullptr;

    // An absent heap node is a deterministic no-op.  Under a valid reference
    // prefix it must also be absent from the reference trajectory.
    if (node == nullptr || node->is_empty) {
        if (ref_node != nullptr && !ref_node->is_empty) {
            Rcpp::stop("Invalid conditional-SMC reference: incompatible node prefix");
        }
        return;
    }
    if (ref_node == nullptr || ref_node->is_empty) {
        Rcpp::stop("Invalid conditional-SMC reference: missing reference node");
    }

    arma::mat region = node->region;
    arma::uvec idx = node->idx;
    arma::mat x = pts.rows(idx);
    const int d = x.n_cols;
    const int n = idx.n_elem;

    std::vector<std::vector<double> > valid_cuts =
        find_valid_cuts(x, region, false);
    int total_cuts = 0;
    for (int axis = 0; axis < d; ++axis)
        total_cuts += static_cast<int>(valid_cuts[axis].size());

    // PPT_one_step_ahead treats a node with no valid split as a forced stop,
    // with proposal and target increment both equal to one.
    if (total_cuts == 0) {
        if (ref_node->S == 1) {
            Rcpp::stop("Invalid conditional-SMC reference: split is no longer valid");
        }
        node->S = 0;
        node->J = -1;
        node->L = NA_REAL;
        node->is_leaf = true;
        return;
    }

    const double log_q = PPT_log_transition_prob(node, ref_node, pts);
    if (!std::isfinite(log_q)) {
        Rcpp::stop("Invalid conditional-SMC reference action under the proposal law");
    }

    double rho_d = rho * std::pow(1.0 + node->depth, -eta);
    if (rho_d < 1e-12)       rho_d = 1e-12;
    if (rho_d > 1.0 - 1e-12) rho_d = 1.0 - 1e-12;

    if (ref_node->S == 0) {
        node->S = 0;
        node->J = -1;
        node->L = NA_REAL;
        node->is_leaf = true;
        log_inc = std::log(1.0 - rho_d) - log_q;
        return;
    }
    if (ref_node->S != 1) {
        Rcpp::stop("Invalid conditional-SMC reference action");
    }

    const int J = ref_node->J;
    const double L = ref_node->L;
    if (J < 0 || J >= d || !std::isfinite(L)) {
        Rcpp::stop("Invalid conditional-SMC reference split");
    }

    int cut_index = -1;
    int n_valid_J = 0;
    for (size_t k = 0; k < valid_cuts[J].size(); ++k) {
        arma::mat candidate_left = region;
        candidate_left(J, 1) = valid_cuts[J][k];
        const int candidate_nL = arma::accu(in_region_nd(x, candidate_left));
        const int candidate_nR = n - candidate_nL;
        if (candidate_nL < min_leaf_n || candidate_nR < min_leaf_n) continue;
        if (std::abs(valid_cuts[J][k] - L) < 1e-10)
            cut_index = static_cast<int>(k);
        ++n_valid_J;
    }
    if (cut_index < 0 || n_valid_J == 0) {
        Rcpp::stop("Invalid conditional-SMC reference cut");
    }

    arma::mat regionL = region;
    arma::mat regionR = region;
    regionL(J, 1) = L;
    regionR(J, 0) = L;
    const arma::uvec inside_L = in_region_nd(x, regionL);
    const int nL = arma::accu(inside_L);
    const int nR = n - nL;
    if (nL < min_leaf_n || nR < min_leaf_n) {
        Rcpp::stop("Invalid conditional-SMC reference child count");
    }

    const double area = arma::prod(region.col(1) - region.col(0));
    const double areaL = arma::prod(regionL.col(1) - regionL.col(0));
    const double areaR = arma::prod(regionR.col(1) - regionR.col(0));
    const double logBF =
        PPT_base_mloglik(nL, areaL, a, b) +
        PPT_base_mloglik(nR, areaR, a, b) -
        PPT_base_mloglik(n, area, a, b);

    // Axis selection is uniform after normalization, and cuts are uniform
    // within the selected axis, exactly as in PPT_one_step_ahead().
    const double log_prior =
        std::log(rho_d) - std::log(static_cast<double>(d)) -
        std::log(static_cast<double>(n_valid_J));
    log_inc = log_prior + logBF - log_q;
    split_node(i, pts, J, L, min_leaf_n);
}





double PPT::PPT_log_transition_prob(const TreeNode* parent_node,
                                   const TreeNode* ref_node,
                                   const arma::mat& pts)
{
    if (!parent_node || !ref_node) return -arma::datum::inf;

    // Get the same information as in your proposal step
    arma::mat region = parent_node->region;
    arma::uvec idx = parent_node->idx;
    arma::mat x = pts.rows(idx);
    int d = x.n_cols;
    int n = idx.n_elem;

    // ---- Recompute the axis priors and logliks as in proposal ----
    arma::vec lamvec = lam * arma::vec(d, arma::fill::ones);
    // Match PPT_one_step_ahead exactly.  In particular, replay must remain a
    // probability law even if the stored scalar axis weight is not 1 / d.
    lamvec /= arma::sum(lamvec);
    arma::vec log_axis_prior = arma::log(lamvec);

    // depth-dependent split prior, matching PPT_one_step_ahead
    double rho_d = rho * std::pow(1.0 + parent_node->depth, -eta);
    if (rho_d < 1e-12)       rho_d = 1e-12;
    if (rho_d > 1.0 - 1e-12) rho_d = 1.0 - 1e-12;

    double area = arma::prod(region.col(1) - region.col(0));
    double loglik_nosplit;
    // loglik_parent = PPT_base_mloglik(n, area, a, b);
    loglik_nosplit = PPT_base_mloglik(n, area, a, b);
    std::vector<std::vector<double>> valid_cuts = find_valid_cuts(x, region, false);

    // 3) Count total cuts
    // int total_cuts = 0;
    // for (int j = 0; j < d; ++j) total_cuts += valid_cuts[j].size();

    // 4) Compute per-axis split marginal loglikelihood
    std::vector<arma::vec> loglik_splits(d);
    for (int axis = 0; axis < d; ++axis) {
        int K = valid_cuts[axis].size();
        arma::vec loglik_per_cut(K); 
        loglik_per_cut.fill(-arma::datum::inf);
        for (int j = 0; j < K; ++j) {
            double cut = valid_cuts[axis][j];
            arma::mat regionL = region; regionL(axis, 1) = cut;
            arma::mat regionR = region; regionR(axis, 0) = cut;
            arma::uvec inside_L = in_region_nd(x, regionL);
            int nL = arma::accu(inside_L);
            int nR = n - nL;
            if (nL < min_leaf_n || nR < min_leaf_n) continue;
            double areaL = arma::prod(regionL.col(1) - regionL.col(0));
            double areaR = arma::prod(regionR.col(1) - regionR.col(0));
            double mL = PPT_base_mloglik(nL, areaL, a, b);
            double mR = PPT_base_mloglik(nR, areaR, a, b);
            // loglik_per_cut(j) = loglik_nosplit - loglik_parent + mL + mR;
            loglik_per_cut(j) = mL + mR;
        }
        loglik_splits[axis] = loglik_per_cut;
    }

    // 5) Collapse to per-axis marginal loglikelihood
    arma::vec loglik_per_axis(d); 
    loglik_per_axis.fill(-arma::datum::inf);
    for (int axis = 0; axis < d; ++axis) {
        arma::vec v = loglik_splits[axis];
        arma::uvec ok = arma::find_finite(v);
        if (ok.n_elem == 0) continue;
        double log_loc_prior = -std::log(ok.n_elem);
        arma::vec vals = v.elem(ok);
        loglik_per_axis(axis) = log_sum_exp(log_loc_prior + vals);
    }

    // 6) Posterior for S (depth-dependent split prior rho_d)
    arma::vec log_post_axis = log_axis_prior + loglik_per_axis;
    double log_post_split = std::log(rho_d) + log_sum_exp(log_post_axis);
    double log_post_nosplit = std::log(1.0 - rho_d) + loglik_nosplit;
    arma::vec log_post(2);
    log_post(0) = log_post_nosplit;
    log_post(1) = log_post_split;
    double Phi = log_sum_exp(log_post);// marginal likelihood
    arma::vec post_S = arma::exp(log_post - Phi);
    if (!post_S.is_finite() || arma::accu(post_S) == 0) { post_S = {1, 0}; }

    // Now: evaluate the **log-probability that S, J, L match the reference node**
    double logprob = 0.0;
    if (ref_node->S == 0) {
        // No split: take probability for S=0
        logprob = log_post(0) - Phi;
    } else if (ref_node->S == 1) {
        // Split: need axis J and cut L
        int J = ref_node->J;
        double L = ref_node->L;
        // (a) log p(S=1)
        logprob = log_post(1) - Phi;

        // (b) axis posterior (as in proposal)
        arma::vec norm_log_post_axis = log_axis_prior + loglik_per_axis - log_sum_exp(log_post_axis);
        arma::vec post_axis = arma::exp(norm_log_post_axis);
        if(J < 0 || J >= d) return -arma::datum::inf;
        logprob += std::log(post_axis[J]);

        // (c) cutpoint posterior
        arma::vec prob_J = loglik_splits[J];
        arma::uvec valid_id = arma::find_finite(prob_J);
        if (valid_id.n_elem == 0) return -arma::datum::inf;
        // Find which valid_id index matches L
        int match_idx = -1;
        for (unsigned int k = 0; k < valid_id.n_elem; ++k) {
            if (std::abs(valid_cuts[J][valid_id[k]] - L) < 1e-10) {
                match_idx = k; break;
            }
        }
        if (match_idx < 0) return -arma::datum::inf; // L not found among valid cuts
        arma::vec log_loc_prior(valid_id.n_elem);
        log_loc_prior.fill(-std::log(valid_id.n_elem));
        arma::vec log_post_loc = log_loc_prior + prob_J.elem(valid_id);
        // Normalize the same log masses used in the numerator.  Previously
        // the uniform 1 / K cut prior appeared only in the numerator, so the
        // replayed conditional cut probabilities summed to 1 / K rather than
        // one and biased the PGAS ancestor weights whenever K differed.
        double norm_post_loc = log_sum_exp(log_post_loc);
        arma::vec post_loc = arma::exp(log_post_loc - norm_post_loc);

        logprob += std::log(post_loc[match_idx]);
    } else {
        // S is not valid
        return -arma::datum::inf;
    }

    return logprob;
}


void PPT::PPT_draw_lambda()
{

    // Draw from the posterior matching PPT_base_mloglik:
    // lambda_leaf | x ~ Ga(a+n, b+area), with the b=0 branch interpreted as
    // the historical improper-prior limit.
    double temp = 0.0;
    for (auto* node : nodes) {
        if (node && node->is_leaf && !node->is_empty) {
            int n = node->idx.n_elem;        // number of points in this leaf
            double area = arma::prod(node->region.col(1) - node->region.col(0));
            node->lambda = R::rgamma(a + n, 1.0 / (b + area));
            temp += area * node->lambda;
        }
    }
    this->total_intensity = temp;

    return;

}

arma::vec PPT::predict_lambda(const arma::mat& XX) {
    int M = XX.n_rows;
    int d = XX.n_cols;
    arma::vec lambda_out(M, arma::fill::zeros);

    for (int i = 0; i < M; ++i) {
        arma::rowvec pt = XX.row(i);
        int idx = 0; // root
        while (idx < (int)nodes.size() && nodes[idx] && !nodes[idx]->is_leaf) {
            int axis = nodes[idx]->J;
            double cut = nodes[idx]->L;
            // Match the training assignment in in_region_nd(): the left child is
            // [.,cut) on every axis except the last, which is closed [.,cut].
            bool go_left = (axis < d - 1) ? (pt(axis) < cut) : (pt(axis) <= cut);
            idx = go_left ? (2*idx + 1) : (2*idx + 2);
        }
        if (idx < (int)nodes.size() && nodes[idx])
            lambda_out(i) = nodes[idx]->lambda;
        else
            lambda_out(i) = NA_REAL;
    }
    return lambda_out;
}

double PPT::PPT_get_lppd(const arma::vec& new_lambda){

    double lppd = - this->total_intensity + arma::sum(arma::log(new_lambda));
    return lppd;

}

double PPT::get_area(int id){
  if(id < 0 || id >= (int)nodes.size() || nodes[id] == nullptr){
    return 0;
  }
  return arma::prod(nodes[id]->region.col(1) - nodes[id]->region.col(0));
}

int PPT::get_counts(int id){
  if(id < 0 || id >= (int)nodes.size() || nodes[id] == nullptr){
    return 0;
  }
  return nodes[id]->idx.n_elem;
}
