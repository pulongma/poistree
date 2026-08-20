// ===========================================================================
// MCMCtree.cpp -- RcppArmadillo reversible-jump MCMC for the Poisson Process
// Tree model (comparison with the SMC). Leaf intensities lambda_l ~ Ga(a,b) are
// integrated out, so the chain moves over TREE space with target
//     pi(T|x) propto prior(T) * m(T),  m(T) = prod_leaf Q0(A_leaf),
//     Q0(A) = Gamma(n+a)/area^{n+a}    (Jeffreys a=0.5,b=0; general b>0 supported).
// Moves: GROW / PRUNE / CHANGE (Green 1995; Chipman-George-McCulloch 1998).
// Acceptance ratios validated against exact enumeration (see verification/).
//
//   Rcpp::sourceCpp("MCMCtree.cpp")        # standalone, or build as part of PPTree
// ===========================================================================
#ifndef _USE_Armadillo
#define _USE_Armadillo
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#endif
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
using namespace Rcpp;

// ---- leaf marginal Q0 (log) and depth-dependent split prior ----------------
static inline double mlogQ0(int n, double area, double a, double b) {
  if (b <= 0.0) return R::lgammafn(n + a) - (n + a) * std::log(area);      // Jeffreys/improper limit
  return R::lgammafn(a + n) - R::lgammafn(a) + a*std::log(b) - (a + n)*std::log(b + area);
}
static inline double mrho(int depth, double alpha, double eta) {
  double r = alpha * std::pow(1.0 + depth, -eta);
  if (r < 1e-12) r = 1e-12; if (r > 1.0 - 1e-12) r = 1.0 - 1e-12; return r;
}
static inline int runif_int(int k) { return (int)(R::unif_rand() * k); }     // 0..k-1

// ---- pointer-based tree node (so GROW/PRUNE can add/remove nodes) -----------
struct MNode {
  arma::mat box;       // d x 2 : [lower, upper] per axis
  arma::uvec idx;      // indices of points in this node
  int depth;
  MNode* L; MNode* R; int axis; double cut;
  MNode(const arma::mat& b, const arma::uvec& i, int d)
    : box(b), idx(i), depth(d), L(nullptr), R(nullptr), axis(-1), cut(NA_REAL) {}
  ~MNode() { if (L) delete L; if (R) delete R; }
  bool   leaf() const { return L == nullptr; }
  double area() const { return arma::prod(box.col(1) - box.col(0)); }
  int    n()    const { return (int)idx.n_elem; }
};
static void collect_leaves(MNode* nd, std::vector<MNode*>& out) {
  if (nd->leaf()) out.push_back(nd); else { collect_leaves(nd->L, out); collect_leaves(nd->R, out); }
}
static void collect_prunable(MNode* nd, std::vector<MNode*>& out) {           // internal, both children leaves
  if (!nd->leaf()) {
    if (nd->L->leaf() && nd->R->leaf()) out.push_back(nd);
    collect_prunable(nd->L, out); collect_prunable(nd->R, out);
  }
}
static int maximum_leaf_depth(MNode* root) {
  std::vector<MNode*> leaves;
  collect_leaves(root, leaves);
  int depth = 0;
  for (size_t k = 0; k < leaves.size(); ++k)
    depth = std::max(depth, leaves[k]->depth);
  return depth;
}

// ---- valid quantile cuts at a node (filtered by min_leaf & a small buffer) --
static inline double qtype1(const arma::vec& s, double p) {                   // type-1 quantile (matches PPTree)
  int n = s.n_elem; double h = (n - 1) * p + 1.0; int i = std::max(1, (int)std::floor(h)); return s(i - 1);
}
static std::vector< std::vector<double> >
mvalid_cuts(MNode* nd, const arma::mat& pts, int min_leaf, int cut_grid_n) {
  int d = nd->box.n_rows, n = nd->n();
  std::vector< std::vector<double> > out(d);
  if (n < 2 * min_leaf) return out;
  for (int j = 0; j < d; ++j) {
    arma::vec colj = pts.col(j); arma::vec s = arma::sort(colj.elem(nd->idx));
    double lo = nd->box(j,0), hi = nd->box(j,1); double last = arma::datum::nan;
    for (int q = 0; q < cut_grid_n; ++q) {
      double p = 0.05 + 0.90 * q / (double)(cut_grid_n - 1);
      double c = qtype1(s, p);
      if (c == last) continue; last = c;
      if (c > lo + 1e-3 && c < hi - 1e-3) {
        int nl = (int)(std::lower_bound(s.begin(), s.end(), c) - s.begin());  // # points < c
        if (nl >= min_leaf && n - nl >= min_leaf) out[j].push_back(c);
      }
    }
  }
  return out;
}
static inline bool has_valid(const std::vector< std::vector<double> >& vc) {
  for (size_t j = 0; j < vc.size(); ++j) if (!vc[j].empty()) return true; return false;
}
static inline bool msplittable(MNode* nd, const arma::mat& pts, int max_depth, int min_leaf, int cg) {
  return nd->leaf() && nd->depth < max_depth && has_valid(mvalid_cuts(nd, pts, min_leaf, cg));
}
static inline double mSfac(MNode* nd, const arma::mat& pts, int max_depth, int min_leaf, int cg,
                           double alpha, double eta) {
  if (nd->depth < max_depth && has_valid(mvalid_cuts(nd, pts, min_leaf, cg)))
    return 1.0 - mrho(nd->depth, alpha, eta);
  return 1.0;
}
static void make_children(MNode* nd, int j, double c, const arma::mat& pts, MNode*& Lo, MNode*& Ro) {
  arma::vec colj = pts.col(j); arma::vec xj = colj.elem(nd->idx);
  arma::uvec lm = arma::find(xj < c), rm = arma::find(xj >= c);
  arma::mat bl = nd->box; bl(j,1) = c; arma::mat br = nd->box; br(j,0) = c;
  Lo = new MNode(bl, nd->idx.elem(lm), nd->depth + 1);
  Ro = new MNode(br, nd->idx.elem(rm), nd->depth + 1);
}
// Draw leaf intensities lambda_l|x ~ Ga(a+n, b+area), evaluate them at the
// requested locations, and retain the exact point-process compensator and
// in-sample log likelihood for the same intensity draw.
static arma::vec predict_grid(MNode* root, const arma::mat& grid,
                              double a, double b, double& integral,
                              double& poisson_loglik) {
  int np = grid.n_rows, d = grid.n_cols; arma::vec out(np, arma::fill::zeros);
  std::vector<MNode*> lv; collect_leaves(root, lv);
  integral = 0.0;
  poisson_loglik = 0.0;
  for (size_t k = 0; k < lv.size(); ++k) {
    MNode* lf = lv[k]; double lam = R::rgamma(a + lf->n(), 1.0 / (b + lf->area()));
    integral += lam * lf->area();
    if (lf->n() > 0)
      poisson_loglik += lf->n() * std::log(std::max(lam, std::numeric_limits<double>::min()));
    arma::uvec inside(np, arma::fill::ones);                  // membership over all axes
    for (int j = 0; j < d; ++j) {
      arma::vec gj = grid.col(j);                             // materialise column (proven pattern)
      inside %= ((gj >= lf->box(j,0)) % (gj <= lf->box(j,1)));
    }
    arma::uvec in_idx = arma::find(inside);
    out.elem(in_idx).fill(lam);
  }
  poisson_loglik -= integral;
  return out;
}
static double tree_loglik(MNode* root, double a, double b) {
  std::vector<MNode*> lv; collect_leaves(root, lv); double s = 0;
  for (size_t k = 0; k < lv.size(); ++k) s += mlogQ0(lv[k]->n(), lv[k]->area(), a, b);
  return s;
}

// ===========================================================================
// [[Rcpp::export]]
Rcpp::List PPT_fit_MCMC(const arma::mat& pts, const arma::mat& grid, const arma::mat& region,
                        int niter = 4000, int burnin = 1000, int max_depth = 8, int min_leaf_n = 5,
                        int cut_grid_n = 30, double a = 0.5, double b = 0.0,
                        double alpha = 0.95, double eta = 2.0, int n_pred = 300) {
  int N = pts.n_rows, d = pts.n_cols, np = grid.n_rows;
  MNode* root = new MNode(region, arma::regspace<arma::uvec>(0, N - 1), 0);
  arma::ivec nl_trace(niter, arma::fill::zeros);
  arma::ivec md_trace(niter, arma::fill::zeros);
  arma::vec  ll_trace(niter, arma::fill::zeros);
  int pred_every = std::max(1, niter / std::max(1, n_pred));
  std::vector<arma::vec> preds;
  std::vector<double> integrals;
  std::vector<double> poisson_loglik;
  long acc_g=0, acc_p=0, acc_c=0, tot_g=0, tot_p=0, tot_c=0;

  for (int it = 0; it < niter + burnin; ++it) {
    double u = R::unif_rand();

    if (u < 1.0/3.0) {                                   // ---- GROW ----
      tot_g++;
      std::vector<MNode*> lv; collect_leaves(root, lv);
      std::vector<MNode*> B;
      for (size_t k=0;k<lv.size();++k) if (msplittable(lv[k],pts,max_depth,min_leaf_n,cut_grid_n)) B.push_back(lv[k]);
      if (!B.empty()) {
        MNode* ell = B[runif_int(B.size())]; int nB = B.size();
        std::vector< std::vector<double> > vc = mvalid_cuts(ell, pts, min_leaf_n, cut_grid_n);
        int j = runif_int(d);
        if (!vc[j].empty()) {
          double c = vc[j][runif_int(vc[j].size())];
          MNode *L, *R; make_children(ell, j, c, pts, L, R);
          double dm = mlogQ0(L->n(),L->area(),a,b) + mlogQ0(R->n(),R->area(),a,b) - mlogQ0(ell->n(),ell->area(),a,b);
          double SL = mSfac(L,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta);
          double SR = mSfac(R,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta);
          double rd = mrho(ell->depth, alpha, eta);
          double dp = std::log(rd) + std::log(SL) + std::log(SR) - std::log(1.0 - rd);
          ell->L=L; ell->R=R; ell->axis=j; ell->cut=c;                       // tentatively attach
          std::vector<MNode*> pr; collect_prunable(root, pr); int Wp = pr.size();
          double logA = dm + dp + std::log((double)nB) - std::log((double)Wp);
          if (std::log(R::unif_rand()) < logA) acc_g++;
          else { delete ell->L; delete ell->R; ell->L=nullptr; ell->R=nullptr; ell->axis=-1; }
        }
      }

    } else if (u < 2.0/3.0) {                            // ---- PRUNE ----
      tot_p++;
      std::vector<MNode*> pr; collect_prunable(root, pr);
      if (!pr.empty()) {
        MNode* v = pr[runif_int(pr.size())]; int nW = pr.size();
        double dm = mlogQ0(v->L->n(),v->L->area(),a,b) + mlogQ0(v->R->n(),v->R->area(),a,b) - mlogQ0(v->n(),v->area(),a,b);
        double SL = mSfac(v->L,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta);
        double SR = mSfac(v->R,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta);
        double rd = mrho(v->depth, alpha, eta);
        double dp = std::log(rd) + std::log(SL) + std::log(SR) - std::log(1.0 - rd);
        std::vector<MNode*> lv; collect_leaves(root, lv); int Bnow=0;         // splittable leaves now
        for (size_t k=0;k<lv.size();++k) if (msplittable(lv[k],pts,max_depth,min_leaf_n,cut_grid_n)) Bnow++;
        int sL = msplittable(v->L,pts,max_depth,min_leaf_n,cut_grid_n)?1:0;
        int sR = msplittable(v->R,pts,max_depth,min_leaf_n,cut_grid_n)?1:0;
        int sV = (v->depth<max_depth && has_valid(mvalid_cuts(v,pts,min_leaf_n,cut_grid_n)))?1:0;
        int Bp = Bnow - sL - sR + sV;                                         // splittable leaves after prune
        double logA = -dm - dp + std::log((double)nW) - std::log((double)Bp);
        if (std::log(R::unif_rand()) < logA) { delete v->L; delete v->R; v->L=nullptr; v->R=nullptr; v->axis=-1; acc_p++; }
      }

    } else {                                             // ---- CHANGE (terminal split) ----
      tot_c++;
      std::vector<MNode*> pr; collect_prunable(root, pr);
      if (!pr.empty()) {
        MNode* v = pr[runif_int(pr.size())];
        std::vector< std::vector<double> > vc = mvalid_cuts(v, pts, min_leaf_n, cut_grid_n);
        int j = runif_int(d);
        if (!vc[j].empty()) {
          double c = vc[j][runif_int(vc[j].size())];
          MNode *nL, *nR; make_children(v, j, c, pts, nL, nR);
          double dm = (mlogQ0(nL->n(),nL->area(),a,b)+mlogQ0(nR->n(),nR->area(),a,b))
                    - (mlogQ0(v->L->n(),v->L->area(),a,b)+mlogQ0(v->R->n(),v->R->area(),a,b));
          double dS = std::log(mSfac(nL,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta))
                    + std::log(mSfac(nR,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta))
                    - std::log(mSfac(v->L,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta))
                    - std::log(mSfac(v->R,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta));
          if (std::log(R::unif_rand()) < dm + dS) { delete v->L; delete v->R; v->L=nL; v->R=nR; v->axis=j; v->cut=c; acc_c++; }
          else { delete nL; delete nR; }
        }
      }
    }

    if (it >= burnin) {
      int t = it - burnin;
      std::vector<MNode*> lv; collect_leaves(root, lv);
      nl_trace(t) = (int)lv.size();
      md_trace(t) = maximum_leaf_depth(root);
      ll_trace(t) = tree_loglik(root, a, b);
      if (t % pred_every == 0) {
        double integral = 0.0, pp_loglik = 0.0;
        preds.push_back(predict_grid(root, grid, a, b, integral, pp_loglik));
        integrals.push_back(integral);
        poisson_loglik.push_back(pp_loglik);
      }
    }
    if ((it & 1023) == 0) Rcpp::checkUserInterrupt();
  }

  int K = preds.size(); arma::mat draws(np, K);
  for (int k = 0; k < K; ++k) draws.col(k) = preds[k];
  arma::vec integral_draws(K), poisson_loglik_draws(K);
  for (int k = 0; k < K; ++k) {
    integral_draws(k) = integrals[k];
    poisson_loglik_draws(k) = poisson_loglik[k];
  }
  arma::vec lam_mean = arma::mean(draws, 1), lam_med = arma::median(draws, 1);
  arma::vec lam_lo = arma::quantile(draws, arma::vec({0.025}), 1);
  arma::vec lam_hi = arma::quantile(draws, arma::vec({0.975}), 1);
  delete root;

  return Rcpp::List::create(
    Rcpp::_["lambda"] = Rcpp::List::create(Rcpp::_["mean"]=lam_mean, Rcpp::_["median"]=lam_med,
              Rcpp::_["lower95"]=lam_lo, Rcpp::_["upper95"]=lam_hi, Rcpp::_["draws"]=draws),
    Rcpp::_["nleaves"]=nl_trace, Rcpp::_["max_depth"]=md_trace,
    Rcpp::_["loglik"]=ll_trace,
    Rcpp::_["poisson_loglik"]=poisson_loglik_draws,
    Rcpp::_["integrated_intensity"]=integral_draws,
    Rcpp::_["accept"]=Rcpp::NumericVector::create(Rcpp::_["grow"]=(double)acc_g/std::max(1L,tot_g),
              Rcpp::_["prune"]=(double)acc_p/std::max(1L,tot_p), Rcpp::_["change"]=(double)acc_c/std::max(1L,tot_c)),
    Rcpp::_["niter"]=niter, Rcpp::_["burnin"]=burnin);
}
