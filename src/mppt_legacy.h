// =====================================================================
// MPPTree.cpp — Multidimensional Gamma-Markov (Smoothed Additive) Poisson
//               Process Tree, Rcpp/RcppArmadillo.  p-DIMENSIONAL extension
//               of MPTreeSmooth.cpp (SA-PPT).
// ---------------------------------------------------------------------
// Additive multiscale intensity over a box D = region (d-dimensional):
//     lambda(x) = sum_{v in path(x)} lambda_v,
// with a prior on the layer SCALES xi.  TWO priors are supported via `model`:
//   model = 1 ("markov") : Gamma-MARKOV coupling (the proposed dependent prior;
//                          smoothing across layers / partition boundaries)
//       xi_root ~ Gamma(a_xi, b_xi);  xi_v | xi_pa ~ Gamma(tau, tau/xi_pa)
//   model = 0 ("ind")    : INDEPENDENT scales (no dependence structure; the
//                          tau -> 0 special case of the proposed prior)
//       xi_v ~ Gamma(a_xi, b_xi)  i.i.d.  (tau unused)
//   both:  lambda_v | xi_v ~ Gamma(kappa, kappa/xi_v);  lambda(x)=sum_path lambda_v
// Under "ind" the lambda_v are independent across nodes; under "markov" siblings
// borrow strength through the shared parent scale.
// Collapsing lambda_v gives the same increment marginal as the additive
// model with a = kappa, b_v = kappa/xi_v (node-dependent):
//     g~(m,A;xi,kappa) = (kappa/xi)^kappa/Gamma(kappa)
//                        * Gamma(kappa+m)/(A+kappa/xi)^{kappa+m}.
//
// TREE PRIOR = Bayesian CART, exactly as in the PPTree package
// (PPT.cpp / PPTmodel.cpp), with
//   * random STOPPING   : P(split | depth d) = rho(d) = alpha (1+d)^{-eta}
//   * random SELECTION  : split axis  J ~ Uniform{1..d}         (weight 1/d)
//   * random SPLITTING  : cut point   L ~ Uniform{valid cuts on axis J}
// GROW proposes (J,L) FROM the prior (J = floor(d*U); L uniform over valid
// cuts on J), so the axis/cut prior and proposal cancel in the RJ ratio
// (same device as MCMCtree.cpp).  Children that cannot split contribute a
// no-split factor of 1 (RJ-consistent form of PPTree's "forced leaf when no
// valid cut").  Membership / routing: a point goes left iff x[J] < cut,
// matching PPT::predict_lambda.  d = 1 with region = [0,1] reproduces
// MPTreeSmooth.cpp.
//
// Self-contained (static helpers, no symbol clash with PPT.cpp / utils.cpp);
// region is d x 2, data / grid / test are n x d — matches PPT_fit_MCMC /
// PPT_fit_SMC, so it drops into the PPTree package after testing.
// Compile: Rcpp::sourceCpp("MPPTree.cpp")
// See ../sa_ppt_formulation.md and GMTree/README.md.
// [COMMON: Chipman-George-McCulloch 1998 (Bayesian CART); Kolaczyk 1999
//  (multiscale Poisson); Zhou-Carin 2015 (Gamma-Markov); Green 1995 (RJMCMC)].
// =====================================================================
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>
using namespace Rcpp;

// A soft-path entry is ignored by the legacy hard-box sampler.  MPPSTree.cpp
// uses it for its optional depth-relative recursive gate: each child records
// the split axis, cut, width of the parent on that axis, and side.
struct SoftGate {
  int axis;
  double cut;
  double parent_width;
  int side;          // -1 = left, +1 = right
};

// ---- node: axis-aligned box, index set, scale, split rule, label count -----
//   leaf  <=>  axis < 0  (cut = NA).  Children of id are 2*id and 2*id+1.
struct Node {
  arma::mat  box;   // d x 2 : [lower, upper] per axis
  arma::uvec idx;   // indices of the data points falling in this box
  double cut;       // split location on axis `axis` (NA for a leaf)
  double xi;        // this node's layer SCALE
  double area;      // prod(box.col(1) - box.col(0))  (cached)
  int axis;         // split axis J (>=0 internal, -1 leaf)
  int depth;
  int m;            // Polya-urn label count  m_v = #{ c_i = v }
  std::vector<SoftGate> soft_path; // root-to-node gates (soft-depth mode only)
};
typedef std::unordered_map<int,Node> Tree;

static inline double box_area(const arma::mat& box){ return arma::prod(box.col(1)-box.col(0)); }

// ---- collapsed increment marginal  log g~(m,A; xi, kappa) ------------------
static inline double lg(int m,double A,double kap,double xi){
  double b = kap/xi;                                   // node-dependent rate
  return kap*std::log(b) - R::lgammafn(kap) + R::lgammafn(kap+m) - (kap+m)*std::log(A+b);
}
// ---- depth-dependent split prior  rho(d) = alpha (1+d)^{-eta}  (clamped) ----
static inline double rho_(int d,double al,double eta){ double r=al*std::pow(1.0+d,-eta); return r>0.999999?0.999999:r; }
static inline double lrho (int d,double al,double eta){ return std::log(rho_(d,al,eta)); }
static inline double l1mrho(int d,double al,double eta){ return std::log(1.0-rho_(d,al,eta)); }
static inline int    runif_int(int k){ return (int)(unif_rand()*k); }        // 0..k-1
// ---- scale prior log-density  log Gamma(xi; shape, rate=shape/mean) --------
static inline double lgam_dens(double xi,double shape,double mean){
  double rate=shape/mean; return shape*std::log(rate)-R::lgammafn(shape)+(shape-1)*std::log(xi)-rate*xi; }
// ---- type-1 quantile (matches PPTree utils::quantile_type1 / MCMCtree) -----
static inline double qtype1(const arma::vec& s,double p){
  int n=s.n_elem; double h=(n-1)*p+1.0; int i=std::max(1,(int)std::floor(h)); return s(i-1); }

// ---- candidate cuts on ONE axis of a box (random-splitting support set) -----
//   mode: 0 = data midpoints, 1 = quantile grid (PPTree default), 2 = uniform.
static std::vector<double> axis_cuts(const arma::mat& pts,const arma::uvec& idx,
    const arma::mat& box,int axis,int mode,int ncand,int nmin){
  std::vector<double> cuts; double lo=box(axis,0), hi=box(axis,1), buf=1e-3;
  int sz=idx.n_elem, nm=nmin<1?1:nmin; if(sz<2*nm) return cuts;
  arma::vec xj=pts.col(axis); arma::vec s=arma::sort(xj.elem(idx));
  if(mode==2){ int M=ncand>30?ncand:30; double dx=(hi-lo)/(M+1); cuts.reserve(M);
    for(int j=1;j<=M;j++){ double c=lo+j*dx;
      int nl=(int)(std::lower_bound(s.begin(),s.end(),c)-s.begin());
      if(nl>=nm && sz-nl>=nm) cuts.push_back(c); }
    return cuts; }
  if(mode==0){                                          // data: midpoints, min-leaf filtered
    for(int j=nm-1;j<=sz-1-nm;j++){ if(!(s[j]<s[j+1])) continue;
      double c=0.5*(s[j]+s[j+1]);
      if(c>lo+buf && c<hi-buf){
        int nl=(int)(std::lower_bound(s.begin(),s.end(),c)-s.begin());
        if(nl>=nm && sz-nl>=nm) cuts.push_back(c); } }
    std::sort(cuts.begin(),cuts.end()); cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    return cuts;
  }
  int Q=ncand<2?2:ncand;                                // quantile grid on [0.05,0.95]
  for(int q=0;q<Q;q++){ double p=0.05+0.90*q/(double)(Q-1); double c=qtype1(s,p);
    if(c>lo+buf && c<hi-buf){ int nl=(int)(std::lower_bound(s.begin(),s.end(),c)-s.begin());
      if(nl>=nm && sz-nl>=nm) cuts.push_back(c); } }
  std::sort(cuts.begin(),cuts.end()); cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
  return cuts;
}
// ---- can this box be split? (depth, min-leaf, at least one valid cut) -------
static bool can_split(const arma::mat& box,const arma::uvec& idx,int depth,
    const arma::mat& pts,int Dmax,int nmin,int mode,int ncand){
  if(depth>=Dmax) return false; int d=box.n_rows, nm=nmin<1?1:nmin;
  if((int)idx.n_elem<2*nm) return false;
  for(int j=0;j<d;j++) if(!axis_cuts(pts,idx,box,j,mode,ncand,nmin).empty()) return true;
  return false;
}
// no-split prior factor for a leaf: log(1-rho_d) if splittable, else log(1)=0
static inline double lSfac(const arma::mat& box,const arma::uvec& idx,int depth,
    const arma::mat& pts,int Dmax,int nmin,int mode,int ncand,double al,double eta){
  return can_split(box,idx,depth,pts,Dmax,nmin,mode,ncand) ? l1mrho(depth,al,eta) : 0.0;
}

// ---- route a point to its leaf, collecting the root->leaf node ids ----------
static void path_ids(const Tree&T,const arma::rowvec& pt,std::vector<int>&p){
  p.clear(); int i=1; p.push_back(1);
  while(T.at(i).axis>=0){ int ax=T.at(i).axis; i=(pt(ax)<T.at(i).cut)?2*i:2*i+1; p.push_back(i); }
}
// ---- growable leaves / prunable internal nodes ------------------------------
static void growable(const Tree&T,const arma::mat& pts,int Dmax,int nmin,int mode,int ncand,std::vector<int>&G){
  G.clear(); for(const auto&kv:T){ const Node&nd=kv.second;
    if(nd.axis<0 && can_split(nd.box,nd.idx,nd.depth,pts,Dmax,nmin,mode,ncand)) G.push_back(kv.first); } }
static void prunable(const Tree&T,std::vector<int>&P){
  P.clear(); for(const auto&kv:T){ int id=kv.first; if(kv.second.axis>=0){
    auto iL=T.find(2*id),iR=T.find(2*id+1);
    if(iL!=T.end()&&iR!=T.end()&&iL->second.axis<0&&iR->second.axis<0) P.push_back(id); } } }

// ---- Polya-urn label sweep: P(c_i=v) prop (kappa+m_v)/(kappa/xi_v + A_v) ----
static void label_sweep(Tree&T,std::vector<int>&c,const arma::mat& pts,double kap){
  int n=pts.n_rows; std::vector<int> p; std::vector<double> w;
  for(int i=0;i<n;i++){
    arma::rowvec pt=pts.row(i); path_ids(T,pt,p); int old=c[i]; T[old].m-=1;
    w.resize(p.size()); double tot=0;
    for(size_t k=0;k<p.size();k++){ Node&nd=T[p[k]]; double wv=(kap+nd.m)/(kap/nd.xi+nd.area); w[k]=wv; tot+=wv; }
    double u=unif_rand()*tot, cum=0; int pick=p.back();
    for(size_t k=0;k<p.size();k++){ cum+=w[k]; if(u<=cum){ pick=p[k]; break; } }
    c[i]=pick; T[pick].m+=1;
  }
}

// ---- xi Metropolis sweep (one node at a time, RW on log xi) -----------------
//   model=1 (markov): target prop  p(xi_v|xi_pa) g~ prod_ch p(xi_B|xi_v)
//   model=0 (ind)   : target prop  Gamma(xi_v;a_xi,b_xi) g~   (node-local)
static int xi_sweep(Tree&T,double kap,double tau,double a_xi,double b_xi,double sd_xi,int model){
  int acc=0; double xi_mean = b_xi>0? a_xi/b_xi : 1.0;
  for(auto&kv:T){ int id=kv.first; Node&nd=kv.second; double A=nd.area;
    int L=2*id,R=2*id+1; bool internal=(nd.axis>=0);
    auto lt=[&](double xi)->double{
      double s = lg(nd.m, A, kap, xi);                            // own increment marginal
      if(model==1){                                               // Gamma-Markov coupling
        if(id==1) s += lgam_dens(xi, a_xi, xi_mean);              // root: Gamma(a_xi,b_xi)
        else      s += lgam_dens(xi, tau, T.at(id/2).xi);         // xi_v | xi_pa ~ Gamma(tau, mean=xi_pa)
        if(internal){ s += lgam_dens(T.at(L).xi, tau, xi) + lgam_dens(T.at(R).xi, tau, xi); }
      } else {                                                    // independent scales
        s += lgam_dens(xi, a_xi, xi_mean);                        // xi_v ~ iid Gamma(a_xi,b_xi)
      }
      return s;
    };
    double cur=nd.xi, prop=std::exp(std::log(cur)+R::rnorm(0.0,sd_xi));
    double la = lt(prop)-lt(cur) + (std::log(prop)-std::log(cur));  // + log-Jacobian of the log-RW
    if(std::log(unif_rand())<la){ nd.xi=prop; acc++; }
  }
  return acc;
}

// ---- kappa Metropolis (RW on log kappa); prior kappa~Gamma(a_k,b_k) ---------
static int kappa_update(Tree&T,double&kap,double a_k,double b_k,double sd_k){
  double cur=kap, prop=std::exp(std::log(cur)+R::rnorm(0.0,sd_k)); if(prop<=0) return 0;
  double dl=0.0;
  for(const auto&kv:T){ const Node&nd=kv.second; dl += lg(nd.m,nd.area,prop,nd.xi)-lg(nd.m,nd.area,cur,nd.xi); }
  dl += (a_k-1.0)*(std::log(prop)-std::log(cur)) - b_k*(prop-cur) + (std::log(prop)-std::log(cur));
  if(std::log(unif_rand())<dl){ kap=prop; return 1; } return 0;
}
// ---- tau Metropolis (RW on log tau); prior tau~Gamma(a_t,b_t) ---------------
static int tau_update(Tree&T,double&tau,double a_t,double b_t,double sd_t){
  double cur=tau, prop=std::exp(std::log(cur)+R::rnorm(0.0,sd_t)); if(prop<=0) return 0;
  double dl=0.0;
  for(const auto&kv:T){ int id=kv.first; if(id==1) continue; double xipa=T.at(id/2).xi;
    dl += lgam_dens(kv.second.xi,prop,xipa)-lgam_dens(kv.second.xi,cur,xipa); }
  dl += (a_t-1.0)*(std::log(prop)-std::log(cur)) - b_t*(prop-cur) + (std::log(prop)-std::log(cur));
  if(std::log(unif_rand())<dl){ tau=prop; return 1; } return 0;
}

// ---- GROW / PRUNE (rates collapsed; a=kappa, b=kappa/xi) --------------------
//   Bayesian-CART prior: axis J = floor(d*U) ~ Uniform{1..d} (SELECTION),
//   cut ~ Uniform{valid cuts on J} (SPLITTING); both cancel prior/proposal.
static int grow_prune(Tree&T,std::vector<int>&c,const arma::mat& pts,double kap,double tau,
                      double al,double eta,int Dmax,int nmin,int mode,int ncand,
                      double a_xi,double b_xi,int model,int&which_move){
  int d=pts.n_cols;
  std::vector<int> G,P; growable(T,pts,Dmax,nmin,mode,ncand,G); prunable(T,P);
  bool do_grow=unif_rand()<0.5; which_move=do_grow?0:1;
  if(do_grow){                                                      // ---- GROW ----
    if(G.empty()) return 0;
    int l=G[runif_int(G.size())]; int dep=T[l].depth; int L=2*l,Rn=2*l+1;
    int J=runif_int(d);                                             // random SELECTION
    std::vector<double> cuts=axis_cuts(pts,T[l].idx,T[l].box,J,mode,ncand,nmin);
    if(cuts.empty()) return 0;                                      // no valid cut on J -> null move
    double s=cuts[runif_int(cuts.size())];                          // random SPLITTING
    double xil=T[l].xi, Al=T[l].area; int mold=T[l].m;
    // child boxes / index sets (left iff x[J] < s), and label redistribution
    arma::vec colJ=pts.col(J);
    std::vector<arma::uword> Lidx,Ridx; std::vector<std::pair<int,int> > reass;
    int ml=0,mL=0,mR=0;
    for(arma::uword t=0;t<T[l].idx.n_elem;t++){ arma::uword pi=T[l].idx[t]; bool left=colJ[pi]<s;
      if(left) Lidx.push_back(pi); else Ridx.push_back(pi);
      if(c[pi]==l){ if(unif_rand()<0.5){ reass.push_back(std::make_pair((int)pi,l)); ml++; }
                    else if(left){ reass.push_back(std::make_pair((int)pi,L)); mL++; }
                    else         { reass.push_back(std::make_pair((int)pi,Rn)); mR++; } }
    }
    int Ssz=ml+mL+mR;                                               // = mold
    arma::uvec Lu(Lidx.size()), Ru(Ridx.size());
    for(size_t t=0;t<Lidx.size();t++) Lu[t]=Lidx[t];
    for(size_t t=0;t<Ridx.size();t++) Ru[t]=Ridx[t];
    arma::mat boxL=T[l].box; boxL(J,1)=s; arma::mat boxR=T[l].box; boxR(J,0)=s;
    double AL=box_area(boxL), AR=box_area(boxR);
    // propose child scales FROM THE PRIOR (cancels in the RJ ratio):
    //   markov -> Gamma(tau, mean=xi_l);   ind -> Gamma(a_xi, b_xi)
    double xiL,xiR;
    if(model==1){ xiL=R::rgamma(tau, xil/tau); xiR=R::rgamma(tau, xil/tau); }
    else        { double sc=b_xi>0?1.0/b_xi:1.0; xiL=R::rgamma(a_xi, sc); xiR=R::rgamma(a_xi, sc); }
    // tentatively attach children so prunable() sees the grown tree
    Node nl; nl.box=boxL; nl.idx=Lu; nl.cut=NA_REAL; nl.xi=xiL; nl.area=AL; nl.axis=-1; nl.depth=dep+1; nl.m=0;
    Node nr; nr.box=boxR; nr.idx=Ru; nr.cut=NA_REAL; nr.xi=xiR; nr.area=AR; nr.axis=-1; nr.depth=dep+1; nr.m=0;
    T[L]=nl; T[Rn]=nr; T[l].axis=J; T[l].cut=s;
    std::vector<int> Pn; prunable(T,Pn);
    double logSL=lSfac(boxL,Lu,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
    double logSR=lSfac(boxR,Ru,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
    double logA = lrho(dep,al,eta) - l1mrho(dep,al,eta) + logSL + logSR
                + lg(ml,Al,kap,xil) + lg(mL,AL,kap,xiL) + lg(mR,AR,kap,xiR) - lg(mold,Al,kap,xil)
                + std::log((double)G.size()) - std::log((double)Pn.size()) + Ssz*std::log(2.0);
    if(std::log(unif_rand())<logA){
      T[l].m=ml; T[L].m=mL; T[Rn].m=mR;
      for(size_t j=0;j<reass.size();j++) c[reass[j].first]=reass[j].second;
      return 1;
    } else { T.erase(L); T.erase(Rn); T[l].axis=-1; T[l].cut=NA_REAL; }
  } else {                                                          // ---- PRUNE ----
    if(P.empty()) return 0;
    int v=P[runif_int(P.size())]; int dep=T[v].depth; int L=2*v,Rn=2*v+1;
    int mv=T[v].m,mL=T[L].m,mR=T[Rn].m,M=mv+mL+mR;
    double Av=T[v].area,AL=T[L].area,AR=T[Rn].area, xiv=T[v].xi,xiL=T[L].xi,xiR=T[Rn].xi;
    double logSL=lSfac(T[L].box,T[L].idx,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
    double logSR=lSfac(T[Rn].box,T[Rn].idx,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
    Node oldL=T[L],oldR=T[Rn]; int Jold=T[v].axis; double oldcut=T[v].cut;
    T.erase(L); T.erase(Rn); T[v].axis=-1; T[v].cut=NA_REAL;
    std::vector<int> Gp; growable(T,pts,Dmax,nmin,mode,ncand,Gp);
    if(Gp.empty()){ T[L]=oldL; T[Rn]=oldR; T[v].axis=Jold; T[v].cut=oldcut; return 0; }
    double logA = -( lrho(dep,al,eta) - l1mrho(dep,al,eta) + logSL + logSR )
                + lg(M,Av,kap,xiv) - lg(mv,Av,kap,xiv) - lg(mL,AL,kap,xiL) - lg(mR,AR,kap,xiR)
                + std::log((double)P.size()) - std::log((double)Gp.size()) - M*std::log(2.0);
    if(std::log(unif_rand())<logA){
      for(size_t idx=0;idx<c.size();idx++) if(c[idx]==L||c[idx]==Rn) c[idx]=v;
      T[v].m=M;
      return 1;
    } else { T[L]=oldL; T[Rn]=oldR; T[v].axis=Jold; T[v].cut=oldcut; }
  }
  return 0;
}

// ---- CHANGE-CUT (repartition a prunable node; anchor scales fixed) ----------
static int change_cut(Tree&T,std::vector<int>&c,const arma::mat& pts,double kap,
                      double al,double eta,int Dmax,int nmin,int mode,int ncand){
  int d=pts.n_cols;
  std::vector<int> P; prunable(T,P); if(P.empty()) return 0;
  int v=P[runif_int(P.size())]; int dep=T[v].depth; int L=2*v,Rn=2*v+1;
  int J2=runif_int(d);
  std::vector<double> cuts=axis_cuts(pts,T[v].idx,T[v].box,J2,mode,ncand,nmin);
  if(cuts.empty()) return 0;
  double s2=cuts[runif_int(cuts.size())];
  double Av=T[v].area, xiv=T[v].xi, xiL=T[L].xi, xiR=T[Rn].xi;
  double oldv = lg(T[v].m,Av,kap,xiv) + lg(T[L].m,T[L].area,kap,xiL) + lg(T[Rn].m,T[Rn].area,kap,xiR);
  double oldS = lSfac(T[L].box,T[L].idx,dep+1,pts,Dmax,nmin,mode,ncand,al,eta)
              + lSfac(T[Rn].box,T[Rn].idx,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
  arma::vec colJ=pts.col(J2);
  std::vector<arma::uword> Lidx,Ridx; std::vector<std::pair<int,int> > reass;
  int nv=0,nL=0,nR=0;
  for(arma::uword t=0;t<T[v].idx.n_elem;t++){ arma::uword pi=T[v].idx[t]; bool left=colJ[pi]<s2;
    if(left) Lidx.push_back(pi); else Ridx.push_back(pi);
    int ci=c[pi];
    if(ci==v||ci==L||ci==Rn){ if(unif_rand()<0.5){ reass.push_back(std::make_pair((int)pi,v)); nv++; }
                              else if(left){ reass.push_back(std::make_pair((int)pi,L)); nL++; }
                              else         { reass.push_back(std::make_pair((int)pi,Rn)); nR++; } }
  }
  arma::uvec Lu(Lidx.size()), Ru(Ridx.size());
  for(size_t t=0;t<Lidx.size();t++) Lu[t]=Lidx[t];
  for(size_t t=0;t<Ridx.size();t++) Ru[t]=Ridx[t];
  arma::mat boxL=T[v].box; boxL(J2,1)=s2; arma::mat boxR=T[v].box; boxR(J2,0)=s2;
  double AL2=box_area(boxL), AR2=box_area(boxR);
  double newv = lg(nv,Av,kap,xiv) + lg(nL,AL2,kap,xiL) + lg(nR,AR2,kap,xiR);
  double newS = lSfac(boxL,Lu,dep+1,pts,Dmax,nmin,mode,ncand,al,eta)
              + lSfac(boxR,Ru,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
  if(std::log(unif_rand())<(newv-oldv)+(newS-oldS)){
    T[v].axis=J2; T[v].cut=s2;
    T[L].box=boxL; T[L].area=AL2; T[L].idx=Lu; T[L].m=nL;
    T[Rn].box=boxR;T[Rn].area=AR2;T[Rn].idx=Ru;T[Rn].m=nR;
    T[v].m=nv;
    for(size_t j=0;j<reass.size();j++) c[reass[j].first]=reass[j].second;
    return 1;
  }
  return 0;
}

static double logsumexp(const arma::vec&v){ if(v.n_elem==0) return -arma::datum::inf; double m=v.max(); return m+std::log(arma::accu(arma::exp(v-m))); }

// ---- posterior summaries (mean/median/2.5%/97.5%) from a draws matrix -------
static void summarize(const arma::mat&D,arma::vec&meanv,arma::vec&medv,arma::vec&lo,arma::vec&hi){
  arma::uword ng=D.n_cols; int S=D.n_rows; meanv.zeros(ng); medv.zeros(ng); lo.zeros(ng); hi.zeros(ng);
  if(S<=0) return; meanv=arma::mean(D,0).t(); medv=arma::median(D,0).t();
  arma::vec P={0.025,0.975}; for(arma::uword gq=0; gq<ng; gq++){ arma::vec q=arma::quantile(D.col(gq),P); lo[gq]=q(0); hi[gq]=q(1); }
}

// ---- build the root node covering `region` (d x 2) --------------------------
static Node make_root(const arma::mat& region,int n,double xi0){
  Node root; root.box=region; root.idx=arma::regspace<arma::uvec>(0,n-1);
  root.cut=NA_REAL; root.xi=xi0>0?xi0:1.0; root.area=box_area(region); root.axis=-1; root.depth=0; root.m=n;
  return root;
}

// Implemented after the soft-tree geometry in irjmcmc_moves.h.  Forward
// declarations let the common hard-chain driver preserve its update schedule
// while selecting the informed move kernel when requested.
static int irj_hard_grow_prune(Tree&tree,std::vector<int>&labels,
    const arma::mat&points,double kappa,double tau,double alpha,double eta,
    int max_depth,int min_leaf_n,int cut_mode,int cut_candidates,
    double a_xi,double b_xi,int model,int&which_move,
    bool conditional_labels);
static int irj_hard_change_cut(Tree&tree,std::vector<int>&labels,
    const arma::mat&points,double kappa,double alpha,double eta,
    int max_depth,int min_leaf_n,int cut_mode,int cut_candidates,
    bool conditional_labels);

// ---- run ONE chain into Draws rows [row0,row0+ns); updates kappa,tau if on ---
static int run_chain(const arma::mat& pts,const arma::mat& grid,const arma::mat& xt,const arma::mat& region,
                     double kappa0,double tau0,double a_xi,double b_xi,double sd_xi,
                     double al,double eta,int Dmax,int nmin,int iters,int burn,int thin,int nmove,int ncc,int cut_mode,int ncand,
                     int update_hyper,double a_k,double b_k,double sd_k,double a_t,double b_t,double sd_t,int model,
                     arma::mat&Draws,arma::vec&loglik,arma::vec&loglik_test,
                     arma::vec&integrated_intensity,int row0,
                     double&nleaf_mean,double&maxdepth_mean,
                     double&kap_mean,double&tau_mean,arma::vec&tree_accept,
                     std::vector<arma::mat>&state_nodes,
                     int proposal_mode=0){
  int n=pts.n_rows; arma::uword ng=grid.n_rows; int nt=xt.n_rows;
  double kappa=kappa0, tau=tau0, xi0=a_xi/std::max(b_xi,1e-12);
  Tree T; T[1]=make_root(region,n,xi0);
  std::vector<int> c(n,1), p;
  double nleaf_sum=0.0, maxdepth_sum=0.0, kap_sum=0.0, tau_sum=0.0;
  double grow_acc=0.0,grow_tot=0.0,prune_acc=0.0,prune_tot=0.0;
  double change_acc=0.0,change_tot=0.0;
  int si=0;
  for(int it=0; it<iters; it++){
    label_sweep(T,c,pts,kappa);
    xi_sweep(T,kappa,tau,a_xi,b_xi,sd_xi,model);
    if(update_hyper){ kappa_update(T,kappa,a_k,b_k,sd_k); if(model==1) tau_update(T,tau,a_t,b_t,sd_t); }
    for(int r=0;r<nmove;r++){
      int which=-1;
      int accepted=proposal_mode==0
        ? grow_prune(T,c,pts,kappa,tau,al,eta,Dmax,nmin,cut_mode,
                     ncand,a_xi,b_xi,model,which)
        : irj_hard_grow_prune(
            T,c,pts,kappa,tau,al,eta,Dmax,nmin,cut_mode,ncand,a_xi,b_xi,
            model,which,proposal_mode==3
          );
      if(which==0){ grow_acc+=accepted; grow_tot+=1.0; }
      else { prune_acc+=accepted; prune_tot+=1.0; }
    }
    for(int r=0;r<ncc;r++){
      change_acc+=proposal_mode==0
        ? change_cut(T,c,pts,kappa,al,eta,Dmax,nmin,cut_mode,ncand)
        : irj_hard_change_cut(
            T,c,pts,kappa,al,eta,Dmax,nmin,cut_mode,ncand,
            proposal_mode==3
          );
      change_tot+=1.0;
    }
    if(it>=burn && (it-burn)%thin==0){
      int row=row0+si;
      std::unordered_map<int,double> lam; lam.reserve(T.size()); double comp=0.0;
      for(const auto&kv:T){ double lv=R::rgamma(kappa+kv.second.m,1.0/(kappa/kv.second.xi+kv.second.area)); lam[kv.first]=lv; comp+=lv*kv.second.area; }
      integrated_intensity[row]=comp;
      for(arma::uword gq=0; gq<ng; gq++){ arma::rowvec pt=grid.row(gq); path_ids(T,pt,p); double sVal=0.0; for(size_t k=0;k<p.size();k++) sVal+=lam[p[k]]; Draws(row,gq)=sVal; }
      double ll=-comp; for(int i=0;i<n;i++){ arma::rowvec pt=pts.row(i); path_ids(T,pt,p); double sVal=0.0; for(size_t k=0;k<p.size();k++) sVal+=lam[p[k]]; ll+=std::log(sVal);} loglik[row]=ll;
      if(nt>0){ double llt=-comp; for(int j=0;j<nt;j++){ arma::rowvec pt=xt.row(j); path_ids(T,pt,p); double sVal=0.0; for(size_t k=0;k<p.size();k++) sVal+=lam[p[k]]; llt+=std::log(sVal);} loglik_test[row]=llt; }
      int nl=0,md=0;
      for(const auto&kv:T){
        if(kv.second.axis<0) nl++;
        if(kv.second.depth>md) md=kv.second.depth;
      }
      // serialize the generative state of this retained draw: one row per
      // node, columns (heap id, axis, cut, lambda, xi, m); every node of
      // the additive multiscale tree carries a rate.
      {
        arma::mat st(T.size(),6);
        int rr=0;
        for(const auto&kv:T){
          const Node&nd=kv.second;
          st(rr,0)=(double)kv.first;
          st(rr,1)=(double)nd.axis;
          st(rr,2)=nd.axis>=0?nd.cut:NA_REAL;
          st(rr,3)=lam[kv.first];
          st(rr,4)=nd.xi;
          st(rr,5)=(double)nd.m;
          rr++;
        }
        state_nodes[row]=st;
      }
      nleaf_sum+=nl; maxdepth_sum+=md; kap_sum+=kappa; tau_sum+=tau; si++;
    }
    if((it & 1023)==0) Rcpp::checkUserInterrupt();
  }
  nleaf_mean=nleaf_sum/std::max(1,si);
  maxdepth_mean=maxdepth_sum/std::max(1,si);
  kap_mean=kap_sum/std::max(1,si);
  tau_mean=tau_sum/std::max(1,si);
  tree_accept.set_size(3);
  tree_accept[0]=grow_acc/std::max(1.0,grow_tot);
  tree_accept[1]=prune_acc/std::max(1.0,prune_tot);
  tree_accept[2]=change_acc/std::max(1.0,change_tot);
  return si;
}

// =====================================================================
// mpptree_chain — a SINGLE p-dim SA-PPT chain (posterior mean/median + 95%
// band + log-lik/logpred + posterior-mean kappa,tau).  update_hyper=1 turns
// on the kappa/tau Metropolis updates.
//   X  : n x d data     grid : ng x d query    Xtest : nt x d (0 rows to skip)
//   region : d x 2 bounding box
//   hp = c(kappa, tau, a_xi, b_xi, sd_xi, a_k, b_k, sd_k, a_t, b_t, sd_t)
// =====================================================================
// [[Rcpp::export]]
List mpptree_chain(arma::mat X, arma::mat grid, arma::mat Xtest, arma::mat region, arma::vec hp, int update_hyper,
                   int model, double al, double eta, int Dmax, int nmin,
                   int iters, int burn, int thin, int nmove, int ncc, int cut_mode, int ncand){
  GetRNGstate();
  double kappa=hp[0],tau=hp[1],a_xi=hp[2],b_xi=hp[3],sd_xi=hp[4],a_k=hp[5],b_k=hp[6],sd_k=hp[7],a_t=hp[8],b_t=hp[9],sd_t=hp[10];
  arma::uword ng=grid.n_rows; int nt=Xtest.n_rows;
  int nsave=0; for(int it=burn; it<iters; it++) if((it-burn)%thin==0) nsave++;
  arma::mat Draws(nsave,ng,arma::fill::zeros); arma::vec loglik(nsave,arma::fill::zeros), loglik_test(nt>0?nsave:0,arma::fill::zeros);
  arma::vec integrated_intensity(nsave,arma::fill::zeros);
  double nlm,mdm,km,tm; arma::vec tree_accept;
  std::vector<arma::mat> state_nodes(nsave);
  int S=run_chain(X,grid,Xtest,region,kappa,tau,a_xi,b_xi,sd_xi,al,eta,Dmax,nmin,iters,burn,thin,nmove,ncc,cut_mode,ncand,
                  update_hyper,a_k,b_k,sd_k,a_t,b_t,sd_t,model,Draws,loglik,
                  loglik_test,integrated_intensity,0,nlm,mdm,km,tm,tree_accept,
                  state_nodes);
  PutRNGstate();
  arma::vec meanv,medv,lo,hi; summarize(Draws,meanv,medv,lo,hi);
  double logpred=(nt>0&&S>0)?(logsumexp(loglik_test)-std::log((double)S)):NA_REAL;
  double loglik_mean=(S>0)?arma::mean(loglik):NA_REAL;
  List sn(S);
  for(int s=0;s<S;s++) sn[s]=state_nodes[s];
  List out=List::create(_["mean"]=meanv,_["median"]=medv,_["lower95"]=lo,_["upper95"]=hi,
                      _["loglik_mean"]=loglik_mean,_["logpred"]=logpred,
                      _["integrated_intensity"]=integrated_intensity,
                      _["kappa_mean"]=km,_["tau_mean"]=tm,
                      _["tree_accept"]=tree_accept,
                      _["ns"]=S,_["nleaf"]=nlm,_["mean_max_depth"]=mdm);
  out["state_nodes"]=sn;
  return out;
}

// =====================================================================
// mpptree_multi — run `chains` chains IN C++, pool ALL recorded draws
// transiently, return EXACT pooled mean/median/95% band + pooled logpred +
// posterior-mean kappa,tau.  One R set.seed() reproduces the run.
// =====================================================================
// [[Rcpp::export]]
List mpptree_multi(arma::mat X, arma::mat grid, arma::mat Xtest, arma::mat region, arma::vec hp, int update_hyper,
                   int model, double al, double eta, int Dmax, int nmin,
                   int iters, int burn, int thin, int nmove, int ncc, int cut_mode, int ncand,
                   int chains, int verbose){
  GetRNGstate();
  double kappa=hp[0],tau=hp[1],a_xi=hp[2],b_xi=hp[3],sd_xi=hp[4],a_k=hp[5],b_k=hp[6],sd_k=hp[7],a_t=hp[8],b_t=hp[9],sd_t=hp[10];
  arma::uword ng=grid.n_rows; int nt=Xtest.n_rows;
  int nsave=0; for(int it=burn; it<iters; it++) if((it-burn)%thin==0) nsave++;
  int total=nsave*chains;
  arma::mat Draws(total,ng,arma::fill::zeros); arma::vec loglik(total,arma::fill::zeros), loglik_test(nt>0?total:0,arma::fill::zeros);
  arma::vec integrated_intensity(total,arma::fill::zeros);
  arma::vec nleaf_ch(chains,arma::fill::zeros), maxdepth_ch(chains,arma::fill::zeros),
    kap_ch(chains,arma::fill::zeros), tau_ch(chains,arma::fill::zeros); int row=0;
  arma::mat tree_accept_ch(chains,3,arma::fill::zeros);
  std::vector<arma::mat> state_nodes(total);
  if(verbose) Rcpp::Rcout<<(model==1?"MPPT":"MPPT-Ind")<<" [RJ-MCMC]: "
                          <<chains<<" chains, "<<nsave<<" draws/chain\n";
  for(int k=0;k<chains;k++){ double nlm,mdm,km,tm; arma::vec tree_accept;
    int got=run_chain(X,grid,Xtest,region,kappa,tau,a_xi,b_xi,sd_xi,al,eta,Dmax,nmin,iters,burn,thin,nmove,ncc,cut_mode,ncand,
                      update_hyper,a_k,b_k,sd_k,a_t,b_t,sd_t,model,Draws,loglik,
                      loglik_test,integrated_intensity,row,nlm,mdm,km,tm,
                      tree_accept,state_nodes);
    nleaf_ch[k]=nlm; maxdepth_ch[k]=mdm; kap_ch[k]=km; tau_ch[k]=tm; row+=got;
    tree_accept_ch.row(k)=tree_accept.t();
    if(verbose){
      Rcpp::Rcout<<"  chain "<<k+1<<"/"<<chains<<" done; leaves="<<nlm
                 <<", max_depth="<<mdm<<", kappa="<<km;
      if(model==1) Rcpp::Rcout<<", tau="<<tm;
      Rcpp::Rcout<<"\n";
    }
    Rcpp::checkUserInterrupt();
  }
  PutRNGstate();
  int S=row; arma::vec meanv,medv,lo,hi; summarize(Draws,meanv,medv,lo,hi);
  double logpred=(nt>0&&S>0)?(logsumexp(loglik_test)-std::log((double)S)):NA_REAL;
  double loglik_mean=(S>0)?arma::mean(loglik):NA_REAL;
  List sn(S);
  for(int s=0;s<S;s++) sn[s]=state_nodes[s];
  List out=List::create(_["mean"]=meanv,_["median"]=medv,_["lower95"]=lo,_["upper95"]=hi,
                      _["loglik_mean"]=loglik_mean,_["logpred"]=logpred,
                      _["integrated_intensity"]=integrated_intensity,
                      _["kappa_mean"]=arma::mean(kap_ch),_["tau_mean"]=arma::mean(tau_ch),
                      _["mean_leaves"]=arma::mean(nleaf_ch),
                      _["mean_max_depth"]=arma::mean(maxdepth_ch),
                      _["tree_accept"]=arma::mean(tree_accept_ch,0).t(),
                      _["ndraws"]=S);
  out["state_nodes"]=sn;
  return out;
}

// =====================================================================
// mpptree_diag — ONE chain returning per-draw TRACES for convergence checks:
//   nleaf, kappa, tau, logdens (training log-lik), and lambda at monitor
//   points `mon` (nm x d) — dimension-invariant functionals comparable across
//   chains with different trees.  Also returns MH acceptance rates
//   (xi / kappa / tau).  [COMMON: Gelman-Rubin R-hat & ESS on scalar
//   functionals; RW-MH tuned to ~0.234 optimal acceptance, Roberts et al. 1997]
// =====================================================================
// [[Rcpp::export]]
List mpptree_diag(arma::mat X, arma::mat mon, arma::mat region, arma::vec hp, int update_hyper,
                  int model, double al, double eta, int Dmax, int nmin,
                  int iters, int burn, int thin, int nmove, int ncc, int cut_mode, int ncand){
  GetRNGstate();
  double kappa=hp[0],tau=hp[1],a_xi=hp[2],b_xi=hp[3],sd_xi=hp[4],a_k=hp[5],b_k=hp[6],sd_k=hp[7],a_t=hp[8],b_t=hp[9],sd_t=hp[10];
  int n=X.n_rows; arma::uword nm=mon.n_rows;
  double xi0=a_xi/std::max(b_xi,1e-12);
  Tree T; T[1]=make_root(region,n,xi0);
  std::vector<int> c(n,1), p;
  int nsave=0; for(int it=burn; it<iters; it++) if((it-burn)%thin==0) nsave++;
  arma::vec tr_nleaf(nsave,arma::fill::zeros), tr_kappa(nsave,arma::fill::zeros),
            tr_tau(nsave,arma::fill::zeros), tr_ldens(nsave,arma::fill::zeros);
  arma::mat tr_mon(nsave, nm, arma::fill::zeros);
  double xi_acc=0, xi_tot=0, k_acc=0, k_tot=0, t_acc=0, t_tot=0;
  double grow_acc=0,grow_tot=0,prune_acc=0,prune_tot=0,change_acc=0,change_tot=0;
  int si=0;
  for(int it=0; it<iters; it++){
    label_sweep(T,c,X,kappa);
    xi_acc += xi_sweep(T,kappa,tau,a_xi,b_xi,sd_xi,model); xi_tot += (double)T.size();
    if(update_hyper){ k_acc+=kappa_update(T,kappa,a_k,b_k,sd_k); k_tot+=1.0;
                      if(model==1){ t_acc+=tau_update(T,tau,a_t,b_t,sd_t); t_tot+=1.0; } }
    for(int r=0;r<nmove;r++){
      int which=-1;
      int accepted=grow_prune(T,c,X,kappa,tau,al,eta,Dmax,nmin,cut_mode,
                              ncand,a_xi,b_xi,model,which);
      if(which==0){ grow_acc+=accepted; grow_tot++; }
      else { prune_acc+=accepted; prune_tot++; }
    }
    for(int r=0;r<ncc;r++){
      change_acc+=change_cut(T,c,X,kappa,al,eta,Dmax,nmin,cut_mode,ncand);
      change_tot++;
    }
    if(it>=burn && (it-burn)%thin==0){
      std::unordered_map<int,double> lam; lam.reserve(T.size()); double comp=0.0;
      for(const auto&kv:T){ double lv=R::rgamma(kappa+kv.second.m,1.0/(kappa/kv.second.xi+kv.second.area)); lam[kv.first]=lv; comp+=lv*kv.second.area; }
      for(arma::uword j=0;j<nm;j++){ arma::rowvec pt=mon.row(j); path_ids(T,pt,p); double sVal=0.0; for(size_t k=0;k<p.size();k++) sVal+=lam[p[k]]; tr_mon(si,j)=sVal; }
      double ll=-comp; for(int i=0;i<n;i++){ arma::rowvec pt=X.row(i); path_ids(T,pt,p); double sVal=0.0; for(size_t k=0;k<p.size();k++) sVal+=lam[p[k]]; ll+=std::log(sVal); }
      int nl=0; for(const auto&kv:T) if(kv.second.axis<0) nl++;
      tr_nleaf[si]=nl; tr_kappa[si]=kappa; tr_tau[si]=tau; tr_ldens[si]=ll; si++;
    }
    if((it & 1023)==0) Rcpp::checkUserInterrupt();
  }
  PutRNGstate();
  arma::vec acc(3);
  acc[0]=xi_acc/std::max(1.0,xi_tot); acc[1]=k_acc/std::max(1.0,k_tot); acc[2]=t_acc/std::max(1.0,t_tot);
  arma::vec tree_accept(3);
  tree_accept[0]=grow_acc/std::max(1.0,grow_tot);
  tree_accept[1]=prune_acc/std::max(1.0,prune_tot);
  tree_accept[2]=change_acc/std::max(1.0,change_tot);
  return List::create(_["nleaf"]=tr_nleaf,_["kappa"]=tr_kappa,_["tau"]=tr_tau,_["logdens"]=tr_ldens,
                      _["mon"]=tr_mon,_["accept"]=acc,
                      _["tree_accept"]=tree_accept,_["ns"]=si);
}
