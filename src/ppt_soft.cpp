// ============================================================================
// PPSTree.cpp -- Soft Poisson Process Tree (terminal-leaf PPT) using
// RcppArmadillo and reversible-jump MCMC.
//
// Only terminal leaves carry intensities,
//
//   lambda(x) = sum_{ell in leaves(T)} lambda_ell phi_ell(x),
//
// where the default recursive logistic gates give nonnegative C-infinity leaf
// bases satisfying sum_ell phi_ell(x)=1.  At a split c on coordinate j,
//
//   G_c(x)=logit^{-1}{gate_j (x_j-c)/(b_j-a_j)},
//   phi_left=phi_parent*(1-G_c), phi_right=phi_parent*G_c.
//
// Thus children sum exactly to their parent.  With one shared slope per
// dimension, the leaf exposure H_ell=int phi_ell is a product of analytic
// one-dimensional rational integrals.  Optional node-relative logistic gates
// instead use each split parent's width, with adaptive one-dimensional
// integration for varying slopes along a path.  The former compact-cubic gate remains
// available as an opt-in compatibility family.  With lambda_ell~Gamma(a,b),
// latent leaf labels give
//
//   lambda_ell | ... ~ Gamma(a+m_ell, b+H_ell).
//
// Grow/prune/change-cut moves propose affected labels from q_left:q_right.
// Their spatial factors cancel the corresponding label proposal probability,
// leaving Gamma-Poisson block, CART-prior, and move-count ratios.
// The wrapper defaults to one shared gate.  Optionally each dimension has its
// own gate_j, with a (possibly lower-truncated) Gamma(shape,rate) prior and a
// coordinate-wise log-RW MH update.
//
// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
// [[Rcpp::plugins(cpp17)]]
// ============================================================================
#include <RcppArmadillo.h>
#include "tree_limits.h"
#include "mcmc_progress.h"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include "soft_logistic.h"
using namespace Rcpp;

struct PPSTGate {
  int axis;
  double cut;
  double parent_width;
  int side;                         // -1 left, +1 right
};

struct PPSTNode {
  arma::mat box;
  arma::uvec idx;                   // hard routing indices, used for cut support
  std::vector<PPSTGate> path;       // recursive soft path
  double cut;
  int axis;
  int depth;
  int m;                            // soft leaf-label count
};
typedef std::unordered_map<int,PPSTNode> PPSTree;

static inline int ppst_runif_int(int k){
  return (int)(R::unif_rand()*k);
}
static arma::vec ppst_expand_positive(const arma::vec&x,arma::uword d,
    const char*name,bool allow_zero=false){
  arma::vec out;
  if(x.n_elem==1) out=arma::vec(d,arma::fill::value(x[0]));
  else if(x.n_elem==d) out=x;
  else stop(std::string(name)+" must have length 1 or dimension d");
  if(!out.is_finite()||
     (allow_zero?arma::any(out<0.0):arma::any(out<=0.0)))
    stop(std::string(name)+(allow_zero?" must be finite and nonnegative":
                                      " must be finite and positive"));
  return out;
}
static inline double ppst_rho(int depth,double alpha,double eta){
  double r=alpha*std::pow(1.0+depth,-eta);
  return std::min(1.0-1e-12,std::max(1e-12,r));
}
static inline double ppst_log_g(int m,double H,double a,double b){
  return R::lgammafn(a+m)-R::lgammafn(a)+a*std::log(b)
       -(a+m)*std::log(b+H);
}
static inline double ppst_logsumexp(const std::vector<double>&x){
  double mx=-std::numeric_limits<double>::infinity();
  for(double z:x) if(z>mx) mx=z;
  if(!std::isfinite(mx)) return mx;
  double s=0.0; for(double z:x) s+=std::exp(z-mx);
  return mx+std::log(s);
}
// ---- candidate cuts and CART support ---------------------------------------
static inline double ppst_qtype1(const arma::vec&s,double p){
  int n=s.n_elem; double h=(n-1)*p+1.0;
  int i=std::max(1,(int)std::floor(h)); return s(i-1);
}

static std::vector<double> ppst_axis_cuts(const arma::mat&pts,
    const arma::uvec&idx,const arma::mat&box,int axis,int mode,int ncand,
    int nmin){
  std::vector<double> cuts;
  double lo=box(axis,0),hi=box(axis,1),buf=1e-3;
  int n=idx.n_elem,nm=std::max(1,nmin);
  if(n<2*nm) return cuts;
  arma::vec col=pts.col(axis),s=arma::sort(col.elem(idx));
  if(mode==2){
    int M=std::max(30,ncand); double dx=(hi-lo)/(M+1.0);
    cuts.reserve(M);
    for(int j=1;j<=M;j++){
      double c=lo+j*dx;
      int nl=(int)(std::lower_bound(s.begin(),s.end(),c)-s.begin());
      if(nl>=nm&&n-nl>=nm) cuts.push_back(c);
    }
    return cuts;
  }
  if(mode==0){
    for(int j=nm-1;j<=n-1-nm;j++){
      if(!(s[j]<s[j+1])) continue;
      double c=0.5*(s[j]+s[j+1]);
      if(c>lo+buf&&c<hi-buf){
        int nl=(int)(std::lower_bound(s.begin(),s.end(),c)-s.begin());
        if(nl>=nm&&n-nl>=nm) cuts.push_back(c);
      }
    }
  }else{
    int Q=std::max(2,ncand);
    for(int q=0;q<Q;q++){
      double p=0.05+0.90*q/(double)(Q-1),c=ppst_qtype1(s,p);
      if(c>lo+buf&&c<hi-buf){
        int nl=(int)(std::lower_bound(s.begin(),s.end(),c)-s.begin());
        if(nl>=nm&&n-nl>=nm) cuts.push_back(c);
      }
    }
  }
  std::sort(cuts.begin(),cuts.end());
  cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
  return cuts;
}

#include "ppt_soft_cuts_cache.h"

static bool ppst_can_split(const PPSTNode&nd,const arma::mat&pts,int Dmax,
    int nmin,int mode,int ncand,PPSTCutsCache*cuts_cache=nullptr){
  if(nd.axis>=0||nd.depth>=Dmax) return false;
  if((int)nd.idx.n_elem<2*std::max(1,nmin)) return false;
  if(cuts_cache) return cuts_cache->any(nd);
  for(arma::uword j=0;j<nd.box.n_rows;j++)
    if(!ppst_axis_cuts(pts,nd.idx,nd.box,j,mode,ncand,nmin).empty())
      return true;
  return false;
}

static double ppst_stop_factor(const PPSTNode&nd,const arma::mat&pts,int Dmax,
    int nmin,int mode,int ncand,double alpha,double eta,
    PPSTCutsCache*cuts_cache=nullptr){
  if(cuts_cache){
    const bool supported=nd.depth<Dmax &&
      (int)nd.idx.n_elem>=2*std::max(1,nmin) && cuts_cache->any(nd);
    return supported ? 1.0-ppst_rho(nd.depth,alpha,eta) : 1.0;
  }
  PPSTNode tmp=nd; tmp.axis=-1;
  return ppst_can_split(tmp,pts,Dmax,nmin,mode,ncand)
    ? 1.0-ppst_rho(nd.depth,alpha,eta) : 1.0;
}

static void ppst_leaves(const PPSTree&T,std::vector<int>&out){
  out.clear(); out.reserve(T.size());
  for(const auto&kv:T) if(kv.second.axis<0) out.push_back(kv.first);
}
static void ppst_growable(const PPSTree&T,const arma::mat&pts,int Dmax,
    int nmin,int mode,int ncand,std::vector<int>&out,
    PPSTCutsCache*cuts_cache=nullptr){
  out.clear();
  for(const auto&kv:T)
    if(ppst_can_split(kv.second,pts,Dmax,nmin,mode,ncand,cuts_cache))
      out.push_back(kv.first);
}
static void ppst_prunable(const PPSTree&T,std::vector<int>&out){
  out.clear();
  for(const auto&kv:T) if(kv.second.axis>=0){
    int id=kv.first,L=2*id,R=2*id+1;
    auto il=T.find(L),ir=T.find(R);
    if(il!=T.end()&&ir!=T.end()&&il->second.axis<0&&ir->second.axis<0)
      out.push_back(id);
  }
}

static void ppst_split_indices(const PPSTNode&nd,const arma::mat&pts,int axis,
    double cut,arma::uvec&li,arma::uvec&ri){
  arma::vec col=pts.col(axis);
  std::vector<arma::uword>L,R;
  for(arma::uword k=0;k<nd.idx.n_elem;k++){
    arma::uword i=nd.idx[k];
    if(col[i]<cut) L.push_back(i); else R.push_back(i);
  }
  li.set_size(L.size()); ri.set_size(R.size());
  for(size_t k=0;k<L.size();k++) li[k]=L[k];
  for(size_t k=0;k<R.size();k++) ri[k]=R[k];
}

static PPSTNode ppst_child(const PPSTNode&parent,const arma::uvec&idx,int axis,
    double cut,int side){
  PPSTNode ch;
  ch.box=parent.box;
  if(side<0) ch.box(axis,1)=cut; else ch.box(axis,0)=cut;
  ch.idx=idx; ch.path=parent.path;
  PPSTGate g={axis,cut,parent.box(axis,1)-parent.box(axis,0),side};
  ch.path.push_back(g);
  ch.cut=NA_REAL; ch.axis=-1; ch.depth=parent.depth+1; ch.m=0;
  return ch;
}

// ---- recursive logistic gate and analytic exposure -------------------------
static double ppst_log_phi_logistic(const PPSTNode&nd,const arma::rowvec&x,
    const arma::mat&region,const arma::vec&gate,bool node_gate=false){
  double ans=0.0;
  for(const PPSTGate&g:nd.path){
    double width=node_gate ? g.parent_width : region(g.axis,1)-region(g.axis,0);
    double z=gate[g.axis]*(x[g.axis]-g.cut)/width;
    ans+=g.side<0 ? pst_logistic_log_right(-z)
                  : pst_logistic_log_right(z);
  }
  return ans;
}

static double ppst_exposure_logistic(const PPSTNode&nd,
    const arma::mat&region,const arma::vec&gate){
  double H=1.0;
  for(arma::uword j=0;j<region.n_rows;j++){
    std::vector<double>cuts;
    std::vector<int>sides;
    for(const PPSTGate&g:nd.path) if(g.axis==(int)j){
      cuts.push_back(g.cut);
      sides.push_back(g.side);
    }
    H*=pst_logistic_path_axis_integral(
      cuts,sides,region(j,0),region(j,1),gate[j]
    );
  }
  return std::max(H,1e-300);
}

// A node-relative gate has the same dimension-specific parameter but a
// different physical slope at every split. Integrate over the full root
// domain: logistic leaf bases extend outside their hard-routing boxes.
static double ppst_exposure_logistic_node(const PPSTNode&nd,
    const arma::mat&region,const arma::vec&gate){
  double H=1.0;
  for(arma::uword j=0;j<region.n_rows;j++){
    std::vector<double> cuts,widths;
    std::vector<int> sides;
    for(const PPSTGate&g:nd.path) if(g.axis==(int)j){
      cuts.push_back(g.cut);
      sides.push_back(g.side);
      widths.push_back(g.parent_width);
    }
    H*=pst_logistic_node_path_axis_integral(
      cuts,sides,widths,region(j,0),region(j,1),gate[j]
    );
  }
  return std::max(H,1e-300);
}

// ---- opt-in compact-cubic gate and exact exposure --------------------------
static inline double ppst_compact_gate_value(
    const PPSTGate&g,double x,double gate){
  double h=g.parent_width/gate;
  double t=(x-(g.cut-h))/(2.0*h);
  double r=t<=0.0?0.0:(t>=1.0?1.0:t*t*(3.0-2.0*t));
  return g.side<0?1.0-r:r;
}
static double ppst_log_phi_compact(const PPSTNode&nd,const arma::rowvec&x,
    const arma::vec&gate){
  double ans=0.0;
  for(const PPSTGate&g:nd.path){
    double q=ppst_compact_gate_value(g,x[g.axis],gate[g.axis]);
    if(q<=0.0) return -std::numeric_limits<double>::infinity();
    ans+=std::log(q);
  }
  return ans;
}
static std::vector<double> ppst_poly_multiply(const std::vector<double>&a,
    const std::vector<double>&b){
  std::vector<double> out(a.size()+b.size()-1,0.0);
  for(size_t i=0;i<a.size();i++)
    for(size_t j=0;j<b.size();j++) out[i+j]+=a[i]*b[j];
  return out;
}
static double ppst_axis_integral(const PPSTNode&nd,int axis,double dom_lo,
    double dom_hi,double gate){
  std::vector<const PPSTGate*> gates;
  std::vector<double> breaks={dom_lo,dom_hi};
  for(const PPSTGate&g:nd.path) if(g.axis==axis){
    gates.push_back(&g); double h=g.parent_width/gate;
    if(g.cut-h>dom_lo&&g.cut-h<dom_hi) breaks.push_back(g.cut-h);
    if(g.cut+h>dom_lo&&g.cut+h<dom_hi) breaks.push_back(g.cut+h);
  }
  if(gates.empty()) return dom_hi-dom_lo;
  std::sort(breaks.begin(),breaks.end());
  std::vector<double> uniq;
  for(double z:breaks)
    if(uniq.empty()||std::abs(z-uniq.back())>1e-13*(1.0+std::abs(z)))
      uniq.push_back(z);
  double total=0.0;
  for(size_t s=0;s+1<uniq.size();s++){
    double lo=uniq[s],hi=uniq[s+1],dx=hi-lo;
    if(dx<=0.0) continue;
    std::vector<double> poly(1,1.0); bool zero=false;
    for(const PPSTGate*gp:gates){
      const PPSTGate&g=*gp; double h=g.parent_width/gate;
      double tl=g.cut-h,tu=g.cut+h;
      if(hi<=tl){ if(g.side>0){zero=true;break;} else continue; }
      if(lo>=tu){ if(g.side<0){zero=true;break;} else continue; }
      double t0=(lo-tl)/(2.0*h),q=dx/(2.0*h);
      double t02=t0*t0,t03=t02*t0,q2=q*q,q3=q2*q;
      std::vector<double> right(4);
      right[0]=3.0*t02-2.0*t03;
      right[1]=6.0*t0*q-6.0*t02*q;
      right[2]=3.0*q2-6.0*t0*q2;
      right[3]=-2.0*q3;
      std::vector<double> f=right;
      if(g.side<0){
        f[0]=1.0-f[0];
        for(size_t k=1;k<f.size();k++) f[k]=-f[k];
      }
      poly=ppst_poly_multiply(poly,f);
    }
    if(zero) continue;
    double val=0.0;
    for(size_t k=0;k<poly.size();k++) val+=poly[k]/(double)(k+1);
    total+=dx*val;
  }
  if(total<0.0&&total>-1e-12*(dom_hi-dom_lo)) total=0.0;
  return std::max(total,1e-300);
}
static double ppst_exposure_compact(const PPSTNode&nd,const arma::mat&region,
    const arma::vec&gate){
  double H=1.0;
  for(arma::uword j=0;j<region.n_rows;j++)
    H*=ppst_axis_integral(nd,j,region(j,0),region(j,1),gate[j]);
  return std::max(H,1e-300);
}

static double ppst_log_phi(const PPSTNode&nd,const arma::rowvec&x,
    const arma::mat&region,const arma::vec&gate,int gate_family){
  return gate_family==1
    ? ppst_log_phi_compact(nd,x,gate)
    : ppst_log_phi_logistic(nd,x,region,gate,gate_family==2);
}

static inline double ppst_phi(const PPSTNode&nd,const arma::rowvec&x,
    const arma::mat&region,const arma::vec&gate,int gate_family){
  double z=ppst_log_phi(nd,x,region,gate,gate_family);
  return z < -745.0 ? 0.0 : std::exp(z);
}

static double ppst_exposure(const PPSTNode&nd,const arma::mat&region,
    const arma::vec&gate,int gate_family){
  if(gate_family==2) return ppst_exposure_logistic_node(nd,region,gate);
  return gate_family==0
    ? ppst_exposure_logistic(nd,region,gate)
    : ppst_exposure_compact(nd,region,gate);
}

#include "ppt_soft_geometry_cache.h"

// [[Rcpp::export]]
List ppstree_geometry(IntegerVector axis,NumericVector cut,
    NumericVector parent_width,IntegerVector side,arma::mat points,
    arma::mat region,arma::vec gate,int gate_scale=0){
  if(gate_scale!=0&&gate_scale!=1) stop("gate_scale must be 0 (root) or 1 (node)");
  const int family=gate_scale==1?2:0;
  int K=axis.size();
  gate=ppst_expand_positive(gate,region.n_rows,"gate");
  if(cut.size()!=K||parent_width.size()!=K||side.size()!=K)
    stop("path vectors must have equal lengths");
  PPSTNode nd; nd.box=region;
  for(int k=0;k<K;k++){
    if(axis[k]<1||axis[k]>(int)region.n_rows||parent_width[k]<=0.0||
       (side[k]!=-1&&side[k]!=1)) stop("invalid soft path");
    PPSTGate g={(int)axis[k]-1,cut[k],parent_width[k],side[k]};
    nd.path.push_back(g);
    if(side[k]<0) nd.box(axis[k]-1,1)=cut[k];
    else nd.box(axis[k]-1,0)=cut[k];
  }
  arma::vec phi(points.n_rows),lp(points.n_rows);
  for(arma::uword i=0;i<points.n_rows;i++){
    lp[i]=ppst_log_phi(nd,points.row(i),region,gate,family);
    phi[i]=lp[i] < -745.0 ? 0.0 : std::exp(lp[i]);
  }
  return List::create(_["phi"]=phi,_["log_phi"]=lp,
                      _["H"]=ppst_exposure(nd,region,gate,family));
}

// ---- leaf labels and gate update -------------------------------------------
static void ppst_label_sweep(PPSTree&T,std::vector<int>&labels,
    const arma::mat&pts,const arma::mat&region,double a,double b,
    const arma::vec&gate,int gate_family,PPSTGeometryCache*cache=nullptr){
  std::vector<int> leaves; ppst_leaves(T,leaves);
  std::vector<double> H(leaves.size()),lw(leaves.size());
  std::vector<const arma::vec*> basis(leaves.size(),nullptr);
  for(size_t k=0;k<leaves.size();k++){
    H[k]=ppst_cached_exposure(T.at(leaves[k]),region,gate,gate_family,cache);
    if(cache) basis[k]=&cache->training(T.at(leaves[k]));
  }
  for(arma::uword i=0;i<pts.n_rows;i++){
    T.at(labels[i]).m--;
    for(size_t k=0;k<leaves.size();k++){
      const PPSTNode&nd=T.at(leaves[k]);
      lw[k]=(cache?(*basis[k])[i]:ppst_log_phi(nd,pts.row(i),region,gate,gate_family))
           +std::log(a+nd.m)-std::log(b+H[k]);
    }
    double den=ppst_logsumexp(lw),u=R::unif_rand(),cum=0.0;
    int pick=leaves.back();
    for(size_t k=0;k<leaves.size();k++){
      cum+=std::exp(lw[k]-den);
      if(u<=cum){pick=leaves[k];break;}
    }
    labels[i]=pick; T.at(pick).m++;
  }
}

static double ppst_gate_target(const PPSTree&T,const std::vector<int>&labels,
    const arma::mat&pts,const arma::mat&region,double a,double b,
    const arma::vec&gate,const arma::vec&a_gate,const arma::vec&b_gate,
    const arma::vec&gate_min,bool gate_shared,int gate_family,
    PPSTGeometryCache*cache=nullptr){
  if(!gate.is_finite()||arma::any(gate<=gate_min)||arma::any(gate<=0.0))
    return -std::numeric_limits<double>::infinity();
  double out=0.0;
  arma::uword prior_n=gate_shared?1:gate.n_elem;
  for(arma::uword j=0;j<prior_n;j++)
    out+=(a_gate[j]-1.0)*std::log(gate[j])-b_gate[j]*gate[j];
  std::unordered_map<int,const arma::vec*> log_basis;
  for(const auto&kv:T) if(kv.second.axis<0){
    out+=ppst_log_g(
      kv.second.m,ppst_cached_exposure(kv.second,region,gate,gate_family,cache),a,b
    );
    if(cache) log_basis[kv.first]=&cache->training(kv.second);
  }
  for(arma::uword i=0;i<pts.n_rows;i++)
    out+=cache?(*log_basis.at(labels[i]))[i]:ppst_log_phi(
      T.at(labels[i]),pts.row(i),region,gate,gate_family
    );
  return out;
}
static int ppst_gate_update(const PPSTree&T,const std::vector<int>&labels,
    const arma::mat&pts,const arma::mat&region,double a,double b,
    arma::vec&gate,const arma::vec&a_gate,const arma::vec&b_gate,
    const arma::vec&sd_gate,const arma::vec&gate_min,bool gate_shared,
    int gate_family,int&which,PPSTGeometryCache*cache=nullptr){
  // `which` is the coordinate updated (ignored for a shared gate)
  if(gate_shared) which=0;
  arma::vec prop=gate;
  double cur=gate[which];
  double proposed=std::exp(std::log(cur)+R::rnorm(0.0,sd_gate[which]));
  if(gate_shared) prop.fill(proposed); else prop[which]=proposed;
  std::unique_ptr<PPSTGeometryCache> proposed_cache;
  if(cache) proposed_cache.reset(new PPSTGeometryCache(*cache,prop));
  double la=ppst_gate_target(T,labels,pts,region,a,b,prop,a_gate,b_gate,
                            gate_min,gate_shared,gate_family,proposed_cache.get())
           -ppst_gate_target(T,labels,pts,region,a,b,gate,a_gate,b_gate,
                             gate_min,gate_shared,gate_family,cache)
           +std::log(prop[which])-std::log(cur);
  if(std::log(R::unif_rand())<la){
    gate=prop;if(cache) *cache=std::move(*proposed_cache);return 1;
  }
  return 0;
}

#include "ppt_soft_pcg.h"

// ---- reversible tree moves -------------------------------------------------
static void ppst_child_log_memberships(const PPSTGate&split,double x,
    double parent_log_phi,const arma::mat&region,const arma::vec&gate,
    int gate_family,double&ll,double&lr){
  if(gate_family==1){
    PPSTGate right=split;right.side=1;
    const double q=ppst_compact_gate_value(right,x,gate[split.axis]);
    const double l=1.0-q;
    ll=l<=0.0?-std::numeric_limits<double>::infinity():parent_log_phi+std::log(l);
    lr=q<=0.0?-std::numeric_limits<double>::infinity():parent_log_phi+std::log(q);
  }else{
    const double width=gate_family==2?split.parent_width:
      region(split.axis,1)-region(split.axis,0);
    const double z=gate[split.axis]*(x-split.cut)/width;
    const double common=std::log1p(std::exp(-std::abs(z)));
    ll=parent_log_phi+(-std::max(z,0.0)-common);
    lr=parent_log_phi+(-std::max(-z,0.0)-common);
  }
}

static double ppst_side_probability(double ll,double lr){
  double den=std::max(ll,lr);
  return std::exp(ll-den)/
    (std::exp(ll-den)+std::exp(lr-den));
}

static int ppst_draw_side(const PPSTNode&left,const PPSTNode&right,
    const arma::mat&pts,arma::uword i,const arma::mat&region,
    const arma::vec&gate,int gate_family,const arma::vec*parent_log_phi){
  double ll,lr;
  if(parent_log_phi){
    const PPSTGate&split=right.path.back();
    ppst_child_log_memberships(split,pts(i,split.axis),(*parent_log_phi)[i],
      region,gate,gate_family,ll,lr);
  }else{
    ll=ppst_log_phi(left,pts.row(i),region,gate,gate_family);
    lr=ppst_log_phi(right,pts.row(i),region,gate,gate_family);
  }
  // Keep the parent terms and the original probability arithmetic: cancelling
  // their common factor changes floating-point rounding for extreme gates.
  const double p_left=ppst_side_probability(ll,lr);
  return R::unif_rand()<p_left ? -1 : 1;
}

static int ppst_grow_prune(PPSTree&T,std::vector<int>&labels,
    const arma::mat&pts,const arma::mat&region,double a,double b,
    const arma::vec&gate,int gate_family,
    double alpha,double eta,int Dmax,int nmin,int mode,int ncand,int&which_move,
    PPSTGeometryCache*cache=nullptr,PPSTCutsCache*cuts_cache=nullptr){
  std::vector<int>G,P;
  ppst_growable(T,pts,Dmax,nmin,mode,ncand,G,cuts_cache);
  ppst_prunable(T,P);
  if(R::unif_rand()<0.5){
    which_move=0;
    if(G.empty()) return 0;
    int v=G[ppst_runif_int(G.size())],L=2*v,R=2*v+1;
    PPSTNode old=T.at(v); int axis=ppst_runif_int(pts.n_cols);
    std::vector<double> cuts=ppst_cached_axis_cuts(old,pts,axis,mode,ncand,nmin,cuts_cache);
    if(cuts.empty()) return 0;
    double cut=cuts[ppst_runif_int(cuts.size())];
    arma::uvec li,ri; ppst_split_indices(old,pts,axis,cut,li,ri);
    PPSTNode nl=ppst_child(old,li,axis,cut,-1);
    PPSTNode nr=ppst_child(old,ri,axis,cut,1);
    std::vector<std::pair<int,int> > reass;
    int mL=0,mR=0;
    const arma::vec*parent_log_phi=nullptr;
    for(size_t i=0;i<labels.size();i++) if(labels[i]==v){
      // Do not populate a training basis for a proposal with no affected labels.
      if(cache&&!parent_log_phi) parent_log_phi=&cache->training(old);
      int dest=ppst_draw_side(
        nl,nr,pts,i,region,gate,gate_family,parent_log_phi
      )<0?L:R;
      reass.push_back(std::make_pair((int)i,dest));
      if(dest==L) mL++; else mR++;
    }
    nl.m=mL; nr.m=mR;
    T[v].axis=axis; T[v].cut=cut; T[v].m=0; T[L]=nl; T[R]=nr;
    std::vector<int>Pnew; ppst_prunable(T,Pnew);
    double Hp=ppst_cached_exposure(old,region,gate,gate_family,cache);
    double HL=ppst_cached_exposure(nl,region,gate,gate_family,cache);
    double HR=ppst_cached_exposure(nr,region,gate,gate_family,cache);
    double dp=std::log(ppst_rho(old.depth,alpha,eta))
      -std::log(1.0-ppst_rho(old.depth,alpha,eta))
      +std::log(ppst_stop_factor(nl,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache))
      +std::log(ppst_stop_factor(nr,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache));
    double logA=ppst_log_g(mL,HL,a,b)+ppst_log_g(mR,HR,a,b)
      -ppst_log_g(old.m,Hp,a,b)+dp
      +std::log((double)G.size())-std::log((double)Pnew.size());
    if(std::log(R::unif_rand())<logA){
      for(const auto&z:reass) labels[z.first]=z.second;
      return 1;
    }
    T.erase(L); T.erase(R); T[v]=old;
    return 0;
  }

  which_move=1;
  if(P.empty()) return 0;
  int v=P[ppst_runif_int(P.size())],L=2*v,R=2*v+1;
  PPSTNode oldp=T.at(v),oldL=T.at(L),oldR=T.at(R);
  int M=oldL.m+oldR.m;
  double Hp=ppst_cached_exposure(oldp,region,gate,gate_family,cache);
  double HL=ppst_cached_exposure(oldL,region,gate,gate_family,cache);
  double HR=ppst_cached_exposure(oldR,region,gate,gate_family,cache);
  double dp=std::log(ppst_rho(oldp.depth,alpha,eta))
    -std::log(1.0-ppst_rho(oldp.depth,alpha,eta))
    +std::log(ppst_stop_factor(oldL,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache))
    +std::log(ppst_stop_factor(oldR,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache));
  T.erase(L); T.erase(R); T[v].axis=-1; T[v].cut=NA_REAL; T[v].m=M;
  std::vector<int>Gnew; ppst_growable(T,pts,Dmax,nmin,mode,ncand,Gnew,cuts_cache);
  if(Gnew.empty()){T[v]=oldp;T[L]=oldL;T[R]=oldR;return 0;}
  double logA=ppst_log_g(M,Hp,a,b)-ppst_log_g(oldL.m,HL,a,b)
    -ppst_log_g(oldR.m,HR,a,b)-dp
    +std::log((double)P.size())-std::log((double)Gnew.size());
  if(std::log(R::unif_rand())<logA){
    for(size_t i=0;i<labels.size();i++)
      if(labels[i]==L||labels[i]==R) labels[i]=v;
    return 1;
  }
  T[v]=oldp; T[L]=oldL; T[R]=oldR;
  return 0;
}

static int ppst_change_cut(PPSTree&T,std::vector<int>&labels,
    const arma::mat&pts,const arma::mat&region,double a,double b,
    const arma::vec&gate,int gate_family,
    double alpha,double eta,int Dmax,int nmin,int mode,int ncand,
    PPSTGeometryCache*cache=nullptr,PPSTCutsCache*cuts_cache=nullptr){
  std::vector<int>P; ppst_prunable(T,P); if(P.empty()) return 0;
  int v=P[ppst_runif_int(P.size())],L=2*v,R=2*v+1;
  PPSTNode parent=T.at(v),oldL=T.at(L),oldR=T.at(R);
  int axis=ppst_runif_int(pts.n_cols);
  std::vector<double>cuts=ppst_cached_axis_cuts(parent,pts,axis,mode,ncand,nmin,cuts_cache);
  if(cuts.empty()) return 0;
  double cut=cuts[ppst_runif_int(cuts.size())];
  arma::uvec li,ri; ppst_split_indices(parent,pts,axis,cut,li,ri);
  PPSTNode nl=ppst_child(parent,li,axis,cut,-1);
  PPSTNode nr=ppst_child(parent,ri,axis,cut,1);
  std::vector<std::pair<int,int> >reass;
  int mL=0,mR=0;
  const arma::vec*parent_log_phi=nullptr;
  for(size_t i=0;i<labels.size();i++) if(labels[i]==L||labels[i]==R){
    if(cache&&!parent_log_phi) parent_log_phi=&cache->training(parent);
    int dest=ppst_draw_side(
      nl,nr,pts,i,region,gate,gate_family,parent_log_phi
    )<0?L:R;
    reass.push_back(std::make_pair((int)i,dest));
    if(dest==L)mL++;else mR++;
  }
  nl.m=mL; nr.m=mR;
  double HoL=ppst_cached_exposure(oldL,region,gate,gate_family,cache);
  double HoR=ppst_cached_exposure(oldR,region,gate,gate_family,cache);
  double HnL=ppst_cached_exposure(nl,region,gate,gate_family,cache);
  double HnR=ppst_cached_exposure(nr,region,gate,gate_family,cache);
  double oldg=ppst_log_g(oldL.m,HoL,a,b)+ppst_log_g(oldR.m,HoR,a,b);
  double newg=ppst_log_g(mL,HnL,a,b)+ppst_log_g(mR,HnR,a,b);
  double oldS=std::log(ppst_stop_factor(oldL,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache))
             +std::log(ppst_stop_factor(oldR,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache));
  double newS=std::log(ppst_stop_factor(nl,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache))
             +std::log(ppst_stop_factor(nr,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache));
  if(std::log(R::unif_rand())<(newg-oldg)+(newS-oldS)){
    T[v].axis=axis;T[v].cut=cut;T[L]=nl;T[R]=nr;
    for(const auto&z:reass) labels[z.first]=z.second;
    return 1;
  }
  return 0;
}

#include "ppt_soft_informed.h"
#include "ppt_soft_surrogate.h"

// ---- posterior draws and chain ---------------------------------------------
static double ppst_eval_intensity(const PPSTree&T,
    const std::unordered_map<int,double>&lam,const arma::rowvec&x,
    const arma::mat&region,const arma::vec&gate,int gate_family){
  double out=0.0;
  for(const auto&kv:lam)
    out+=kv.second*ppst_phi(
      T.at(kv.first),x,region,gate,gate_family
    );
  return std::max(out,1e-300);
}

static int ppst_run_chain(const arma::mat&pts,const arma::mat&grid,
    const arma::mat&xt,const arma::mat&region,double a,double b,
    const arma::vec&gate0,const arma::vec&a_gate,const arma::vec&b_gate,
    const arma::vec&sd_gate,const arma::vec&gate_min,bool gate_shared,
    int gate_family,double alpha,double eta,int Dmax,int nmin,int iters,
    int burn,int thin,
    int nmove,int ncc,int mode,int ncand,int update_gate,arma::mat&draws,
    arma::vec&loglik,arma::vec&loglik_test,arma::vec&integrated_intensity,
    int row0,double&mean_leaves,
    double&mean_max_depth,arma::vec&mean_gate,arma::vec&accept,
    arma::vec&gate_accept,
    std::vector<arma::mat>&state_nodes,arma::mat&state_gate,bool informed,bool verbose,
    bool pcg,double ram_target,double ram_decay,int ram_adapt,PPSTPCGStats&pcg_stats,
    bool cache_geometry,bool cache_cuts,double tau,double eps){
  int n=pts.n_rows,si=0; double nls=0.0,mds=0.0;
  // informed=TRUE means the exact neighborhood sampler for RJ-MCMC and the
  // hard-surrogate proposals for PCG.
  const bool exact_informed=informed&&!pcg,surrogate_informed=informed&&pcg;
  arma::vec gate=gate0,gs(gate0.n_elem,arma::fill::zeros);
  if(pcg) pcg_stats.initialize(sd_gate,gate_shared);
  PPSTGeometryCache geometry(pts,region,gate,gate_family,&grid,cache_geometry);
  PPSTGeometryCache*cache=cache_geometry?&geometry:nullptr;
  PPSTCutsCache cuts(pts,nmin,mode,ncand);
  PPSTCutsCache*cuts_cache=cache_cuts&&!exact_informed?&cuts:nullptr;
  if(surrogate_informed&&!cuts_cache) stop("informed pcg requires cache_cuts");
  PPSTSurrogate surrogate(pts,region,cuts,a,b,alpha,eta,tau,eps);
  PPSTree T; PPSTNode root;
  root.box=region; root.idx=arma::regspace<arma::uvec>(0,n-1);
  root.cut=NA_REAL;root.axis=-1;root.depth=0;root.m=n;T[1]=root;
  std::vector<int>labels(n,1);
  double ag=0,tg=0,ap=0,tp=0,ac=0,tc=0;
  arma::vec aga(gate.n_elem,arma::fill::zeros);
  arma::vec tga(gate.n_elem,arma::fill::zeros);
  PPSTIContext informed_context(pts,region,gate,a,b,alpha,eta,
                                Dmax,nmin,mode,ncand,gate_family);
  PPTMCMCProgress progress(iters, verbose);
  for(int it=0;it<iters;it++){
    if(pcg){
      ppst_pcg_block(T,labels,pts,region,a,b,gate,a_gate,b_gate,gate_min,
        gate_shared,gate_family,update_gate!=0,pcg_stats,cache);
    }else{
    ppst_label_sweep(T,labels,pts,region,a,b,gate,gate_family,cache);
    if(update_gate){
      // one Metropolis proposal per coordinate (systematic scan), or one
      // proposal for a shared gate
      int nup=gate_shared?1:(int)gate.n_elem;
      for(int which=0;which<nup;which++){
        int accepted=ppst_gate_update(T,labels,pts,region,a,b,gate,a_gate,b_gate,
                                      sd_gate,gate_min,gate_shared,gate_family,
                                      which,cache);
        aga[which]+=accepted; tga[which]++;
      }
    }
    }
    if(exact_informed && update_gate) informed_context.clear_gate_cache();
    PPSTINeighborhood informed_gp,informed_change;
    bool valid_gp=false,valid_change=false;
    for(int r=0;r<nmove;r++){
      int move=-1,ok;
      if(exact_informed) ok=ppsti_step(T,labels,informed_context,0,
                                informed_gp,valid_gp,move);
      else if(surrogate_informed) ok=ppsts_grow_prune(
        T,labels,pts,region,a,b,gate,gate_family,alpha,eta,
        Dmax,nmin,mode,ncand,move,cache,cuts_cache,surrogate
      );
      else ok=ppst_grow_prune(
        T,labels,pts,region,a,b,gate,gate_family,alpha,eta,
        Dmax,nmin,mode,ncand,move,cache,cuts_cache
      );
      if(cache) cache->trim(T);
      if(cuts_cache) cuts_cache->trim(T);
      if(move==0){ag+=ok;tg++;}else if(move==1){ap+=ok;tp++;}
    }
    for(int r=0;r<ncc;r++){
      if(exact_informed){
        int move=-1;
        ac+=ppsti_step(T,labels,informed_context,1,
                      informed_change,valid_change,move);
      }else if(surrogate_informed) ac+=ppsts_change_cut(
        T,labels,pts,region,a,b,gate,gate_family,alpha,eta,Dmax,nmin,
        mode,ncand,cache,cuts_cache,surrogate
      );
      else ac+=ppst_change_cut(
        T,labels,pts,region,a,b,gate,gate_family,alpha,eta,Dmax,nmin,
        mode,ncand,cache,cuts_cache
      );
      if(cache) cache->trim(T);
      if(cuts_cache) cuts_cache->trim(T);
      tc++;
    }
    if(surrogate_informed) surrogate.trim(T);
    if(pcg) ppst_pcg_adapt(pcg_stats,it,ram_target,ram_decay,ram_adapt);
    if(it>=burn&&(it-burn)%thin==0){
      int row=row0+si; std::vector<int>leaves; ppst_leaves(T,leaves);
      std::unordered_map<int,double>lam; double comp=0.0;
      for(int id:leaves){
        double H=ppst_cached_exposure(T.at(id),region,gate,gate_family,cache);
        double lv=R::rgamma(a+T.at(id).m,1.0/(b+H));
        lam[id]=lv;comp+=lv*H;
      }
      integrated_intensity[row]=comp;
      if(cache){
        draws.row(row)=cache->intensity(T,lam,grid).t();
        const arma::vec fitted=cache->intensity(T,lam,pts);
        double ll=-comp;
        for(int i=0;i<n;i++) ll+=std::log(fitted[i]);
        loglik[row]=ll;
        if(xt.n_rows>0){
          const arma::vec predicted=cache->intensity(T,lam,xt);
          double lt=-comp;
          for(arma::uword i=0;i<xt.n_rows;i++) lt+=std::log(predicted[i]);
          loglik_test[row]=lt;
        }
      }else{
        for(arma::uword j=0;j<grid.n_rows;j++)
          draws(row,j)=ppst_eval_intensity(T,lam,grid.row(j),region,gate,gate_family);
        double ll=-comp;
        for(int i=0;i<n;i++)
          ll+=std::log(ppst_eval_intensity(T,lam,pts.row(i),region,gate,gate_family));
        loglik[row]=ll;
        if(xt.n_rows>0){
          double lt=-comp;
          for(arma::uword i=0;i<xt.n_rows;i++)
            lt+=std::log(ppst_eval_intensity(T,lam,xt.row(i),region,gate,gate_family));
          loglik_test[row]=lt;
        }
      }
      int md=0;
      for(const auto&kv:T) if(kv.second.depth>md) md=kv.second.depth;
      // serialize the generative state of this retained draw: one row per
      // node, columns (heap id, axis, cut, lambda, xi, m).  Internal nodes
      // carry no rate in the terminal-leaf model (lambda = 0).
      {
        arma::mat st(T.size(),6);
        int rr=0;
        for(const auto&kv:T){
          const PPSTNode&nd=kv.second;
          st(rr,0)=(double)kv.first;
          st(rr,1)=(double)nd.axis;
          st(rr,2)=nd.axis>=0?nd.cut:NA_REAL;
          st(rr,3)=nd.axis<0?lam[kv.first]:0.0;
          st(rr,4)=NA_REAL;
          st(rr,5)=(double)nd.m;
          rr++;
        }
        state_nodes[row]=st;
        state_gate.row(row)=gate.t();
      }
      nls+=leaves.size();mds+=md;gs+=gate;si++;
    }
    progress.increment();
  }
  mean_leaves=nls/std::max(1,si);
  mean_max_depth=mds/std::max(1,si);
  mean_gate=gs/std::max(1,si);
  accept.set_size(3);
  accept[0]=ag/std::max(1.0,tg);accept[1]=ap/std::max(1.0,tp);
  accept[2]=ac/std::max(1.0,tc);
  gate_accept.set_size(gate.n_elem);
  if(pcg){
    gate_accept.fill(pcg_stats.acceptance());
  }else if(gate_shared){
    gate_accept.fill(aga[0]/std::max(1.0,tga[0]));
  }else{
    for(arma::uword j=0;j<gate.n_elem;j++)
      gate_accept[j]=aga[j]/std::max(1.0,tga[j]);
  }
  return si;
}

static void ppst_summarize(const arma::mat&D,arma::vec&mn,arma::vec&md,
    arma::vec&lo,arma::vec&hi){
  mn=arma::mean(D,0).t();md=arma::median(D,0).t();
  lo.set_size(D.n_cols);hi.set_size(D.n_cols);
  arma::vec p={0.025,0.975};
  for(arma::uword j=0;j<D.n_cols;j++){
    arma::vec q=arma::quantile(D.col(j),p);lo[j]=q[0];hi[j]=q[1];
  }
}
static double ppst_lse_vec(const arma::vec&x){
  double mx=x.max();return mx+std::log(arma::sum(arma::exp(x-mx)));
}

// [[Rcpp::export]]
List ppstree_multi(arma::mat X,arma::mat grid,arma::mat Xtest,arma::mat region,
    double a,double b,arma::vec gate,arma::vec a_gate,arma::vec b_gate,
    arma::vec sd_gate,arma::vec gate_min,int gate_shared,double alpha,
    double eta,double Dmax,int nmin,int iters,int burn,int thin,int nmove,int ncc,
    int cut_mode,int ncand,int update_gate,int gate_family,int chains,
    int verbose,bool informed=false,bool pcg=false,double ram_target=0.234,
    double ram_decay=0.7,int ram_adapt=0,bool cache_geometry=true,
    bool cache_cuts=true,double proposal_temperature=0.5,
    double proposal_defensive=0.1){
  const int depth = ppt_checked_depth(Dmax);
  ppst_pcg_controls(pcg,informed,ram_target,ram_decay,ram_adapt,burn,
                    proposal_temperature,proposal_defensive);
  if(X.n_rows==0||X.n_cols==0) stop("X must be a non-empty matrix");
  if(region.n_rows!=X.n_cols||region.n_cols!=2||
     grid.n_cols!=X.n_cols||Xtest.n_cols!=X.n_cols)
    stop("X, grid, Xtest, and region have incompatible dimensions");
  if(chains<1||iters<=burn||burn<0||thin<1||depth<0||nmin<1||
     nmove<0||ncc<0||ncand<2)
    stop("invalid MCMC or tree controls");
  if(gate_family<0||gate_family>2)
    stop("gate_family must be 0 (root logistic), 1 (compact), or 2 (node logistic)");
  if(a<=0.0||b<=0.0) stop("a and b must be positive for the soft PPT");
  gate=ppst_expand_positive(gate,X.n_cols,"gate");
  a_gate=ppst_expand_positive(a_gate,X.n_cols,"a_gate");
  b_gate=ppst_expand_positive(b_gate,X.n_cols,"b_gate");
  sd_gate=ppst_expand_positive(sd_gate,X.n_cols,"sd_gate");
  gate_min=ppst_expand_positive(gate_min,X.n_cols,"gate_min",true);
  if(arma::any(gate<=gate_min)) stop("every gate must exceed gate_min");
  if(pcg&&gate_shared&&arma::any(gate!=gate[0]))
    stop("shared pcg gates must have the same initial value");
  int ns=0;for(int it=burn;it<iters;it++)if((it-burn)%thin==0)ns++;
  int total=ns*chains,row=0;
  arma::mat D(total,grid.n_rows,arma::fill::zeros);
  arma::vec ll(total,arma::fill::zeros);
  arma::vec llt(Xtest.n_rows>0?total:0,arma::fill::zeros);
  arma::vec integrated_intensity(total,arma::fill::zeros);
  arma::vec leaves(chains),maxdepth(chains);arma::mat gates(chains,X.n_cols);
  arma::mat acc(chains,3),gacc(chains,X.n_cols);
  std::vector<arma::mat> state_nodes(total);
  arma::mat state_gate(total,X.n_cols,arma::fill::zeros);
  List ram_covariance(chains),ram_factor(chains);
  IntegerVector ram_updates(chains),ram_failures(chains);
  arma::vec joint_accept(chains,arma::fill::zeros),joint_alpha(chains,arma::fill::zeros);
  if(verbose) Rcpp::Rcerr<<(pcg?"S-PPT [PCG with RAM]: ":
    (informed?"S-PPT [informed RJ-MCMC]: ":"S-PPT [RJ-MCMC]: "))<<chains
                         <<" chains, "<<ns<<" draws/chain\n";
  for(int k=0;k<chains;k++){
    if(verbose) Rcpp::Rcerr<<"  chain "<<k+1<<"/"<<chains<<"\n";
    arma::vec ak,gak,gm;double nl,md;PPSTPCGStats pcg_stats;
    int got=ppst_run_chain(X,grid,Xtest,region,a,b,gate,a_gate,b_gate,sd_gate,
      gate_min,gate_shared,gate_family,alpha,eta,depth,nmin,iters,burn,thin,
      nmove,ncc,cut_mode,ncand,update_gate,D,ll,llt,integrated_intensity,
      row,nl,md,gm,ak,gak,state_nodes,state_gate,informed,verbose != 0,
      pcg,ram_target,ram_decay,ram_adapt,pcg_stats,cache_geometry,cache_cuts,
      proposal_temperature,proposal_defensive);
    row+=got;leaves[k]=nl;maxdepth[k]=md;gates.row(k)=gm.t();acc.row(k)=ak.t();
    gacc.row(k)=gak.t();
    if(pcg){
      ram_factor[k]=pcg_stats.factor;
      ram_covariance[k]=arma::mat(pcg_stats.factor*pcg_stats.factor.t());
      ram_updates[k]=pcg_stats.updates;ram_failures[k]=pcg_stats.failures;
      joint_accept[k]=pcg_stats.acceptance();joint_alpha[k]=pcg_stats.mean_alpha();
    }
    if(verbose) Rcpp::Rcerr<<"  chain "<<k+1<<"/"<<chains<<" done; leaves="
      <<nl<<", max_depth="<<md<<", gate="<<arma::mean(gm)
      <<", gate_accept="<<arma::mean(gak)<<"\n";
  }
  arma::vec mn,md,lo,hi;ppst_summarize(D,mn,md,lo,hi);
  double lp=(Xtest.n_rows>0&&row>0)
    ?ppst_lse_vec(llt)-std::log((double)row):NA_REAL;
  List sn(row);
  for(int s=0;s<row;s++) sn[s]=state_nodes[s];
  List out=List::create(
    _["mean"]=mn,_["median"]=md,_["lower95"]=lo,_["upper95"]=hi,
    _["draws"]=D.t(),
    _["loglik_mean"]=row>0?arma::mean(ll):NA_REAL,_["logpred"]=lp,
    _["loglik_draws"]=ll,
    _["integrated_intensity"]=integrated_intensity,
    _["chain_mean_leaves"]=leaves,
    _["chain_mean_max_depth"]=maxdepth,
    _["mean_leaves"]=arma::mean(leaves),
    _["mean_max_depth"]=arma::mean(maxdepth),
    _["kappa_mean"]=NA_REAL,_["tau_mean"]=NA_REAL,
    _["gate_mean"]=arma::mean(gates,0).t(),
    _["accept"]=arma::mean(acc,0).t(),
    _["gate_accept"]=arma::mean(gacc,0).t(),_["ndraws"]=row);
  out["state_nodes"]=sn;
  out["state_gate"]=state_gate.rows(0,std::max(0,row-1));
  if(pcg){
    out["ram_covariance"]=ram_covariance;out["ram_factor"]=ram_factor;
    out["ram_updates"]=ram_updates;out["ram_failures"]=ram_failures;
    out["gate_joint_accept"]=arma::mean(joint_accept);
    out["chain_gate_joint_accept"]=joint_accept;
    out["chain_gate_joint_accept_prob"]=joint_alpha;
  }
  return out;
}

// [[Rcpp::export]]
List ppstree_diag(arma::mat X,arma::mat mon,arma::mat region,
    double a,double b,arma::vec gate,arma::vec a_gate,arma::vec b_gate,
    arma::vec sd_gate,arma::vec gate_min,int gate_shared,double alpha,
    double eta,double Dmax,int nmin,int iters,int burn,int thin,int nmove,int ncc,
    int cut_mode,int ncand,int update_gate,int gate_family,bool informed=false,
    bool verbose=false,bool pcg=false,double ram_target=0.234,
    double ram_decay=0.7,int ram_adapt=0,bool cache_geometry=true,
    bool cache_cuts=true,double proposal_temperature=0.5,
    double proposal_defensive=0.1){
  const int depth = ppt_checked_depth(Dmax);
  ppst_pcg_controls(pcg,informed,ram_target,ram_decay,ram_adapt,burn,
                    proposal_temperature,proposal_defensive);
  if(X.n_rows==0||X.n_cols==0) stop("X must be a non-empty matrix");
  if(region.n_rows!=X.n_cols||region.n_cols!=2||mon.n_cols!=X.n_cols)
    stop("X, mon, and region have incompatible dimensions");
  if(iters<=burn||burn<0||thin<1||depth<0||nmin<1||
     nmove<0||ncc<0||ncand<2)
    stop("invalid MCMC or tree controls");
  if(gate_family<0||gate_family>2)
    stop("gate_family must be 0 (root logistic), 1 (compact), or 2 (node logistic)");
  if(a<=0.0||b<=0.0) stop("a and b must be positive for the soft PPT");
  gate=ppst_expand_positive(gate,X.n_cols,"gate");
  a_gate=ppst_expand_positive(a_gate,X.n_cols,"a_gate");
  b_gate=ppst_expand_positive(b_gate,X.n_cols,"b_gate");
  sd_gate=ppst_expand_positive(sd_gate,X.n_cols,"sd_gate");
  gate_min=ppst_expand_positive(gate_min,X.n_cols,"gate_min",true);
  if(arma::any(gate<=gate_min)) stop("every gate must exceed gate_min");
  if(pcg&&gate_shared&&arma::any(gate!=gate[0]))
    stop("shared pcg gates must have the same initial value");

  int n=X.n_rows,ns=0,si=0;
  for(int it=burn;it<iters;it++)
    if((it-burn)%thin==0) ns++;

  PPSTGeometryCache geometry(X,region,gate,gate_family,&mon,cache_geometry);
  PPSTGeometryCache*cache=cache_geometry?&geometry:nullptr;
  PPSTCutsCache cuts(X,nmin,cut_mode,ncand);
  const bool exact_informed=informed&&!pcg,surrogate_informed=informed&&pcg;
  PPSTCutsCache*cuts_cache=cache_cuts&&!exact_informed?&cuts:nullptr;
  if(surrogate_informed&&!cuts_cache) stop("informed pcg requires cache_cuts");
  PPSTSurrogate surrogate(X,region,cuts,a,b,alpha,eta,proposal_temperature,proposal_defensive);
  PPSTree T; PPSTNode root;
  root.box=region; root.idx=arma::regspace<arma::uvec>(0,n-1);
  root.cut=NA_REAL;root.axis=-1;root.depth=0;root.m=n;T[1]=root;
  std::vector<int>labels(n,1);

  PPSTPCGStats pcg_stats;
  if(pcg) pcg_stats.initialize(sd_gate,gate_shared);
  arma::vec tr_nleaf(ns,arma::fill::zeros);
  arma::mat tr_gate(ns,gate.n_elem,arma::fill::zeros);
  arma::vec tr_logdens(ns,arma::fill::zeros);
  arma::mat tr_mon(ns,mon.n_rows,arma::fill::zeros);
  double ag=0,tg=0,ap=0,tp=0,ac=0,tc=0;
  arma::vec aga(gate.n_elem,arma::fill::zeros);
  arma::vec tga(gate.n_elem,arma::fill::zeros);
  PPSTIContext informed_context(X,region,gate,a,b,alpha,eta,
                                depth,nmin,cut_mode,ncand,gate_family);

  PPTMCMCProgress progress(iters, verbose);
  for(int it=0;it<iters;it++){
    if(pcg){
      ppst_pcg_block(T,labels,X,region,a,b,gate,a_gate,b_gate,gate_min,
        gate_shared,gate_family,update_gate!=0,pcg_stats,cache);
    }else{
    ppst_label_sweep(T,labels,X,region,a,b,gate,gate_family,cache);
    if(update_gate){
      // Match the fitting backend: one proposal per dimension, or one
      // proposal for a shared gate.
      int nup=gate_shared?1:(int)gate.n_elem;
      for(int which=0;which<nup;which++){
        int accepted=ppst_gate_update(T,labels,X,region,a,b,gate,a_gate,b_gate,
                                      sd_gate,gate_min,gate_shared,gate_family,
                                      which,cache);
        aga[which]+=accepted; tga[which]++;
      }
    }
    }
    if(exact_informed && update_gate) informed_context.clear_gate_cache();
    PPSTINeighborhood informed_gp,informed_change;
    bool valid_gp=false,valid_change=false;
    for(int r=0;r<nmove;r++){
      int move=-1,ok;
      if(exact_informed) ok=ppsti_step(T,labels,informed_context,0,
                                informed_gp,valid_gp,move);
      else if(surrogate_informed) ok=ppsts_grow_prune(
        T,labels,X,region,a,b,gate,gate_family,alpha,eta,
        depth,nmin,cut_mode,ncand,move,cache,cuts_cache,surrogate
      );
      else ok=ppst_grow_prune(
        T,labels,X,region,a,b,gate,gate_family,alpha,eta,
        depth,nmin,cut_mode,ncand,move,cache,cuts_cache
      );
      if(cache) cache->trim(T);
      if(cuts_cache) cuts_cache->trim(T);
      if(move==0){ag+=ok;tg++;}
      else if(move==1){ap+=ok;tp++;}
    }
    for(int r=0;r<ncc;r++){
      if(exact_informed){
        int move=-1;
        ac+=ppsti_step(T,labels,informed_context,1,
                      informed_change,valid_change,move);
      }else if(surrogate_informed) ac+=ppsts_change_cut(
        T,labels,X,region,a,b,gate,gate_family,alpha,eta,depth,nmin,
        cut_mode,ncand,cache,cuts_cache,surrogate
      );
      else ac+=ppst_change_cut(
        T,labels,X,region,a,b,gate,gate_family,alpha,eta,depth,nmin,
        cut_mode,ncand,cache,cuts_cache
      );
      if(cache) cache->trim(T);
      if(cuts_cache) cuts_cache->trim(T);
      tc++;
    }

    if(surrogate_informed) surrogate.trim(T);
    if(pcg) ppst_pcg_adapt(pcg_stats,it,ram_target,ram_decay,ram_adapt);
    if(it>=burn&&(it-burn)%thin==0){
      std::vector<int>leaves;
      ppst_leaves(T,leaves);
      std::unordered_map<int,double>lam;
      double compensator=0.0;
      for(int id:leaves){
        double H=ppst_cached_exposure(T.at(id),region,gate,gate_family,cache);
        double lv=R::rgamma(a+T.at(id).m,1.0/(b+H));
        lam[id]=lv;
        compensator+=lv*H;
      }
      double ll=-compensator;
      if(cache){
        tr_mon.row(si)=cache->intensity(T,lam,mon).t();
        const arma::vec fitted=cache->intensity(T,lam,X);
        for(int i=0;i<n;i++) ll+=std::log(fitted[i]);
      }else{
        for(arma::uword j=0;j<mon.n_rows;j++)
          tr_mon(si,j)=ppst_eval_intensity(T,lam,mon.row(j),region,gate,gate_family);
        for(int i=0;i<n;i++)
          ll+=std::log(ppst_eval_intensity(T,lam,X.row(i),region,gate,gate_family));
      }
      tr_nleaf[si]=leaves.size();
      tr_gate.row(si)=gate.t();
      tr_logdens[si]=ll;
      si++;
    }
    progress.increment();
  }

  arma::vec accept(3);
  accept[0]=ag/std::max(1.0,tg);
  accept[1]=ap/std::max(1.0,tp);
  accept[2]=ac/std::max(1.0,tc);
  arma::vec gate_accept(gate.n_elem);
  if(pcg){
    gate_accept.fill(pcg_stats.acceptance());
  }else if(gate_shared){
    gate_accept.fill(aga[0]/std::max(1.0,tga[0]));
  }else{
    for(arma::uword j=0;j<gate.n_elem;j++)
      gate_accept[j]=aga[j]/std::max(1.0,tga[j]);
  }
  List out=List::create(
    _["nleaf"]=tr_nleaf,_["gate"]=tr_gate,_["logdens"]=tr_logdens,
    _["mon"]=tr_mon,_["accept"]=accept,_["gate_accept"]=gate_accept,_["ns"]=si);
  if(pcg){
    out["ram_covariance"]=arma::mat(pcg_stats.factor*pcg_stats.factor.t());
    out["ram_factor"]=pcg_stats.factor;
    out["ram_updates"]=pcg_stats.updates;out["ram_failures"]=pcg_stats.failures;
    out["gate_joint_accept"]=pcg_stats.acceptance();
    out["gate_joint_accept_prob"]=pcg_stats.mean_alpha();
  }
  return out;
}

// Internal exact-transition inspector for small labeled states. The fitting
// backend never enumerates labels; this diagnostic deliberately does so to
// check its Poisson-binomial action masses against independent finite targets.
// [[Rcpp::export]]
List ppstree_informed_transition(arma::mat X,arma::mat region,
    arma::mat splits,IntegerVector labels,arma::vec gate,double a,double b,
    double alpha,double eta,double Dmax,int nmin,int cut_mode,int ncand,
    int gate_family,int kind){
  const int depth = ppt_checked_depth(Dmax);
  if(X.n_rows==0||X.n_rows>8||X.n_cols==0)
    stop("informed transition inspection requires 1 to 8 observations");
  if(region.n_rows!=X.n_cols||region.n_cols!=2||splits.n_cols!=3||
     (arma::uword)labels.size()!=X.n_rows)
    stop("incompatible transition-inspection dimensions");
  if(a<=0||b<=0||depth<0||nmin<1||ncand<2||
     (gate_family<0||gate_family>2)||(kind!=0&&kind!=1))
    stop("invalid transition-inspection controls");
  gate=ppst_expand_positive(gate,X.n_cols,"gate");
  PPSTree tree; PPSTNode root;
  root.box=region;root.idx=arma::regspace<arma::uvec>(0,X.n_rows-1);
  root.axis=-1;root.cut=NA_REAL;root.depth=0;root.m=0;tree[1]=root;
  PPSTIContext context(X,region,gate,a,b,alpha,eta,depth,nmin,cut_mode,ncand,gate_family);
  arma::uvec order=arma::sort_index(splits.col(0));
  for(arma::uword k:order){
    if(!std::isfinite(splits(k,0))||!std::isfinite(splits(k,1))||
       splits(k,0)>std::numeric_limits<int>::max()||splits(k,0)<1||
       splits(k,1)<0||splits(k,1)>=(double)X.n_cols)
      stop("invalid split index in informed transition inspection");
    int id=(int)splits(k,0),axis=(int)splits(k,1); double cut=splits(k,2);
    if(splits(k,0)!=id||splits(k,1)!=axis||
       !std::isfinite(cut)||tree.find(id)==tree.end()||tree.at(id).axis>=0)
      stop("invalid split topology in informed transition inspection");
    PPSTNode parent=tree.at(id);
    auto support=context.cuts(parent);
    if(std::find(support->axis[axis].begin(),support->axis[axis].end(),cut)==
       support->axis[axis].end()) stop("split outside candidate support");
    arma::uvec li,ri;ppst_split_indices(parent,X,axis,cut,li,ri);
    tree[id].axis=axis;tree[id].cut=cut;
    tree[2*id]=ppst_child(parent,li,axis,cut,-1);
    tree[2*id+1]=ppst_child(parent,ri,axis,cut,1);
  }
  std::vector<int> current_labels(labels.begin(),labels.end());
  for(size_t i=0;i<current_labels.size();i++){
    auto it=tree.find(current_labels[i]);
    if(it==tree.end()||it->second.axis>=0||
       !std::isfinite(ppst_log_phi(it->second,X.row(i),region,gate,gate_family)))
      stop("labels must select terminal leaves with positive basis values");
    it->second.m++;
  }
  PPSTINeighborhood neighborhood=ppsti_neighborhood(tree,current_labels,context,kind);
  std::vector<List> neighbors;
  for(const PPSTIAction&act:neighborhood.actions){
    int m=(int)act.affected.size(),count=act.move==1?1:(1<<m);
    for(int mask=0;mask<count;mask++){
      PPSTree next=tree;std::vector<int> next_labels=current_labels;
      double logalloc=0.0;int nleft=0;
      if(act.move!=1) for(int k=0;k<m;k++){
        bool left=(mask&(1<<k))!=0;int i=act.affected[k];
        next_labels[i]=2*act.id+(left?0:1);
        logalloc+=left?act.logleft[k]:act.logright[k];
        if(left) nleft++;
      }
      if(!std::isfinite(logalloc)) continue;
      ppsti_apply(act,next,next_labels,nleft);
      std::vector<int> ids;
      for(const auto&kv:next) if(kv.second.axis>=0) ids.push_back(kv.first);
      std::sort(ids.begin(),ids.end());arma::mat next_splits(ids.size(),3);
      for(size_t k=0;k<ids.size();k++){
        next_splits(k,0)=ids[k];next_splits(k,1)=next.at(ids[k]).axis;
        next_splits(k,2)=next.at(ids[k]).cut;
      }
      double logq=act.logq+logalloc;
      double ratio=act.logratio[act.move==1?0:nleft];
      double logeta=logq+0.5*ratio;
      neighbors.push_back(List::create(_["splits"]=next_splits,
        _["labels"]=next_labels,_["move"]=act.move,_["log_q0"]=logq,
        _["log_ratio"]=ratio,_["log_eta"]=logeta,
        _["probability"]=std::exp(logeta-neighborhood.logZ)));
    }
  }
  List result(neighbors.size());
  for(size_t k=0;k<neighbors.size();k++) result[k]=neighbors[k];
  return List::create(_["log_normalizer"]=neighborhood.logZ,_["neighbors"]=result);
}

// Internal deterministic checks for the RAM equation and allocation-collapsed
// target. These helpers are exported to the package namespace, not its API.
// [[Rcpp::export]]
List ppstree_ram_inspect(arma::mat factor,arma::vec direction,
    double accept_prob,double target=0.234,double decay=0.7,int iteration=1){
  if(iteration<1||!std::isfinite(decay)||decay<=0.5||decay>1.0)
    stop("iteration must be positive and decay must be in (0.5, 1]");
  if(factor.n_rows==0||factor.n_rows!=factor.n_cols||
     direction.n_elem!=factor.n_rows||!factor.is_finite()||
     arma::any(factor.diag()<=0.0)||
     (factor.n_rows>1&&arma::accu(arma::abs(arma::trimatu(factor,1)))!=0.0))
    stop("factor must be a finite lower triangular matrix with positive diagonal");
  const double step=std::min(1.0,(double)factor.n_rows*std::pow((double)iteration,-decay));
  const bool success=ppst_ram_update(factor,direction,accept_prob,target,step);
  return List::create(_["factor"]=factor,
    _["covariance"]=arma::mat(factor*factor.t()),_["eta"]=step,_["success"]=success);
}

// [[Rcpp::export]]
List ppstree_pcg_inspect(arma::mat X,arma::mat region,arma::mat splits,
    arma::vec gate,arma::vec lambda,arma::vec a_gate,arma::vec b_gate,
    arma::vec gate_min,int gate_shared,int gate_family,bool cache_geometry=true){
  if(X.n_rows==0||X.n_cols==0||region.n_rows!=X.n_cols||region.n_cols!=2||
     splits.n_cols!=3||!X.is_finite()||!region.is_finite()||
     arma::any(region.col(1)<=region.col(0))||(gate_family<0||gate_family>2))
    stop("invalid pcg inspection geometry");
  gate=ppst_expand_positive(gate,X.n_cols,"gate");
  a_gate=ppst_expand_positive(a_gate,X.n_cols,"a_gate");
  b_gate=ppst_expand_positive(b_gate,X.n_cols,"b_gate");
  gate_min=ppst_expand_positive(gate_min,X.n_cols,"gate_min",true);
  if(arma::any(gate<=gate_min)||(gate_shared&&arma::any(gate!=gate[0])))
    stop("invalid pcg inspection gates");
  PPSTree tree;PPSTNode root;
  root.box=region;root.idx=arma::regspace<arma::uvec>(0,X.n_rows-1);
  root.axis=-1;root.cut=NA_REAL;root.depth=0;root.m=0;tree[1]=root;
  const arma::uvec order=arma::sort_index(splits.col(0));
  for(arma::uword k:order){
    if(!std::isfinite(splits(k,0))||!std::isfinite(splits(k,1))||
       splits(k,0)<1||splits(k,0)>(std::numeric_limits<int>::max()-1)/2||
       splits(k,1)<0||splits(k,1)>=(double)X.n_cols)
      stop("invalid pcg inspection split indices");
    const int id=(int)splits(k,0),axis=(int)splits(k,1);
    const double cut=splits(k,2);
    if(splits(k,0)!=id||splits(k,1)!=axis||tree.find(id)==tree.end()||
       tree.at(id).axis>=0||!std::isfinite(cut)||
       cut<=tree.at(id).box(axis,0)||cut>=tree.at(id).box(axis,1))
      stop("invalid pcg inspection split topology");
    PPSTNode parent=tree.at(id);arma::uvec left,right;
    ppst_split_indices(parent,X,axis,cut,left,right);
    tree[id].axis=axis;tree[id].cut=cut;
    tree[2*id]=ppst_child(parent,left,axis,cut,-1);
    tree[2*id+1]=ppst_child(parent,right,axis,cut,1);
  }
  std::vector<int> leaves;ppst_leaves(tree,leaves);std::sort(leaves.begin(),leaves.end());
  if(lambda.n_elem!=leaves.size()||!lambda.is_finite()||arma::any(lambda<=0.0))
    stop("lambda must be positive, with one value per leaf in ascending node order");
  PPSTGeometryCache geometry(X,region,gate,gate_family,nullptr,cache_geometry);
  PPSTGeometryCache*cache=cache_geometry?&geometry:nullptr;
  const arma::vec log_rate=arma::log(lambda);
  const PPSTPCGEvaluation value=ppst_pcg_evaluate(tree,leaves,log_rate,X,region,
    gate,a_gate,b_gate,gate_min,gate_shared,gate_family,cache);
  if(!std::isfinite(value.target)) stop("non-finite pcg inspection target");
  arma::mat probability(X.n_rows,leaves.size()),phi(X.n_rows,leaves.size());
  arma::vec H(leaves.size());std::vector<double> lw(leaves.size());
  for(size_t k=0;k<leaves.size();k++){
    H[k]=ppst_cached_exposure(tree.at(leaves[k]),region,gate,gate_family,cache);
    for(arma::uword i=0;i<X.n_rows;i++)
      phi(i,k)=std::exp(value.log_weights(i,k)-log_rate[k]);
  }
  for(arma::uword i=0;i<X.n_rows;i++){
    for(size_t k=0;k<leaves.size();k++) lw[k]=value.log_weights(i,k);
    const double den=value.log_normalizers[i];
    for(size_t k=0;k<leaves.size();k++) probability(i,k)=std::exp(lw[k]-den);
  }
  const int p=gate_shared?1:(int)gate.n_elem;
  return List::create(_["log_target"]=value.target,
    _["log_scale_target"]=value.target+arma::sum(arma::log(gate.head(p))),
    _["allocation_prob"]=probability,_["log_normalizers"]=value.log_normalizers,
    _["phi"]=phi,_["exposure"]=H,
    _["leaf_ids"]=leaves);
}

// Internal test/diagnostic helper: report ordinary cut support without sampling.
// [[Rcpp::export]]
List ppstree_cuts_inspect(arma::mat X,arma::mat region,arma::mat splits,
    double Dmax,int nmin,int cut_mode,int ncand,bool cache_cuts=true){
  const int depth=ppt_checked_depth(Dmax);
  if(X.n_rows==0||X.n_cols==0||region.n_rows!=X.n_cols||region.n_cols!=2||
     splits.n_cols!=3||!X.is_finite()||!region.is_finite()||
     arma::any(region.col(1)<=region.col(0))||nmin<1||ncand<2||
     cut_mode<0||cut_mode>2)
    stop("invalid cut inspection controls or geometry");
  PPSTree tree;PPSTNode root;
  root.box=region;root.idx=arma::regspace<arma::uvec>(0,X.n_rows-1);
  root.axis=-1;root.cut=NA_REAL;root.depth=0;root.m=0;tree[1]=root;
  const arma::uvec order=arma::sort_index(splits.col(0));
  for(arma::uword k:order){
    if(!std::isfinite(splits(k,0))||!std::isfinite(splits(k,1))||
       splits(k,0)<1||splits(k,0)>(std::numeric_limits<int>::max()-1)/2||
       splits(k,1)<0||splits(k,1)>=(double)X.n_cols)
      stop("invalid cut inspection split indices");
    const int id=(int)splits(k,0),axis=(int)splits(k,1);
    const double cut=splits(k,2);
    if(splits(k,0)!=id||splits(k,1)!=axis||tree.find(id)==tree.end()||
       tree.at(id).axis>=0||!std::isfinite(cut)||
       cut<=tree.at(id).box(axis,0)||cut>=tree.at(id).box(axis,1))
      stop("invalid cut inspection split topology");
    PPSTNode parent=tree.at(id);arma::uvec left,right;
    ppst_split_indices(parent,X,axis,cut,left,right);
    tree[id].axis=axis;tree[id].cut=cut;
    tree[2*id]=ppst_child(parent,left,axis,cut,-1);
    tree[2*id+1]=ppst_child(parent,right,axis,cut,1);
  }
  std::vector<int> ids;
  for(const auto&kv:tree) ids.push_back(kv.first);
  std::sort(ids.begin(),ids.end());
  PPSTCutsCache cuts(X,nmin,cut_mode,ncand);
  PPSTCutsCache*cache=cache_cuts?&cuts:nullptr;
  List out_cuts(ids.size());LogicalVector can_split(ids.size());
  // Repeat identical queries to exercise cache hits as well as cold values.
  for(int pass=0;pass<2;pass++) for(size_t k=0;k<ids.size();k++){
    const PPSTNode&node=tree.at(ids[k]);
    can_split[k]=ppst_can_split(node,X,depth,nmin,cut_mode,ncand,cache);
    List axis(X.n_cols);
    for(arma::uword j=0;j<X.n_cols;j++)
      axis[j]=ppst_cached_axis_cuts(node,X,j,cut_mode,ncand,nmin,cache);
    out_cuts[k]=axis;
  }
  cuts.trim(tree);
  return List::create(_["node_ids"]=ids,_["cuts"]=out_cuts,
    _["can_split"]=can_split,_["cache_entries"]=(double)cuts.size());
}

// Internal deterministic check for proposal labels, including zero support.
// parent_path columns are zero-based axis, cut, parent width, and side (-1/+1).
// [[Rcpp::export]]
List ppstree_side_inspect(arma::mat X,arma::mat region,arma::mat parent_path,
    int axis,double cut,double parent_width,arma::vec gate,int gate_family,
    bool cache_geometry=true){
  if(X.n_cols==0||region.n_rows!=X.n_cols||region.n_cols!=2||
     parent_path.n_cols!=4||!X.is_finite()||!region.is_finite()||
     !parent_path.is_finite()||arma::any(region.col(1)<=region.col(0))||
     axis<0||axis>=(int)X.n_cols||!std::isfinite(cut)||
     !std::isfinite(parent_width)||parent_width<=0.0||
     gate_family<0||gate_family>2)
    stop("invalid proposal-side inspection geometry");
  gate=ppst_expand_positive(gate,X.n_cols,"gate");
  PPSTNode parent;
  parent.box=region;parent.axis=-1;parent.cut=NA_REAL;
  parent.depth=parent_path.n_rows;parent.m=0;
  for(arma::uword k=0;k<parent_path.n_rows;k++){
    if(parent_path(k,0)<0.0||parent_path(k,0)>=(double)X.n_cols||
       parent_path(k,0)!=std::floor(parent_path(k,0))||parent_path(k,2)<=0.0||
       (parent_path(k,3)!=-1.0&&parent_path(k,3)!=1.0))
      stop("invalid proposal-side inspection parent path");
    PPSTGate g={(int)parent_path(k,0),parent_path(k,1),parent_path(k,2),
      (int)parent_path(k,3)};
    parent.path.push_back(g);
  }
  PPSTNode left=parent,right=parent;
  PPSTGate split={axis,cut,parent_width,1};right.path.push_back(split);
  split.side=-1;left.path.push_back(split);
  PPSTGeometryCache geometry(X,region,gate,gate_family,nullptr,cache_geometry);
  const arma::vec*prefix=cache_geometry?&geometry.training(parent):nullptr;
  arma::mat logs(X.n_rows,2);arma::vec probability(X.n_rows);
  for(arma::uword i=0;i<X.n_rows;i++){
    if(prefix) ppst_child_log_memberships(split,X(i,axis),(*prefix)[i],
      region,gate,gate_family,logs(i,0),logs(i,1));
    else{
      logs(i,0)=ppst_log_phi(left,X.row(i),region,gate,gate_family);
      logs(i,1)=ppst_log_phi(right,X.row(i),region,gate,gate_family);
    }
    probability[i]=ppst_side_probability(logs(i,0),logs(i,1));
  }
  return List::create(_["log_children"]=logs,_["p_left"]=probability);
}

// Internal exact-row-map check for the quadrature adapter. Core has no map.
// [[Rcpp::export]]
IntegerVector ppstree_training_background_inspect(arma::mat X,arma::mat region,
    arma::vec gate,int gate_family){
  if(X.n_cols==0||region.n_rows!=X.n_cols||region.n_cols!=2||
     !X.is_finite()||!region.is_finite()||
     arma::any(region.col(1)<=region.col(0))||gate_family<0||gate_family>2)
    stop("invalid training/background inspection geometry");
  gate=ppst_expand_positive(gate,X.n_cols,"gate");
  PPSTGeometryCache geometry(X,region,gate,gate_family);
  IntegerVector rows(X.n_rows,NA_INTEGER);
#ifdef POISTREE_SPATIAL_QUADRATURE_H
  const arma::uvec*map=geometry.training_background_rows();
  if(map) for(arma::uword i=0;i<X.n_rows;i++)
    if((*map)[i]<qpp_background.n_rows) rows[i]=(*map)[i]+1;
#endif
  return rows;
}
