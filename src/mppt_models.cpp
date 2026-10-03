// =====================================================================
// MPPSTree.cpp -- Multidimensional Poisson Process Soft Tree (MPPSTree)
//
// Non-destructive extension of MPPTree.cpp.  Eight model options:
//   model = 0 : MPPTree-Ind       (legacy hard boxes, independent xi)
//   model = 1 : MPPTree-Markov    (legacy hard boxes, Gamma-Markov xi)
//   model = 2 : MPPTree-Soft-Global      (full-domain logistic gates,
//                                         Gamma-Markov xi)
//   model = 3 : MPPTree-Soft-Global-Ind  (full-domain logistic gates,
//                                         independent xi)
//   model = 4 : MPPTree-Soft      (depth-relative soft gates, Gamma-Markov xi)
//   model = 5 : MPPTree-Soft-Ind  (depth-relative soft gates, independent xi)
//   model = 6 : MPPTree-Soft      (recursive logistic gates, Gamma-Markov xi)
//   model = 7 : MPPTree-Soft-Ind  (recursive logistic gates, independent xi)
//
// The soft model replaces 1{x in D_v} by a nonnegative separable basis
// phi_v(x).  On coordinate j, for global interval [a_j,b_j] and a node
// interval [l_j,u_j],
//
//   B_[l,u](x) = G_l(x) - G_u(x),
//   G_a(x)=1, G_b(x)=0, G_c(x)=logit^{-1}((x-c)/h_j),
//   h_j=(b_j-a_j)/gate_j.
//
// Thus phi_v(x)=prod_j B_[l_j,u_j](x_j), children sum exactly to their
// parent, and H_v=int_D phi_v is a product of analytic 1-D integrals.
// The unified R wrapper uses the recursive C-infinity logistic family (models
// 6 and 7) by default.  At cut c on coordinate j it uses
//   right(x)=logit^{-1}{gate_j*(x_j-c)/(b_j-a_j)}, left(x)=1-right(x).
// One slope is shared by all splits on a dimension, so path exposures reduce
// to analytic rational integrals.  By default one dimensionless inverse
// bandwidth is shared across coordinates. The dimension-specific option uses
// gate_j ~ Gamma(a_gate_j,b_gate_j), in the shape-rate convention.
//
// The opt-in compatibility option uses a recursive compact cubic gate at
// split s:
//   h_s = width(parent along split axis) /
//         {gate_j * (1 + parent_depth)^gate_depth},
//   phi_left=phi_parent*q_s, phi_right=phi_parent*(1-q_s).
// The cubic transition is nonnegative and compact on [s-h_s,s+h_s].  Products
// along a path remain separable by coordinate.  Their one-dimensional
// exposures are integrated exactly as piecewise polynomials, so arbitrary-d
// exposure remains a product of analytic 1-D integrals.  A recursive gate is
// necessary: putting different logistic bandwidths directly into
// G_lower-G_upper can make that difference negative.
//
// Given layer labels, conjugacy is unchanged:
//   lambda_v | ... ~ Gamma(kappa+m_v, kappa/xi_v + H_v).
// The spatial factor prod_i phi_{c_i}(x_i) is included in every soft
// label, tree, and gate update.
//
// Compile: Rcpp::sourceCpp("MPPSTree.cpp")
// =====================================================================
// [[Rcpp::depends(RcppArmadillo)]]
// Indirection prevents sourceCpp from compiling the legacy translation unit
// separately; the legacy implementation is included once into MPPSTree.cpp.
#define MPPSTREE_LEGACY_SOURCE "mppt_legacy.h"
#include MPPSTREE_LEGACY_SOURCE
#include "soft_logistic.h"

// ---- stable scalar helpers -------------------------------------------------
static inline double mpps_softplus(double z){
  return std::max(z,0.0)+std::log1p(std::exp(-std::abs(z)));
}
static inline double mpps_logsigmoid(double z){ return -mpps_softplus(-z); }
// log(1-exp(x)), x <= 0
static inline double mpps_log1mexp(double x){
  if(x>=0.0) return -arma::datum::inf;
  return x < -0.6931471805599453 ? std::log1p(-std::exp(x)) : std::log(-std::expm1(x));
}
static inline double mpps_lse2(double a,double b){
  double m=std::max(a,b); return m+std::log(std::exp(a-m)+std::exp(b-m));
}
static inline bool mpps_at_boundary(double x,double b,double width){
  return std::abs(x-b)<=1e-12*(1.0+std::abs(width));
}
static arma::vec mpps_expand_positive(const arma::vec&x,arma::uword d,
    const char*name){
  arma::vec out;
  if(x.n_elem==1) out=arma::vec(d,arma::fill::value(x[0]));
  else if(x.n_elem==d) out=x;
  else stop(std::string(name)+" must have length 1 or dimension d");
  if(!out.is_finite()||arma::any(out<=0.0))
    stop(std::string(name)+" must be finite and positive");
  return out;
}

// ---- soft box basis and exact exposure ------------------------------------
static double mpps_log_interval_basis(double x,double lo,double hi,
    double dom_lo,double dom_hi,double gate){
  double width=dom_hi-dom_lo, h=width/gate;
  bool at_lo=mpps_at_boundary(lo,dom_lo,width);
  bool at_hi=mpps_at_boundary(hi,dom_hi,width);
  if(at_lo && at_hi) return 0.0;                       // root interval
  if(at_lo){                                           // 1-sigma((x-hi)/h)
    return mpps_logsigmoid((hi-x)/h);
  }
  if(at_hi){                                           // sigma((x-lo)/h)
    return mpps_logsigmoid((x-lo)/h);
  }
  double zl=(x-lo)/h, zu=(x-hi)/h;                    // zl > zu
  // sigma(zl)-sigma(zu)
  return mpps_logsigmoid(zl)+mpps_logsigmoid(-zu)+mpps_log1mexp(zu-zl);
}

static double mpps_log_phi_global_node(const Node&nd,const arma::rowvec&x,
    const arma::mat&region,const arma::vec&gate){
  double ans=0.0;
  for(arma::uword j=0;j<region.n_rows;j++)
    ans += mpps_log_interval_basis(x[j],nd.box(j,0),nd.box(j,1),
                                  region(j,0),region(j,1),gate[j]);
  return ans;
}

static double mpps_integrated_boundary(double cut,double dom_lo,double dom_hi,
    double gate){
  double width=dom_hi-dom_lo;
  if(mpps_at_boundary(cut,dom_lo,width)) return width;
  if(mpps_at_boundary(cut,dom_hi,width)) return 0.0;
  double h=width/gate;
  return h*(mpps_softplus((dom_hi-cut)/h)-mpps_softplus((dom_lo-cut)/h));
}

static double mpps_exposure_global_node(const Node&nd,const arma::mat&region,
    const arma::vec&gate){
  double H=1.0;
  for(arma::uword j=0;j<region.n_rows;j++){
    double z=mpps_integrated_boundary(nd.box(j,0),region(j,0),region(j,1),gate[j])
            -mpps_integrated_boundary(nd.box(j,1),region(j,0),region(j,1),gate[j]);
    H *= std::max(z,1e-300);
  }
  return std::max(H,1e-300);
}

// ---- depth-relative recursive compact-cubic gate ---------------------------
// right(t)=3t^2-2t^3 on t in (0,1), with t=(x-(cut-h))/(2h);
// left=1-right.  Both are in [0,1] and sum to one exactly.
static inline double mpps_depth_gate_value(const SoftGate&g,double x,double gate){
  double h=g.parent_width/gate;
  double t=(x-(g.cut-h))/(2.0*h);
  double right=t<=0.0?0.0:(t>=1.0?1.0:t*t*(3.0-2.0*t));
  return g.side<0 ? 1.0-right : right;
}

static double mpps_log_phi_depth_node(const Node&nd,const arma::rowvec&x,
    const arma::vec&gate){
  double ans=0.0;
  for(const SoftGate&g:nd.soft_path){
    double q=mpps_depth_gate_value(g,x[g.axis],gate[g.axis]);
    if(q<=0.0) return -arma::datum::inf;
    ans+=std::log(q);
  }
  return ans;
}

static std::vector<double> mpps_poly_multiply(const std::vector<double>&a,
    const std::vector<double>&b){
  std::vector<double> out(a.size()+b.size()-1,0.0);
  for(size_t i=0;i<a.size();i++)
    for(size_t j=0;j<b.size();j++) out[i+j]+=a[i]*b[j];
  return out;
}

// Exact integral of the product of all recursive gates on one coordinate.
// Breakpoints make every factor polynomial on each segment.  Coefficients use
// u=(x-segment_lower)/segment_width in [0,1] for numerical stability.
static double mpps_depth_axis_integral(const Node&nd,int axis,double dom_lo,
    double dom_hi,double gate){
  std::vector<const SoftGate*> gates;
  std::vector<double> breaks={dom_lo,dom_hi};
  for(const SoftGate&g:nd.soft_path) if(g.axis==axis){
    gates.push_back(&g);
    double h=g.parent_width/gate;
    if(g.cut-h>dom_lo&&g.cut-h<dom_hi) breaks.push_back(g.cut-h);
    if(g.cut+h>dom_lo&&g.cut+h<dom_hi) breaks.push_back(g.cut+h);
  }
  if(gates.empty()) return dom_hi-dom_lo;
  std::sort(breaks.begin(),breaks.end());
  std::vector<double> uniq; uniq.reserve(breaks.size());
  for(double z:breaks)
    if(uniq.empty()||std::abs(z-uniq.back())>1e-13*(1.0+std::abs(z)))
      uniq.push_back(z);

  double total=0.0;
  for(size_t s=0;s+1<uniq.size();s++){
    double lo=uniq[s],hi=uniq[s+1],dx=hi-lo;
    if(dx<=0.0) continue;
    std::vector<double> poly(1,1.0);
    bool zero=false;
    for(const SoftGate*gp:gates){
      const SoftGate&g=*gp; double h=g.parent_width/gate;
      double tr_lo=g.cut-h,tr_hi=g.cut+h;
      if(hi<=tr_lo){
        if(g.side>0){ zero=true; break; } // right gate is zero
        continue;                         // left gate is one
      }
      if(lo>=tr_hi){
        if(g.side<0){ zero=true; break; } // left gate is zero
        continue;                         // right gate is one
      }
      // This whole segment lies inside the cubic transition.
      double t0=(lo-tr_lo)/(2.0*h),q=dx/(2.0*h);
      double t02=t0*t0,t03=t02*t0,q2=q*q,q3=q2*q;
      std::vector<double> right(4);
      right[0]=3.0*t02-2.0*t03;
      right[1]=6.0*t0*q-6.0*t02*q;
      right[2]=3.0*q2-6.0*t0*q2;
      right[3]=-2.0*q3;
      std::vector<double> factor=right;
      if(g.side<0){
        factor[0]=1.0-factor[0];
        for(size_t k=1;k<factor.size();k++) factor[k]=-factor[k];
      }
      poly=mpps_poly_multiply(poly,factor);
    }
    if(zero) continue;
    double integral=0.0;
    for(size_t k=0;k<poly.size();k++) integral+=poly[k]/(double)(k+1);
    total+=dx*integral;
  }
  if(total<0.0&&total>-1e-12*(dom_hi-dom_lo)) total=0.0;
  return std::max(total,1e-300);
}

static double mpps_exposure_depth_node(const Node&nd,const arma::mat&region,
    const arma::vec&gate){
  double H=1.0;
  for(arma::uword j=0;j<region.n_rows;j++)
    H*=mpps_depth_axis_integral(nd,(int)j,region(j,0),region(j,1),gate[j]);
  return std::max(H,1e-300);
}

// ---- recursive C-infinity logistic gate and analytic exposure --------------
static double mpps_log_phi_logistic_node(const Node&nd,
    const arma::rowvec&x,const arma::mat&region,const arma::vec&gate){
  double ans=0.0;
  for(const SoftGate&g:nd.soft_path){
    double width=region(g.axis,1)-region(g.axis,0);
    double z=gate[g.axis]*(x[g.axis]-g.cut)/width;
    ans+=g.side<0 ? pst_logistic_log_right(-z)
                  : pst_logistic_log_right(z);
  }
  return ans;
}

static double mpps_exposure_logistic_node(const Node&nd,
    const arma::mat&region,const arma::vec&gate){
  double H=1.0;
  for(arma::uword j=0;j<region.n_rows;j++){
    std::vector<double>cuts;
    std::vector<int>sides;
    for(const SoftGate&g:nd.soft_path) if(g.axis==(int)j){
      cuts.push_back(g.cut);
      sides.push_back(g.side);
    }
    H*=pst_logistic_path_axis_integral(
      cuts,sides,region(j,0),region(j,1),gate[j]
    );
  }
  return std::max(H,1e-300);
}

static double mpps_log_phi_node(const Node&nd,const arma::rowvec&x,
    const arma::mat&region,const arma::vec&gate,int gate_mode=0){
  if(gate_mode==1) return mpps_log_phi_depth_node(nd,x,gate);
  if(gate_mode==2) return mpps_log_phi_logistic_node(nd,x,region,gate);
  return mpps_log_phi_global_node(nd,x,region,gate);
}

static inline double mpps_phi_node(const Node&nd,const arma::rowvec&x,
    const arma::mat&region,const arma::vec&gate,int gate_mode=0){
  double z=mpps_log_phi_node(nd,x,region,gate,gate_mode);
  return z < -745.0 ? 0.0 : std::exp(z);
}

static double mpps_exposure_node(const Node&nd,const arma::mat&region,
    const arma::vec&gate,int gate_mode=0){
  if(gate_mode==1) return mpps_exposure_depth_node(nd,region,gate);
  if(gate_mode==2) return mpps_exposure_logistic_node(nd,region,gate);
  return mpps_exposure_global_node(nd,region,gate);
}

// [[Rcpp::export]]
List mppstree_soft_global_geometry(arma::mat box,arma::mat points,arma::mat region,
    arma::vec gate){
  gate=mpps_expand_positive(gate,region.n_rows,"gate");
  Node nd; nd.box=box;
  arma::vec phi(points.n_rows),logphi(points.n_rows);
  for(arma::uword i=0;i<points.n_rows;i++){
    logphi[i]=mpps_log_phi_global_node(nd,points.row(i),region,gate);
    phi[i]=logphi[i] < -745.0 ? 0.0 : std::exp(logphi[i]);
  }
  return List::create(_["phi"]=phi,_["log_phi"]=logphi,
                      _["H"]=mpps_exposure_global_node(nd,region,gate));
}

// Geometry checker for a recursive depth-relative path.  axis and side are
// 1-based and {-1,+1}, respectively, at the R interface.
// [[Rcpp::export]]
List mppstree_soft_geometry(IntegerVector axis,NumericVector cut,
    NumericVector parent_width,IntegerVector side,arma::mat points,
    arma::mat region,arma::vec gate,double gate_depth=0.0){
  int K=axis.size();
  gate=mpps_expand_positive(gate,region.n_rows,"gate");
  if(gate_depth<0.0) stop("gate_depth must be nonnegative");
  if(cut.size()!=K||parent_width.size()!=K||side.size()!=K)
    stop("axis, cut, parent_width, and side must have equal lengths");
  Node nd;
  for(int k=0;k<K;k++){
    if(axis[k]<1||axis[k]>(int)region.n_rows||parent_width[k]<=0.0||
       (side[k]!=-1&&side[k]!=1))
      stop("invalid recursive gate path");
    double effective_width=parent_width[k]/std::pow(1.0+(double)k,gate_depth);
    SoftGate g={(int)axis[k]-1,cut[k],effective_width,side[k]};
    nd.soft_path.push_back(g);
  }
  arma::vec phi(points.n_rows),logphi(points.n_rows);
  for(arma::uword i=0;i<points.n_rows;i++){
    logphi[i]=mpps_log_phi_depth_node(nd,points.row(i),gate);
    phi[i]=logphi[i] < -745.0 ? 0.0 : std::exp(logphi[i]);
  }
  return List::create(_["phi"]=phi,_["log_phi"]=logphi,
                      _["H"]=mpps_exposure_depth_node(nd,region,gate));
}

// Geometry checker for the default recursive logistic path.
// [[Rcpp::export]]
List mppstree_logistic_geometry(IntegerVector axis,NumericVector cut,
    IntegerVector side,arma::mat points,arma::mat region,arma::vec gate){
  int K=axis.size();
  gate=mpps_expand_positive(gate,region.n_rows,"gate");
  if(cut.size()!=K||side.size()!=K)
    stop("axis, cut, and side must have equal lengths");
  Node nd;
  for(int k=0;k<K;k++){
    if(axis[k]<1||axis[k]>(int)region.n_rows||
       (side[k]!=-1&&side[k]!=1))
      stop("invalid recursive logistic path");
    SoftGate g={(int)axis[k]-1,cut[k],
                region(axis[k]-1,1)-region(axis[k]-1,0),side[k]};
    nd.soft_path.push_back(g);
  }
  arma::vec phi(points.n_rows),logphi(points.n_rows);
  for(arma::uword i=0;i<points.n_rows;i++){
    logphi[i]=mpps_log_phi_logistic_node(nd,points.row(i),region,gate);
    phi[i]=logphi[i] < -745.0 ? 0.0 : std::exp(logphi[i]);
  }
  return List::create(_["phi"]=phi,_["log_phi"]=logphi,
                      _["H"]=mpps_exposure_logistic_node(nd,region,gate));
}

// ---- soft label sweep over every node -------------------------------------
static void mpps_label_sweep(Tree&T,std::vector<int>&c,const arma::mat&pts,
    const arma::mat&region,double kap,const arma::vec&gate,int gate_mode){
  std::vector<int> ids; ids.reserve(T.size());
  for(const auto&kv:T) ids.push_back(kv.first);
  std::vector<double> lw(ids.size()),exposure(ids.size());
  // T and gate are fixed throughout this Gibbs sweep.  In particular, the
  // depth-relative piecewise-polynomial H_v must not be recomputed n times.
  for(size_t k=0;k<ids.size();k++)
    exposure[k]=mpps_exposure_node(T.at(ids[k]),region,gate,gate_mode);
  for(arma::uword i=0;i<pts.n_rows;i++){
    int old=c[i]; T[old].m-=1; double mx=-arma::datum::inf;
    for(size_t k=0;k<ids.size();k++){
      Node&nd=T[ids[k]]; double H=exposure[k];
      lw[k]=mpps_log_phi_node(nd,pts.row(i),region,gate,gate_mode)
           +std::log(kap+nd.m)-std::log(kap/nd.xi+H);
      if(lw[k]>mx) mx=lw[k];
    }
    double tot=0.0; for(double z:lw) tot+=std::exp(z-mx);
    double u=unif_rand()*tot,cum=0.0; int pick=ids.back();
    for(size_t k=0;k<ids.size();k++){
      cum+=std::exp(lw[k]-mx); if(u<=cum){ pick=ids[k]; break; }
    }
    c[i]=pick; T[pick].m+=1;
  }
}

// ---- soft xi, kappa, and dimension-specific gate MH updates ---------------
static int mpps_xi_sweep(Tree&T,double kap,double tau,double a_xi,double b_xi,
    double sd_xi,const arma::mat&region,const arma::vec&gate,
    int scale_model,int gate_mode){
  int acc=0; double xi_mean=b_xi>0?a_xi/b_xi:1.0;
  for(auto&kv:T){
    int id=kv.first; Node&nd=kv.second; int L=2*id,Rn=2*id+1;
    bool internal=(nd.axis>=0); double H=mpps_exposure_node(nd,region,gate,gate_mode);
    auto lt=[&](double xi)->double{
      double s=lg(nd.m,H,kap,xi);
      if(scale_model==1){
        if(id==1) s+=lgam_dens(xi,a_xi,xi_mean);
        else s+=lgam_dens(xi,tau,T.at(id/2).xi);
        if(internal)
          s+=lgam_dens(T.at(L).xi,tau,xi)+lgam_dens(T.at(Rn).xi,tau,xi);
      }else{
        s+=lgam_dens(xi,a_xi,xi_mean);
      }
      return s;
    };
    double cur=nd.xi,prop=std::exp(std::log(cur)+R::rnorm(0.0,sd_xi));
    double la=lt(prop)-lt(cur)+std::log(prop)-std::log(cur);
    if(std::log(unif_rand())<la){ nd.xi=prop; acc++; }
  }
  return acc;
}

static int mpps_kappa_update(Tree&T,double&kap,double a_k,double b_k,double sd_k,
    const arma::mat&region,const arma::vec&gate,int gate_mode){
  double cur=kap,prop=std::exp(std::log(cur)+R::rnorm(0.0,sd_k));
  if(prop<=0.0) return 0;
  double dl=0.0;
  for(const auto&kv:T){
    double H=mpps_exposure_node(kv.second,region,gate,gate_mode);
    dl+=lg(kv.second.m,H,prop,kv.second.xi)-lg(kv.second.m,H,cur,kv.second.xi);
  }
  dl+=(a_k-1.0)*(std::log(prop)-std::log(cur))-b_k*(prop-cur)
     +std::log(prop)-std::log(cur);
  if(std::log(unif_rand())<dl){ kap=prop; return 1; }
  return 0;
}

static double mpps_gate_target(const Tree&T,const std::vector<int>&c,
    const arma::mat&pts,const arma::mat&region,double kap,const arma::vec&gate,
    const arma::vec&a_gate,const arma::vec&b_gate,int gate_mode,
    bool gate_shared){
  if(!gate.is_finite()||arma::any(gate<=0.0)) return -arma::datum::inf;
  double out=0.0;
  arma::uword prior_n=gate_shared?1:gate.n_elem;
  for(arma::uword j=0;j<prior_n;j++)
    out+=(a_gate[j]-1.0)*std::log(gate[j])-b_gate[j]*gate[j];
  for(const auto&kv:T){
    double H=mpps_exposure_node(kv.second,region,gate,gate_mode);
    out+=lg(kv.second.m,H,kap,kv.second.xi);
  }
  for(arma::uword i=0;i<pts.n_rows;i++)
    out+=mpps_log_phi_node(T.at(c[i]),pts.row(i),region,gate,gate_mode);
  return out;
}

static int mpps_gate_update(const Tree&T,const std::vector<int>&c,
    const arma::mat&pts,const arma::mat&region,double kap,arma::vec&gate,
    const arma::vec&a_gate,const arma::vec&b_gate,const arma::vec&sd_gate,
    int gate_mode,bool gate_shared,int&which){
  which=gate_shared?0:runif_int(gate.n_elem);
  arma::vec prop=gate;
  double cur=gate[which];
  double proposed=std::exp(std::log(cur)+R::rnorm(0.0,sd_gate[which]));
  if(gate_shared) prop.fill(proposed); else prop[which]=proposed;
  double la=mpps_gate_target(T,c,pts,region,kap,prop,a_gate,b_gate,gate_mode,
                            gate_shared)
           -mpps_gate_target(T,c,pts,region,kap,gate,a_gate,b_gate,gate_mode,
                             gate_shared)
           +std::log(prop[which])-std::log(cur);
  if(std::log(unif_rand())<la){ gate=prop; return 1; }
  return 0;
}

// ---- reversible label proposal used inside soft tree moves ----------------
static double mpps_alloc_logq(int dest,int parent,int L,int Rn,const arma::rowvec&pt,
    const Node&np,const Node&nl,const Node&nr,const arma::mat&region,
    const arma::vec&gate,int gate_mode){
  if(dest==parent) return std::log(0.5);
  double a=mpps_log_phi_node(nl,pt,region,gate,gate_mode);
  double b=mpps_log_phi_node(nr,pt,region,gate,gate_mode);
  double den=mpps_lse2(a,b);
  if(dest==L) return std::log(0.5)+a-den;
  if(dest==Rn) return std::log(0.5)+b-den;
  return -arma::datum::inf;
}

static int mpps_draw_alloc(int parent,int L,int Rn,const arma::rowvec&pt,
    const Node&np,const Node&nl,const Node&nr,const arma::mat&region,
    const arma::vec&gate,int gate_mode,double&logq){
  if(unif_rand()<0.5){ logq=std::log(0.5); return parent; }
  double a=mpps_log_phi_node(nl,pt,region,gate,gate_mode);
  double b=mpps_log_phi_node(nr,pt,region,gate,gate_mode);
  double den=mpps_lse2(a,b),pL=std::exp(a-den);
  if(unif_rand()<pL){ logq=std::log(0.5)+a-den; return L; }
  logq=std::log(0.5)+b-den; return Rn;
}

// ---- soft GROW / PRUNE -----------------------------------------------------
static int mpps_grow_prune(Tree&T,std::vector<int>&c,const arma::mat&pts,
    const arma::mat&region,double kap,double tau,double a_xi,double b_xi,
    int scale_model,const arma::vec&gate,int gate_mode,double gate_depth,
    double al,double eta,int Dmax,
    int nmin,int mode,int ncand,int&which_move){
  int d=pts.n_cols; std::vector<int> G,P;
  growable(T,pts,Dmax,nmin,mode,ncand,G); prunable(T,P);
  bool do_grow=unif_rand()<0.5; which_move=do_grow?0:1;
  if(do_grow){
    if(G.empty()) return 0;
    int v=G[runif_int(G.size())],dep=T[v].depth,L=2*v,Rn=2*v+1;
    int J=runif_int(d);
    std::vector<double> cuts=axis_cuts(pts,T[v].idx,T[v].box,J,mode,ncand,nmin);
    if(cuts.empty()) return 0;
    double cut=cuts[runif_int(cuts.size())],xip=T[v].xi;
    arma::vec col=pts.col(J); std::vector<arma::uword> li,ri;
    for(arma::uword t=0;t<T[v].idx.n_elem;t++){
      arma::uword ii=T[v].idx[t]; if(col[ii]<cut) li.push_back(ii); else ri.push_back(ii);
    }
    arma::uvec Lu(li.size()),Ru(ri.size());
    for(size_t t=0;t<li.size();t++) Lu[t]=li[t];
    for(size_t t=0;t<ri.size();t++) Ru[t]=ri[t];
    arma::mat boxL=T[v].box,boxR=T[v].box; boxL(J,1)=cut; boxR(J,0)=cut;
    Node nl,nr;
    double xiL=scale_model==1
      ? R::rgamma(tau,xip/tau)
          : R::rgamma(a_xi,1.0/b_xi);
    double xiR=scale_model==1
      ? R::rgamma(tau,xip/tau)
          : R::rgamma(a_xi,1.0/b_xi);
    nl.box=boxL; nl.idx=Lu; nl.cut=NA_REAL; nl.xi=xiL;
    nl.area=box_area(boxL); nl.axis=-1; nl.depth=dep+1; nl.m=0;
    nr.box=boxR; nr.idx=Ru; nr.cut=NA_REAL; nr.xi=xiR;
    nr.area=box_area(boxR); nr.axis=-1; nr.depth=dep+1; nr.m=0;
    nl.soft_path=T[v].soft_path; nr.soft_path=T[v].soft_path;
    if(gate_mode!=0){
      double pw=T[v].box(J,1)-T[v].box(J,0);
      if(gate_mode==1)
        pw/=std::pow(1.0+(double)dep,gate_depth);
      SoftGate gl={J,cut,pw,-1},gr={J,cut,pw,1};
      nl.soft_path.push_back(gl); nr.soft_path.push_back(gr);
    }
    int mold=T[v].m; Node oldp=T[v];
    T[L]=nl; T[Rn]=nr; T[v].axis=J; T[v].cut=cut;
    std::vector<int> Pn; prunable(T,Pn);
    int mv=0,mL=0,mR=0; double logq=0.0,dsp=0.0;
    std::vector<std::pair<int,int> > reass;
    for(arma::uword i=0;i<pts.n_rows;i++) if(c[i]==v){
      double lqi; int z=mpps_draw_alloc(v,L,Rn,pts.row(i),T[v],T[L],T[Rn],
                                        region,gate,gate_mode,lqi);
      logq+=lqi;
      dsp+=mpps_log_phi_node(T.at(z),pts.row(i),region,gate,gate_mode)
          -mpps_log_phi_node(oldp,pts.row(i),region,gate,gate_mode);
      reass.push_back(std::make_pair((int)i,z));
      if(z==v) mv++; else if(z==L) mL++; else mR++;
    }
    double Hp=mpps_exposure_node(T[v],region,gate,gate_mode);
    double HL=mpps_exposure_node(T[L],region,gate,gate_mode);
    double HR=mpps_exposure_node(T[Rn],region,gate,gate_mode);
    double logSL=lSfac(boxL,Lu,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
    double logSR=lSfac(boxR,Ru,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
    double logA=lrho(dep,al,eta)-l1mrho(dep,al,eta)+logSL+logSR
      +lg(mv,Hp,kap,xip)+lg(mL,HL,kap,T[L].xi)+lg(mR,HR,kap,T[Rn].xi)
      -lg(mold,Hp,kap,xip)+dsp
      +std::log((double)G.size())-std::log((double)Pn.size())-logq;
    if(std::log(unif_rand())<logA){
      T[v].m=mv; T[L].m=mL; T[Rn].m=mR;
      for(const auto&z:reass) c[z.first]=z.second;
      return 1;
    }else{
      T.erase(L); T.erase(Rn); T[v]=oldp;
    }
  }else{
    if(P.empty()) return 0;
    int v=P[runif_int(P.size())],dep=T[v].depth,L=2*v,Rn=2*v+1;
    Node oldp=T[v],oldL=T[L],oldR=T[Rn];
    int M=oldp.m+oldL.m+oldR.m;
    double logqrev=0.0,dsp=0.0;
    for(arma::uword i=0;i<pts.n_rows;i++) if(c[i]==v||c[i]==L||c[i]==Rn){
      logqrev+=mpps_alloc_logq(c[i],v,L,Rn,pts.row(i),oldp,oldL,oldR,
                              region,gate,gate_mode);
      dsp+=mpps_log_phi_node(oldp,pts.row(i),region,gate,gate_mode)
          -mpps_log_phi_node(T.at(c[i]),pts.row(i),region,gate,gate_mode);
    }
    double Hp=mpps_exposure_node(oldp,region,gate,gate_mode);
    double HL=mpps_exposure_node(oldL,region,gate,gate_mode);
    double HR=mpps_exposure_node(oldR,region,gate,gate_mode);
    double logSL=lSfac(oldL.box,oldL.idx,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
    double logSR=lSfac(oldR.box,oldR.idx,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
    T.erase(L); T.erase(Rn); T[v].axis=-1; T[v].cut=NA_REAL;
    std::vector<int> Gp; growable(T,pts,Dmax,nmin,mode,ncand,Gp);
    if(Gp.empty()){ T[v]=oldp; T[L]=oldL; T[Rn]=oldR; return 0; }
    double logA=-(lrho(dep,al,eta)-l1mrho(dep,al,eta)+logSL+logSR)
      +lg(M,Hp,kap,oldp.xi)-lg(oldp.m,Hp,kap,oldp.xi)
      -lg(oldL.m,HL,kap,oldL.xi)-lg(oldR.m,HR,kap,oldR.xi)+dsp
      +std::log((double)P.size())-std::log((double)Gp.size())+logqrev;
    if(std::log(unif_rand())<logA){
      for(size_t i=0;i<c.size();i++) if(c[i]==L||c[i]==Rn) c[i]=v;
      T[v].m=M;
      return 1;
    }else{
      T[v]=oldp; T[L]=oldL; T[Rn]=oldR;
    }
  }
  return 0;
}

// ---- soft CHANGE-CUT for a prunable node ----------------------------------
static int mpps_change_cut(Tree&T,std::vector<int>&c,const arma::mat&pts,
    const arma::mat&region,double kap,const arma::vec&gate,int gate_mode,
    double gate_depth,
    double al,double eta,
    int Dmax,int nmin,int mode,int ncand){
  int d=pts.n_cols; std::vector<int>P; prunable(T,P); if(P.empty()) return 0;
  int v=P[runif_int(P.size())],dep=T[v].depth,L=2*v,Rn=2*v+1,J=runif_int(d);
  std::vector<double> cuts=axis_cuts(pts,T[v].idx,T[v].box,J,mode,ncand,nmin);
  if(cuts.empty()) return 0;
  double cut=cuts[runif_int(cuts.size())];
  Node np=T[v],oldL=T[L],oldR=T[Rn],newL=oldL,newR=oldR;
  arma::vec col=pts.col(J); std::vector<arma::uword>li,ri;
  for(arma::uword t=0;t<T[v].idx.n_elem;t++){
    arma::uword ii=T[v].idx[t]; if(col[ii]<cut) li.push_back(ii); else ri.push_back(ii);
  }
  arma::uvec Lu(li.size()),Ru(ri.size());
  for(size_t t=0;t<li.size();t++) Lu[t]=li[t];
  for(size_t t=0;t<ri.size();t++) Ru[t]=ri[t];
  newL.box=np.box; newL.box(J,1)=cut; newL.idx=Lu; newL.area=box_area(newL.box); newL.m=0;
  newR.box=np.box; newR.box(J,0)=cut; newR.idx=Ru; newR.area=box_area(newR.box); newR.m=0;
  newL.soft_path=np.soft_path; newR.soft_path=np.soft_path;
  if(gate_mode!=0){
    double pw=np.box(J,1)-np.box(J,0);
    if(gate_mode==1)
      pw/=std::pow(1.0+(double)dep,gate_depth);
    SoftGate gl={J,cut,pw,-1},gr={J,cut,pw,1};
    newL.soft_path.push_back(gl); newR.soft_path.push_back(gr);
  }
  int mv=0,mL=0,mR=0; double logqf=0.0,logqr=0.0,dsp=0.0;
  std::vector<std::pair<int,int> >reass;
  for(arma::uword i=0;i<pts.n_rows;i++) if(c[i]==v||c[i]==L||c[i]==Rn){
    int old=c[i]; double lq;
    int z=mpps_draw_alloc(v,L,Rn,pts.row(i),np,newL,newR,region,gate,gate_mode,lq);
    logqf+=lq;
    logqr+=mpps_alloc_logq(old,v,L,Rn,pts.row(i),np,oldL,oldR,region,gate,gate_mode);
    const Node&nn=(z==v?np:(z==L?newL:newR));
    const Node&on=(old==v?np:(old==L?oldL:oldR));
    dsp+=mpps_log_phi_node(nn,pts.row(i),region,gate,gate_mode)
        -mpps_log_phi_node(on,pts.row(i),region,gate,gate_mode);
    reass.push_back(std::make_pair((int)i,z));
    if(z==v) mv++; else if(z==L) mL++; else mR++;
  }
  double Hp=mpps_exposure_node(np,region,gate,gate_mode);
  double HoL=mpps_exposure_node(oldL,region,gate,gate_mode);
  double HoR=mpps_exposure_node(oldR,region,gate,gate_mode);
  double HnL=mpps_exposure_node(newL,region,gate,gate_mode);
  double HnR=mpps_exposure_node(newR,region,gate,gate_mode);
  double oldg=lg(np.m,Hp,kap,np.xi)+lg(oldL.m,HoL,kap,oldL.xi)+lg(oldR.m,HoR,kap,oldR.xi);
  double newg=lg(mv,Hp,kap,np.xi)+lg(mL,HnL,kap,newL.xi)+lg(mR,HnR,kap,newR.xi);
  double oldS=lSfac(oldL.box,oldL.idx,dep+1,pts,Dmax,nmin,mode,ncand,al,eta)
             +lSfac(oldR.box,oldR.idx,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
  double newS=lSfac(newL.box,newL.idx,dep+1,pts,Dmax,nmin,mode,ncand,al,eta)
             +lSfac(newR.box,newR.idx,dep+1,pts,Dmax,nmin,mode,ncand,al,eta);
  if(std::log(unif_rand())<(newg-oldg)+(newS-oldS)+dsp+logqr-logqf){
    T[v].axis=J; T[v].cut=cut; T[v].m=mv;
    T[L]=newL; T[L].m=mL; T[Rn]=newR; T[Rn].m=mR;
    for(const auto&z:reass) c[z.first]=z.second;
    return 1;
  }
  return 0;
}

#include "irjmcmc_moves.h"

// ---- intensity draw and log likelihood under the soft basis ---------------
static void mpps_draw_lambda(const Tree&T,const arma::mat&region,double kap,
    const arma::vec&gate,int gate_mode,
    std::unordered_map<int,double>&lam,double&comp){
  lam.clear(); lam.reserve(T.size()); comp=0.0;
  for(const auto&kv:T){
    double H=mpps_exposure_node(kv.second,region,gate,gate_mode);
    double lv=R::rgamma(kap+kv.second.m,1.0/(kap/kv.second.xi+H));
    lam[kv.first]=lv; comp+=lv*H;
  }
}

static double mpps_intensity(const Tree&T,const std::unordered_map<int,double>&lam,
    const arma::rowvec&x,const arma::mat&region,const arma::vec&gate,
    int gate_mode){
  double out=0.0;
  for(const auto&kv:T)
    out+=lam.at(kv.first)*mpps_phi_node(kv.second,x,region,gate,gate_mode);
  return std::max(out,1e-300);
}

struct MPPSGateConfig {
  arma::vec gate;
  arma::vec a;
  arma::vec b;
  arma::vec sd;
  double depth;
  bool shared;
};

// Legacy hp has 16 entries and denotes one shared gate.  The dimension-
// specific layout is [11 scale entries, gate[d], a[d], b[d], sd[d], depth].
static MPPSGateConfig mpps_parse_gate_config(const arma::vec&hp,arma::uword d){
  MPPSGateConfig z;
  if(hp.n_elem==16){
    z.gate=arma::vec(d,arma::fill::value(hp[11]));
    z.a=arma::vec(d,arma::fill::value(hp[12]));
    z.b=arma::vec(d,arma::fill::value(hp[13]));
    z.sd=arma::vec(d,arma::fill::value(hp[14]));
    z.depth=hp[15]; z.shared=true;
  }else{
    arma::uword need=12+4*d;
    if(hp.n_elem!=need)
      stop("dimension-specific hp must have 12 + 4*d entries");
    arma::uword k=11;
    z.gate=hp.subvec(k,k+d-1); k+=d;
    z.a=hp.subvec(k,k+d-1); k+=d;
    z.b=hp.subvec(k,k+d-1); k+=d;
    z.sd=hp.subvec(k,k+d-1); k+=d;
    z.depth=hp[k]; z.shared=false;
  }
  z.gate=mpps_expand_positive(z.gate,d,"gate");
  z.a=mpps_expand_positive(z.a,d,"a_gate");
  z.b=mpps_expand_positive(z.b,d,"b_gate");
  z.sd=mpps_expand_positive(z.sd,d,"sd_gate");
  if(!std::isfinite(z.depth)||z.depth<0.0)
    stop("gate_depth must be finite and nonnegative");
  return z;
}

// ---- one MPPTree-Soft chain ------------------------------------------------
static int mpps_run_soft_chain(const arma::mat&pts,const arma::mat&grid,
    const arma::mat&xt,const arma::mat&region,double kappa0,double tau0,
    double a_xi,double b_xi,double sd_xi,const arma::vec&gate0,
    const arma::vec&a_gate,const arma::vec&b_gate,const arma::vec&sd_gate,
    bool gate_shared,double gate_depth,double al,double eta,
    int Dmax,int nmin,
    int iters,int burn,int thin,int nmove,int ncc,int cut_mode,int ncand,
    int update_hyper,int update_gate,int scale_model,int gate_mode,double a_k,double b_k,double sd_k,double a_t,
    double b_t,double sd_t,arma::mat&Draws,arma::vec&loglik,
    arma::vec&loglik_test,arma::vec&integrated_intensity,int row0,
    double&nleaf_mean,double&kap_mean,
    double&tau_mean,arma::vec&gate_mean,arma::vec&gate_accept,
    double&maxdepth_mean,arma::vec&tree_accept,
    std::vector<arma::mat>&state_nodes,arma::mat&state_gate,
    int proposal_mode=0){
  int n=pts.n_rows,nt=xt.n_rows; arma::uword ng=grid.n_rows;
  double kap=kappa0,tau=tau0,xi0=a_xi/std::max(b_xi,1e-12);
  arma::vec gate=gate0;
  Tree T; T[1]=make_root(region,n,xi0); std::vector<int>c(n,1);
  double nls=0.0,ds=0.0,ks=0.0,ts=0.0; int si=0;
  arma::vec gs(gate.n_elem,arma::fill::zeros);
  arma::vec gacc(gate.n_elem,arma::fill::zeros);
  arma::vec gtot(gate.n_elem,arma::fill::zeros);
  double grow_acc=0.0,grow_tot=0.0,prune_acc=0.0,prune_tot=0.0;
  double change_acc=0.0,change_tot=0.0;
  for(int it=0;it<iters;it++){
    mpps_label_sweep(T,c,pts,region,kap,gate,gate_mode);
    mpps_xi_sweep(T,kap,tau,a_xi,b_xi,sd_xi,region,gate,scale_model,gate_mode);
    if(update_hyper){
      mpps_kappa_update(T,kap,a_k,b_k,sd_k,region,gate,gate_mode);
      if(scale_model==1) tau_update(T,tau,a_t,b_t,sd_t);
    }
    if(update_gate){
      int which=0;
      int accepted=mpps_gate_update(T,c,pts,region,kap,gate,a_gate,b_gate,
                                    sd_gate,gate_mode,gate_shared,which);
      gacc[which]+=accepted; gtot[which]+=1.0;
    }
    for(int r=0;r<nmove;r++){
      int which=-1;
      int accepted=proposal_mode==0
        ? mpps_grow_prune(
            T,c,pts,region,kap,tau,a_xi,b_xi,scale_model,gate,gate_mode,
            gate_depth,al,eta,Dmax,nmin,cut_mode,ncand,which
          )
        : irj_soft_grow_prune(
            T,c,pts,region,kap,tau,a_xi,b_xi,scale_model,gate,gate_mode,
            gate_depth,al,eta,Dmax,nmin,cut_mode,ncand,which,
            proposal_mode==3
          );
      if(which==0){ grow_acc+=accepted; grow_tot+=1.0; }
      else { prune_acc+=accepted; prune_tot+=1.0; }
    }
    for(int r=0;r<ncc;r++){
      change_acc+=proposal_mode==0
        ? mpps_change_cut(
            T,c,pts,region,kap,gate,gate_mode,gate_depth,al,eta,Dmax,nmin,
            cut_mode,ncand
          )
        : irj_soft_change_cut(
            T,c,pts,region,kap,gate,gate_mode,gate_depth,al,eta,Dmax,nmin,
            cut_mode,ncand,proposal_mode==3
          );
      change_tot+=1.0;
    }
    if(it>=burn&&(it-burn)%thin==0){
      int row=row0+si; std::unordered_map<int,double>lam; double comp;
      mpps_draw_lambda(T,region,kap,gate,gate_mode,lam,comp);
      integrated_intensity[row]=comp;
      for(arma::uword j=0;j<ng;j++)
        Draws(row,j)=mpps_intensity(T,lam,grid.row(j),region,gate,gate_mode);
      double ll=-comp;
      for(int i=0;i<n;i++)
        ll+=std::log(mpps_intensity(T,lam,pts.row(i),region,gate,gate_mode));
      loglik[row]=ll;
      if(nt>0){
        double llt=-comp;
        for(int i=0;i<nt;i++)
          llt+=std::log(mpps_intensity(T,lam,xt.row(i),region,gate,gate_mode));
        loglik_test[row]=llt;
      }
      int nl=0,md=0;
      for(const auto&kv:T){
        if(kv.second.axis<0) nl++;
        if(kv.second.depth>md) md=kv.second.depth;
      }
      // serialize the generative state of this retained draw: one row per
      // node, columns (heap id, axis, cut, lambda, xi, m); every node of
      // the additive multiscale tree carries a rate, and the gate vector
      // of the draw is stored alongside.
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
        state_gate.row(row)=gate.t();
      }
      nls+=nl; ds+=md; ks+=kap; ts+=tau; gs+=gate; si++;
    }
    if((it&1023)==0) Rcpp::checkUserInterrupt();
  }
  nleaf_mean=nls/std::max(1,si); kap_mean=ks/std::max(1,si);
  tau_mean=scale_model==1 ? ts/std::max(1,si) : NA_REAL;
  gate_mean=gs/std::max(1,si);
  gate_accept.set_size(gate.n_elem);
  if(gate_shared){
    double a=gacc[0]/std::max(1.0,gtot[0]);
    gate_accept.fill(a);
  }else{
    for(arma::uword j=0;j<gate.n_elem;j++)
      gate_accept[j]=gacc[j]/std::max(1.0,gtot[j]);
  }
  maxdepth_mean=ds/std::max(1,si);
  tree_accept.set_size(3);
  tree_accept[0]=grow_acc/std::max(1.0,grow_tot);
  tree_accept[1]=prune_acc/std::max(1.0,prune_tot);
  tree_accept[2]=change_acc/std::max(1.0,change_tot);
  return si;
}

static int mpps_run_chain(const arma::mat&pts,const arma::mat&grid,const arma::mat&xt,
    const arma::mat&region,const arma::vec&hp,int update_hyper,int update_gate,
    int model,
    double al,double eta,int Dmax,int nmin,int iters,int burn,int thin,int nmove,
    int ncc,int cut_mode,int ncand,arma::mat&Draws,arma::vec&loglik,
    arma::vec&loglik_test,arma::vec&integrated_intensity,int row0,
    double&nl,double&km,double&tm,arma::vec&gm,
    arma::vec&ga,double&dm,arma::vec&tree_accept,
    std::vector<arma::mat>&state_nodes,arma::mat&state_gate,
    int proposal_mode=0){
  if(model<2){
    gm=arma::vec(region.n_rows,arma::fill::value(NA_REAL));
    ga=arma::vec(region.n_rows,arma::fill::value(NA_REAL));
    return run_chain(pts,grid,xt,region,hp[0],hp[1],hp[2],hp[3],hp[4],
      al,eta,Dmax,nmin,iters,burn,thin,nmove,ncc,cut_mode,ncand,update_hyper,
      hp[5],hp[6],hp[7],hp[8],hp[9],hp[10],model,Draws,loglik,loglik_test,
      integrated_intensity,row0,nl,dm,km,tm,tree_accept,state_nodes,
      proposal_mode);
  }
  int scale_model=(model==2||model==4||model==6)?1:0;
  int gate_mode=(model==4||model==5)?1:
                ((model==6||model==7)?2:0);
  MPPSGateConfig g=mpps_parse_gate_config(hp,region.n_rows);
  return mpps_run_soft_chain(pts,grid,xt,region,hp[0],hp[1],hp[2],hp[3],hp[4],
    g.gate,g.a,g.b,g.sd,g.shared,g.depth,al,eta,Dmax,nmin,iters,burn,thin,nmove,ncc,
    cut_mode,ncand,update_hyper,update_gate,scale_model,gate_mode,hp[5],hp[6],
    hp[7],hp[8],hp[9],hp[10],Draws,loglik,loglik_test,
    integrated_intensity,row0,nl,km,tm,gm,ga,dm,tree_accept,
    state_nodes,state_gate,proposal_mode);
}

// [[Rcpp::export]]
List mppstree_chain(arma::mat X,arma::mat grid,arma::mat Xtest,arma::mat region,
    arma::vec hp,int update_hyper,int update_gate,int model,double al,double eta,
    int Dmax,int nmin,
    int iters,int burn,int thin,int nmove,int ncc,int cut_mode,int ncand){
  if(hp.n_elem<11) stop("hp is too short; use MPPSTree_fit.R");
  if(model<0||model>7)
    stop("model must be 0 (ind), 1 (markov), 2 (soft_global), "
         "3 (soft_global_ind), 4 (compact), 5 (compact_ind), "
         "6 (soft), or 7 (soft_ind)");
  GetRNGstate();
  int ns=0; for(int it=burn;it<iters;it++) if((it-burn)%thin==0) ns++;
  arma::mat D(ns,grid.n_rows,arma::fill::zeros);
  arma::vec ll(ns,arma::fill::zeros),llt(Xtest.n_rows>0?ns:0,arma::fill::zeros);
  arma::vec integrated_intensity(ns,arma::fill::zeros);
  double nl,km,tm,dm; arma::vec gm,ga,tree_accept;
  std::vector<arma::mat> state_nodes(ns);
  arma::mat state_gate(ns,X.n_cols,arma::fill::zeros);
  int got=mpps_run_chain(X,grid,Xtest,region,hp,update_hyper,update_gate,model,
    al,eta,Dmax,nmin,iters,burn,thin,nmove,ncc,cut_mode,ncand,D,ll,llt,
    integrated_intensity,0,nl,km,tm,gm,ga,dm,tree_accept,
    state_nodes,state_gate);
  PutRNGstate();
  arma::vec mn,md,lo,hi; summarize(D,mn,md,lo,hi);
  double lp=(Xtest.n_rows>0&&got>0)?logsumexp(llt)-std::log((double)got):NA_REAL;
  List sn(got);
  for(int s=0;s<got;s++) sn[s]=state_nodes[s];
  List out=List::create(_["mean"]=mn,_["median"]=md,_["lower95"]=lo,_["upper95"]=hi,
    _["loglik_mean"]=got>0?arma::mean(ll):NA_REAL,_["logpred"]=lp,
    _["integrated_intensity"]=integrated_intensity,
    _["kappa_mean"]=km,
    _["tau_mean"]=(model==0||model==3||model==5||model==7)?NA_REAL:tm,
    _["gate_mean"]=gm,
    _["gate_accept"]=ga,_["tree_accept"]=tree_accept,
    _["mean_leaves"]=nl,_["mean_max_depth"]=dm,
    _["ndraws"]=got);
  out["state_nodes"]=sn;
  out["state_gate"]=model>=2
    ?state_gate.rows(0,std::max(0,got-1))
    :arma::mat(std::max(1,got),X.n_cols,arma::fill::value(NA_REAL));
  return out;
}

// [[Rcpp::export]]
List mppstree_multi(arma::mat X,arma::mat grid,arma::mat Xtest,arma::mat region,
    arma::vec hp,int update_hyper,int update_gate,int model,double al,double eta,
    int Dmax,int nmin,
    int iters,int burn,int thin,int nmove,int ncc,int cut_mode,int ncand,
    int chains,int verbose,int proposal_mode=0){
  if(hp.n_elem<11) stop("hp is too short; use MPPSTree_fit.R");
  if(model<0||model>7)
    stop("model must be 0 (ind), 1 (markov), 2 (soft_global), "
         "3 (soft_global_ind), 4 (compact), 5 (compact_ind), "
         "6 (soft), or 7 (soft_ind)");
  if(proposal_mode!=0&&proposal_mode!=3)
    stop("proposal_mode must be 0 (RJ-MCMC) or 3 (combined iRJ-MCMC)");
  GetRNGstate();
  int ns=0; for(int it=burn;it<iters;it++) if((it-burn)%thin==0) ns++;
  int total=ns*chains,row=0;
  arma::mat D(total,grid.n_rows,arma::fill::zeros);
  arma::vec ll(total,arma::fill::zeros),llt(Xtest.n_rows>0?total:0,arma::fill::zeros);
  arma::vec integrated_intensity(total,arma::fill::zeros);
  arma::vec nlv(chains),dmv(chains),kv(chains),tv(chains);
  arma::mat gv(chains,X.n_cols,arma::fill::value(NA_REAL));
  arma::mat gav(chains,X.n_cols,arma::fill::value(NA_REAL));
  arma::mat tav(chains,3,arma::fill::zeros);
  std::vector<arma::mat> state_nodes(total);
  arma::mat state_gate(total,X.n_cols,arma::fill::zeros);
  const char*name=(model==5||model==7||model==3)?"S-MPPT":
                  ((model==4||model==6||model==2)?"S-MPPT-Markov":
                  (model==1?"MPPT-Markov":"MPPT"));
  const char*algorithm=proposal_mode==0?"RJ-MCMC":"iRJ-MCMC";
  if(verbose) Rcpp::Rcout<<name<<" ["<<algorithm<<"]: "<<chains
                         <<" chains, "<<ns<<" draws/chain\n";
  for(int k=0;k<chains;k++){
    double nl,km,tm,dm; arma::vec gm,ga,tree_accept;
    int got=mpps_run_chain(X,grid,Xtest,region,hp,update_hyper,update_gate,model,
      al,eta,Dmax,nmin,iters,burn,thin,nmove,ncc,cut_mode,ncand,D,ll,llt,
      integrated_intensity,row,nl,km,tm,gm,ga,dm,tree_accept,
      state_nodes,state_gate,proposal_mode);
    row+=got; nlv[k]=nl; dmv[k]=dm; kv[k]=km; tv[k]=tm;
    gv.row(k)=gm.t(); gav.row(k)=ga.t();
    tav.row(k)=tree_accept.t();
    if(verbose){
      Rcpp::Rcout<<"  chain "<<k+1<<"/"<<chains<<" done; leaves="<<nl
                 <<", max_depth="<<dm
                 <<", kappa="<<km;
      if(model==1||model==2||model==4||model==6)
        Rcpp::Rcout<<", tau="<<tm;
      if(model>=2) Rcpp::Rcout<<", gate="<<arma::mean(gm)
                              <<", gate_accept="<<arma::mean(ga);
      Rcpp::Rcout<<"\n";
    }
    Rcpp::checkUserInterrupt();
  }
  PutRNGstate();
  arma::vec mn,md,lo,hi; summarize(D,mn,md,lo,hi);
  double lp=(Xtest.n_rows>0&&row>0)?logsumexp(llt)-std::log((double)row):NA_REAL;
  List out=List::create(_["mean"]=mn,_["median"]=md,_["lower95"]=lo,_["upper95"]=hi,
    _["draws"]=D.t(),
    _["loglik_mean"]=row>0?arma::mean(ll):NA_REAL,_["logpred"]=lp,
    _["loglik_draws"]=ll,
    _["integrated_intensity"]=integrated_intensity,
    _["chain_mean_leaves"]=nlv,
    _["chain_mean_max_depth"]=dmv,
    _["kappa_mean"]=arma::mean(kv),
    _["tau_mean"]=(model==0||model==3||model==5||model==7)
      ?NA_REAL:arma::mean(tv),
    _["gate_mean"]=model>=2?arma::mean(gv,0).t():
      arma::vec(X.n_cols,arma::fill::value(NA_REAL)),
    _["gate_accept"]=model>=2?arma::mean(gav,0).t():
      arma::vec(X.n_cols,arma::fill::value(NA_REAL)),
    _["tree_accept"]=arma::mean(tav,0).t(),
    _["mean_leaves"]=arma::mean(nlv),
    _["mean_max_depth"]=arma::mean(dmv),
    _["ndraws"]=row);
  List sn(row);
  for(int s=0;s<row;s++) sn[s]=state_nodes[s];
  out["state_nodes"]=sn;
  out["state_gate"]=model>=2
    ?state_gate.rows(0,std::max(0,row-1))
    :arma::mat(std::max(1,row),X.n_cols,arma::fill::value(NA_REAL));
  return out;
}

// [[Rcpp::export]]
List mppstree_diag(arma::mat X,arma::mat mon,arma::mat region,arma::vec hp,
    int update_hyper,int update_gate,int model,double al,double eta,int Dmax,
    int nmin,int iters,int burn,int thin,int nmove,int ncc,int cut_mode,
    int ncand){
  if(model<0||model>7)
    stop("model must be 0 (ind), 1 (markov), 2 (soft_global), "
         "3 (soft_global_ind), 4 (compact), 5 (compact_ind), "
         "6 (soft), or 7 (soft_ind)");
  if(model<2){
    arma::vec oldhp=hp.subvec(0,10);
    return mpptree_diag(X,mon,region,oldhp,update_hyper,model,al,eta,Dmax,nmin,
      iters,burn,thin,nmove,ncc,cut_mode,ncand);
  }
  if(hp.n_elem<11) stop("hp is too short; use MPPSTree_fit.R");
  GetRNGstate();
  int n=X.n_rows,ns=0; for(int it=burn;it<iters;it++) if((it-burn)%thin==0) ns++;
  int scale_model=(model==2||model==4||model==6)?1:0;
  int gate_mode=(model==4||model==5)?1:
                ((model==6||model==7)?2:0);
  MPPSGateConfig gc=mpps_parse_gate_config(hp,region.n_rows);
  double kap=hp[0],tau=hp[1],xi0=hp[2]/std::max(hp[3],1e-12);
  arma::vec gate=gc.gate;
  Tree T; T[1]=make_root(region,n,xi0); std::vector<int>c(n,1);
  arma::vec trn(ns),trk(ns),trt(ns),trl(ns);
  arma::mat trg(ns,gate.n_elem,arma::fill::zeros);
  arma::mat trm(ns,mon.n_rows,arma::fill::zeros);
  double xa=0,xtot=0,ka=0,kt=0,ta=0,tt=0; int si=0;
  double grow_acc=0,grow_tot=0,prune_acc=0,prune_tot=0;
  double change_acc=0,change_tot=0;
  arma::vec ga(gate.n_elem,arma::fill::zeros);
  arma::vec gt(gate.n_elem,arma::fill::zeros);
  for(int it=0;it<iters;it++){
    mpps_label_sweep(T,c,X,region,kap,gate,gate_mode);
    xa+=mpps_xi_sweep(T,kap,tau,hp[2],hp[3],hp[4],region,gate,scale_model,
                      gate_mode);
    xtot+=T.size();
    if(update_hyper){
      ka+=mpps_kappa_update(T,kap,hp[5],hp[6],hp[7],region,gate,gate_mode); kt++;
      if(scale_model==1){ ta+=tau_update(T,tau,hp[8],hp[9],hp[10]); tt++; }
    }
    if(update_gate){
      int which=0;
      int accepted=mpps_gate_update(T,c,X,region,kap,gate,gc.a,gc.b,gc.sd,
                                    gate_mode,gc.shared,which);
      ga[which]+=accepted; gt[which]++;
    }
    for(int r=0;r<nmove;r++){
      int which=-1;
      int accepted=mpps_grow_prune(
        T,c,X,region,kap,tau,hp[2],hp[3],scale_model,gate,gate_mode,
        gc.depth,al,eta,Dmax,nmin,cut_mode,ncand,which
      );
      if(which==0){ grow_acc+=accepted; grow_tot++; }
      else { prune_acc+=accepted; prune_tot++; }
    }
    for(int r=0;r<ncc;r++){
      change_acc+=mpps_change_cut(
        T,c,X,region,kap,gate,gate_mode,gc.depth,al,eta,Dmax,nmin,cut_mode,
        ncand
      );
      change_tot++;
    }
    if(it>=burn&&(it-burn)%thin==0){
      std::unordered_map<int,double>lam; double comp;
      mpps_draw_lambda(T,region,kap,gate,gate_mode,lam,comp);
      for(arma::uword j=0;j<mon.n_rows;j++)
        trm(si,j)=mpps_intensity(T,lam,mon.row(j),region,gate,gate_mode);
      double ll=-comp;
      for(int i=0;i<n;i++)
        ll+=std::log(mpps_intensity(T,lam,X.row(i),region,gate,gate_mode));
      int nl=0; for(const auto&kv:T) if(kv.second.axis<0) nl++;
      trn[si]=nl; trk[si]=kap; trt[si]=scale_model==1?tau:NA_REAL;
      trg.row(si)=gate.t(); trl[si]=ll; si++;
    }
    if((it&1023)==0) Rcpp::checkUserInterrupt();
  }
  PutRNGstate();
  arma::vec acc(3); acc[0]=xa/std::max(1.0,xtot); acc[1]=ka/std::max(1.0,kt);
  acc[2]=scale_model==1?ta/std::max(1.0,tt):NA_REAL;
  arma::vec gate_acc(gate.n_elem);
  if(gc.shared){
    gate_acc.fill(ga[0]/std::max(1.0,gt[0]));
  }else{
    for(arma::uword j=0;j<gate.n_elem;j++)
      gate_acc[j]=ga[j]/std::max(1.0,gt[j]);
  }
  arma::vec tree_accept(3);
  tree_accept[0]=grow_acc/std::max(1.0,grow_tot);
  tree_accept[1]=prune_acc/std::max(1.0,prune_tot);
  tree_accept[2]=change_acc/std::max(1.0,change_tot);
  return List::create(_["nleaf"]=trn,_["kappa"]=trk,_["tau"]=trt,_["gate"]=trg,
    _["logdens"]=trl,_["mon"]=trm,_["accept"]=acc,
    _["gate_accept"]=gate_acc,_["tree_accept"]=tree_accept,_["ns"]=si);
}

// Deterministic depth-one guardrail for the informed grow/prune algebra.  It
// compares the log Metropolis ratio used by the kernel with the same ratio
// assembled independently as log pi(split) + log q(prune) minus
// log pi(root) - log q(grow).  This internal export is exercised by package
// tests for hard, compact-soft, and logistic-soft independent-scale models.
// [[Rcpp::export]]
List mppstree_irj_balance_check(arma::mat X,arma::mat region,int soft,
    arma::vec gate,int gate_mode,double kappa,double a_xi,double b_xi,
    double alpha,double eta,int cut_mode,int ncand){
  if(X.n_rows<4||X.n_cols!=region.n_rows)
    stop("balance check requires at least four points with matching dimension");
  if(!(kappa>0.0&&a_xi>0.0&&b_xi>0.0))
    stop("balance-check scale parameters must be positive");
  if(soft){
    if(gate_mode!=1&&gate_mode!=2)
      stop("soft balance check requires compact or logistic gating");
    gate=mpps_expand_positive(gate,region.n_rows,"gate");
  }else{
    gate=arma::vec(region.n_rows,arma::fill::ones);
    gate_mode=0;
  }

  int n=X.n_rows,d=X.n_cols,min_leaf_n=1,max_depth=1;
  double xi_parent=a_xi/b_xi;
  Tree root_tree;
  root_tree[1]=make_root(region,n,xi_parent);
  int axis=-1;
  std::vector<double> cuts;
  for(int j=0;j<d;j++){
    cuts=axis_cuts(X,root_tree[1].idx,root_tree[1].box,j,cut_mode,ncand,
                   min_leaf_n);
    if(!cuts.empty()){ axis=j; break; }
  }
  if(axis<0) stop("balance-check data have no valid depth-one cut");
  int chosen=(int)cuts.size()/2;
  double cut=cuts[chosen];
  double xi_left=0.75*xi_parent,xi_right=1.25*xi_parent;

  Node parent=root_tree[1],left_node,right_node;
  if(soft){
    irj_soft_children(parent,axis,cut,gate_mode,0.0,left_node,right_node);
  }else{
    left_node=Node(); right_node=Node();
    left_node.box=parent.box; left_node.box(axis,1)=cut;
    right_node.box=parent.box; right_node.box(axis,0)=cut;
    left_node.depth=1; right_node.depth=1;
    left_node.axis=-1; right_node.axis=-1;
  }
  left_node.cut=NA_REAL; right_node.cut=NA_REAL;
  left_node.xi=xi_left; right_node.xi=xi_right;
  left_node.area=box_area(left_node.box);
  right_node.area=box_area(right_node.box);
  std::vector<arma::uword> left_index,right_index;
  for(int i=0;i<n;i++){
    if(X(i,axis)<cut) left_index.push_back(i);
    else right_index.push_back(i);
  }
  left_node.idx=arma::uvec(left_index.size());
  right_node.idx=arma::uvec(right_index.size());
  for(size_t k=0;k<left_index.size();k++) left_node.idx[k]=left_index[k];
  for(size_t k=0;k<right_index.size();k++) right_node.idx[k]=right_index[k];

  double A_parent=soft
    ?mpps_exposure_node(parent,region,gate,gate_mode):parent.area;
  double A_left=soft
    ?mpps_exposure_node(left_node,region,gate,gate_mode):left_node.area;
  double A_right=soft
    ?mpps_exposure_node(right_node,region,gate,gate_mode):right_node.area;
  std::vector<int> affected(n);
  for(int i=0;i<n;i++) affected[i]=i;
  std::vector<double> cut_logp;
  if(soft){
    irj_soft_cut_logp(parent,affected,X,region,kappa,gate,gate_mode,0.0,
      axis,xi_parent,xi_parent,cuts,cut_logp);
  }else{
    irj_hard_cut_logp(parent,affected,X,axis,kappa,xi_parent,xi_parent,
                      cuts,cut_logp);
  }

  int m_parent=0,m_left=0,m_right=0;
  double log_q_labels=0.0,spatial_difference=0.0;
  for(int i=0;i<n;i++){
    double log_phi_parent=0.0,log_phi_left=0.0,log_phi_right=0.0;
    if(soft){
      log_phi_parent=mpps_log_phi_node(parent,X.row(i),region,gate,gate_mode);
      log_phi_left=mpps_log_phi_node(left_node,X.row(i),region,gate,gate_mode);
      log_phi_right=mpps_log_phi_node(right_node,X.row(i),region,gate,gate_mode);
    }
    double lp_parent=log_phi_parent+std::log(kappa+m_parent)-
      std::log(kappa/xi_parent+A_parent);
    double lp_left=log_phi_left+std::log(kappa+m_left)-
      std::log(kappa/xi_left+A_left);
    double lp_right=log_phi_right+std::log(kappa+m_right)-
      std::log(kappa/xi_right+A_right);
    int label;
    if(soft){
      // Compact gates may be exactly zero away from their transition. Choose
      // a deterministic allocation that remains in the proposal support.
      if(i%3==1&&std::isfinite(lp_left)) label=2;
      else if(i%3==2&&std::isfinite(lp_right)) label=3;
      else label=1;
    }else{
      label=(i%3==0)?1:(X(i,axis)<cut?2:3);
      // A hard child on the opposite side has zero support.
      if(X(i,axis)<cut) lp_right=-arma::datum::inf;
      else lp_left=-arma::datum::inf;
    }
    double denominator=irj_lse3(lp_parent,lp_left,lp_right);
    if(label==1){
      log_q_labels+=lp_parent-denominator; m_parent++;
    }else if(label==2){
      log_q_labels+=lp_left-denominator; m_left++;
      spatial_difference+=log_phi_left-log_phi_parent;
    }else{
      log_q_labels+=lp_right-denominator; m_right++;
      spatial_difference+=log_phi_right-log_phi_parent;
    }
  }

  double left_stop=lSfac(left_node.box,left_node.idx,1,X,max_depth,
    min_leaf_n,cut_mode,ncand,alpha,eta);
  double right_stop=lSfac(right_node.box,right_node.idx,1,X,max_depth,
    min_leaf_n,cut_mode,ncand,alpha,eta);
  double tree_difference=lrho(0,alpha,eta)-l1mrho(0,alpha,eta)+
    left_stop+right_stop;
  double collapsed_difference=
    lg(m_parent,A_parent,kappa,xi_parent)+
    lg(m_left,A_left,kappa,xi_left)+
    lg(m_right,A_right,kappa,xi_right)-
    lg(n,A_parent,kappa,xi_parent);
  double kernel_grow=tree_difference+collapsed_difference+
    spatial_difference-log_q_labels-
    std::log((double)cuts.size())-cut_logp[chosen];
  double kernel_prune=-tree_difference-collapsed_difference-
    spatial_difference+log_q_labels+
    std::log((double)cuts.size())+cut_logp[chosen];

  double xi_mean=a_xi/b_xi;
  double log_prior_parent=lgam_dens(xi_parent,a_xi,xi_mean);
  double log_prior_children=lgam_dens(xi_left,a_xi,xi_mean)+
                            lgam_dens(xi_right,a_xi,xi_mean);
  double log_target_root=l1mrho(0,alpha,eta)+log_prior_parent+
    lg(n,A_parent,kappa,xi_parent);
  double log_target_split=lrho(0,alpha,eta)-std::log((double)d)-
    std::log((double)cuts.size())+left_stop+right_stop+
    log_prior_parent+log_prior_children+
    lg(m_parent,A_parent,kappa,xi_parent)+
    lg(m_left,A_left,kappa,xi_left)+
    lg(m_right,A_right,kappa,xi_right)+spatial_difference;
  double log_q_grow=std::log(0.5)-std::log((double)d)+
    cut_logp[chosen]+log_prior_children+log_q_labels;
  double log_q_prune=std::log(0.5);
  double direct=log_target_split+log_q_prune-log_target_root-log_q_grow;

  return List::create(
    _["kernel_grow_log_ratio"]=kernel_grow,
    _["kernel_prune_log_ratio"]=kernel_prune,
    _["direct_log_ratio"]=direct,
    _["grow_direct_error"]=kernel_grow-direct,
    _["grow_prune_cycle_error"]=kernel_grow+kernel_prune,
    _["cut_probability"]=std::exp(cut_logp[chosen]),
    _["label_probability"]=std::exp(log_q_labels),
    _["axis"]=axis+1,_["cut"]=cut,
    _["n_candidates"]=(int)cuts.size()
  );
}
