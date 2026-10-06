#ifndef POISTREE_PPT_SOFT_INFORMED_H
#define POISTREE_PPT_SOFT_INFORMED_H

// Exact locally balanced MH on the collapsed augmented state (tree, labels).
// The ordinary grow/prune and change proposal laws include Bernoulli label
// proposals. Their spatial factors cancel from the MH ratio, leaving a ratio
// that depends on the proposed left count, not on which observations go left.
// Sum these balanced masses by a Poisson-binomial recursion, then sample a
// count and a conditional Bernoulli subset. No allocation space is enumerated
// by the fitting algorithm.
#include <memory>
#include <sstream>
#include <iomanip>
#include <string>

static inline double ppsti_ladd(double a,double b){
  if(a==-std::numeric_limits<double>::infinity()) return b;
  if(b==-std::numeric_limits<double>::infinity()) return a;
  double hi=std::max(a,b),lo=std::min(a,b);
  return hi+std::log1p(std::exp(lo-hi));
}
static int ppsti_sample_log(const std::vector<double>&lp){
  double den=ppst_logsumexp(lp);
  if(!std::isfinite(den)) return -1;
  double u=R::unif_rand(),cum=0.0;
  int last=-1;
  for(size_t k=0;k<lp.size();k++) if(std::isfinite(lp[k])){
    last=(int)k; cum+=std::exp(lp[k]-den);
    if(u<cum) return (int)k;
  }
  return last;
}

struct PPSTICuts { std::vector<std::vector<double> > axis; bool any=false; };
struct PPSTIRouting { std::vector<double> left,right; };

// These bounded caches contain only geometry and gate-dependent basis terms.
// Cut support survives label/gate changes; exposures/routing are cleared after
// a gate sweep. Neighborhood masses are separately invalidated after every
// label/gate sweep, and reused after rejected tree proposals.
struct PPSTIContext {
  const arma::mat &pts,&region;
  const arma::vec &gate;
  double a,b,alpha,eta;
  int Dmax,nmin,mode,ncand,family;
  std::unordered_map<std::string,std::shared_ptr<PPSTICuts> > cuts_cache;
  std::unordered_map<std::string,double> exposure_cache;
  std::unordered_map<std::string,std::shared_ptr<PPSTIRouting> > routing_cache;
  PPSTIContext(const arma::mat&x,const arma::mat&r,const arma::vec&g,
      double aa,double bb,double al,double et,int dep,int nm,int md,int nc,int f)
    :pts(x),region(r),gate(g),a(aa),b(bb),alpha(al),eta(et),Dmax(dep),
     nmin(nm),mode(md),ncand(nc),family(f){}
  std::string key(const PPSTNode&nd) const {
    std::ostringstream s; s<<std::setprecision(17);
    for(const PPSTGate&g:nd.path)
      s<<g.axis<<':'<<g.cut<<':'<<g.parent_width<<':'<<g.side<<';';
    return s.str();
  }
  std::shared_ptr<PPSTICuts> cuts(const PPSTNode&nd){
    std::string k=key(nd);
    auto it=cuts_cache.find(k); if(it!=cuts_cache.end()) return it->second;
    std::shared_ptr<PPSTICuts> out(new PPSTICuts);
    out->axis.resize(pts.n_cols);
    if(nd.depth<Dmax && (int)nd.idx.n_elem>=2*std::max(1,nmin)){
      for(arma::uword j=0;j<pts.n_cols;j++){
        out->axis[j]=ppst_axis_cuts(pts,nd.idx,nd.box,j,mode,ncand,nmin);
        if(!out->axis[j].empty()) out->any=true;
      }
    }
    if(cuts_cache.size()>=512) cuts_cache.clear();
    cuts_cache[k]=out; return out;
  }
  double stop(const PPSTNode&nd){
    return cuts(nd)->any ? std::log1p(-ppst_rho(nd.depth,alpha,eta)) : 0.0;
  }
  double exposure(const PPSTNode&nd){
    std::string k=key(nd);
    auto it=exposure_cache.find(k); if(it!=exposure_cache.end()) return it->second;
    double H=ppst_exposure(nd,region,gate,family);
    if(exposure_cache.size()>=4096) exposure_cache.clear();
    exposure_cache[k]=H; return H;
  }
  std::shared_ptr<PPSTIRouting> routing(const PPSTNode&parent,int axis,double cut){
    double width=parent.box(axis,1)-parent.box(axis,0);
    std::ostringstream ss; ss<<std::setprecision(17)<<axis<<':'<<cut<<':'<<width;
    std::string k=ss.str();
    auto it=routing_cache.find(k); if(it!=routing_cache.end()) return it->second;
    std::shared_ptr<PPSTIRouting> out(new PPSTIRouting);
    out->left.resize(pts.n_rows); out->right.resize(pts.n_rows);
    PPSTGate split={axis,cut,width,1};
    for(arma::uword i=0;i<pts.n_rows;i++){
      if(family!=1){
        const double scale=family==2 ? width : region(axis,1)-region(axis,0);
        double z=gate[axis]*(pts(i,axis)-cut)/scale;
        out->left[i]=pst_logistic_log_right(-z);
        out->right[i]=pst_logistic_log_right(z);
      }else{
        double p=ppst_compact_gate_value(split,pts(i,axis),gate[axis]);
        out->left[i]=std::log1p(-p); out->right[i]=std::log(p);
      }
    }
    if(routing_cache.size()>=128) routing_cache.clear();
    routing_cache[k]=out; return out;
  }
  void clear_gate_cache(){ exposure_cache.clear(); routing_cache.clear(); }
};

struct PPSTIAction {
  int move=0,id=1,axis=-1;
  double cut=NA_REAL,logq=0.0,logmass=0.0;
  PPSTNode left,right;
  std::vector<int> affected;
  std::vector<double> logleft,logright,logcount,logratio;
};
struct PPSTINeighborhood {
  std::vector<PPSTIAction> actions;
  double logZ=-std::numeric_limits<double>::infinity();
};

// Distribution of the left-label count on a specified range, in log space.
// O(m^2) arithmetic and O(m) memory, including deterministic compact gates.
static std::vector<double> ppsti_counts(const PPSTIAction&act,int begin,int end){
  const double neg=-std::numeric_limits<double>::infinity();
  std::vector<double> out(end-begin+1,neg); out[0]=0.0;
  int used=0;
  for(int z=begin;z<end;z++){
    double ll=act.logleft[z],lr=act.logright[z];
    ++used;
    if((used&255)==0) Rcpp::checkUserInterrupt();
    for(int k=used;k>=0;k--){
      double stay=k<used ? out[k]+lr : neg;
      double add=k>0 ? out[k-1]+ll : neg;
      out[k]=ppsti_ladd(stay,add);
    }
  }
  return out;
}

static bool ppsti_is_cherry(const PPSTree&T,int id){
  auto p=T.find(id),l=T.find(2*id),r=T.find(2*id+1);
  return p!=T.end()&&p->second.axis>=0&&l!=T.end()&&r!=T.end()&&
    l->second.axis<0&&r->second.axis<0;
}

static PPSTINeighborhood ppsti_neighborhood(const PPSTree&T,
    const std::vector<int>&labels,PPSTIContext&ctx,int kind){
  PPSTINeighborhood out;
  std::vector<int> G,P;
  for(const auto&kv:T){
    if(kv.second.axis<0 && ctx.cuts(kv.second)->any) G.push_back(kv.first);
    if(ppsti_is_cherry(T,kv.first)) P.push_back(kv.first);
  }
  // A stable action order also makes diagnostic comparisons independent of
  // unordered_map insertion order and cache hits.
  std::sort(G.begin(),G.end()); std::sort(P.begin(),P.end());
  std::vector<int> nodes=kind==0 ? G : P;
  for(int id:nodes){
    Rcpp::checkUserInterrupt();
    const PPSTNode&parent=T.at(id);
    auto cuts=ctx.cuts(parent);
    double oldg;
    double oldstop=0.0;
    if(kind==0) oldg=ppst_log_g(parent.m,ctx.exposure(parent),ctx.a,ctx.b);
    else{
      const PPSTNode&l=T.at(2*id); const PPSTNode&r=T.at(2*id+1);
      oldg=ppst_log_g(l.m,ctx.exposure(l),ctx.a,ctx.b)+
           ppst_log_g(r.m,ctx.exposure(r),ctx.a,ctx.b);
      oldstop=ctx.stop(l)+ctx.stop(r);
    }
    std::vector<int> affected;
    for(size_t i=0;i<labels.size();i++)
      if(kind==0 ? labels[i]==id : (labels[i]==2*id||labels[i]==2*id+1))
        affected.push_back((int)i);
    const int m=(int)affected.size();
    for(arma::uword j=0;j<ctx.pts.n_cols;j++){
      for(double cut:cuts->axis[j]){
        Rcpp::checkUserInterrupt();
        PPSTIAction act; act.move=kind==0?0:2; act.id=id; act.axis=j;
        act.cut=cut; act.affected=affected;
        auto route=ctx.routing(parent,j,cut);
        act.logleft.resize(m);act.logright.resize(m);
        for(int z=0;z<m;z++){
          act.logleft[z]=route->left[affected[z]];
          act.logright[z]=route->right[affected[z]];
        }
        arma::uvec li,ri; ppst_split_indices(parent,ctx.pts,j,cut,li,ri);
        act.left=ppst_child(parent,li,j,cut,-1);
        act.right=ppst_child(parent,ri,j,cut,1);
        double HL=ctx.exposure(act.left),HR=ctx.exposure(act.right);
        double change=ctx.stop(act.left)+ctx.stop(act.right)-oldstop;
        if(kind==0){
          int Pnew=(int)P.size()+1-(id>1 && ppsti_is_cherry(T,id/2)?1:0);
          double rho=ppst_rho(parent.depth,ctx.alpha,ctx.eta);
          change+=std::log(rho)-std::log1p(-rho)+
                  std::log((double)G.size())-std::log((double)Pnew);
          act.logq=std::log(0.5)-std::log((double)G.size());
        }else act.logq=-std::log((double)P.size());
        act.logq-=std::log((double)ctx.pts.n_cols)+std::log((double)cuts->axis[j].size());
        act.logcount=ppsti_counts(act,0,m);
        act.logratio.resize(m+1);
        std::vector<double> mass(m+1);
        for(int k=0;k<=m;k++){
          act.logratio[k]=ppst_log_g(k,HL,ctx.a,ctx.b)+
            ppst_log_g(m-k,HR,ctx.a,ctx.b)-oldg+change;
          mass[k]=act.logcount[k]+0.5*act.logratio[k];
        }
        act.logmass=act.logq+ppst_logsumexp(mass);
        out.actions.push_back(std::move(act));
      }
    }
  }
  if(kind==0) for(int id:P){
    const PPSTNode&parent=T.at(id),&l=T.at(2*id),&r=T.at(2*id+1);
    int Gnew=(int)G.size()-(ctx.cuts(l)->any?1:0)-(ctx.cuts(r)->any?1:0)+1;
    double rho=ppst_rho(parent.depth,ctx.alpha,ctx.eta);
    double ratio=ppst_log_g(l.m+r.m,ctx.exposure(parent),ctx.a,ctx.b)-
      ppst_log_g(l.m,ctx.exposure(l),ctx.a,ctx.b)-
      ppst_log_g(r.m,ctx.exposure(r),ctx.a,ctx.b)-
      std::log(rho)+std::log1p(-rho)-ctx.stop(l)-ctx.stop(r)+
      std::log((double)P.size())-std::log((double)Gnew);
    PPSTIAction act; act.move=1; act.id=id;
    act.logq=std::log(0.5)-std::log((double)P.size());
    act.logratio.push_back(ratio); act.logmass=act.logq+0.5*ratio;
    out.actions.push_back(std::move(act));
  }
  std::vector<double> masses; masses.reserve(out.actions.size());
  for(const auto&act:out.actions) masses.push_back(act.logmass);
  out.logZ=ppst_logsumexp(masses);
  return out;
}

// Draw Bernoulli allocations conditional on their total, using a binary
// divide-and-conquer recursion. Unlike a full forward/backward table this
// needs O(m) working memory; total arithmetic remains O(m^2).
static void ppsti_subset(const PPSTIAction&act,int begin,int end,int count,
    std::vector<int>&labels){
  int m=end-begin;
  if(count==0 || count==m){
    int dest=2*act.id+(count==0?1:0);
    for(int i=begin;i<end;i++) labels[act.affected[i]]=dest;
    return;
  }
  int mid=begin+m/2;
  auto left=ppsti_counts(act,begin,mid),right=ppsti_counts(act,mid,end);
  std::vector<double> lp(left.size(),-std::numeric_limits<double>::infinity());
  for(int k=0;k<(int)left.size();k++)
    if(count-k>=0 && count-k<(int)right.size()) lp[k]=left[k]+right[count-k];
  int nleft=ppsti_sample_log(lp);
  if(nleft<0) stop("Invalid conditional allocation distribution in informed MH");
  ppsti_subset(act,begin,mid,nleft,labels);
  ppsti_subset(act,mid,end,count-nleft,labels);
}

static void ppsti_apply(const PPSTIAction&act,PPSTree&T,std::vector<int>&labels,
    int nleft){
  int id=act.id,L=2*id,R=2*id+1;
  if(act.move==1){
    int m=T.at(L).m+T.at(R).m;
    T.erase(L);T.erase(R); T.at(id).axis=-1; T.at(id).cut=NA_REAL;T.at(id).m=m;
    for(int&z:labels) if(z==L||z==R) z=id;
  }else{
    T.at(id).axis=act.axis;T.at(id).cut=act.cut;T.at(id).m=0;
    T[L]=act.left;T[R]=act.right;
    T.at(L).m=nleft; T.at(R).m=(int)act.affected.size()-nleft;
  }
}

static int ppsti_step(PPSTree&T,std::vector<int>&labels,PPSTIContext&ctx,
    int kind,PPSTINeighborhood&current,bool&valid,int&which){
  which=-1;
  if(!valid){current=ppsti_neighborhood(T,labels,ctx,kind);valid=true;}
  if(!std::isfinite(current.logZ)) return 0;
  std::vector<double> mass;
  for(const auto&act:current.actions) mass.push_back(act.logmass);
  int pick=ppsti_sample_log(mass);
  if(pick<0) return 0;
  const PPSTIAction&act=current.actions[pick]; which=act.move;
  PPSTree proposal=T; std::vector<int> proposed_labels=labels;
  int k=0;
  if(act.move!=1){
    std::vector<double> counts(act.logratio.size());
    for(size_t j=0;j<counts.size();j++) counts[j]=act.logcount[j]+0.5*act.logratio[j];
    k=ppsti_sample_log(counts);
    if(k<0) stop("Invalid count distribution in informed MH");
    ppsti_subset(act,0,(int)act.affected.size(),k,proposed_labels);
  }
  ppsti_apply(act,proposal,proposed_labels,k);
  PPSTINeighborhood next=ppsti_neighborhood(proposal,proposed_labels,ctx,kind);
  if(!std::isfinite(next.logZ)) stop("Missing reverse neighborhood in informed MH");
  if(std::log(R::unif_rand())<current.logZ-next.logZ){
    T.swap(proposal); labels.swap(proposed_labels); current=std::move(next);
    return 1;
  }
  return 0;
}

#endif
