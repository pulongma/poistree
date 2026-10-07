#ifndef POISTREE_PPT_SOFT_SURROGATE_H
#define POISTREE_PPT_SOFT_SURROGATE_H

struct PPSTSurrogateEntry {
  std::vector<int> axis,per_axis;
  std::vector<double> cut,score;
  double logS=-std::numeric_limits<double>::infinity();
};

class PPSTSurrogate {
  const arma::mat &pts_;
  PPSTCutsCache&cuts_;
  const double a_,b_,alpha_,eta_,tau_,eps_;
  std::unordered_map<std::string,PPSTSurrogateEntry> entries_;

  bool quadrature() const {
#ifdef POISTREE_SPATIAL_QUADRATURE_H
    return qpp_active;
#else
    return false;
#endif
  }
  void box_rows(const arma::mat&box,std::vector<int>&rows,double&total) const {
    rows.clear(); total=0.0;
#ifdef POISTREE_SPATIAL_QUADRATURE_H
    if(quadrature()){
      for(arma::uword i=0;i<qpp_weights.n_elem;i++){
        bool inside=true;
        for(arma::uword j=0;j<box.n_rows&&inside;j++){
          const double x=qpp_background(i,j);
          const bool last=box(j,1)==qpp_region(j,1);
          inside=x>=box(j,0)&&(last?x<=box(j,1):x<box(j,1));
        }
        if(inside){rows.push_back((int)i);total+=qpp_weights[i];}
      }
      return;
    }
#endif
    total=arma::prod(box.col(1)-box.col(0));
  }
public:
  PPSTSurrogate(const arma::mat&pts,PPSTCutsCache&cuts,
      double a,double b,double alpha,double eta,double tau,double eps)
    :pts_(pts),cuts_(cuts),a_(a),b_(b),alpha_(alpha),eta_(eta),
     tau_(tau),eps_(eps){}
  double tau() const {return tau_;}
  double eps() const {return eps_;}

  PPSTSurrogateEntry& entry(const PPSTNode&nd){
    const std::string k=PPSTCutsCache::key(nd);
    auto found=entries_.find(k);
    if(found!=entries_.end()) return found->second;
    PPSTSurrogateEntry e;
    const arma::uword d=pts_.n_cols;
    e.per_axis.assign(d,0);
    const double rho=ppst_rho(nd.depth,alpha_,eta_);
    const double base=std::log(rho)-std::log1p(-rho)-std::log((double)d);
    std::vector<int> rows; double area;
    box_rows(nd.box,rows,area);
    const double parent=ppst_log_g((int)nd.idx.n_elem,area,a_,b_);
    for(arma::uword j=0;j<d;j++){
      const std::vector<double>&cuts=cuts_.axis(nd,j);
      if(cuts.empty()) continue;
      e.per_axis[j]=(int)cuts.size();
      const arma::vec col=pts_.col(j),x=arma::sort(col.elem(nd.idx));
      std::vector<double> xs,prefix;
#ifdef POISTREE_SPATIAL_QUADRATURE_H
      if(quadrature()){
        std::sort(rows.begin(),rows.end(),[&](int p,int q){
          return qpp_background(p,j)<qpp_background(q,j);});
        xs.resize(rows.size());prefix.resize(rows.size()+1);prefix[0]=0.0;
        for(size_t r=0;r<rows.size();r++){
          xs[r]=qpp_background(rows[r],j);prefix[r+1]=prefix[r]+qpp_weights[rows[r]];
        }
      }
#endif
      const double lo=nd.box(j,0),hi=nd.box(j,1);
      for(double c:cuts){
        const int nL=(int)(std::lower_bound(x.begin(),x.end(),c)-x.begin());
        const int nR=(int)x.n_elem-nL;
        double AL;
        if(quadrature()) AL=prefix[std::lower_bound(xs.begin(),xs.end(),c)-xs.begin()];
        else AL=area*(c-lo)/(hi-lo);
        const double AR=std::max(0.0,area-AL);
        e.axis.push_back((int)j);e.cut.push_back(c);
        e.score.push_back(ppst_log_g(nL,AL,a_,b_)+ppst_log_g(nR,AR,a_,b_)-parent
                          +base-std::log((double)cuts.size()));
      }
    }
    std::vector<double> t(e.score.size());
    for(size_t k=0;k<t.size();k++) t[k]=tau_*e.score[k];
    e.logS=ppst_logsumexp(t);
    return entries_[k]=std::move(e);
  }
  // Return the index of the current split in its candidate table.
  int current(const PPSTNode&nd){
    const PPSTSurrogateEntry&e=entry(nd);
    for(size_t k=0;k<e.cut.size();k++)
      if(e.axis[k]==nd.axis&&e.cut[k]==nd.cut) return (int)k;
    stop("surrogate: current split is not a candidate cut");
    return -1;
  }
  void trim(const PPSTree&tree){
    std::unordered_set<std::string> active;
    for(const auto&kv:tree) active.insert(PPSTCutsCache::key(kv.second));
    for(auto it=entries_.begin();it!=entries_.end();){
      if(active.find(it->first)==active.end()) it=entries_.erase(it);
      else ++it;
    }
  }
};

// Return a log proposal probability for the informed-uniform mixture.

static inline double ppsts_log_q(double informed,double logW,long n,double eps){
  return ppsti_ladd(std::log1p(-eps)+informed-logW,std::log(eps)-std::log((double)n));
}

struct PPSTSGrowSet {
  std::vector<int> nodes; std::vector<const PPSTSurrogateEntry*> entry;
  double logW; long n;
};

static PPSTSGrowSet ppsts_grow_set(const PPSTree&T,const std::vector<int>&G,
    PPSTSurrogate&sur){
  PPSTSGrowSet out; out.nodes=G; out.n=0; out.entry.resize(G.size());
  std::vector<double> logS(G.size());
  for(size_t v=0;v<G.size();v++){
    out.entry[v]=&sur.entry(T.at(G[v]));
    logS[v]=out.entry[v]->logS; out.n+=(long)out.entry[v]->cut.size();
  }
  out.logW=ppst_logsumexp(logS);
  return out;
}

struct PPSTSPruneSet {
  std::vector<int> nodes,current; std::vector<const PPSTSurrogateEntry*> entry;
  double logW;
};

static PPSTSPruneSet ppsts_prune_set(const PPSTree&T,const std::vector<int>&P,
    PPSTSurrogate&sur){
  PPSTSPruneSet out; out.nodes=P; out.current.resize(P.size()); out.entry.resize(P.size());
  std::vector<double> w(P.size());
  for(size_t v=0;v<P.size();v++){
    const PPSTNode&nd=T.at(P[v]);
    out.current[v]=sur.current(nd);
    out.entry[v]=&sur.entry(nd);
    w[v]=-sur.tau()*out.entry[v]->score[out.current[v]];
  }
  out.logW=ppst_logsumexp(w);
  return out;
}

static double ppsts_log_grow_q(const PPSTSurrogateEntry&e,int k,const PPSTSGrowSet&set,
    double tau,double eps){
  return ppsts_log_q(tau*e.score[k],set.logW,set.n,eps);
}
static double ppsts_log_prune_q(const PPSTSurrogateEntry&e,int k,const PPSTSPruneSet&set,
    double tau,double eps){
  return ppsts_log_q(-tau*e.score[k],set.logW,(long)set.nodes.size(),eps);
}
static double ppsts_log_change_q(const PPSTSurrogateEntry&e,int k,double tau,double eps){
  return ppsts_log_q(tau*e.score[k],e.logS,(long)e.cut.size(),eps);
}

// Draw a candidate index from a node's proposal mixture.
static int ppsts_draw_in_node(const PPSTSurrogateEntry&e,double tau,double eps){
  if(R::unif_rand()<eps) return ppst_runif_int((int)e.cut.size());
  std::vector<double> t(e.score.size());
  for(size_t k=0;k<t.size();k++) t[k]=tau*e.score[k];
  const int k=ppsti_sample_log(t);
  if(k<0) stop("surrogate: invalid candidate distribution");
  return k;
}

static int ppsts_grow_prune(PPSTree&T,std::vector<int>&labels,
    const arma::mat&pts,const arma::mat&region,double a,double b,
    const arma::vec&gate,int gate_family,double alpha,double eta,int Dmax,
    int nmin,int mode,int ncand,int&which_move,PPSTGeometryCache*cache,
    PPSTCutsCache*cuts_cache,PPSTSurrogate&sur){
  const double tau=sur.tau(),eps=sur.eps(),d=pts.n_cols;
  std::vector<int>G,P;
  ppst_growable(T,pts,Dmax,nmin,mode,ncand,G,cuts_cache);
  ppst_prunable(T,P);
  if(R::unif_rand()<0.5){
    which_move=0;
    if(G.empty()) return 0;
    const PPSTSGrowSet set=ppsts_grow_set(T,G,sur);
    std::vector<double> w(G.size());
    const bool uniform=R::unif_rand()<eps;
    for(size_t u=0;u<G.size();u++)
      w[u]=uniform?std::log((double)set.entry[u]->cut.size()):set.entry[u]->logS;
    const int pick=ppsti_sample_log(w),v=G[pick];
    const PPSTSurrogateEntry&e=*set.entry[pick];
    int k;
    if(uniform) k=ppst_runif_int((int)e.cut.size());
    else{
      std::vector<double> t(e.score.size());
      for(size_t z=0;z<t.size();z++) t[z]=tau*e.score[z];
      k=ppsti_sample_log(t);
    }
    const double log_qf=ppsts_log_grow_q(e,k,set,tau,eps);
    const int axis=e.axis[k],L=2*v,R=2*v+1;
    const double cut=e.cut[k];
    PPSTNode old=T.at(v);
    arma::uvec li,ri; ppst_split_indices(old,pts,axis,cut,li,ri);
    PPSTNode nl=ppst_child(old,li,axis,cut,-1);
    PPSTNode nr=ppst_child(old,ri,axis,cut,1);
    std::vector<std::pair<int,int> > reass;
    int mL=0,mR=0;
    const arma::vec*parent_log_phi=nullptr;
    for(size_t i=0;i<labels.size();i++) if(labels[i]==v){
      if(cache&&!parent_log_phi) parent_log_phi=&cache->training(old);
      int dest=ppst_draw_side(nl,nr,pts,i,region,gate,gate_family,parent_log_phi)<0?L:R;
      reass.push_back(std::make_pair((int)i,dest));
      if(dest==L) mL++; else mR++;
    }
    nl.m=mL; nr.m=mR;
    T[v].axis=axis; T[v].cut=cut; T[v].m=0; T[L]=nl; T[R]=nr;
    std::vector<int>Pnew; ppst_prunable(T,Pnew);
    const PPSTSPruneSet reverse=ppsts_prune_set(T,Pnew,sur);
    const double log_qr=ppsts_log_prune_q(e,k,reverse,tau,eps);
    const double Hp=ppst_cached_exposure(old,region,gate,gate_family,cache);
    const double HL=ppst_cached_exposure(nl,region,gate,gate_family,cache);
    const double HR=ppst_cached_exposure(nr,region,gate,gate_family,cache);
    const double dp=std::log(ppst_rho(old.depth,alpha,eta))
      -std::log(1.0-ppst_rho(old.depth,alpha,eta))
      +std::log(ppst_stop_factor(nl,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache))
      +std::log(ppst_stop_factor(nr,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache));
    const double logA=ppst_log_g(mL,HL,a,b)+ppst_log_g(mR,HR,a,b)
      -ppst_log_g(old.m,Hp,a,b)+dp-std::log(d)-std::log((double)e.per_axis[axis])
      -log_qf+log_qr;
    if(std::log(R::unif_rand())<logA){
      for(const auto&z:reass) labels[z.first]=z.second;
      return 1;
    }
    T.erase(L); T.erase(R); T[v]=old;
    return 0;
  }

  which_move=1;
  if(P.empty()) return 0;
  const PPSTSPruneSet set=ppsts_prune_set(T,P,sur);
  int u;
  if(R::unif_rand()<eps) u=ppst_runif_int((int)P.size());
  else{
    std::vector<double> w(P.size());
    for(size_t z=0;z<P.size();z++) w[z]=-tau*set.entry[z]->score[set.current[z]];
    u=ppsti_sample_log(w);
  }
  const int v=P[u],k=set.current[u],L=2*v,R=2*v+1;
  const PPSTSurrogateEntry&e=*set.entry[u];
  const double log_qf=ppsts_log_prune_q(e,k,set,tau,eps);
  PPSTNode oldp=T.at(v),oldL=T.at(L),oldR=T.at(R);
  const int M=oldL.m+oldR.m,axis=oldp.axis;
  const double Hp=ppst_cached_exposure(oldp,region,gate,gate_family,cache);
  const double HL=ppst_cached_exposure(oldL,region,gate,gate_family,cache);
  const double HR=ppst_cached_exposure(oldR,region,gate,gate_family,cache);
  const double dp=std::log(ppst_rho(oldp.depth,alpha,eta))
    -std::log(1.0-ppst_rho(oldp.depth,alpha,eta))
    +std::log(ppst_stop_factor(oldL,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache))
    +std::log(ppst_stop_factor(oldR,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache));
  T.erase(L); T.erase(R); T[v].axis=-1; T[v].cut=NA_REAL; T[v].m=M;
  std::vector<int>Gnew; ppst_growable(T,pts,Dmax,nmin,mode,ncand,Gnew,cuts_cache);
  if(Gnew.empty()){T[v]=oldp;T[L]=oldL;T[R]=oldR;return 0;}
  const PPSTSGrowSet reverse=ppsts_grow_set(T,Gnew,sur);
  const double log_qr=ppsts_log_grow_q(e,k,reverse,tau,eps);
  const int K=e.per_axis[axis];
  const double logA=ppst_log_g(M,Hp,a,b)-ppst_log_g(oldL.m,HL,a,b)
    -ppst_log_g(oldR.m,HR,a,b)-dp+std::log(d)+std::log((double)K)
    -log_qf+log_qr;
  if(std::log(R::unif_rand())<logA){
    for(size_t i=0;i<labels.size();i++)
      if(labels[i]==L||labels[i]==R) labels[i]=v;
    return 1;
  }
  T[v]=oldp; T[L]=oldL; T[R]=oldR;
  return 0;
}

static int ppsts_change_cut(PPSTree&T,std::vector<int>&labels,
    const arma::mat&pts,const arma::mat&region,double a,double b,
    const arma::vec&gate,int gate_family,double alpha,double eta,int Dmax,
    int nmin,int mode,int ncand,PPSTGeometryCache*cache,PPSTCutsCache*cuts_cache,
    PPSTSurrogate&sur){
  std::vector<int>P; ppst_prunable(T,P); if(P.empty()) return 0;
  const int v=P[ppst_runif_int(P.size())],L=2*v,R=2*v+1;
  PPSTNode parent=T.at(v),oldL=T.at(L),oldR=T.at(R);
  const PPSTSurrogateEntry&e=sur.entry(parent);
  const int k_old=sur.current(parent),k_new=ppsts_draw_in_node(e,sur.tau(),sur.eps());
  const int axis=e.axis[k_new];
  const double cut=e.cut[k_new];
  arma::uvec li,ri; ppst_split_indices(parent,pts,axis,cut,li,ri);
  PPSTNode nl=ppst_child(parent,li,axis,cut,-1);
  PPSTNode nr=ppst_child(parent,ri,axis,cut,1);
  std::vector<std::pair<int,int> >reass;
  int mL=0,mR=0;
  const arma::vec*parent_log_phi=nullptr;
  for(size_t i=0;i<labels.size();i++) if(labels[i]==L||labels[i]==R){
    if(cache&&!parent_log_phi) parent_log_phi=&cache->training(parent);
    int dest=ppst_draw_side(nl,nr,pts,i,region,gate,gate_family,parent_log_phi)<0?L:R;
    reass.push_back(std::make_pair((int)i,dest));
    if(dest==L)mL++;else mR++;
  }
  nl.m=mL; nr.m=mR;
  const double HoL=ppst_cached_exposure(oldL,region,gate,gate_family,cache);
  const double HoR=ppst_cached_exposure(oldR,region,gate,gate_family,cache);
  const double HnL=ppst_cached_exposure(nl,region,gate,gate_family,cache);
  const double HnR=ppst_cached_exposure(nr,region,gate,gate_family,cache);
  const double oldg=ppst_log_g(oldL.m,HoL,a,b)+ppst_log_g(oldR.m,HoR,a,b);
  const double newg=ppst_log_g(mL,HnL,a,b)+ppst_log_g(mR,HnR,a,b);
  const double oldS=std::log(ppst_stop_factor(oldL,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache))
                   +std::log(ppst_stop_factor(oldR,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache));
  const double newS=std::log(ppst_stop_factor(nl,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache))
                   +std::log(ppst_stop_factor(nr,pts,Dmax,nmin,mode,ncand,alpha,eta,cuts_cache));
  const double logA=(newg-oldg)+(newS-oldS)
    +std::log((double)e.per_axis[parent.axis])-std::log((double)e.per_axis[axis])
    +ppsts_log_change_q(e,k_old,sur.tau(),sur.eps())
    -ppsts_log_change_q(e,k_new,sur.tau(),sur.eps());
  if(std::log(R::unif_rand())<logA){
    T[v].axis=axis;T[v].cut=cut;T[L]=nl;T[R]=nr;
    for(const auto&z:reass) labels[z.first]=z.second;
    return 1;
  }
  return 0;
}

#endif
