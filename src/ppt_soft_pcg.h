#ifndef POISTREE_PPT_SOFT_PCG_H
#define POISTREE_PPT_SOFT_PCG_H

struct PPSTPCGStats {
  arma::mat factor;
  arma::vec direction;
  double last_alpha=0.0;
  bool pending=false;
  int updates=0,failures=0,proposals=0,accepted=0;
  double alpha_sum=0.0;
  void initialize(const arma::vec&sd,bool shared){
    const arma::uword p=shared?1:sd.n_elem;
    factor=arma::diagmat(sd.head(p));
    updates=failures=proposals=accepted=0;alpha_sum=0.0;pending=false;
  }
  double acceptance() const {return accepted/std::max(1.0,(double)proposals);}
  double mean_alpha() const {return alpha_sum/std::max(1.0,(double)proposals);}
};

static void ppst_pcg_controls(bool pcg,bool informed,double target,
    double decay,int adapt,int burn,double tau,double eps){
  if(!pcg) return;
  if(informed&&(!std::isfinite(tau)||tau<=0.0||tau>1.0))
    stop("proposal_temperature must be in (0, 1]");
  if(informed&&(!std::isfinite(eps)||eps<0.0||eps>=1.0))
    stop("proposal_defensive must be in [0, 1)");
  if(!std::isfinite(target)||target<=0.0||target>=1.0)
    stop("ram_target must be strictly between zero and one");
  if(!std::isfinite(decay)||decay<=0.5||decay>1.0)
    stop("ram_decay must be in (0.5, 1]");
  if(adapt<0||adapt>burn) stop("ram_adapt must be between zero and burn");
}

// Apply Vihola's robust adaptive Metropolis update to the Cholesky factor.

static bool ppst_ram_update(arma::mat&factor,const arma::vec&u,
    double alpha,double target,double step){
  const arma::uword p=factor.n_rows;
  if(p==0||factor.n_cols!=p||u.n_elem!=p||!factor.is_finite()||
     !u.is_finite()||!std::isfinite(alpha)||alpha<0.0||alpha>1.0||
     !std::isfinite(target)||target<=0.0||target>=1.0||
     !std::isfinite(step)||step<=0.0||step>1.0) return false;
  const double norm=arma::norm(u,2);
  if(!(norm>0.0)||!std::isfinite(norm)) return false;
  const double beta=step*(alpha-target);
  arma::mat next=factor;
  arma::vec x=std::sqrt(std::abs(beta))*(factor*(u/norm));
  if(!x.is_finite()) return false;
  const double sign=beta<0.0?-1.0:1.0;
  for(arma::uword k=0;k<p;k++){
    const double diag=next(k,k),ratio=x[k]/diag;
    if(!(diag>0.0)||!std::isfinite(ratio)) return false;
    if(sign<0.0 && !(std::abs(ratio)<1.0)) return false;
    const double r=sign>0.0 ? std::hypot(diag,x[k])
      : diag*std::sqrt((1.0-ratio)*(1.0+ratio));
    const double c=r/diag;
    if(!(r>0.0)||!std::isfinite(r)||!std::isfinite(c)) return false;
    next(k,k)=r;
    for(arma::uword j=k+1;j<p;j++){
      next(j,k)=(next(j,k)+sign*ratio*x[j])/c;
      x[j]=c*x[j]-ratio*next(j,k);
    }
  }
  if(!next.is_finite()) return false;
  factor=std::move(next);
  return true;
}

// Draw a Gamma rate and return its logarithm.
static double ppst_pcg_log_rate(double shape,double rate){
  double value;
  if(shape<1.0)
    value=std::log(R::rgamma(shape+1.0,1.0))+
      std::log(R::runif(0.0,1.0))/shape-std::log(rate);
  else value=std::log(R::rgamma(shape,1.0))-std::log(rate);
  if(!std::isfinite(value)) stop("non-finite temporary log intensity in pcg");
  return value;
}

struct PPSTPCGEvaluation {
  double target=-std::numeric_limits<double>::infinity();
  arma::mat log_weights;
  arma::vec log_normalizers;
};

// Evaluate the allocation-collapsed target with the selected gate family.

static PPSTPCGEvaluation ppst_pcg_evaluate(const PPSTree&T,
    const std::vector<int>&leaves,const arma::vec&log_rate,
    const arma::mat&pts,const arma::mat&region,const arma::vec&gate,
    const arma::vec&a_gate,const arma::vec&b_gate,const arma::vec&gate_min,
    bool shared,int family,PPSTGeometryCache*cache=nullptr){
  PPSTPCGEvaluation out;
  if(!gate.is_finite()||arma::any(gate<=gate_min)||arma::any(gate<=0.0))
    return out;
  double target=0.0;
  const arma::uword p=shared?1:gate.n_elem;
  for(arma::uword j=0;j<p;j++)
    target+=(a_gate[j]-1.0)*std::log(gate[j])-b_gate[j]*gate[j];
  out.log_weights.set_size(pts.n_rows,leaves.size());
  for(size_t k=0;k<leaves.size();k++){
    const PPSTNode&nd=T.at(leaves[k]);
    const double H=ppst_cached_exposure(nd,region,gate,family,cache);
    if(!std::isfinite(H)||H<0.0) return out;
    if(H>0.0) target-=std::exp(log_rate[k]+std::log(H));
    if(cache){
      const arma::vec&basis=cache->training(nd);
      for(arma::uword i=0;i<pts.n_rows;i++) out.log_weights(i,k)=log_rate[k]+basis[i];
    }else{
      for(arma::uword i=0;i<pts.n_rows;i++)
        out.log_weights(i,k)=log_rate[k]+ppst_log_phi(nd,pts.row(i),region,gate,family);
    }
  }
  out.log_normalizers.set_size(pts.n_rows);
  std::vector<double> lw(leaves.size());
  for(arma::uword i=0;i<pts.n_rows;i++){
    for(size_t k=0;k<leaves.size();k++) lw[k]=out.log_weights(i,k);
    const double den=ppst_logsumexp(lw);
    if(!std::isfinite(den)) return out;
    out.log_normalizers[i]=den;
    target+=den;
  }
  if(std::isfinite(target)) out.target=target;
  return out;
}

static void ppst_pcg_restore(PPSTree&T,std::vector<int>&labels,
    const std::vector<int>&leaves,const PPSTPCGEvaluation&value){
  for(int id:leaves) T.at(id).m=0;
  std::vector<double> lw(leaves.size());
  for(size_t i=0;i<labels.size();i++){
    for(size_t k=0;k<leaves.size();k++) lw[k]=value.log_weights(i,k);
    const double den=value.log_normalizers[i],u=R::runif(0.0,1.0);
    if(!std::isfinite(den)) stop("invalid allocation probabilities in pcg");
    double cumulative=0.0;int pick=-1;
    for(size_t k=0;k<leaves.size();k++){
      if(!std::isfinite(lw[k])) continue;
      pick=leaves[k];
      cumulative+=std::exp(lw[k]-den);
      if(u<cumulative) break;
    }
    if(pick<0) stop("empty allocation support in pcg");
    labels[i]=pick;T.at(pick).m++;
  }
}

static void ppst_pcg_block(PPSTree&T,std::vector<int>&labels,
    const arma::mat&pts,const arma::mat&region,double a,double b,
    arma::vec&gate,const arma::vec&a_gate,const arma::vec&b_gate,
    const arma::vec&gate_min,bool shared,int family,bool update_gate,
    PPSTPCGStats&stats,PPSTGeometryCache*cache=nullptr){
  stats.pending=false;
  std::vector<int> leaves;ppst_leaves(T,leaves);
  arma::vec log_rate(leaves.size());
  for(size_t k=0;k<leaves.size();k++){
    const PPSTNode&nd=T.at(leaves[k]);
    log_rate[k]=ppst_pcg_log_rate(a+nd.m,b+ppst_cached_exposure(nd,region,gate,family,cache));
  }
  PPSTPCGEvaluation current=ppst_pcg_evaluate(T,leaves,log_rate,pts,region,
    gate,a_gate,b_gate,gate_min,shared,family,cache);
  if(!std::isfinite(current.target)) stop("non-finite current gate target in pcg");
  if(update_gate){
    const arma::uword p=stats.factor.n_rows;
    arma::vec u(p);
    for(arma::uword j=0;j<p;j++) u[j]=R::rnorm(0.0,1.0);
    const arma::vec theta=arma::log(gate.head(p));
    const arma::vec proposed_theta=theta+stats.factor*u;
    arma::vec proposal=gate;
    if(shared) proposal.fill(std::exp(proposed_theta[0]));
    else proposal=arma::exp(proposed_theta);
    std::unique_ptr<PPSTGeometryCache> proposed_cache;
    if(cache) proposed_cache.reset(new PPSTGeometryCache(*cache,proposal));
    PPSTPCGEvaluation proposed=ppst_pcg_evaluate(T,leaves,log_rate,pts,region,
      proposal,a_gate,b_gate,gate_min,shared,family,proposed_cache.get());
    const double log_ratio=proposed.target-current.target+
      arma::sum(proposed_theta-theta);
    double alpha=0.0;
    if(std::isfinite(proposed.target)&&!std::isnan(log_ratio))
      alpha=log_ratio>=0.0?1.0:std::exp(log_ratio);
    stats.proposals++;stats.alpha_sum+=alpha;
    if(R::runif(0.0,1.0)<alpha){
      gate=proposal;current=std::move(proposed);stats.accepted++;
      if(cache) *cache=std::move(*proposed_cache);
    }

    stats.direction=std::move(u);stats.last_alpha=alpha;stats.pending=true;
  }

  ppst_pcg_restore(T,labels,leaves,current);
}

static void ppst_pcg_adapt(PPSTPCGStats&stats,int iteration,
    double ram_target,double ram_decay,int ram_adapt){
  if(stats.pending&&iteration<ram_adapt){
    const double step=std::min(1.0,(double)stats.factor.n_rows*
      std::pow(iteration+1.0,-ram_decay));
    if(ppst_ram_update(stats.factor,stats.direction,stats.last_alpha,ram_target,step))
      stats.updates++;
    else stats.failures++;
  }
  stats.pending=false;
}

#endif
