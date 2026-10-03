#ifndef POISTREE_IRJMCMC_MOVES_H
#define POISTREE_IRJMCMC_MOVES_H

// Cut-informed reversible-jump moves for MPPT and S-MPPT.  The proposal
// upgrades are deliberately separated from the surrounding Gibbs/MH update
// schedule.  The production iRJ kernel combines informed cuts with sequential
// conditional allocation proposals.  Every non-uniform proposal mass is
// included in the forward/reverse MH ratio.

static inline double irj_lse2(double x,double y){
  double m=std::max(x,y);
  if(!std::isfinite(m)) return m;
  return m+std::log(std::exp(x-m)+std::exp(y-m));
}

static inline double irj_lse3(double x,double y,double z){
  double m=std::max(x,std::max(y,z));
  if(!std::isfinite(m)) return m;
  return m+std::log(std::exp(x-m)+std::exp(y-m)+std::exp(z-m));
}

static void irj_normalize(std::vector<double>&logp){
  if(logp.empty()) return;
  double mx=-std::numeric_limits<double>::infinity();
  int n_pos_inf=0;
  for(double z:logp){
    if(z==std::numeric_limits<double>::infinity()) n_pos_inf++;
    else if(std::isfinite(z)) mx=std::max(mx,z);
  }
  if(n_pos_inf>0){
    double value=-std::log((double)n_pos_inf);
    for(double&z:logp)
      z=(z==std::numeric_limits<double>::infinity())
        ?value:-std::numeric_limits<double>::infinity();
    return;
  }
  if(!std::isfinite(mx)){
    double value=-std::log((double)logp.size());
    std::fill(logp.begin(),logp.end(),value);
    return;
  }
  double total=0.0;
  for(double&z:logp){
    if(std::isfinite(z)){ z-=mx; total+=std::exp(z); }
    else z=-std::numeric_limits<double>::infinity();
  }
  if(!std::isfinite(total)||total<=0.0){
    double value=-std::log((double)logp.size());
    std::fill(logp.begin(),logp.end(),value);
    return;
  }
  double log_total=std::log(total);
  for(double&z:logp) if(std::isfinite(z)) z-=log_total;
}

static int irj_pick(const std::vector<double>&logp){
  double u=unif_rand(),cum=0.0;
  for(size_t k=0;k<logp.size();k++){
    cum+=std::exp(logp[k]);
    if(u<=cum) return (int)k;
  }
  return (int)logp.size()-1;
}

static int irj_find_cut(const std::vector<double>&cuts,double cut){
  for(size_t k=0;k<cuts.size();k++) if(cuts[k]==cut) return (int)k;
  return -1;
}

// ---- hard MPPT -------------------------------------------------------------

static void irj_hard_cut_logp(const Node&node,const std::vector<int>&affected,
    const arma::mat&points,int axis,double kappa,double xi_left,
    double xi_right,const std::vector<double>&cuts,
    std::vector<double>&logp){
  std::vector<double> values;
  values.reserve(affected.size());
  for(int i:affected) values.push_back(points(i,axis));
  std::sort(values.begin(),values.end());
  double lower=node.box(axis,0),width=node.box(axis,1)-lower;
  logp.assign(cuts.size(),0.0);
  for(size_t k=0;k<cuts.size();k++){
    int n_left=(int)(std::lower_bound(values.begin(),values.end(),cuts[k])-
                     values.begin());
    int n_right=(int)values.size()-n_left;
    double A_left=node.area*(cuts[k]-lower)/width;
    double A_right=node.area-A_left;
    logp[k]=lg(n_left,A_left,kappa,xi_left)+
            lg(n_right,A_right,kappa,xi_right);
  }
  irj_normalize(logp);
}

static int irj_hard_grow_prune(Tree&tree,std::vector<int>&labels,
    const arma::mat&points,double kappa,double tau,double alpha,double eta,
    int max_depth,int min_leaf_n,int cut_mode,int cut_candidates,
    double a_xi,double b_xi,int model,int&which_move,
    bool conditional_labels){
  int d=points.n_cols,n=points.n_rows;
  std::vector<int> growable_nodes,prunable_nodes;
  growable(tree,points,max_depth,min_leaf_n,cut_mode,cut_candidates,
           growable_nodes);
  prunable(tree,prunable_nodes);
  bool do_grow=unif_rand()<0.5;
  which_move=do_grow?0:1;

  if(do_grow){
    if(growable_nodes.empty()) return 0;
    int parent=growable_nodes[runif_int(growable_nodes.size())];
    int depth=tree[parent].depth,left=2*parent,right=left+1;
    int axis=runif_int(d);
    std::vector<double> cuts=axis_cuts(
      points,tree[parent].idx,tree[parent].box,axis,cut_mode,
      cut_candidates,min_leaf_n
    );
    if(cuts.empty()) return 0;
    double xi_parent=tree[parent].xi,A_parent=tree[parent].area;
    int old_count=tree[parent].m;
    std::vector<int> affected;
    for(int i=0;i<n;i++) if(labels[i]==parent) affected.push_back(i);
    std::vector<double> logp;
    irj_hard_cut_logp(tree[parent],affected,points,axis,kappa,xi_parent,
                      xi_parent,cuts,logp);
    int chosen=irj_pick(logp);
    double cut=cuts[chosen];

    arma::vec coordinate=points.col(axis);
    std::vector<arma::uword> left_index,right_index;
    for(arma::uword k=0;k<tree[parent].idx.n_elem;k++){
      arma::uword i=tree[parent].idx[k];
      if(coordinate[i]<cut) left_index.push_back(i);
      else right_index.push_back(i);
    }
    arma::uvec left_u(left_index.size()),right_u(right_index.size());
    for(size_t k=0;k<left_index.size();k++) left_u[k]=left_index[k];
    for(size_t k=0;k<right_index.size();k++) right_u[k]=right_index[k];
    arma::mat left_box=tree[parent].box; left_box(axis,1)=cut;
    arma::mat right_box=tree[parent].box; right_box(axis,0)=cut;
    double A_left=box_area(left_box),A_right=box_area(right_box);
    double xi_left,xi_right;
    if(model==1){
      xi_left=R::rgamma(tau,xi_parent/tau);
      xi_right=R::rgamma(tau,xi_parent/tau);
    }else{
      double scale=b_xi>0?1.0/b_xi:1.0;
      xi_left=R::rgamma(a_xi,scale);
      xi_right=R::rgamma(a_xi,scale);
    }

    std::vector<std::pair<int,int> > reassignment;
    int m_parent=0,m_left=0,m_right=0;
    double log_q=0.0;
    for(int i:affected){
      bool goes_left=coordinate[i]<cut;
      double lp_parent=0.0,lp_child=0.0;
      if(conditional_labels){
        lp_parent=std::log(kappa+m_parent)-
                  std::log(kappa/xi_parent+A_parent);
        lp_child=goes_left
          ? std::log(kappa+m_left)-std::log(kappa/xi_left+A_left)
          : std::log(kappa+m_right)-std::log(kappa/xi_right+A_right);
      }
      double denominator=irj_lse2(lp_parent,lp_child);
      if(std::log(unif_rand())<lp_parent-denominator){
        reassignment.push_back(std::make_pair(i,parent));
        log_q+=lp_parent-denominator;
        m_parent++;
      }else{
        log_q+=lp_child-denominator;
        if(goes_left){
          reassignment.push_back(std::make_pair(i,left)); m_left++;
        }else{
          reassignment.push_back(std::make_pair(i,right)); m_right++;
        }
      }
    }

    Node left_node;
    left_node.box=left_box; left_node.idx=left_u; left_node.cut=NA_REAL;
    left_node.xi=xi_left; left_node.area=A_left; left_node.axis=-1;
    left_node.depth=depth+1; left_node.m=0;
    Node right_node;
    right_node.box=right_box; right_node.idx=right_u; right_node.cut=NA_REAL;
    right_node.xi=xi_right; right_node.area=A_right; right_node.axis=-1;
    right_node.depth=depth+1; right_node.m=0;
    tree[left]=left_node; tree[right]=right_node;
    tree[parent].axis=axis; tree[parent].cut=cut;
    std::vector<int> new_prunable;
    prunable(tree,new_prunable);
    double left_stop=lSfac(left_box,left_u,depth+1,points,max_depth,
      min_leaf_n,cut_mode,cut_candidates,alpha,eta);
    double right_stop=lSfac(right_box,right_u,depth+1,points,max_depth,
      min_leaf_n,cut_mode,cut_candidates,alpha,eta);
    double log_accept=lrho(depth,alpha,eta)-l1mrho(depth,alpha,eta)+
      left_stop+right_stop+
      lg(m_parent,A_parent,kappa,xi_parent)+
      lg(m_left,A_left,kappa,xi_left)+
      lg(m_right,A_right,kappa,xi_right)-
      lg(old_count,A_parent,kappa,xi_parent)+
      std::log((double)growable_nodes.size())-
      std::log((double)new_prunable.size())-
      log_q-std::log((double)cuts.size())-logp[chosen];
    if(std::log(unif_rand())<log_accept){
      tree[parent].m=m_parent; tree[left].m=m_left; tree[right].m=m_right;
      for(const auto&entry:reassignment) labels[entry.first]=entry.second;
      return 1;
    }
    tree.erase(left); tree.erase(right);
    tree[parent].axis=-1; tree[parent].cut=NA_REAL;
    return 0;
  }

  if(prunable_nodes.empty()) return 0;
  int parent=prunable_nodes[runif_int(prunable_nodes.size())];
  int depth=tree[parent].depth,left=2*parent,right=left+1;
  int m_parent=tree[parent].m,m_left=tree[left].m,m_right=tree[right].m;
  int merged_count=m_parent+m_left+m_right;
  double A_parent=tree[parent].area,A_left=tree[left].area,A_right=tree[right].area;
  double xi_parent=tree[parent].xi,xi_left=tree[left].xi,xi_right=tree[right].xi;
  std::vector<double> cuts=axis_cuts(
    points,tree[parent].idx,tree[parent].box,tree[parent].axis,cut_mode,
    cut_candidates,min_leaf_n
  );
  int old_index=irj_find_cut(cuts,tree[parent].cut);
  if(old_index<0) return 0;
  std::vector<int> affected;
  for(int i=0;i<n;i++)
    if(labels[i]==parent||labels[i]==left||labels[i]==right)
      affected.push_back(i);
  std::vector<double> logp;
  // This is the reverse grow proposal, which scores prospective children
  // using the parent scale before their scales are drawn from the prior.
  irj_hard_cut_logp(tree[parent],affected,points,tree[parent].axis,kappa,
                    xi_parent,xi_parent,cuts,logp);
  arma::vec coordinate=points.col(tree[parent].axis);
  int replay_parent=0,replay_left=0,replay_right=0;
  double reverse_log_q=0.0;
  for(int i:affected){
    bool goes_left=coordinate[i]<tree[parent].cut;
    double lp_parent=0.0,lp_child=0.0;
    if(conditional_labels){
      lp_parent=std::log(kappa+replay_parent)-
                std::log(kappa/xi_parent+A_parent);
      lp_child=goes_left
        ? std::log(kappa+replay_left)-std::log(kappa/xi_left+A_left)
        : std::log(kappa+replay_right)-std::log(kappa/xi_right+A_right);
    }
    double denominator=irj_lse2(lp_parent,lp_child);
    if(labels[i]==parent){
      reverse_log_q+=lp_parent-denominator; replay_parent++;
    }else{
      reverse_log_q+=lp_child-denominator;
      if(goes_left) replay_left++; else replay_right++;
    }
  }
  double left_stop=lSfac(tree[left].box,tree[left].idx,depth+1,points,
    max_depth,min_leaf_n,cut_mode,cut_candidates,alpha,eta);
  double right_stop=lSfac(tree[right].box,tree[right].idx,depth+1,points,
    max_depth,min_leaf_n,cut_mode,cut_candidates,alpha,eta);
  Node old_left=tree[left],old_right=tree[right];
  int old_axis=tree[parent].axis; double old_cut=tree[parent].cut;
  tree.erase(left); tree.erase(right);
  tree[parent].axis=-1; tree[parent].cut=NA_REAL;
  std::vector<int> new_growable;
  growable(tree,points,max_depth,min_leaf_n,cut_mode,cut_candidates,new_growable);
  if(new_growable.empty()){
    tree[left]=old_left; tree[right]=old_right;
    tree[parent].axis=old_axis; tree[parent].cut=old_cut;
    return 0;
  }
  double log_accept=-(lrho(depth,alpha,eta)-l1mrho(depth,alpha,eta)+
    left_stop+right_stop)+
    lg(merged_count,A_parent,kappa,xi_parent)-
    lg(m_parent,A_parent,kappa,xi_parent)-
    lg(m_left,A_left,kappa,xi_left)-lg(m_right,A_right,kappa,xi_right)+
    std::log((double)prunable_nodes.size())-
    std::log((double)new_growable.size())+reverse_log_q+
    std::log((double)cuts.size())+logp[old_index];
  if(std::log(unif_rand())<log_accept){
    for(size_t i=0;i<labels.size();i++)
      if(labels[i]==left||labels[i]==right) labels[i]=parent;
    tree[parent].m=merged_count;
    return 1;
  }
  tree[left]=old_left; tree[right]=old_right;
  tree[parent].axis=old_axis; tree[parent].cut=old_cut;
  return 0;
}

static int irj_hard_change_cut(Tree&tree,std::vector<int>&labels,
    const arma::mat&points,double kappa,double alpha,double eta,
    int max_depth,int min_leaf_n,int cut_mode,int cut_candidates,
    bool conditional_labels){
  int d=points.n_cols,n=points.n_rows;
  std::vector<int> prunable_nodes;
  prunable(tree,prunable_nodes);
  if(prunable_nodes.empty()) return 0;
  int parent=prunable_nodes[runif_int(prunable_nodes.size())];
  int depth=tree[parent].depth,left=2*parent,right=left+1;
  double A_parent=tree[parent].area,xi_parent=tree[parent].xi;
  double xi_left=tree[left].xi,xi_right=tree[right].xi;
  std::vector<int> affected;
  for(int i=0;i<n;i++)
    if(labels[i]==parent||labels[i]==left||labels[i]==right)
      affected.push_back(i);
  int new_axis=runif_int(d);
  std::vector<double> new_cuts=axis_cuts(
    points,tree[parent].idx,tree[parent].box,new_axis,cut_mode,
    cut_candidates,min_leaf_n
  );
  if(new_cuts.empty()) return 0;
  std::vector<double> new_logp;
  irj_hard_cut_logp(tree[parent],affected,points,new_axis,kappa,xi_left,
                    xi_right,new_cuts,new_logp);
  int new_index=irj_pick(new_logp);
  double new_cut=new_cuts[new_index];
  std::vector<double> old_cuts=axis_cuts(
    points,tree[parent].idx,tree[parent].box,tree[parent].axis,cut_mode,
    cut_candidates,min_leaf_n
  );
  int old_index=irj_find_cut(old_cuts,tree[parent].cut);
  if(old_index<0) return 0;
  std::vector<double> old_logp;
  irj_hard_cut_logp(tree[parent],affected,points,tree[parent].axis,kappa,
                    xi_left,xi_right,old_cuts,old_logp);

  arma::vec new_coordinate=points.col(new_axis);
  arma::vec old_coordinate=points.col(tree[parent].axis);
  arma::mat new_left_box=tree[parent].box; new_left_box(new_axis,1)=new_cut;
  arma::mat new_right_box=tree[parent].box; new_right_box(new_axis,0)=new_cut;
  double new_A_left=box_area(new_left_box),new_A_right=box_area(new_right_box);
  double old_A_left=tree[left].area,old_A_right=tree[right].area;
  double old_cut=tree[parent].cut;
  std::vector<arma::uword> left_index,right_index;
  for(arma::uword k=0;k<tree[parent].idx.n_elem;k++){
    arma::uword i=tree[parent].idx[k];
    if(new_coordinate[i]<new_cut) left_index.push_back(i);
    else right_index.push_back(i);
  }
  arma::uvec left_u(left_index.size()),right_u(right_index.size());
  for(size_t k=0;k<left_index.size();k++) left_u[k]=left_index[k];
  for(size_t k=0;k<right_index.size();k++) right_u[k]=right_index[k];

  std::vector<std::pair<int,int> > reassignment;
  int new_parent=0,new_left=0,new_right=0;
  int old_parent=0,old_left=0,old_right=0;
  double forward_log_q=0.0,reverse_log_q=0.0;
  for(int i:affected){
    bool goes_new_left=new_coordinate[i]<new_cut;
    double lp_parent=0.0,lp_child=0.0;
    if(conditional_labels){
      lp_parent=std::log(kappa+new_parent)-
                std::log(kappa/xi_parent+A_parent);
      lp_child=goes_new_left
        ? std::log(kappa+new_left)-std::log(kappa/xi_left+new_A_left)
        : std::log(kappa+new_right)-std::log(kappa/xi_right+new_A_right);
    }
    double denominator=irj_lse2(lp_parent,lp_child);
    if(std::log(unif_rand())<lp_parent-denominator){
      reassignment.push_back(std::make_pair(i,parent));
      forward_log_q+=lp_parent-denominator; new_parent++;
    }else{
      forward_log_q+=lp_child-denominator;
      if(goes_new_left){
        reassignment.push_back(std::make_pair(i,left)); new_left++;
      }else{
        reassignment.push_back(std::make_pair(i,right)); new_right++;
      }
    }

    bool goes_old_left=old_coordinate[i]<old_cut;
    double old_lp_parent=0.0,old_lp_child=0.0;
    if(conditional_labels){
      old_lp_parent=std::log(kappa+old_parent)-
                    std::log(kappa/xi_parent+A_parent);
      old_lp_child=goes_old_left
        ? std::log(kappa+old_left)-std::log(kappa/xi_left+old_A_left)
        : std::log(kappa+old_right)-std::log(kappa/xi_right+old_A_right);
    }
    double old_denominator=irj_lse2(old_lp_parent,old_lp_child);
    if(labels[i]==parent){
      reverse_log_q+=old_lp_parent-old_denominator; old_parent++;
    }else{
      reverse_log_q+=old_lp_child-old_denominator;
      if(goes_old_left) old_left++; else old_right++;
    }
  }

  double old_target=lg(tree[parent].m,A_parent,kappa,xi_parent)+
    lg(tree[left].m,old_A_left,kappa,xi_left)+
    lg(tree[right].m,old_A_right,kappa,xi_right);
  double new_target=lg(new_parent,A_parent,kappa,xi_parent)+
    lg(new_left,new_A_left,kappa,xi_left)+
    lg(new_right,new_A_right,kappa,xi_right);
  double old_stop=lSfac(tree[left].box,tree[left].idx,depth+1,points,
    max_depth,min_leaf_n,cut_mode,cut_candidates,alpha,eta)+
    lSfac(tree[right].box,tree[right].idx,depth+1,points,max_depth,
      min_leaf_n,cut_mode,cut_candidates,alpha,eta);
  double new_stop=lSfac(new_left_box,left_u,depth+1,points,max_depth,
    min_leaf_n,cut_mode,cut_candidates,alpha,eta)+
    lSfac(new_right_box,right_u,depth+1,points,max_depth,min_leaf_n,
      cut_mode,cut_candidates,alpha,eta);
  double log_accept=new_target-old_target+new_stop-old_stop+
    reverse_log_q-forward_log_q-
    std::log((double)new_cuts.size())-new_logp[new_index]+
    std::log((double)old_cuts.size())+old_logp[old_index];
  if(std::log(unif_rand())<log_accept){
    tree[parent].axis=new_axis; tree[parent].cut=new_cut;
    tree[parent].m=new_parent;
    tree[left].box=new_left_box; tree[left].area=new_A_left;
    tree[left].idx=left_u; tree[left].m=new_left;
    tree[right].box=new_right_box; tree[right].area=new_A_right;
    tree[right].idx=right_u; tree[right].m=new_right;
    for(const auto&entry:reassignment) labels[entry.first]=entry.second;
    return 1;
  }
  return 0;
}

// ---- soft S-MPPT -----------------------------------------------------------

static void irj_soft_children(const Node&parent,int axis,double cut,
    int gate_mode,double gate_depth,Node&left,Node&right){
  left=Node(); right=Node();
  left.box=parent.box; left.box(axis,1)=cut;
  right.box=parent.box; right.box(axis,0)=cut;
  left.soft_path=parent.soft_path; right.soft_path=parent.soft_path;
  left.depth=parent.depth+1; right.depth=parent.depth+1;
  left.axis=-1; right.axis=-1; left.m=0; right.m=0;
  if(gate_mode!=0){
    double parent_width=parent.box(axis,1)-parent.box(axis,0);
    if(gate_mode==1)
      parent_width/=std::pow(1.0+(double)parent.depth,gate_depth);
    SoftGate left_gate={axis,cut,parent_width,-1};
    SoftGate right_gate={axis,cut,parent_width,1};
    left.soft_path.push_back(left_gate);
    right.soft_path.push_back(right_gate);
  }
}

static void irj_soft_cut_logp(const Node&parent,
    const std::vector<int>&affected,const arma::mat&points,
    const arma::mat&region,double kappa,const arma::vec&gate,int gate_mode,
    double gate_depth,int axis,double xi_left,double xi_right,
    const std::vector<double>&cuts,std::vector<double>&logp){
  std::vector<double> values;
  values.reserve(affected.size());
  for(int i:affected) values.push_back(points(i,axis));
  std::sort(values.begin(),values.end());
  logp.assign(cuts.size(),0.0);
  for(size_t k=0;k<cuts.size();k++){
    int n_left=(int)(std::lower_bound(values.begin(),values.end(),cuts[k])-
                     values.begin());
    int n_right=(int)values.size()-n_left;
    Node left,right;
    irj_soft_children(parent,axis,cuts[k],gate_mode,gate_depth,left,right);
    // Use the same closed-form exposure routine as the target density. Its
    // logistic implementation has a deterministic adaptive-quadrature
    // fallback for ill-conditioned partial fractions; no box-volume surrogate
    // enters the production proposal score.
    double A_left=mpps_exposure_node(left,region,gate,gate_mode);
    double A_right=mpps_exposure_node(right,region,gate,gate_mode);
    logp[k]=lg(n_left,A_left,kappa,xi_left)+
            lg(n_right,A_right,kappa,xi_right);
  }
  irj_normalize(logp);
}

static int irj_soft_grow_prune(Tree&tree,std::vector<int>&labels,
    const arma::mat&points,const arma::mat&region,double kappa,double tau,
    double a_xi,double b_xi,int scale_model,const arma::vec&gate,
    int gate_mode,double gate_depth,double alpha,double eta,int max_depth,
    int min_leaf_n,int cut_mode,int cut_candidates,int&which_move,
    bool conditional_labels){
  int d=points.n_cols,n=points.n_rows;
  std::vector<int> growable_nodes,prunable_nodes;
  growable(tree,points,max_depth,min_leaf_n,cut_mode,cut_candidates,
           growable_nodes);
  prunable(tree,prunable_nodes);
  bool do_grow=unif_rand()<0.5;
  which_move=do_grow?0:1;

  if(do_grow){
    if(growable_nodes.empty()) return 0;
    int parent=growable_nodes[runif_int(growable_nodes.size())];
    int depth=tree[parent].depth,left=2*parent,right=left+1;
    int axis=runif_int(d);
    std::vector<double> cuts=axis_cuts(
      points,tree[parent].idx,tree[parent].box,axis,cut_mode,
      cut_candidates,min_leaf_n
    );
    if(cuts.empty()) return 0;
    double xi_parent=tree[parent].xi;
    std::vector<int> affected;
    for(int i=0;i<n;i++) if(labels[i]==parent) affected.push_back(i);
    std::vector<double> logp;
    irj_soft_cut_logp(tree[parent],affected,points,region,kappa,gate,
      gate_mode,gate_depth,axis,xi_parent,xi_parent,cuts,logp);
    int chosen=irj_pick(logp);
    double cut=cuts[chosen];
    arma::vec coordinate=points.col(axis);
    std::vector<arma::uword> left_index,right_index;
    for(arma::uword k=0;k<tree[parent].idx.n_elem;k++){
      arma::uword i=tree[parent].idx[k];
      if(coordinate[i]<cut) left_index.push_back(i);
      else right_index.push_back(i);
    }
    arma::uvec left_u(left_index.size()),right_u(right_index.size());
    for(size_t k=0;k<left_index.size();k++) left_u[k]=left_index[k];
    for(size_t k=0;k<right_index.size();k++) right_u[k]=right_index[k];
    Node left_node,right_node;
    irj_soft_children(tree[parent],axis,cut,gate_mode,gate_depth,
                      left_node,right_node);
    left_node.idx=left_u; right_node.idx=right_u;
    left_node.area=box_area(left_node.box);
    right_node.area=box_area(right_node.box);
    left_node.cut=NA_REAL; right_node.cut=NA_REAL;
    double xi_left=scale_model==1
      ? R::rgamma(tau,xi_parent/tau)
      : R::rgamma(a_xi,1.0/b_xi);
    double xi_right=scale_model==1
      ? R::rgamma(tau,xi_parent/tau)
      : R::rgamma(a_xi,1.0/b_xi);
    left_node.xi=xi_left; right_node.xi=xi_right;
    int old_count=tree[parent].m;
    Node old_parent=tree[parent];
    double A_parent=mpps_exposure_node(old_parent,region,gate,gate_mode);
    double A_left=mpps_exposure_node(left_node,region,gate,gate_mode);
    double A_right=mpps_exposure_node(right_node,region,gate,gate_mode);
    tree[left]=left_node; tree[right]=right_node;
    tree[parent].axis=axis; tree[parent].cut=cut;

    int m_parent=0,m_left=0,m_right=0;
    double log_q=0.0,spatial_difference=0.0;
    std::vector<std::pair<int,int> > reassignment;
    for(int i:affected){
      double log_phi_parent=mpps_log_phi_node(
        old_parent,points.row(i),region,gate,gate_mode
      );
      double log_phi_left=mpps_log_phi_node(
        tree.at(left),points.row(i),region,gate,gate_mode
      );
      double log_phi_right=mpps_log_phi_node(
        tree.at(right),points.row(i),region,gate,gate_mode
      );
      double lp_parent,lp_left,lp_right;
      if(conditional_labels){
        lp_parent=log_phi_parent+std::log(kappa+m_parent)-
                  std::log(kappa/xi_parent+A_parent);
        lp_left=log_phi_left+std::log(kappa+m_left)-
                std::log(kappa/xi_left+A_left);
        lp_right=log_phi_right+std::log(kappa+m_right)-
                 std::log(kappa/xi_right+A_right);
      }else{
        double child_denominator=irj_lse2(log_phi_left,log_phi_right);
        lp_parent=std::log(0.5);
        lp_left=std::log(0.5)+log_phi_left-child_denominator;
        lp_right=std::log(0.5)+log_phi_right-child_denominator;
      }
      double denominator=irj_lse3(lp_parent,lp_left,lp_right);
      double log_u=std::log(unif_rand());
      if(log_u<lp_parent-denominator){
        reassignment.push_back(std::make_pair(i,parent));
        log_q+=lp_parent-denominator; m_parent++;
      }else if(log_u<irj_lse2(lp_parent,lp_left)-denominator){
        reassignment.push_back(std::make_pair(i,left));
        log_q+=lp_left-denominator; m_left++;
        spatial_difference+=log_phi_left-log_phi_parent;
      }else{
        reassignment.push_back(std::make_pair(i,right));
        log_q+=lp_right-denominator; m_right++;
        spatial_difference+=log_phi_right-log_phi_parent;
      }
    }
    std::vector<int> new_prunable;
    prunable(tree,new_prunable);
    double left_stop=lSfac(left_node.box,left_u,depth+1,points,max_depth,
      min_leaf_n,cut_mode,cut_candidates,alpha,eta);
    double right_stop=lSfac(right_node.box,right_u,depth+1,points,max_depth,
      min_leaf_n,cut_mode,cut_candidates,alpha,eta);
    double log_accept=lrho(depth,alpha,eta)-l1mrho(depth,alpha,eta)+
      left_stop+right_stop+
      lg(m_parent,A_parent,kappa,xi_parent)+
      lg(m_left,A_left,kappa,xi_left)+
      lg(m_right,A_right,kappa,xi_right)-
      lg(old_count,A_parent,kappa,xi_parent)+spatial_difference+
      std::log((double)growable_nodes.size())-
      std::log((double)new_prunable.size())-
      log_q-std::log((double)cuts.size())-logp[chosen];
    if(std::log(unif_rand())<log_accept){
      tree[parent].m=m_parent; tree[left].m=m_left; tree[right].m=m_right;
      for(const auto&entry:reassignment) labels[entry.first]=entry.second;
      return 1;
    }
    tree.erase(left); tree.erase(right); tree[parent]=old_parent;
    return 0;
  }

  if(prunable_nodes.empty()) return 0;
  int parent=prunable_nodes[runif_int(prunable_nodes.size())];
  int depth=tree[parent].depth,left=2*parent,right=left+1;
  Node old_parent=tree[parent],old_left=tree[left],old_right=tree[right];
  int merged_count=old_parent.m+old_left.m+old_right.m;
  std::vector<double> cuts=axis_cuts(
    points,old_parent.idx,old_parent.box,old_parent.axis,cut_mode,
    cut_candidates,min_leaf_n
  );
  int old_index=irj_find_cut(cuts,old_parent.cut);
  if(old_index<0) return 0;
  std::vector<int> affected;
  for(int i=0;i<n;i++)
    if(labels[i]==parent||labels[i]==left||labels[i]==right)
      affected.push_back(i);
  Node parent_leaf=old_parent;
  parent_leaf.axis=-1; parent_leaf.cut=NA_REAL;
  std::vector<double> logp;
  irj_soft_cut_logp(parent_leaf,affected,points,region,kappa,gate,gate_mode,
    gate_depth,old_parent.axis,old_parent.xi,old_parent.xi,cuts,logp);
  double A_parent=mpps_exposure_node(old_parent,region,gate,gate_mode);
  double A_left=mpps_exposure_node(old_left,region,gate,gate_mode);
  double A_right=mpps_exposure_node(old_right,region,gate,gate_mode);
  int replay_parent=0,replay_left=0,replay_right=0;
  double reverse_log_q=0.0,spatial_difference=0.0;
  for(int i:affected){
    double log_phi_parent=mpps_log_phi_node(
      old_parent,points.row(i),region,gate,gate_mode
    );
    double log_phi_left=mpps_log_phi_node(
      old_left,points.row(i),region,gate,gate_mode
    );
    double log_phi_right=mpps_log_phi_node(
      old_right,points.row(i),region,gate,gate_mode
    );
    double lp_parent,lp_left,lp_right;
    if(conditional_labels){
      lp_parent=log_phi_parent+std::log(kappa+replay_parent)-
                std::log(kappa/old_parent.xi+A_parent);
      lp_left=log_phi_left+std::log(kappa+replay_left)-
              std::log(kappa/old_left.xi+A_left);
      lp_right=log_phi_right+std::log(kappa+replay_right)-
               std::log(kappa/old_right.xi+A_right);
    }else{
      double child_denominator=irj_lse2(log_phi_left,log_phi_right);
      lp_parent=std::log(0.5);
      lp_left=std::log(0.5)+log_phi_left-child_denominator;
      lp_right=std::log(0.5)+log_phi_right-child_denominator;
    }
    double denominator=irj_lse3(lp_parent,lp_left,lp_right);
    if(labels[i]==parent){
      reverse_log_q+=lp_parent-denominator; replay_parent++;
    }else if(labels[i]==left){
      reverse_log_q+=lp_left-denominator; replay_left++;
      spatial_difference+=log_phi_parent-log_phi_left;
    }else{
      reverse_log_q+=lp_right-denominator; replay_right++;
      spatial_difference+=log_phi_parent-log_phi_right;
    }
  }
  double left_stop=lSfac(old_left.box,old_left.idx,depth+1,points,max_depth,
    min_leaf_n,cut_mode,cut_candidates,alpha,eta);
  double right_stop=lSfac(old_right.box,old_right.idx,depth+1,points,max_depth,
    min_leaf_n,cut_mode,cut_candidates,alpha,eta);
  tree.erase(left); tree.erase(right);
  tree[parent].axis=-1; tree[parent].cut=NA_REAL;
  std::vector<int> new_growable;
  growable(tree,points,max_depth,min_leaf_n,cut_mode,cut_candidates,new_growable);
  if(new_growable.empty()){
    tree[parent]=old_parent; tree[left]=old_left; tree[right]=old_right;
    return 0;
  }
  double log_accept=-(lrho(depth,alpha,eta)-l1mrho(depth,alpha,eta)+
    left_stop+right_stop)+
    lg(merged_count,A_parent,kappa,old_parent.xi)-
    lg(old_parent.m,A_parent,kappa,old_parent.xi)-
    lg(old_left.m,A_left,kappa,old_left.xi)-
    lg(old_right.m,A_right,kappa,old_right.xi)+spatial_difference+
    std::log((double)prunable_nodes.size())-
    std::log((double)new_growable.size())+reverse_log_q+
    std::log((double)cuts.size())+logp[old_index];
  if(std::log(unif_rand())<log_accept){
    for(size_t i=0;i<labels.size();i++)
      if(labels[i]==left||labels[i]==right) labels[i]=parent;
    tree[parent].m=merged_count;
    return 1;
  }
  tree[parent]=old_parent; tree[left]=old_left; tree[right]=old_right;
  return 0;
}

static int irj_soft_change_cut(Tree&tree,std::vector<int>&labels,
    const arma::mat&points,const arma::mat&region,double kappa,
    const arma::vec&gate,int gate_mode,double gate_depth,double alpha,
    double eta,int max_depth,int min_leaf_n,int cut_mode,int cut_candidates,
    bool conditional_labels){
  int d=points.n_cols,n=points.n_rows;
  std::vector<int> prunable_nodes;
  prunable(tree,prunable_nodes);
  if(prunable_nodes.empty()) return 0;
  int parent=prunable_nodes[runif_int(prunable_nodes.size())];
  int depth=tree[parent].depth,left=2*parent,right=left+1;
  Node parent_leaf=tree[parent];
  parent_leaf.axis=-1; parent_leaf.cut=NA_REAL;
  std::vector<int> affected;
  for(int i=0;i<n;i++)
    if(labels[i]==parent||labels[i]==left||labels[i]==right)
      affected.push_back(i);
  int new_axis=runif_int(d);
  std::vector<double> new_cuts=axis_cuts(
    points,tree[parent].idx,tree[parent].box,new_axis,cut_mode,
    cut_candidates,min_leaf_n
  );
  if(new_cuts.empty()) return 0;
  std::vector<double> new_logp;
  irj_soft_cut_logp(parent_leaf,affected,points,region,kappa,gate,gate_mode,
    gate_depth,new_axis,tree[left].xi,tree[right].xi,new_cuts,new_logp);
  int new_index=irj_pick(new_logp);
  double new_cut=new_cuts[new_index];
  std::vector<double> old_cuts=axis_cuts(
    points,tree[parent].idx,tree[parent].box,tree[parent].axis,cut_mode,
    cut_candidates,min_leaf_n
  );
  int old_index=irj_find_cut(old_cuts,tree[parent].cut);
  if(old_index<0) return 0;
  std::vector<double> old_logp;
  irj_soft_cut_logp(parent_leaf,affected,points,region,kappa,gate,gate_mode,
    gate_depth,tree[parent].axis,tree[left].xi,tree[right].xi,
    old_cuts,old_logp);

  Node old_parent=tree[parent],old_left=tree[left],old_right=tree[right];
  Node new_left,new_right;
  irj_soft_children(old_parent,new_axis,new_cut,gate_mode,gate_depth,
                    new_left,new_right);
  arma::vec coordinate=points.col(new_axis);
  std::vector<arma::uword> left_index,right_index;
  for(arma::uword k=0;k<old_parent.idx.n_elem;k++){
    arma::uword i=old_parent.idx[k];
    if(coordinate[i]<new_cut) left_index.push_back(i);
    else right_index.push_back(i);
  }
  arma::uvec left_u(left_index.size()),right_u(right_index.size());
  for(size_t k=0;k<left_index.size();k++) left_u[k]=left_index[k];
  for(size_t k=0;k<right_index.size();k++) right_u[k]=right_index[k];
  new_left.idx=left_u; new_right.idx=right_u;
  new_left.area=box_area(new_left.box); new_right.area=box_area(new_right.box);
  new_left.cut=NA_REAL; new_right.cut=NA_REAL;
  new_left.xi=old_left.xi; new_right.xi=old_right.xi;
  double A_parent=mpps_exposure_node(old_parent,region,gate,gate_mode);
  double old_A_left=mpps_exposure_node(old_left,region,gate,gate_mode);
  double old_A_right=mpps_exposure_node(old_right,region,gate,gate_mode);
  double new_A_left=mpps_exposure_node(new_left,region,gate,gate_mode);
  double new_A_right=mpps_exposure_node(new_right,region,gate,gate_mode);

  int new_parent_count=0,new_left_count=0,new_right_count=0;
  int old_parent_count=0,old_left_count=0,old_right_count=0;
  double forward_log_q=0.0,reverse_log_q=0.0,spatial_difference=0.0;
  std::vector<std::pair<int,int> > reassignment;
  for(int i:affected){
    double log_phi_parent=mpps_log_phi_node(
      old_parent,points.row(i),region,gate,gate_mode
    );
    double new_log_phi_left=mpps_log_phi_node(
      new_left,points.row(i),region,gate,gate_mode
    );
    double new_log_phi_right=mpps_log_phi_node(
      new_right,points.row(i),region,gate,gate_mode
    );
    double lp_parent,lp_left,lp_right;
    if(conditional_labels){
      lp_parent=log_phi_parent+std::log(kappa+new_parent_count)-
                std::log(kappa/old_parent.xi+A_parent);
      lp_left=new_log_phi_left+std::log(kappa+new_left_count)-
              std::log(kappa/new_left.xi+new_A_left);
      lp_right=new_log_phi_right+std::log(kappa+new_right_count)-
               std::log(kappa/new_right.xi+new_A_right);
    }else{
      double child_denominator=irj_lse2(new_log_phi_left,new_log_phi_right);
      lp_parent=std::log(0.5);
      lp_left=std::log(0.5)+new_log_phi_left-child_denominator;
      lp_right=std::log(0.5)+new_log_phi_right-child_denominator;
    }
    double denominator=irj_lse3(lp_parent,lp_left,lp_right);
    double log_u=std::log(unif_rand());
    int proposed;
    if(log_u<lp_parent-denominator){
      proposed=parent; forward_log_q+=lp_parent-denominator;
      new_parent_count++;
    }else if(log_u<irj_lse2(lp_parent,lp_left)-denominator){
      proposed=left; forward_log_q+=lp_left-denominator;
      new_left_count++;
    }else{
      proposed=right; forward_log_q+=lp_right-denominator;
      new_right_count++;
    }
    reassignment.push_back(std::make_pair(i,proposed));

    double old_log_phi_left=mpps_log_phi_node(
      old_left,points.row(i),region,gate,gate_mode
    );
    double old_log_phi_right=mpps_log_phi_node(
      old_right,points.row(i),region,gate,gate_mode
    );
    double old_lp_parent,old_lp_left,old_lp_right;
    if(conditional_labels){
      old_lp_parent=log_phi_parent+std::log(kappa+old_parent_count)-
                    std::log(kappa/old_parent.xi+A_parent);
      old_lp_left=old_log_phi_left+std::log(kappa+old_left_count)-
                  std::log(kappa/old_left.xi+old_A_left);
      old_lp_right=old_log_phi_right+std::log(kappa+old_right_count)-
                   std::log(kappa/old_right.xi+old_A_right);
    }else{
      double child_denominator=irj_lse2(old_log_phi_left,old_log_phi_right);
      old_lp_parent=std::log(0.5);
      old_lp_left=std::log(0.5)+old_log_phi_left-child_denominator;
      old_lp_right=std::log(0.5)+old_log_phi_right-child_denominator;
    }
    double old_denominator=irj_lse3(
      old_lp_parent,old_lp_left,old_lp_right
    );
    if(labels[i]==parent){
      reverse_log_q+=old_lp_parent-old_denominator; old_parent_count++;
    }else if(labels[i]==left){
      reverse_log_q+=old_lp_left-old_denominator; old_left_count++;
    }else{
      reverse_log_q+=old_lp_right-old_denominator; old_right_count++;
    }
    const Node&new_node=proposed==parent?old_parent:
                         (proposed==left?new_left:new_right);
    const Node&old_node=labels[i]==parent?old_parent:
                         (labels[i]==left?old_left:old_right);
    spatial_difference+=mpps_log_phi_node(
      new_node,points.row(i),region,gate,gate_mode
    )-mpps_log_phi_node(old_node,points.row(i),region,gate,gate_mode);
  }

  double old_target=lg(old_parent.m,A_parent,kappa,old_parent.xi)+
    lg(old_left.m,old_A_left,kappa,old_left.xi)+
    lg(old_right.m,old_A_right,kappa,old_right.xi);
  double new_target=lg(new_parent_count,A_parent,kappa,old_parent.xi)+
    lg(new_left_count,new_A_left,kappa,new_left.xi)+
    lg(new_right_count,new_A_right,kappa,new_right.xi);
  double old_stop=lSfac(old_left.box,old_left.idx,depth+1,points,max_depth,
    min_leaf_n,cut_mode,cut_candidates,alpha,eta)+
    lSfac(old_right.box,old_right.idx,depth+1,points,max_depth,
      min_leaf_n,cut_mode,cut_candidates,alpha,eta);
  double new_stop=lSfac(new_left.box,new_left.idx,depth+1,points,max_depth,
    min_leaf_n,cut_mode,cut_candidates,alpha,eta)+
    lSfac(new_right.box,new_right.idx,depth+1,points,max_depth,
      min_leaf_n,cut_mode,cut_candidates,alpha,eta);
  double log_accept=new_target-old_target+new_stop-old_stop+
    spatial_difference+reverse_log_q-forward_log_q-
    std::log((double)new_cuts.size())-new_logp[new_index]+
    std::log((double)old_cuts.size())+old_logp[old_index];
  if(std::log(unif_rand())<log_accept){
    tree[parent].axis=new_axis; tree[parent].cut=new_cut;
    tree[parent].m=new_parent_count;
    tree[left]=new_left; tree[left].m=new_left_count;
    tree[right]=new_right; tree[right].m=new_right_count;
    for(const auto&entry:reassignment) labels[entry.first]=entry.second;
    return 1;
  }
  return 0;
}

#endif
