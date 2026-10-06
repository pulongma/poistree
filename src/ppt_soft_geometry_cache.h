#ifndef POISTREE_PPT_SOFT_GEOMETRY_CACHE_H
#define POISTREE_PPT_SOFT_GEOMETRY_CACHE_H

#include <string>
#include <unordered_set>

// Geometry is local to one chain and one gate vector.  Full path keys include
// the split parent's width, so a node-relative gate cannot reuse a value from
// a former ancestor geometry.  Counts and intensities are deliberately absent
// from this cache: changing allocations does not change a leaf basis.
struct PPSTGeometryEntry {
  arma::vec training,background;
  bool training_ready=false,background_ready=false,exposure_ready=false;
  double exposure=0.0;
};

class PPSTGeometryCache {
  const arma::mat *training_,*region_,*background_=nullptr,*prediction_=nullptr;
  const arma::vec *weights_=nullptr;
  arma::vec gate_;
  int family_,prediction_role_=-1;
  std::unordered_map<std::string,PPSTGeometryEntry> entries_;

  template<class Value> static void append(std::string&key,const Value&value){
    key.append(reinterpret_cast<const char*>(&value),sizeof(Value));
  }
  static std::string key(const std::vector<PPSTGate>&path){
    std::string out;
    out.reserve(path.size()*(2*sizeof(int)+2*sizeof(double)));
    for(const PPSTGate&g:path){
      append(out,g.axis);append(out,g.cut);append(out,g.parent_width);append(out,g.side);
    }
    return out;
  }
  static bool same_points(const arma::mat&a,const arma::mat&b){
    return a.n_rows==b.n_rows&&a.n_cols==b.n_cols&&
      (a.memptr()==b.memptr()||arma::approx_equal(a,b,"absdiff",0.0));
  }
  arma::vec& values(PPSTGeometryEntry&entry,bool background){
    return background?entry.background:entry.training;
  }
  bool& ready(PPSTGeometryEntry&entry,bool background){
    return background?entry.background_ready:entry.training_ready;
  }
  const arma::vec& basis(const std::vector<PPSTGate>&path,bool background){
    PPSTGeometryEntry&entry=entries_[key(path)];
    if(ready(entry,background)) return values(entry,background);
    const arma::mat&points=background?*background_:*training_;
    if(path.empty()){
      values(entry,background).zeros(points.n_rows);
      ready(entry,background)=true;
      return values(entry,background);
    }
    std::vector<PPSTGate> parent=path;
    const PPSTGate split=parent.back();parent.pop_back();
    const arma::vec&prefix=basis(parent,background);
    std::vector<PPSTGate> sibling=path;sibling.back().side=-split.side;
    PPSTGeometryEntry&other=entries_[key(sibling)];
    arma::vec&left=values(split.side<0?entry:other,background);
    arma::vec&right=values(split.side>0?entry:other,background);
    left.set_size(points.n_rows);right.set_size(points.n_rows);
    const double*coordinate=points.colptr(split.axis);
    if(family_==1){
      PPSTGate right_split=split;right_split.side=1;
      for(arma::uword i=0;i<points.n_rows;i++){
        const double q=ppst_compact_gate_value(right_split,coordinate[i],gate_[split.axis]);
        const double l=1.0-q;
        left[i]=l<=0.0?-std::numeric_limits<double>::infinity():prefix[i]+std::log(l);
        right[i]=q<=0.0?-std::numeric_limits<double>::infinity():prefix[i]+std::log(q);
      }
    }else{
      const double width=family_==2?split.parent_width:
        (*region_)(split.axis,1)-(*region_)(split.axis,0);
      for(arma::uword i=0;i<points.n_rows;i++){
        // Match the original multiplication/division and path addition order.
        // The expensive log1p(exp()) term is shared by the two children.
        const double z=gate_[split.axis]*(coordinate[i]-split.cut)/width;
        const double common=std::log1p(std::exp(-std::abs(z)));
        left[i]=prefix[i]+(-std::max(z,0.0)-common);
        right[i]=prefix[i]+(-std::max(-z,0.0)-common);
      }
    }
    ready(entry,background)=true;ready(other,background)=true;
    return values(entry,background);
  }

public:
  PPSTGeometryCache(const arma::mat&training,const arma::mat&region,
      const arma::vec&gate,int family,const arma::mat*prediction=nullptr):
    training_(&training),region_(&region),prediction_(prediction),gate_(gate),family_(family){
#ifdef POISTREE_SPATIAL_QUADRATURE_H
    if(qpp_active){background_=&qpp_background;weights_=&qpp_weights;}
#endif
    if(prediction_){
      if(background_&&same_points(*prediction_,*background_)) prediction_role_=1;
      else if(same_points(*prediction_,*training_)) prediction_role_=0;
    }
  }
  // A proposed gate vector starts empty; acceptance moves its cache into the
  // chain, whereas rejection destroys only the proposal's geometry.
  PPSTGeometryCache(const PPSTGeometryCache&current,const arma::vec&gate):
    training_(current.training_),region_(current.region_),
    background_(current.background_),prediction_(current.prediction_),
    weights_(current.weights_),gate_(gate),family_(current.family_),
    prediction_role_(current.prediction_role_){}

  const arma::vec& training(const PPSTNode&node){return basis(node.path,false);}
  double exposure(const PPSTNode&node){
    PPSTGeometryEntry&entry=entries_[key(node.path)];
    if(entry.exposure_ready) return entry.exposure;
    double H;
    if(background_){
      const arma::vec&lp=basis(node.path,true);
      H=0.0;
      // Keep the scalar quadrature summation and underflow convention intact.
      for(arma::uword i=0;i<weights_->n_elem;i++)
        H+=(*weights_)[i]*(lp[i]<-745.0?0.0:std::exp(lp[i]));
    }else H=ppst_exposure(node,*region_,gate_,family_);
    entry.exposure=H;entry.exposure_ready=true;
    return H;
  }
  void trim(const PPSTree&tree){
    std::unordered_set<std::string> active;active.reserve(tree.size());
    for(const auto&item:tree) active.insert(key(item.second.path));
    for(auto it=entries_.begin();it!=entries_.end();){
      if(active.find(it->first)==active.end()) it=entries_.erase(it);
      else ++it;
    }
  }
  arma::vec intensity(const PPSTree&tree,
      const std::unordered_map<int,double>&lambda,const arma::mat&points){
    int role=-1;
    if(&points==training_) role=0;
    else if(&points==background_) role=1;
    else if(&points==prediction_) role=prediction_role_;
    if(role<0){
      // Other prediction/test sets are temporary, bounded by that one set.
      PPSTGeometryCache temporary(points,*region_,gate_,family_);
      return temporary.intensity(tree,lambda,points);
    }
    arma::vec out(points.n_rows,arma::fill::zeros);
    // lambda's iteration order is also used by ppst_eval_intensity().
    for(const auto&item:lambda){
      const arma::vec&lp=basis(tree.at(item.first).path,role==1);
      for(arma::uword i=0;i<points.n_rows;i++)
        out[i]+=item.second*(lp[i]<-745.0?0.0:std::exp(lp[i]));
    }
    for(arma::uword i=0;i<out.n_elem;i++) out[i]=std::max(out[i],1e-300);
    return out;
  }
};

static inline double ppst_cached_exposure(const PPSTNode&node,
    const arma::mat&region,const arma::vec&gate,int family,PPSTGeometryCache*cache){
  return cache?cache->exposure(node):ppst_exposure(node,region,gate,family);
}

#endif
