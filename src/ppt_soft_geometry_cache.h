#ifndef POISTREE_PPT_SOFT_GEOMETRY_CACHE_H
#define POISTREE_PPT_SOFT_GEOMETRY_CACHE_H

#include <memory>
#include <string>
#include <unordered_set>

struct PPSTGeometryEntry {
  arma::vec training,background;
  bool training_ready=false,background_ready=false,exposure_ready=false;
  double exposure=0.0;
};

struct PPSTBackgroundRows {
  arma::uvec rows,unmatched;
  arma::uword matched=0;
};

class PPSTGeometryCache {
  const arma::mat *training_,*region_,*background_=nullptr,*prediction_=nullptr;
  const arma::vec *weights_=nullptr;
  arma::vec gate_;
  int family_,prediction_role_=-1;
  std::unordered_map<std::string,PPSTGeometryEntry> entries_;
  std::shared_ptr<const PPSTBackgroundRows> background_rows_;

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
  static std::string row_key(const arma::mat&points,arma::uword row){
    std::string out;
    out.reserve(points.n_cols*sizeof(double));
    for(arma::uword j=0;j<points.n_cols;j++) append(out,points(row,j));
    return out;
  }
  static std::shared_ptr<const PPSTBackgroundRows> map_background_rows(
      const arma::mat&training,const arma::mat&background){
    if(training.n_cols!=background.n_cols) return nullptr;
    auto out=std::make_shared<PPSTBackgroundRows>();
    out->rows.set_size(training.n_rows);
    std::unordered_map<std::string,arma::uword> lookup;
    lookup.reserve(background.n_rows);

    for(arma::uword i=0;i<background.n_rows;i++)
      lookup.emplace(row_key(background,i),i);
    std::vector<arma::uword> unmatched;
    unmatched.reserve(training.n_rows);
    for(arma::uword i=0;i<training.n_rows;i++){
      auto found=lookup.find(row_key(training,i));
      if(found==lookup.end()){
        out->rows[i]=background.n_rows;unmatched.push_back(i);
      }else{
        out->rows[i]=found->second;out->matched++;
      }
    }
    out->unmatched=arma::uvec(unmatched);
    return out;
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
      if(background&&family_!=1) values(entry,background).ones(points.n_rows);
      else values(entry,background).zeros(points.n_rows);
      ready(entry,background)=true;
      return values(entry,background);
    }
    std::vector<PPSTGate> parent=path;
    const PPSTGate split=parent.back();parent.pop_back();
    std::vector<PPSTGate> sibling=path;sibling.back().side=-split.side;
    PPSTGeometryEntry&other=entries_[key(sibling)];
    arma::vec&left=values(split.side<0?entry:other,background);
    arma::vec&right=values(split.side>0?entry:other,background);
    left.set_size(points.n_rows);right.set_size(points.n_rows);
    if(background&&family_!=1){
      const arma::vec&parent_phi=basis(parent,true);
      const double width=family_==2?split.parent_width:
        (*region_)(split.axis,1)-(*region_)(split.axis,0);
      const double*coordinate=points.colptr(split.axis);
      const double g=gate_[split.axis];
      for(arma::uword i=0;i<points.n_rows;i++){
        const double z=g*(coordinate[i]-split.cut)/width;
        const double e=std::exp(-std::abs(z)),big=1.0/(1.0+e),small=e*big;
        left[i]=parent_phi[i]*(z>0.0?small:big);
        right[i]=parent_phi[i]*(z>0.0?big:small);
      }
      ready(entry,true)=true;ready(other,true)=true;
      return values(entry,true);
    }

    const PPSTBackgroundRows*rows=family_==1&&!background&&background_rows_&&
      background_rows_->matched?background_rows_.get():nullptr;
    if(rows){

      basis(path,true);
      const arma::vec&left_background=values(split.side<0?entry:other,true);
      const arma::vec&right_background=values(split.side>0?entry:other,true);
      for(arma::uword i=0;i<points.n_rows;i++){
        const arma::uword q=rows->rows[i];
        if(q<background_->n_rows){
          left[i]=left_background[q];right[i]=right_background[q];
        }
      }
      if(rows->unmatched.n_elem==0){
        ready(entry,background)=true;ready(other,background)=true;
        return values(entry,background);
      }
    }
    const arma::vec&prefix=basis(parent,background);
    const arma::uword count=rows?rows->unmatched.n_elem:points.n_rows;
    const double*coordinate=points.colptr(split.axis);
    if(family_==1){
      PPSTGate right_split=split;right_split.side=1;
      for(arma::uword j=0;j<count;j++){
        const arma::uword i=rows?rows->unmatched[j]:j;
        const double q=ppst_compact_gate_value(right_split,coordinate[i],gate_[split.axis]);
        const double l=1.0-q;
        left[i]=l<=0.0?-std::numeric_limits<double>::infinity():prefix[i]+std::log(l);
        right[i]=q<=0.0?-std::numeric_limits<double>::infinity():prefix[i]+std::log(q);
      }
    }else{
      const double width=family_==2?split.parent_width:
        (*region_)(split.axis,1)-(*region_)(split.axis,0);
      for(arma::uword j=0;j<count;j++){
        const arma::uword i=rows?rows->unmatched[j]:j;

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
      const arma::vec&gate,int family,const arma::mat*prediction=nullptr,
      bool reuse_background=true):
    training_(&training),region_(&region),prediction_(prediction),gate_(gate),family_(family){
#ifdef POISTREE_SPATIAL_QUADRATURE_H
    if(qpp_active){background_=&qpp_background;weights_=&qpp_weights;}
#endif
    if(background_&&reuse_background)
      background_rows_=map_background_rows(training,*background_);
    if(prediction_){
      if(background_&&same_points(*prediction_,*background_)) prediction_role_=1;
      else if(same_points(*prediction_,*training_)) prediction_role_=0;
    }
  }
  // Create an empty geometry cache for a proposed gate vector.

  PPSTGeometryCache(const PPSTGeometryCache&current,const arma::vec&gate):
    training_(current.training_),region_(current.region_),
    background_(current.background_),prediction_(current.prediction_),
    weights_(current.weights_),gate_(gate),family_(current.family_),
    prediction_role_(current.prediction_role_),background_rows_(current.background_rows_){}

  const arma::uvec* training_background_rows() const{
    return background_rows_?&background_rows_->rows:nullptr;
  }
  const arma::vec& training(const PPSTNode&node){return basis(node.path,false);}
  double exposure(const PPSTNode&node){
    PPSTGeometryEntry&entry=entries_[key(node.path)];
    if(entry.exposure_ready) return entry.exposure;
    double H;
    if(background_){
      const arma::vec&lp=basis(node.path,true);
      H=0.0;

      if(family_==1) for(arma::uword i=0;i<weights_->n_elem;i++)
        H+=(*weights_)[i]*(lp[i]<-745.0?0.0:std::exp(lp[i]));
      else for(arma::uword i=0;i<weights_->n_elem;i++) H+=(*weights_)[i]*lp[i];
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

      PPSTGeometryCache temporary(points,*region_,gate_,family_,nullptr,false);
      return temporary.intensity(tree,lambda,points);
    }
    arma::vec out(points.n_rows,arma::fill::zeros);

    for(const auto&item:lambda){
      const arma::vec&lp=basis(tree.at(item.first).path,role==1);
      if(role==1&&family_!=1) for(arma::uword i=0;i<points.n_rows;i++) out[i]+=item.second*lp[i];
      else for(arma::uword i=0;i<points.n_rows;i++)
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
