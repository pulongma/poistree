#ifndef POISTREE_PPT_SOFT_CUTS_CACHE_H
#define POISTREE_PPT_SOFT_CUTS_CACHE_H

#include <string>
#include <unordered_set>

// Candidate cuts depend on hard routing and the fixed data/fit controls, not
// on latent labels or soft gates. Keep each axis lazy so a support query still
// stops at the first nonempty axis, exactly as the direct calculation does.
struct PPSTCutsEntry {
  std::vector<std::vector<double> > axis;
  std::vector<unsigned char> ready;
  int any=-1;                       // -1 unknown, 0 no supported axis, 1 supported
  explicit PPSTCutsEntry(arma::uword d):axis(d),ready(d,0){}
};

class PPSTCutsCache {
  const arma::mat& points_;
  const int nmin_,mode_,ncand_;
  std::unordered_map<std::string,PPSTCutsEntry> entries_;

public:
  static std::string key(const PPSTNode&node){
    std::string out;
    const size_t fields=2*sizeof(int)+2*sizeof(double);
    out.reserve(sizeof(int)+node.path.size()*fields+
                node.box.n_elem*sizeof(double));
    out.append(reinterpret_cast<const char*>(&node.depth),sizeof(node.depth));
    for(const PPSTGate&g:node.path){
      out.append(reinterpret_cast<const char*>(&g.axis),sizeof(g.axis));
      out.append(reinterpret_cast<const char*>(&g.cut),sizeof(g.cut));
      out.append(reinterpret_cast<const char*>(&g.parent_width),sizeof(g.parent_width));
      out.append(reinterpret_cast<const char*>(&g.side),sizeof(g.side));
    }
    out.append(reinterpret_cast<const char*>(node.box.memptr()),
               node.box.n_elem*sizeof(double));
    return out;
  }
private:
  PPSTCutsEntry& entry(const PPSTNode&node){
    return entries_.try_emplace(key(node),points_.n_cols).first->second;
  }
  void fill_axis(PPSTCutsEntry&value,const PPSTNode&node,int axis){
    if(value.ready[axis]) return;
    value.axis[axis]=ppst_axis_cuts(points_,node.idx,node.box,axis,
                                  mode_,ncand_,nmin_);
    value.ready[axis]=1;
    if(!value.axis[axis].empty()) value.any=1;
  }
public:
  PPSTCutsCache(const arma::mat&points,int nmin,int mode,int ncand)
    :points_(points),nmin_(nmin),mode_(mode),ncand_(ncand){}
  const std::vector<double>& axis(const PPSTNode&node,int which){
    PPSTCutsEntry&value=entry(node);
    fill_axis(value,node,which);
    return value.axis[which];
  }
  bool any(const PPSTNode&node){
    PPSTCutsEntry&value=entry(node);
    if(value.any>=0) return value.any!=0;
    for(arma::uword j=0;j<points_.n_cols;j++){
      fill_axis(value,node,j);
      if(value.any==1) return true;
    }
    value.any=0;
    return false;
  }
  // Retain internal nodes too: pruning or changing a split needs their cuts.
  // Each attempt adds only proposed children; rejected/removed paths expire
  // immediately, so memory is bounded by the current tree, not chain length.
  void trim(const PPSTree&tree){
    std::unordered_set<std::string> active;
    active.reserve(tree.size());
    for(const auto&kv:tree) active.insert(key(kv.second));
    for(auto it=entries_.begin();it!=entries_.end();){
      if(active.find(it->first)==active.end()) it=entries_.erase(it);
      else ++it;
    }
  }
  size_t size() const { return entries_.size(); }
};

static std::vector<double> ppst_cached_axis_cuts(const PPSTNode&node,
    const arma::mat&points,int axis,int mode,int ncand,int nmin,
    PPSTCutsCache*cache){
  if(cache) return cache->axis(node,axis);
  return ppst_axis_cuts(points,node.idx,node.box,axis,mode,ncand,nmin);
}

#endif
