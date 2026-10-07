#ifndef POISTREE_SOFT_LOGISTIC_H
#define POISTREE_SOFT_LOGISTIC_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Evaluate a stable logistic right-routing probability.
static inline double pst_logistic_right(double z){
  if(z>=0.0){
    double e=std::exp(-z);
    return 1.0/(1.0+e);
  }
  double e=std::exp(z);
  return e/(1.0+e);
}

static inline double pst_logistic_log_right(double z){
  return -std::max(-z,0.0)-std::log1p(std::exp(-std::abs(z)));
}

static inline double pst_logistic_path_product(
    double x,const std::vector<double>&cuts,const std::vector<int>&sides,
    double dom_lo,double dom_hi,double gate){
  double width=dom_hi-dom_lo,out=1.0;
  for(size_t k=0;k<cuts.size();k++){
    double right=pst_logistic_right(gate*(x-cuts[k])/width);
    out*=sides[k]<0 ? 1.0-right : right;
  }
  return out;
}

// Evaluate a product of logistic path gates stably.

static inline double pst_logistic_path_product_stable(
    double x,const std::vector<double>&cuts,const std::vector<int>&sides,
    double dom_lo,double dom_hi,double gate){
  double width=dom_hi-dom_lo,log_out=0.0;
  for(size_t k=0;k<cuts.size();k++){
    double z=gate*(x-cuts[k])/width;
    log_out+=sides[k]<0 ? pst_logistic_log_right(-z)
                        : pst_logistic_log_right(z);
  }
  if(log_out<std::log(std::numeric_limits<double>::min())) return 0.0;
  return std::exp(log_out);
}

// Integrate one interval using adaptive Simpson quadrature.

static double pst_logistic_adaptive_simpson(
    const std::vector<double>&cuts,const std::vector<int>&sides,
    double dom_lo,double dom_hi,double gate,double left,double right,
    double f_left,double f_mid,double f_right,double whole,
    double abs_tol,int depth){
  double mid=0.5*(left+right);
  double left_mid=0.5*(left+mid),right_mid=0.5*(mid+right);
  double f_left_mid=pst_logistic_path_product_stable(
    left_mid,cuts,sides,dom_lo,dom_hi,gate
  );
  double f_right_mid=pst_logistic_path_product_stable(
    right_mid,cuts,sides,dom_lo,dom_hi,gate
  );
  double left_part=(mid-left)*(f_left+4.0*f_left_mid+f_mid)/6.0;
  double right_part=(right-mid)*(f_mid+4.0*f_right_mid+f_right)/6.0;
  double refined=left_part+right_part;
  double error=refined-whole;
  double tolerance=15.0*(abs_tol+5e-13*std::abs(refined));
  if(depth<=0||std::abs(error)<=tolerance)
    return refined+error/15.0;
  return pst_logistic_adaptive_simpson(
      cuts,sides,dom_lo,dom_hi,gate,left,mid,
      f_left,f_left_mid,f_mid,left_part,0.5*abs_tol,depth-1
    )+
    pst_logistic_adaptive_simpson(
      cuts,sides,dom_lo,dom_hi,gate,mid,right,
      f_mid,f_right_mid,f_right,right_part,0.5*abs_tol,depth-1
    );
}

// Numerically integrate a logistic path as a fallback for analytic integration.

static double pst_logistic_path_axis_integral_numeric(
    const std::vector<double>&cuts,const std::vector<int>&sides,
    double dom_lo,double dom_hi,double gate){
  double width=dom_hi-dom_lo;
  std::vector<double>breaks;
  breaks.reserve(cuts.size()+2);
  breaks.push_back(dom_lo);
  for(size_t k=0;k<cuts.size();k++)
    if(cuts[k]>dom_lo&&cuts[k]<dom_hi) breaks.push_back(cuts[k]);
  breaks.push_back(dom_hi);
  std::sort(breaks.begin(),breaks.end());
  std::vector<double>::iterator last=std::unique(
    breaks.begin(),breaks.end(),
    [width](double x,double y){
      return std::abs(x-y)<=32.0*std::numeric_limits<double>::epsilon()*
        std::max(width,std::max(std::abs(x),std::abs(y)));
    }
  );
  breaks.erase(last,breaks.end());

  double total=0.0;
  const double total_abs_tol=
    8.0*std::numeric_limits<double>::epsilon()*width;
  for(size_t k=1;k<breaks.size();k++){
    double left=breaks[k-1],right=breaks[k];
    if(!(right>left)) continue;
    double mid=0.5*(left+right);
    double f_left=pst_logistic_path_product_stable(
      left,cuts,sides,dom_lo,dom_hi,gate
    );
    double f_mid=pst_logistic_path_product_stable(
      mid,cuts,sides,dom_lo,dom_hi,gate
    );
    double f_right=pst_logistic_path_product_stable(
      right,cuts,sides,dom_lo,dom_hi,gate
    );
    double whole=(right-left)*(f_left+4.0*f_mid+f_right)/6.0;
    double interval_tol=total_abs_tol*(right-left)/width;
    total+=pst_logistic_adaptive_simpson(
      cuts,sides,dom_lo,dom_hi,gate,left,right,
      f_left,f_mid,f_right,whole,interval_tol,24
    );
  }
  return std::max(total,std::numeric_limits<double>::min());
}

// Integrate a product of logistic gates sharing one coordinate slope.

static double pst_logistic_path_axis_integral(
    const std::vector<double>&cuts,const std::vector<int>&sides,
    double dom_lo,double dom_hi,double gate){
  size_t m=cuts.size();
  double width=dom_hi-dom_lo;
  if(m==0) return width;
  if(!(width>0.0)||!(gate>0.0)||!std::isfinite(gate))
    return std::numeric_limits<double>::quiet_NaN();
  if(gate<1e-4||gate>500.0)
    return pst_logistic_path_axis_integral_numeric(
      cuts,sides,dom_lo,dom_hi,gate
    );

  std::vector<long double>a(m);
  int nright=0;
  long double left_constant=1.0L;
  for(size_t k=0;k<m;k++){
    long double ck=(cuts[k]-dom_lo)/width;
    a[k]=std::exp((long double)gate*(ck-0.5L));
    if(sides[k]>0) nright++;
    else left_constant*=a[k];
  }

  for(size_t i=0;i<m;i++) for(size_t j=i+1;j<m;j++){
    long double scale=std::max(std::abs(a[i]),std::abs(a[j]));
    if(std::abs(a[i]-a[j])<=1e-10L*std::max(1.0L,scale))
      return pst_logistic_path_axis_integral_numeric(
        cuts,sides,dom_lo,dom_hi,gate
      );
  }

  long double u0=std::exp(-(long double)gate/2.0L);
  long double u1=std::exp((long double)gate/2.0L);
  long double sum=0.0L;
  long double sum_abs=0.0L;
  if(nright==0){
    sum=std::log(u1/u0);
    sum_abs=std::abs(sum);
  }

  for(size_t i=0;i<m;i++){
    long double denom=1.0L;
    for(size_t j=0;j<m;j++) if(j!=i) denom*=a[j]-a[i];
    if(nright==0) denom*=(-a[i]);
    long double numerator=left_constant;
    if(nright>0){
      long double p=1.0L;
      for(int r=1;r<nright;r++) p*=(-a[i]);
      numerator*=p;
    }
    long double coefficient=numerator/denom;
    long double log_ratio=std::log((u1+a[i])/(u0+a[i]));
    long double term=coefficient*log_ratio;
    if(!std::isfinite(term))
      return pst_logistic_path_axis_integral_numeric(
        cuts,sides,dom_lo,dom_hi,gate
      );
    sum+=term;
    sum_abs+=std::abs(term);
  }

  long double cancellation=sum_abs/
    std::max(std::abs(sum),std::numeric_limits<long double>::min());
  if(!std::isfinite(cancellation)||cancellation>1e5L)
    return pst_logistic_path_axis_integral_numeric(
      cuts,sides,dom_lo,dom_hi,gate
    );

  long double answer=(long double)width*sum/(long double)gate;
  double out=(double)answer;
  if(!std::isfinite(out)||out<=0.0||out>width*(1.0+1e-7))
    return pst_logistic_path_axis_integral_numeric(
      cuts,sides,dom_lo,dom_hi,gate
    );
  return std::min(out,width);
}

// Evaluate a path product with node-relative logistic slopes.

static inline double pst_logistic_node_path_product(
    double x,const std::vector<double>&cuts,const std::vector<int>&sides,
    const std::vector<double>&parent_widths,double gate){
  double log_out=0.0;
  for(size_t k=0;k<cuts.size();k++){
    double z=gate*(x-cuts[k])/parent_widths[k];
    log_out+=sides[k]<0 ? pst_logistic_log_right(-z)
                        : pst_logistic_log_right(z);
  }
  return std::exp(log_out);
}

static double pst_logistic_node_adaptive_simpson(
    const std::vector<double>&cuts,const std::vector<int>&sides,
    const std::vector<double>&parent_widths,double gate,
    double left,double right,double f_left,double f_mid,double f_right,
    double whole,double abs_tol,int depth){
  double mid=left+0.5*(right-left);
  double left_mid=left+0.5*(mid-left),right_mid=mid+0.5*(right-mid);
  if(left_mid<=left||right_mid>=right) return whole;
  double fl=pst_logistic_node_path_product(
    left_mid,cuts,sides,parent_widths,gate);
  double fr=pst_logistic_node_path_product(
    right_mid,cuts,sides,parent_widths,gate);
  double lower=(mid-left)*(f_left+4.0*fl+f_mid)/6.0;
  double upper=(right-mid)*(f_mid+4.0*fr+f_right)/6.0;
  double refined=lower+upper,error=refined-whole;
  if(depth<=0||std::abs(error)<=15.0*(abs_tol+5e-13*std::abs(refined)))
    return std::max(0.0,refined+error/15.0);
  return pst_logistic_node_adaptive_simpson(
      cuts,sides,parent_widths,gate,left,mid,f_left,fl,f_mid,lower,
      0.5*abs_tol,depth-1)+
    pst_logistic_node_adaptive_simpson(
      cuts,sides,parent_widths,gate,mid,right,f_mid,fr,f_right,upper,
      0.5*abs_tol,depth-1);
}

static inline double pst_logistic_node_path_axis_integral(
    const std::vector<double>&cuts,const std::vector<int>&sides,
    const std::vector<double>&parent_widths,
    double dom_lo,double dom_hi,double gate){
  const double width=dom_hi-dom_lo;
  if(cuts.empty()) return width;
  if(!(width>0.0)||!std::isfinite(width)||!(gate>0.0)||
     !std::isfinite(gate)||cuts.size()!=sides.size()||
     cuts.size()!=parent_widths.size())
    return std::numeric_limits<double>::quiet_NaN();
  bool common_width=true;
  for(size_t k=0;k<cuts.size();k++){
    if(!(parent_widths[k]>0.0)||!std::isfinite(parent_widths[k])||
       !std::isfinite(cuts[k])||(sides[k]!=-1&&sides[k]!=1))
      return std::numeric_limits<double>::quiet_NaN();
    common_width=common_width&&parent_widths[k]==parent_widths[0];
  }

  const double effective_gate=gate*(width/parent_widths[0]);
  if(common_width&&effective_gate>=1e-4&&effective_gate<=500.0)
    return pst_logistic_path_axis_integral(
      cuts,sides,dom_lo,dom_hi,effective_gate);

  std::vector<double>breaks;
  breaks.reserve(15*cuts.size()+2);
  breaks.push_back(dom_lo);
  breaks.push_back(dom_hi);
  const double offsets[]={-32,-16,-8,-4,-2,-1,0,1,2,4,8,16,32};
  for(size_t k=0;k<cuts.size();k++){
    const double bandwidth=parent_widths[k]/gate;
    for(double offset:offsets){
      double point=cuts[k]+offset*bandwidth;
      if(point>dom_lo&&point<dom_hi) breaks.push_back(point);
    }
  }
  std::sort(breaks.begin(),breaks.end());
  breaks.erase(std::unique(breaks.begin(),breaks.end()),breaks.end());
  double total=0.0;
  const double total_abs_tol=8.0*std::numeric_limits<double>::epsilon()*width;
  for(size_t k=1;k<breaks.size();k++){
    double left=breaks[k-1],right=breaks[k];
    if(!(right>left)) continue;
    double mid=left+0.5*(right-left);
    double fl=pst_logistic_node_path_product(left,cuts,sides,parent_widths,gate);
    double fm=pst_logistic_node_path_product(mid,cuts,sides,parent_widths,gate);
    double fr=pst_logistic_node_path_product(right,cuts,sides,parent_widths,gate);
    double whole=(right-left)*(fl+4.0*fm+fr)/6.0;
    total+=pst_logistic_node_adaptive_simpson(
      cuts,sides,parent_widths,gate,left,right,fl,fm,fr,whole,
      total_abs_tol*(right-left)/width,24);
  }
  return std::max(std::numeric_limits<double>::min(),std::min(total,width));
}

#endif
