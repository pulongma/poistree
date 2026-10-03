#ifndef POISTREE_SOFT_LOGISTIC_H
#define POISTREE_SOFT_LOGISTIC_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Stable logistic right-routing probability.
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

// Stable evaluation used by the quadrature fallback.  Accumulating on the
// log scale avoids losing a small, but still representable, path product when
// a long path contains several gates close to zero.
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

// Adaptive Simpson integration on one interval.  The caller splits at every
// gate centre first, so steep logistic transitions occur at interval
// boundaries rather than being hidden inside an initially coarse panel.
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

// Deterministic fallback for ill-conditioned partial fractions, nearly
// coincident poles, or extreme gate values.  This is deliberately adaptive:
// the former fixed 256-panel rule was accurate for moderate gates but could
// miss narrow transitions for large gates and deep paths.
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

// Closed-form integral for a product of complementary logistic gates sharing
// one slope on a coordinate.  With t=(x-dom_lo)/(dom_hi-dom_lo),
// u=exp(gate*(t-1/2)), and a_k=exp(gate*(c_k-1/2)), the integrand times dx
// is a proper rational function of u.  Distinct tree cuts give simple poles,
// so partial fractions reduce the exposure to a finite sum of logarithms.
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

  // Nearly repeated poles make the simple-pole representation unstable.
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

  // A finite partial-fraction sum can still be grossly inaccurate when large
  // signed terms cancel.  Long-double arithmetic gives about 18 decimal
  // digits on common platforms; use quadrature before cancellation can erase
  // more than roughly five decimal digits.  This conservative threshold also
  // covers platforms where long double has the same precision as double.
  // This check is what protects deep paths
  // with many splits on the same coordinate at small or moderate gates.
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

#endif
