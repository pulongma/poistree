#ifndef _USE_Armadillo
#define _USE_Armadillo
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#endif


#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <cmath>
#endif


// #ifndef _USE_VECTOR
// #define _USE_VECTOR
// #include <vector>
// #endif


using namespace Rcpp;

#include "PPT.h"



#ifndef _USE_SMCtree
#define _USE_SMCtree

class SMCtree: public PPT {
public: 
	// members
	int P; // number of particles
	int niter; // number of PG iterations
	double resample_thresh; // resampling rate
	double force_mid_cut; 

    std::vector< PPT > particles; // List of particles (each a tree)
    arma::vec weights;     // Particle weights over time
    arma::imat ancestors;  // Ancestor indices
    arma::vec ESS_hist;    // Effective sample size history
    double logZ_hat = 0.0; // SMC log marginal-evidence estimate, relative to the root model (Z/Q0(D))
    arma::vec logZ_inc;    // per-step increment Delta_t (sum = logZ_hat; aligned with ESS_hist)
    arma::vec logZ_run;    // running cumulative log-evidence estimate (cumsum of logZ_inc; last active = logZ_hat)


    // Constructor
    SMCtree(int P_) : P(P_) {
        particles.resize(P);
    }
    SMCtree(int P_, int niter_, double resample_thresh_) : P(P_), niter(niter_), 
    	resample_thresh(resample_thresh_)
    {
        particles.resize(P);
    }

    // Destructor: clean up all TreeNode* in all particles
    ~SMCtree() {clear();}

    // Clear all SMC objects and free memory
    void clear() {
	    for (auto& tree : particles) {
	        tree.clear();   // this calls PPT::clear(), which deletes all TreeNode*'s
	    }
	    particles.clear();    // clears the vector itself
	    weights.reset();      // set to empty
	    ancestors.reset();
	    ESS_hist.reset();
	}
   // member function declarations
    void init(const arma::mat& region_root, int n, int max_depth, int min_leaf_n, 
        double a=.5, double b=0.0,
        double rho=0.95, int cut_grid_n=30
    	){
	    for (int p = 0; p < P; ++p) {
	        particles[p].initialize(region_root, n, max_depth, 
	        	min_leaf_n, a, b, rho, cut_grid_n);
	    
        }
        resample_thresh = 0.5;
    }


    void init_AS(const PPT& ref_tree, const arma::mat& region_root, int n, 
    	int max_depth, int min_leaf_n, 
        double a=.5, double b=0.0,
        double rho=0.95, int cut_grid_n=30
    	) {
        int P = particles.size();
        for (int p = 0; p < P; ++p) {
            if (p == 0) {
                // First particle is a deep copy of the reference tree
                particles[p] = ref_tree.deep_copy();
            } else {
                particles[p].initialize(region_root, n, max_depth, min_leaf_n, 
	        	a, b, rho, cut_grid_n);
            }
        }

    }

    // void run_smc(const arma::mat& pts, int max_depth, int min_leaf_n,
    //              double resample_thresh, double rho, int cut_grid_n,
    //              double lam, bool verbose,
    //              double a, double b, double w0, double u0);
    // void resample(const arma::vec& weights);
    // PPT deep_copy(const PPT& particle); // for single tree
    SMCtree deep_copy() {  // for multiple tree 
        SMCtree copy(particles.size());
        for (size_t i = 0; i < particles.size(); ++i) {
            copy.particles[i] = particles[i].deep_copy();
        }
        return copy;
    }
    // main SMC functions
    int get_MAP_index() const {
    	int idmax = arma::index_max(weights);
        return idmax;
    }
    PPT get_trajectory(const arma::imat& ancestors,
                const std::vector< PPT >& particles);

    // PPT routines 
    void PPT_SMC(const arma::mat& pts,
                   int max_depth,
				   int min_leaf_n, 
				   double a, double b,
				   double rho, double lam, int cut_grid_n,
                   bool verbose = true);


    void PPT_cSMC(PPT*  ref_tree,
                           const arma::mat& pts,
                           int max_depth,
						   int min_leaf_n, 
						   double a, double b, 
						   double rho, double lam, int cut_grid_n,
                           bool verbose = true);
	Rcpp::List PPT_PGAS(const arma::mat& pts, const arma::mat& grid,
						   int niter, int max_depth, 
						   int min_leaf_n, 
						   double a, double b,
						   double rho, double lam, int cut_grid_n,                          
                           bool verbose=true);



    /*****************************************************************************/

};

#endif
