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

using namespace Rcpp;

#include "PPT.h"

#ifndef _USE_SMCtree
#define _USE_SMCtree

class SMCtree: public PPT {
public:

	int P;
	int niter;
	double resample_thresh = 0.5;
	double force_mid_cut = 0.0;

    std::vector< PPT > particles;
    arma::vec weights;
    arma::vec ESS_hist;
    double logZ_hat = 0.0;
    arma::vec logZ_inc;
    arma::vec logZ_run;

    SMCtree(int P_) : P(P_) {
        particles.resize(P);
    }
    SMCtree(int P_, int niter_, double resample_thresh_) : P(P_), niter(niter_),
    	resample_thresh(resample_thresh_)
    {
        particles.resize(P);
    }

    // Delete all particle trees.
    ~SMCtree() {clear();}

    // Clear particle trees, weights, and ESS history.
    void clear() {
	    for (auto& tree : particles) {
	        tree.clear();
	    }
	    particles.clear();
	    weights.reset();
	    ESS_hist.reset();
	}

   void init(const arma::mat& region_root, int n, int max_depth, int min_leaf_n,
        double a=.5, double b=0.0,
        double rho=0.95, int cut_grid_n=50,
        double max_aspect_ratio=std::numeric_limits<double>::infinity()
    	){
        ppt_check_dense_storage(max_depth, P);
	    for (int p = 0; p < P; ++p) {
	        particles[p].initialize(region_root, n, max_depth,
	        	min_leaf_n, a, b, rho, cut_grid_n, 2.0,
                max_aspect_ratio);

        }
    }

    void init_AS(const PPT& ref_tree, const arma::mat& region_root, int n,
        int max_depth, int min_leaf_n,
        double a=.5, double b=0.0,
        double rho=0.95, int cut_grid_n=50,
        double max_aspect_ratio=std::numeric_limits<double>::infinity()
    	) {
        ppt_check_dense_storage(max_depth, P);
        int P = particles.size();
        for (int p = 0; p < P; ++p) {
            if (p == 0) {

                particles[p] = ref_tree.deep_copy();
            } else {
                particles[p].initialize(region_root, n, max_depth, min_leaf_n,
	        	a, b, rho, cut_grid_n, 2.0, max_aspect_ratio);
            }
        }

    }

    SMCtree deep_copy() {
        SMCtree copy(particles.size());
        for (size_t i = 0; i < particles.size(); ++i) {
            copy.particles[i] = particles[i].deep_copy();
        }
        return copy;
    }

    int get_MAP_index() const {
    	int idmax = arma::index_max(weights);
        return idmax;
    }

    void PPT_SMC(const arma::mat& pts,
                   int max_depth,
				   int min_leaf_n,
				   double a, double b,
				   double rho, double lam, int cut_grid_n,
                   double max_aspect_ratio,
                   bool verbose = true);

    void PPT_cSMC(PPT*  ref_tree,
                           const arma::mat& pts,
                           int max_depth,
						   int min_leaf_n,
						   double a, double b,
						   double rho, double lam, int cut_grid_n,
                           double max_aspect_ratio,
                           bool verbose = true);
	Rcpp::List PPT_PGAS(const arma::mat& pts, const arma::mat& grid,
						   int niter, int max_depth,
						   int min_leaf_n,
						   double a, double b,
						   double rho, double lam, int cut_grid_n,
                           double max_aspect_ratio,
                           bool verbose=true);

};

#include "soft_smc.h"

class SoftSMCtree {
public:
    SoftModel M;
    int P;
    bool use_as;
    SoftGateTable G;

    using AxisPathKey = std::vector<std::pair<double, int>>;
    std::vector<std::map<AxisPathKey, double>> axis_integrals;
    arma::vec integral_gate, table_gate;
    double cached_axis_log_integral(const arma::vec& gate,
        const std::vector<SoftGateStep>& path, int j,
        double extra_cut = 0.0, int extra_side = 0);
    double cached_log_exposure(const arma::vec& gate,
        const std::vector<SoftGateStep>& path);
    arma::vec cur_gate;
    std::vector<SoftPathNode> store;
    std::vector<SoftParticleS> particles;
    arma::vec weights;
    double final_ess = 0.0;
    int as_moved = 0;
    int as_stamp = 0;
    int n_expanded = 0;
    int n_resampled = 0;
    bool resample_node = true;
    double ess_threshold = 1.0;
    bool alloc_rates = false;

    struct RootAxisCache {
        bool ready = false, exact = false;
        bool hard = false;
        double gate = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> axisL, axisR, logHL, logHR, logPsi, count_norm;
        std::vector<int> proxy_count;
        std::vector<double> proxy_logHL, proxy_logHR;
        std::shared_ptr<SoftExactAxis> exact_axis;
    };
    std::vector<RootAxisCache> root_axis_cache;
    bool cache_root = true;
    size_t root_axis_builds = 0, root_axis_hits = 0;
    size_t exact_prefix_builds = 0, exact_forced_routes = 0, exact_boundary_draws = 0;
    bool lazy_gates = false;
    std::vector<std::vector<double> > hard_cuts;
    std::vector<std::vector<int> > hard_order;
    std::vector<int> hard_bins;
    struct HardAxisExposure { std::vector<double> logL, logR; };
    std::vector<std::map<AxisPathKey, std::shared_ptr<HardAxisExposure> > > hard_exposures;
    arma::vec hard_exposure_gate;
    size_t hard_bin_builds = 0, hard_exposure_builds = 0, hard_exposure_hits = 0;
    size_t true_exposure_builds = 0;

    SoftSMCtree(const SoftModel& M_, int P_, bool use_as_) : M(M_), P(P_), use_as(use_as_) {
        ppt_check_soft_storage(M.Dmax);
        if (M.n_nodes != ppt_tree_slots(M.Dmax) + 1)
            Rcpp::stop("Invalid soft tree allocation size.");
        lazy_gates = M.hard_proposal;
    }

    void build_gate_table(const arma::vec& gate);
    std::vector<SoftGateStep> path_of(int v) const;
    int make_root();
    int make_child(int parent, int cand, int side, std::vector<int>&& pts);
    std::shared_ptr<SoftExactAxis> exact_cut_axis(SoftPathNode& node, int cand,
                                                double& delta);
    void expand_node(int v);
    void initialize_hard_bins();
    std::shared_ptr<HardAxisExposure> hard_axis_exposures(
        const std::vector<SoftGateStep>& path, int j);
    void expand_hard_node(int v);
    void ensure_true_exposure(int v, int cand);
    void expand_level(int level);
    bool sample_position(int t, SoftRef* ref);
    void resample(SoftRef* ref);
    void export_particle(const SoftParticleS& p, SoftRef& out);

    const SoftNodeP& ref_decision(SoftRef& ref, int h);
    int ref_bit(SoftRef& ref, int i, int h);
    double glued_logQ(SoftRef& ref, int h, const std::vector<int>& pts,
                      const std::vector<SoftGateStep>& path);
    double as_log_weight(const SoftParticleS& p, SoftRef& ref);

    void sweep(SoftRef* ref_in, const arma::vec& gate, SoftRef& ref_out);
    void label_sweep(SoftRef& ref, const arma::vec& gate);
    int gate_update(SoftRef& ref, arma::vec& gate, const arma::vec& a_gate,
                    const arma::vec& b_gate, const arma::vec& sd_gate,
                    const arma::vec& gate_min, bool shared);
    double log_target_gate(const SoftRef& ref, const arma::vec& gate,
                           const arma::vec& a_gate, const arma::vec& b_gate,
                           const arma::vec& gate_min, bool shared);
    void refresh_exposures(SoftRef& ref, const arma::vec& gate);
    Rcpp::List PGAS(const arma::mat& grid, const arma::mat& xtest,
                    arma::vec gate, const arma::vec& a_gate, const arma::vec& b_gate,
                    const arma::vec& sd_gate, const arma::vec& gate_min, bool gate_shared,
                    int niter, int burn, int thin, int label_sweeps, bool update_gate,
                    bool verbose);
};

#endif
