#ifndef _USE_Armadillo
#define _USE_Armadillo
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#endif


#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#endif

#ifndef _USE_Utils
#define _USE_Utils
#include "utils.h"
#endif

#ifndef _USE_PPT
#define _USE_PPT

#include "tree_limits.h"


class PPT {
public:
    std::vector<TreeNode*> nodes; // Array-based tree

    // (in-class initialisers so a default-constructed / freshly copied PPT never
    //  has indeterminate tuning parameters)
    int dim = 0;
    int max_depth = 0;
    int max_nodes = 0;

    // global hyper parameters
    double a = 0.5;
    double b = 0.0;

    // tree tuning parameters
    int min_leaf_n = 1;
    int cut_grid_n = 50;
    double max_aspect_ratio = std::numeric_limits<double>::infinity();
    double rho = 0.5;   // base split prob (alpha) at depth 0
    double lam = 1.0;   // axis selection prob (normalised per node)
    double eta = 2.0;   // depth penalty: P(split | depth d) = rho * (1+d)^{-eta}

    double loglik = 0.0;
    double total_intensity = 0.0;


public:
    // Constructors
    PPT() {}

    // Deep copy constructor from vector<TreeNode*>
    PPT(const std::vector<TreeNode*>& nodes_) {
        nodes.reserve(nodes_.size());
        for (const auto* node : nodes_) {
            if (node) nodes.push_back(node->clone());
            else      nodes.push_back(nullptr);
        }
    }
    // Copy constructor (deep copy)
    PPT(const PPT& other) {
        copy_params(other);   // FIX: scalar tuning params were NOT copied before,
        nodes.reserve(other.nodes.size());  // so resampled/deep-copied particles ran
        for (const auto* node : other.nodes) {   // the proposal with garbage rho /
            if (node) nodes.push_back(node->clone());  // min_leaf_n / cut_grid_n.
            else      nodes.push_back(nullptr);
        }
    }
    // Assignment operator (deep copy)
    PPT& operator=(const PPT& other) {
        if (this != &other) {
            clear();
            copy_params(other);          // FIX: carry the scalar params too
            nodes.reserve(other.nodes.size());
            for (const auto* node : other.nodes) {
                if (node) nodes.push_back(node->clone());
                else      nodes.push_back(nullptr);
            }
        }
        return *this;
    }

    // Copy all scalar hyper-parameters / tuning settings (not the nodes).
    void copy_params(const PPT& other) {
        dim = other.dim; max_depth = other.max_depth; max_nodes = other.max_nodes;
        a = other.a; b = other.b;
        min_leaf_n = other.min_leaf_n; cut_grid_n = other.cut_grid_n;
        max_aspect_ratio = other.max_aspect_ratio;
        rho = other.rho; lam = other.lam; eta = other.eta;
        loglik = other.loglik; total_intensity = other.total_intensity;
    }
    // Destructor
    ~PPT() { clear(); }

    // Clear all nodes and free memory
    void clear() {
        for (auto* node : nodes) {
            delete node;
        }
        nodes.clear();
    }

    // Resize (frees old nodes)
    void resize(size_t n) {
        clear();
        nodes.resize(n, nullptr);
    }

    // Accessors
    size_t size() const { return nodes.size(); }
    TreeNode* operator[](size_t i) { return nodes[i]; }
    const TreeNode* operator[](size_t i) const { return nodes[i]; }
    std::vector<TreeNode*>& get_nodes() { return nodes; }
    const std::vector<TreeNode*>& get_nodes() const { return nodes; }

    // Replace all nodes (deep copy from another vector)
    void set_nodes(const std::vector<TreeNode*>& other_nodes) {
        clear();
        nodes.reserve(other_nodes.size());
        for (const auto* node : other_nodes) {
            if (node) nodes.push_back(node->clone());
            else      nodes.push_back(nullptr);
        }
    }

    // Deep copy method
    PPT deep_copy() const { return PPT(*this); }


    // member functions
    void init_tree_prior(double rho_, double lam_, int min_leaf_n_, int cut_grid_n_){
        this->min_leaf_n = min_leaf_n_;
        this->cut_grid_n = cut_grid_n_;
        this->rho = rho_;
        this->lam = lam_;
    }

    void init_intensity_par(double a_, double b_){
        this->a = a_;
        this->b = b_;     
    }

    void initialize(const arma::mat& region_root, int n, int max_depth_) {

        // Regions are stored as a d x 2 matrix: rows are dimensions and the
        // two columns are lower/upper bounds.
        this->dim = region_root.n_rows;
        this->max_depth = ppt_checked_depth(max_depth_);
        this->a = 0.5;
        this->b = 0;
        this->min_leaf_n = 1;
        this->cut_grid_n = 50;
        this->max_aspect_ratio = std::numeric_limits<double>::infinity();
        this->rho = 0.5;
        this->eta = 2.0;   // depth penalty on the split prior (Chipman et al. 1998)
        this->lam = 1.0/this->dim;
        // Calculate total number of nodes (array-based binary tree)
        max_nodes = ppt_tree_slots(max_depth);
        // Resize and initialize with nullptr
        nodes.assign(max_nodes, nullptr);

        // Root indices (0-based indexing)
        arma::uvec root_idx = arma::regspace<arma::uvec>(0, n-1);

        // Create the root node (is_empty=false, is_leaf=true)
        TreeNode* root = new TreeNode(region_root, root_idx, 0, false, true);
        nodes[0] = root;

    }
    void initialize(const arma::mat& region_root, int n, int max_depth_, int min_leaf_n_,
        double a_=0.5, double b_=0.0,
        double rho_=0.5, int cut_grid_n_=50, double eta_=2.0,
        double max_aspect_ratio_=std::numeric_limits<double>::infinity()) {

        // Regions are stored as a d x 2 matrix: rows are dimensions and the
        // two columns are lower/upper bounds.
        this->dim = region_root.n_rows;
        this->max_depth = ppt_checked_depth(max_depth_);
        this->a = a_;
        this->b = b_;
        this->min_leaf_n = min_leaf_n_;
        this->cut_grid_n = cut_grid_n_;
        this->max_aspect_ratio = max_aspect_ratio_;
        this->rho = rho_;
        this->eta = eta_;   // depth penalty on the split prior
        this->lam = 1.0/this->dim;
        // Calculate total number of nodes (array-based binary tree)
        max_nodes = ppt_tree_slots(max_depth);
        // Resize and initialize with nullptr
        nodes.assign(max_nodes, nullptr);

        // Root indices (0-based indexing)
        arma::uvec root_idx = arma::regspace<arma::uvec>(0, n-1);

        // Create the root node (is_empty=false, is_leaf=true)
        TreeNode* root = new TreeNode(region_root, root_idx, 0, false, true);
        nodes[0] = root;
        // double area = arma::prod(region_root.col(1) - region_root.col(0));
        // double beta = b + area; 
        // nodes[0]->lambda =  R::rgamma(a*area + n, 1.0 / beta);


    }

    /**************************************************************/
    // generic helper
    void split_node(int i, const arma::mat& pts, int axis, double cut, int min_leaf_n);
    std::vector<std::vector<double> > find_valid_cuts(const arma::mat& x, 
      const arma::mat& region, bool force_mid_cut=false);

    Rcpp::List to_R_list();
    double get_area(int id);
    int get_counts(int id);

    /**************************************************************/
    // PP tree routines
    void PPT_one_step_ahead(double& log_inc, int i, const arma::mat& pts);
    // Apply the action stored at node i of a complete reference tree while
    // evaluating the same importance increment as PPT_one_step_ahead().  This
    // is the conditioning operation required by a valid conditional SMC
    // kernel: the reference action is never redrawn.
    void PPT_force_reference_step(double& log_inc, int i,
                                  const PPT& ref_tree,
                                  const arma::mat& pts);
    void PPT_draw_lambda();
    arma::vec predict_lambda(const arma::mat& XX); 
    double PPT_get_lppd(const arma::vec& new_lambda);

    double PPT_base_mloglik(int n, double area, double a=.5, double b=0.0){
        // For b > 0 this is the normalized Ga(a,b) marginal
        //   b^a Gamma(n+a) / {Gamma(a) (b+area)^(n+a)}.
        // The b=0 branch deliberately retains the package's historical
        // improper-prior limit.  In particular a=0.5 reproduces exactly
        // Gamma(n+0.5) / area^(n+0.5).
        if (b > 0.0) {
            return lgamma(n + a) - lgamma(a) + a * log(b) -
                   (n + a) * log(b + area);
        }
        return lgamma(n + a) - (n + a) * log(area);
    }
     

    double PPT_log_transition_prob(const TreeNode* parent_node,
                                   const TreeNode* ref_node,
                                   const arma::mat& pts);

    // Conditional Poisson-process log likelihood for the sampled leaf rates.
    // Collapsed Gamma-Poisson scores are used only in tree-update calculations.
    void get_TreeLoglik(){
        double total = 0.0;
        for (const auto* node : nodes) {
            if (node && node->is_leaf && !node->is_empty) {
                int n = node->idx.n_elem;        // number of points in this leaf
                double area = arma::prod(node->region.col(1) - node->region.col(0));
                if (n > 0) total += n * std::log(node->lambda);
                total -= node->lambda * area;
            }
        }
        this->loglik = total;
        return;
    }
    /**************************************************************/

};



#endif

















