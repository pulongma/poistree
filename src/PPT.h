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
    std::vector<TreeNode*> nodes;

    int dim = 0;
    int max_depth = 0;
    int max_nodes = 0;

    double a = 0.5;
    double b = 0.0;

    int min_leaf_n = 1;
    int cut_grid_n = 50;
    double max_aspect_ratio = std::numeric_limits<double>::infinity();
    double rho = 0.5;
    double lam = 1.0;
    double eta = 2.0;

    double loglik = 0.0;
    double total_intensity = 0.0;

public:

    PPT() {}

    // Deep-copy a vector of tree nodes.
    PPT(const std::vector<TreeNode*>& nodes_) {
        nodes.reserve(nodes_.size());
        for (const auto* node : nodes_) {
            if (node) nodes.push_back(node->clone());
            else      nodes.push_back(nullptr);
        }
    }
    // Deep-copy a tree and its tuning settings.
    PPT(const PPT& other) {
        copy_params(other);
        nodes.reserve(other.nodes.size());
        for (const auto* node : other.nodes) {
            if (node) nodes.push_back(node->clone());
            else      nodes.push_back(nullptr);
        }
    }
    // Replace this tree with a deep copy.
    PPT& operator=(const PPT& other) {
        if (this != &other) {
            clear();
            copy_params(other);
            nodes.reserve(other.nodes.size());
            for (const auto* node : other.nodes) {
                if (node) nodes.push_back(node->clone());
                else      nodes.push_back(nullptr);
            }
        }
        return *this;
    }

    // Copy scalar hyperparameters and tuning settings.
    void copy_params(const PPT& other) {
        dim = other.dim; max_depth = other.max_depth; max_nodes = other.max_nodes;
        a = other.a; b = other.b;
        min_leaf_n = other.min_leaf_n; cut_grid_n = other.cut_grid_n;
        max_aspect_ratio = other.max_aspect_ratio;
        rho = other.rho; lam = other.lam; eta = other.eta;
        loglik = other.loglik; total_intensity = other.total_intensity;
    }

    ~PPT() { clear(); }

    // Delete all nodes and release their storage.
    void clear() {
        for (auto* node : nodes) {
            delete node;
        }
        nodes.clear();
    }

    // Resize node storage, deleting existing nodes.
    void resize(size_t n) {
        clear();
        nodes.resize(n, nullptr);
    }

    size_t size() const { return nodes.size(); }
    TreeNode* operator[](size_t i) { return nodes[i]; }
    const TreeNode* operator[](size_t i) const { return nodes[i]; }
    std::vector<TreeNode*>& get_nodes() { return nodes; }
    const std::vector<TreeNode*>& get_nodes() const { return nodes; }

    // Replace all nodes with deep copies.
    void set_nodes(const std::vector<TreeNode*>& other_nodes) {
        clear();
        nodes.reserve(other_nodes.size());
        for (const auto* node : other_nodes) {
            if (node) nodes.push_back(node->clone());
            else      nodes.push_back(nullptr);
        }
    }

    // Return a deep copy of the tree.
    PPT deep_copy() const { return PPT(*this); }

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

        this->dim = region_root.n_rows;
        this->max_depth = ppt_checked_depth(max_depth_);
        this->a = 0.5;
        this->b = 0;
        this->min_leaf_n = 1;
        this->cut_grid_n = 50;
        this->max_aspect_ratio = std::numeric_limits<double>::infinity();
        this->rho = 0.5;
        this->eta = 2.0;
        this->lam = 1.0/this->dim;

        max_nodes = ppt_tree_slots(max_depth);

        nodes.assign(max_nodes, nullptr);

        arma::uvec root_idx = arma::regspace<arma::uvec>(0, n-1);

        TreeNode* root = new TreeNode(region_root, root_idx, 0, false, true);
        nodes[0] = root;

    }
    void initialize(const arma::mat& region_root, int n, int max_depth_, int min_leaf_n_,
        double a_=0.5, double b_=0.0,
        double rho_=0.5, int cut_grid_n_=50, double eta_=2.0,
        double max_aspect_ratio_=std::numeric_limits<double>::infinity()) {

        this->dim = region_root.n_rows;
        this->max_depth = ppt_checked_depth(max_depth_);
        this->a = a_;
        this->b = b_;
        this->min_leaf_n = min_leaf_n_;
        this->cut_grid_n = cut_grid_n_;
        this->max_aspect_ratio = max_aspect_ratio_;
        this->rho = rho_;
        this->eta = eta_;
        this->lam = 1.0/this->dim;

        max_nodes = ppt_tree_slots(max_depth);

        nodes.assign(max_nodes, nullptr);

        arma::uvec root_idx = arma::regspace<arma::uvec>(0, n-1);

        TreeNode* root = new TreeNode(region_root, root_idx, 0, false, true);
        nodes[0] = root;

    }

    void split_node(int i, const arma::mat& pts, int axis, double cut, int min_leaf_n);
    std::vector<std::vector<double> > find_valid_cuts(const arma::mat& x,
      const arma::mat& region, bool force_mid_cut=false);

    Rcpp::List to_R_list();
    double get_area(int id);
    int get_counts(int id);

    void PPT_one_step_ahead(double& log_inc, int i, const arma::mat& pts);
    // Apply a reference action and evaluate its importance increment.

    void PPT_force_reference_step(double& log_inc, int i,
                                  const PPT& ref_tree,
                                  const arma::mat& pts);
    void PPT_draw_lambda();
    arma::vec predict_lambda(const arma::mat& XX);
    double PPT_get_lppd(const arma::vec& new_lambda);

    double PPT_base_mloglik(int n, double area, double a=.5, double b=0.0){

        if (b > 0.0) {
            return lgamma(n + a) - lgamma(a) + a * log(b) -
                   (n + a) * log(b + area);
        }
        return lgamma(n + a) - (n + a) * log(area);
    }

    double PPT_log_transition_prob(const TreeNode* parent_node,
                                   const TreeNode* ref_node,
                                   const arma::mat& pts);

    // Compute the point-process log likelihood for the sampled leaf rates.

    void get_TreeLoglik(){
        double total = 0.0;
        for (const auto* node : nodes) {
            if (node && node->is_leaf && !node->is_empty) {
                int n = node->idx.n_elem;
                double area = arma::prod(node->region.col(1) - node->region.col(0));
                if (n > 0) total += n * std::log(node->lambda);
                total -= node->lambda * area;
            }
        }
        this->loglik = total;
        return;
    }

};

#endif

