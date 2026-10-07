#include "SMCtree.h"

static SoftModel hard_probe_model(int exact_max = 0) {
  SoftModel M;
  M.X = {{0.05, 0.86}, {0.16, 0.71}, {0.28, 0.62},
         {0.54, 0.40}, {0.77, 0.24}, {0.94, 0.11}};
  M.region = {{0.0, 1.0}, {0.0, 1.0}};
  M.grid = {{0.2, 0.5, 0.8}, {0.3, 0.6}};
  M.a = 0.5; M.b = 0.15; M.rho = 0.9; M.eta = 1.0;
  M.Dmax = 2; M.n_nodes = ppt_tree_slots(M.Dmax) + 1;
  M.exact_max = exact_max; M.defensive = 0.2; M.hard_proposal = true;
  M.init_tables();
  return M;
}

// Inspect a node; path columns are one-based axis, cut, and side (-1/+1).

// [[Rcpp::export]]
Rcpp::List SPPT_hard_node_probe(const arma::mat& x, const arma::mat& region,
    Rcpp::List grids, const arma::vec& gate, const arma::mat& path,
    double temperature = 0.5, double defensive = 0.1, int selected = -1) {
  const int d = x.n_cols, n = x.n_rows;
  if (n < 1 || d < 1 || region.n_rows != (arma::uword)d || region.n_cols != 2 ||
      gate.n_elem != (arma::uword)d || grids.size() != d || path.n_cols != 3 || path.n_rows > 8 ||
      !x.is_finite() || !region.is_finite() || !gate.is_finite() || !path.is_finite())
    Rcpp::stop("invalid diagnostic dimensions or nonfinite input");
  if (!(temperature > 0 && temperature <= 1) || !(defensive > 0 && defensive < 1))
    Rcpp::stop("invalid surrogate controls");
  SoftModel M;
  M.X = x; M.region = region; M.a = 0.5; M.b = 0.15; M.rho = 0.9; M.eta = 1;
  M.Dmax = path.n_rows + 1; M.n_nodes = ppt_tree_slots(M.Dmax) + 1;
  M.exact_max = 0; M.defensive = 0; M.hard_proposal = true;
  M.proposal_temperature = temperature; M.proposal_defensive = defensive;
  M.grid.resize(d);
  for (int j = 0; j < d; ++j) {
    if (!(region(j, 1) > region(j, 0)) || gate[j] <= 0) Rcpp::stop("invalid region or gate");
    Rcpp::NumericVector cuts = grids[j];
    if (cuts.size() == 0) Rcpp::stop("empty grid");
    for (double c : cuts) {
      if (!(c > region(j, 0) && c < region(j, 1))) Rcpp::stop("invalid cut");
      M.grid[j].push_back(c);
    }
  }
  M.init_tables();
  SoftSMCtree smc(M, 2, true);
  smc.cur_gate = gate; smc.build_gate_table(gate);
  int node = smc.make_root();
  std::vector<SoftGateStep> steps;
  for (arma::uword r = 0; r < path.n_rows; ++r) {
    const int j = (int)path(r, 0) - 1, side = path(r, 2);
    if (path(r, 0) != j + 1 || j < 0 || j >= d || path(r, 2) != side || (side != -1 && side != 1))
      Rcpp::stop("invalid path axis or side");
    if (!(path(r, 1) > region(j, 0) && path(r, 1) < region(j, 1))) Rcpp::stop("invalid path cut");
    steps.push_back(SoftGateStep{j, path(r, 1), side});
    SoftPathNode child;
    child.parent = node; child.depth = r + 1;
    child.heap = 2 * smc.store[node].heap + (side > 0);
    child.gate = steps.back(); child.pts = smc.store[node].pts;
    child.axisH = smc.store[node].axisH;
    child.axisH[j] = smc.cached_axis_log_integral(gate, steps, j);
    child.logH = std::accumulate(child.axisH.begin(), child.axisH.end(), 0.0);
    smc.store.push_back(std::move(child)); node = smc.store.size() - 1;
  }
  smc.expand_node(node);
  const size_t before = smc.G.evaluations;
  if (selected >= 0) {
    if (selected >= smc.G.n_cand) Rcpp::stop("selected cut index out of range");
    smc.ensure_true_exposure(node, selected);
    for (int i : smc.store[node].pts) {
      smc.G.left(i, selected); smc.G.right(i, selected); smc.G.r(i, selected);
    }
  }
  const auto& A = smc.store[node];
  return Rcpp::List::create(Rcpp::_["count"] = A.proxy_count,
    Rcpp::_["proxy_logHL"] = A.proxy_logHL, Rcpp::_["proxy_logHR"] = A.proxy_logHR,
    Rcpp::_["logHL"] = A.logHL, Rcpp::_["logHR"] = A.logHR,
    Rcpp::_["logH"] = A.logH, Rcpp::_["logQA"] = A.logQA,
    Rcpp::_["score"] = A.score, Rcpp::_["logprior"] = A.logprior, Rcpp::_["logq"] = A.logq,
    Rcpp::_["axis"] = smc.G.cand_axis, Rcpp::_["cuts"] = smc.G.cand_cut,
    Rcpp::_["gate_evaluations_before"] = (double)before,
    Rcpp::_["gate_evaluations_after"] = (double)smc.G.evaluations,
    Rcpp::_["true_exposure_builds"] = (double)smc.true_exposure_builds,
    Rcpp::_["dense_cells"] = (double)smc.G.loglft.size(),
    Rcpp::_["bin_builds"] = (double)smc.hard_bin_builds);
}

// [[Rcpp::export]]
Rcpp::List SPPT_hard_cache_probe(const arma::mat& gates,
    Rcpp::NumericVector temperature, Rcpp::NumericVector defensive,
    Rcpp::LogicalVector hard, bool cached) {
  const int n = gates.n_rows;
  if (gates.n_cols != 2 || n < 1 || !gates.is_finite() || arma::any(arma::vectorise(gates) <= 0) ||
      temperature.size() != n || defensive.size() != n || hard.size() != n)
    Rcpp::stop("invalid cache diagnostic inputs");
  SoftSMCtree smc(hard_probe_model(), 2, true);
  smc.cache_root = cached;
  Rcpp::List output(n);
  Rcpp::NumericVector builds(n), hits(n), bins(n), true_builds(n);
  for (int k = 0; k < n; ++k) {
    if (!(temperature[k] > 0 && temperature[k] <= 1) || !(defensive[k] > 0 && defensive[k] < 1) ||
        Rcpp::LogicalVector::is_na(hard[k])) Rcpp::stop("invalid controls");
    smc.M.hard_proposal = hard[k]; smc.M.proposal_temperature = temperature[k];
    smc.M.proposal_defensive = defensive[k]; smc.cur_gate = gates.row(k).t();
    smc.build_gate_table(smc.cur_gate); smc.store.clear();
    const int root = smc.make_root(); smc.expand_node(root);
    const auto& A = smc.store[root];
    output[k] = Rcpp::List::create(Rcpp::_["score"] = A.score, Rcpp::_["logq"] = A.logq,
      Rcpp::_["count"] = A.proxy_count, Rcpp::_["proxy_logHL"] = A.proxy_logHL,
      Rcpp::_["proxy_logHR"] = A.proxy_logHR, Rcpp::_["logHL"] = A.logHL,
      Rcpp::_["logHR"] = A.logHR);
    builds[k] = smc.root_axis_builds; hits[k] = smc.root_axis_hits;
    bins[k] = smc.hard_bin_builds; true_builds[k] = smc.true_exposure_builds;
  }
  return Rcpp::List::create(Rcpp::_["output"] = output, Rcpp::_["builds"] = builds,
    Rcpp::_["hits"] = hits, Rcpp::_["bins"] = bins, Rcpp::_["true_builds"] = true_builds);
}

// [[Rcpp::export]]
Rcpp::List SPPT_hard_gate_probe(const arma::mat& gates, Rcpp::List requests,
    bool lazy = true, int capacity = 65536) {
  if (gates.n_cols != 2 || gates.n_rows != (arma::uword)requests.size() ||
      !gates.is_finite() || arma::any(arma::vectorise(gates) <= 0) || capacity < 0)
    Rcpp::stop("invalid gate diagnostic inputs");
  SoftSMCtree smc(hard_probe_model(), 2, true);
  smc.lazy_gates = lazy; smc.G.max_axis_entries = capacity;
  Rcpp::List values(requests.size());
  Rcpp::NumericVector evaluations(requests.size()), entries(requests.size());
  for (int k = 0; k < requests.size(); ++k) {
    smc.build_gate_table(gates.row(k).t());
    Rcpp::IntegerMatrix query = requests[k];
    if (query.ncol() != 2) Rcpp::stop("requests need point/candidate columns");
    Rcpp::NumericMatrix out(query.nrow(), 3);
    for (int r = 0; r < query.nrow(); ++r) {
      const int i = query(r, 0) - 1, c = query(r, 1) - 1;
      if (i < 0 || i >= smc.M.n() || c < 0 || c >= smc.G.n_cand) Rcpp::stop("invalid query");
      out(r, 0) = smc.G.left(i, c); out(r, 1) = smc.G.right(i, c); out(r, 2) = smc.G.r(i, c);
    }
    values[k] = out; evaluations[k] = smc.G.evaluations;
    for (const auto& axis : smc.G.sparse) entries[k] += axis.size();
  }
  return Rcpp::List::create(Rcpp::_["values"] = values,
    Rcpp::_["evaluations"] = evaluations, Rcpp::_["entries"] = entries);
}

// [[Rcpp::export]]
Rcpp::List SPPT_hard_gate_chain(bool hard, int exact_max, bool update_gate,
    bool shared, bool lazy) {
  SoftModel M = hard_probe_model(exact_max); M.hard_proposal = hard;
  SoftSMCtree smc(M, 4, true); smc.lazy_gates = lazy;
  arma::vec gate = {4.0, shared ? 4.0 : 7.0};
  arma::vec a = {3.0, 3.0}, b = {0.5, 0.5}, sd = {0.2, 0.2}, zero(2, arma::fill::zeros);
  const arma::mat grid = {{0.25, 0.25}, {0.75, 0.75}};
  return smc.PGAS(grid, arma::mat(0, 2), gate, a, b, sd, zero,
    shared, 60, 10, 1, 1, update_gate, false);
}

// [[Rcpp::export]]
Rcpp::List SPPT_forced_hard_probe(Rcpp::IntegerVector bits,
    const arma::vec& gate, double temperature = 0.5, double defensive = 0.1) {
  if (bits.size() != 6 || gate.n_elem != 2 || !gate.is_finite() || arma::any(gate <= 0) ||
      !(temperature > 0 && temperature <= 1) || !(defensive > 0 && defensive < 1))
    Rcpp::stop("invalid forced-proposal inputs");
  SoftModel M = hard_probe_model();
  M.proposal_temperature = temperature; M.proposal_defensive = defensive;
  SoftSMCtree smc(M, 2, true); smc.cur_gate = gate; smc.build_gate_table(gate);
  const int root = smc.make_root(); smc.expand_node(root);
  smc.particles.assign(2, SoftParticleS());
  smc.particles[0].rec[1] = SoftRecord{root, -1, -1};
  smc.particles[1].rec[1] = SoftRecord{root, 0, -1};
  SoftRef ref; ref.node.resize(M.n_nodes); ref.z.resize(6); ref.valid = true;

  ref.node[1].S = 1; ref.node[1].J = 0; ref.node[1].L = 0.2;
  for (int i = 0; i < 6; ++i) {
    if (bits[i] != 0 && bits[i] != 1) Rcpp::stop("invalid bit");
    ref.z[i] = bits[i] ? 2 : 3;
  }
  smc.sample_position(1, &ref);
  const auto& A = smc.store[root];
  return Rcpp::List::create(Rcpp::_["logw"] = smc.particles[0].logw,
    Rcpp::_["logq"] = A.logq, Rcpp::_["logprior"] = A.logprior,
    Rcpp::_["logHL"] = A.logHL[0], Rcpp::_["logHR"] = A.logHR[0],
    Rcpp::_["proxy_logHL"] = A.proxy_logHL[0], Rcpp::_["proxy_logHR"] = A.proxy_logHR[0],
    Rcpp::_["logQA"] = A.logQA, Rcpp::_["x"] = M.X,
    Rcpp::_["gate_evaluations"] = (double)smc.G.evaluations,
    Rcpp::_["true_exposure_builds"] = (double)smc.true_exposure_builds);
}
