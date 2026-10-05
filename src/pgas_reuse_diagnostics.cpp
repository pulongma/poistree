// Internal regression entry points for the real PGAS cache and routing paths.
#include "SMCtree.h"

static SoftModel reuse_test_model(int exact_max, double defensive) {
  SoftModel M;
  M.X = {{0.05, 0.86}, {0.16, 0.71}, {0.28, 0.62},
         {0.54, 0.40}, {0.77, 0.24}, {0.94, 0.11}};
  M.region = {{0.0, 1.0}, {0.0, 1.0}};
  M.grid = {{0.2, 0.5, 0.8}, {0.3, 0.6}};
  M.a = 0.5; M.b = 0.15; M.rho = 0.9; M.eta = 1.0;
  M.Dmax = 2; M.n_nodes = ppt_tree_slots(M.Dmax) + 1;
  M.exact_max = exact_max; M.defensive = defensive;
  M.init_tables();
  return M;
}

// [[Rcpp::export]]
Rcpp::List SPPT_root_reuse_probe(const arma::mat& gates, int exact_max, bool cached) {
  if (gates.n_cols != 2 || gates.n_rows == 0 || !gates.is_finite() || arma::any(arma::vectorise(gates) <= 0))
    Rcpp::stop("gates must have two positive finite columns");
  SoftSMCtree smc(reuse_test_model(exact_max, 0.2), 2, true);
  smc.cache_root = cached;
  Rcpp::List outputs(gates.n_rows);
  Rcpp::NumericVector builds(gates.n_rows), hits(gates.n_rows), prefixes(gates.n_rows);
  size_t prefix_builds = 0;
  for (arma::uword step = 0; step < gates.n_rows; ++step) {
    smc.cur_gate = gates.row(step).t();
    smc.build_gate_table(smc.cur_gate);
    smc.store.clear();
    const int root = smc.make_root();
    smc.expand_node(root);
    SoftPathNode& A = smc.store[root];
    if (smc.M.n() <= exact_max)
      for (int j = 0; j < smc.M.d(); ++j)
        if (soft_exact_prepare_prefix(*A.exact_axis.at(j))) ++prefix_builds;
    outputs[step] = Rcpp::List::create(
      Rcpp::_["score"] = A.score, Rcpp::_["logq"] = A.logq,
      Rcpp::_["logHL"] = A.logHL, Rcpp::_["logHR"] = A.logHR,
      Rcpp::_["count_norm"] = A.count_norm);
    builds[step] = smc.root_axis_builds; hits[step] = smc.root_axis_hits;
    prefixes[step] = prefix_builds;
  }
  return Rcpp::List::create(Rcpp::_["output"] = outputs,
    Rcpp::_["builds"] = builds, Rcpp::_["hits"] = hits, Rcpp::_["prefixes"] = prefixes);
}

// [[Rcpp::export]]
Rcpp::List SPPT_root_reuse_chain(bool cached, bool update_gate, bool shared, int exact_max) {
  SoftSMCtree smc(reuse_test_model(exact_max, 0.2), 4, true);
  smc.cache_root = cached;
  arma::vec gate = {4.0, shared ? 4.0 : 7.0};
  arma::vec a = {3.0, 3.0}, b = {0.5, 0.5}, sd = {0.2, 0.2}, zero(2, arma::fill::zeros);
  const arma::mat grid = {{0.25, 0.25}, {0.75, 0.75}};
  return smc.PGAS(grid, arma::mat(0, 2), gate, a, b, sd, zero,
                  shared, 80, 10, 1, 1, update_gate, false);
}

// [[Rcpp::export]]
Rcpp::List SPPT_forced_exact_probe(Rcpp::IntegerVector bits, double defensive) {
  if (bits.size() != 6 || defensive < 0 || defensive >= 1)
    Rcpp::stop("require six bits and 0 <= defensive < 1");
  SoftSMCtree smc(reuse_test_model(150, defensive), 2, true);
  smc.cur_gate = arma::vec({4.0, 7.0});
  smc.build_gate_table(smc.cur_gate);
  const int root = smc.make_root();
  smc.expand_node(root);
  smc.particles.assign(2, SoftParticleS());
  smc.particles[0].rec[1] = SoftRecord{root, -1, -1};
  // Only the reference particle advances, isolating the no-replay path.
  smc.particles[1].rec[1] = SoftRecord{root, 0, -1};
  SoftRef ref;
  ref.node.resize(smc.M.n_nodes); ref.z.resize(6); ref.valid = true;
  ref.node[1].S = 1; ref.node[1].J = 0; ref.node[1].L = 0.5;
  for (int i = 0; i < 6; ++i) {
    if (bits[i] != 0 && bits[i] != 1) Rcpp::stop("bits must be 0 or 1");
    ref.z[i] = bits[i] ? 2 : 3;
  }
  const std::vector<double> score = smc.store[root].score;
  const std::vector<double> logq = smc.store[root].logq;
  const double logQA = smc.store[root].logQA;
  smc.sample_position(1, &ref);
  const auto& particle = smc.particles[0];
  return Rcpp::List::create(
    Rcpp::_["logw"] = particle.logw, Rcpp::_["score"] = score,
    Rcpp::_["logq"] = logq, Rcpp::_["logQA"] = logQA,
    Rcpp::_["left"] = smc.store[particle.rec.at(2).node].pts,
    Rcpp::_["right"] = smc.store[particle.rec.at(3).node].pts,
    Rcpp::_["prefix_builds"] = (double)smc.exact_prefix_builds,
    Rcpp::_["forced_routes"] = (double)smc.exact_forced_routes);
}
