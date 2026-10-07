#ifndef POISTREE_TREE_LIMITS_H
#define POISTREE_TREE_LIMITS_H
#include <Rcpp.h>
#include <cmath>
#include <cstddef>
#include <limits>

inline int ppt_checked_depth(double depth, int minimum = 0) {
  if (!std::isfinite(depth) || depth != std::floor(depth) ||
      depth < minimum || depth > 20.0)
    Rcpp::stop("max_depth must be a finite scalar integer between %d and 20.", minimum);
  return static_cast<int>(depth);
}

inline int ppt_tree_slots(int depth) {
  ppt_checked_depth(depth);
  return static_cast<int>((std::size_t(1) << (depth + 1)) - 1);
}

inline int ppt_tree_steps(int depth) {
  ppt_checked_depth(depth);
  return static_cast<int>((std::size_t(1) << depth) - 1);
}

inline void ppt_check_dense_storage(int depth, double particles) {
  const int slots = ppt_tree_slots(depth);
  if (!std::isfinite(particles) || particles < 1 ||
      particles != std::floor(particles) ||
      particles > std::numeric_limits<int>::max())
    Rcpp::stop("particles must be a positive representable integer.");

  if ((2.0L * particles + 3.0L) * slots > 16777216.0L)
    Rcpp::stop("Dense tree storage limit exceeded; reduce max_depth or particles.");
}

inline void ppt_check_soft_storage(int depth) {
  ppt_checked_depth(depth, 1);

  if (4.0L * (static_cast<long double>(ppt_tree_slots(depth)) + 1.0L) > 4194304.0L)
    Rcpp::stop("Soft PGAS tree storage limit exceeded; reduce max_depth to at most 19.");
}
#endif
