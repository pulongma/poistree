## Standalone exact-target verification for the production iRJ-MCMC kernel.
##
## This is intentionally not a routine testthat test: it compiles the native
## source and runs four long chains.  From the package root, run
##
##   Rscript tests/verification/verify-irjmcmc-depth1.R
##
## Set POISTREE_IRJ_VERIFY_ITER to change the default 100,000 iterations.

suppressPackageStartupMessages(library(Rcpp))

script_argument <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_argument[grepl("^--file=", script_argument)])
script_dir <- if (length(script_file)) {
  dirname(normalizePath(script_file[[1L]]))
} else {
  normalizePath("tests/verification")
}
package_root <- normalizePath(file.path(script_dir, "..", ".."))

cache_dir <- tempfile("poistree-irj-verification-")
dir.create(cache_dir)
Rcpp::sourceCpp(
  file.path(package_root, "src", "mppt_models.cpp"),
  cacheDir = cache_dir,
  rebuild = TRUE,
  verbose = FALSE
)

x <- c(0.12, 0.31, 0.47, 0.60, 0.88)
X <- matrix(x, ncol = 1L)
monitor <- matrix(c(0.10, 0.50, 0.90), ncol = 1L)
region <- matrix(c(0, 1), nrow = 1L)

n <- length(x)
kappa <- 1.2
a_xi <- 2
b_xi <- 1
gate <- 8
alpha <- 0.6
eta <- 1
max_depth <- 1L
min_leaf_n <- 1L
cut_candidates <- 3L

quantile_type1 <- function(z, probability) {
  h <- (length(z) - 1) * probability + 1
  z[max(1L, floor(h))]
}

sorted_x <- sort(x)
cuts <- unique(sort(vapply(
  0.05 + 0.90 * (0:(cut_candidates - 1L)) / (cut_candidates - 1L),
  function(probability) quantile_type1(sorted_x, probability),
  numeric(1)
)))
cuts <- cuts[
  cuts > 1e-3 & cuts < 1 - 1e-3 &
    vapply(cuts, function(cut) {
      sum(x < cut) >= min_leaf_n && sum(x >= cut) >= min_leaf_n
    }, logical(1))
]

collapsed_factor <- function(m, exposure, xi) {
  exp(
    kappa * log(kappa / xi) - lgamma(kappa) + lgamma(kappa + m) -
      (kappa + m) * log(exposure + kappa / xi)
  )
}

## For a node with count m and exposure A, integrate out its independent
## xi ~ Gamma(a_xi, b_xi).  The first vector is
##
##   q(m,A) = integral g(m,A;xi) pi(xi) dxi,
##
## and the second is E(lambda | m,A), averaged over xi under that posterior.
node_integrals <- function(exposure) {
  q <- vapply(0:n, function(m) {
    integrate(
      function(xi) {
        collapsed_factor(m, exposure, xi) *
          dgamma(xi, shape = a_xi, rate = b_xi)
      },
      0, Inf, rel.tol = 1e-9
    )$value
  }, numeric(1))
  mean_lambda <- vapply(0:n, function(m) {
    numerator <- integrate(
      function(xi) {
        (kappa + m) / (kappa / xi + exposure) *
          collapsed_factor(m, exposure, xi) *
          dgamma(xi, shape = a_xi, rate = b_xi)
      },
      0, Inf, rel.tol = 1e-9
    )$value
    numerator / q[[m + 1L]]
  }, numeric(1))
  list(q = q, mean_lambda = mean_lambda)
}

exact_target <- function(soft = FALSE) {
  root <- node_integrals(1)
  root_weight <- (1 - alpha) * root$q[[n + 1L]]
  root_intensity <- rep(root$mean_lambda[[n + 1L]], nrow(monitor))

  split_weights <- numeric(length(cuts))
  split_intensity <- matrix(0, nrow = length(cuts), ncol = nrow(monitor))

  for (k in seq_along(cuts)) {
    cut <- cuts[[k]]
    if (soft) {
      phi_left_x <- plogis(-gate * (x - cut))
      phi_right_x <- 1 - phi_left_x
      phi_left_monitor <- plogis(-gate * (monitor[, 1L] - cut))
      phi_right_monitor <- 1 - phi_left_monitor
      left_exposure <- integrate(
        function(z) plogis(-gate * (z - cut)), 0, 1,
        rel.tol = 1e-11
      )$value
      right_exposure <- 1 - left_exposure
      configurations <- as.matrix(expand.grid(rep(list(0:2), n)))
    } else {
      left_exposure <- cut
      right_exposure <- 1 - cut
      configurations <- as.matrix(expand.grid(rep(list(0:1), n)))
    }

    left <- node_integrals(left_exposure)
    right <- node_integrals(right_exposure)
    configuration_weight <- numeric(nrow(configurations))
    configuration_intensity <- matrix(
      0, nrow = nrow(configurations), ncol = nrow(monitor)
    )

    for (r in seq_len(nrow(configurations))) {
      labels <- configurations[r, ]
      m_parent <- sum(labels == 0L)
      if (soft) {
        m_left <- sum(labels == 1L)
        m_right <- sum(labels == 2L)
        spatial_factor <- prod(ifelse(
          labels == 0L, 1,
          ifelse(labels == 1L, phi_left_x, phi_right_x)
        ))
        configuration_intensity[r, ] <-
          root$mean_lambda[[m_parent + 1L]] +
          left$mean_lambda[[m_left + 1L]] * phi_left_monitor +
          right$mean_lambda[[m_right + 1L]] * phi_right_monitor
      } else {
        m_left <- sum(labels == 1L & x < cut)
        m_right <- sum(labels == 1L & x >= cut)
        spatial_factor <- 1
        configuration_intensity[r, ] <-
          root$mean_lambda[[m_parent + 1L]] +
          ifelse(
            monitor[, 1L] < cut,
            left$mean_lambda[[m_left + 1L]],
            right$mean_lambda[[m_right + 1L]]
          )
      }
      configuration_weight[[r]] <-
        root$q[[m_parent + 1L]] *
        left$q[[m_left + 1L]] *
        right$q[[m_right + 1L]] * spatial_factor
    }

    split_weights[[k]] <- alpha / length(cuts) * sum(configuration_weight)
    normalized <- configuration_weight / sum(configuration_weight)
    split_intensity[k, ] <- colSums(normalized * configuration_intensity)
  }

  normalizer <- root_weight + sum(split_weights)
  list(
    split_probability = sum(split_weights) / normalizer,
    mean_intensity = (
      root_weight * root_intensity +
        colSums(split_weights * split_intensity)
    ) / normalizer
  )
}

iterations <- as.integer(Sys.getenv("POISTREE_IRJ_VERIFY_ITER", "100000"))
burn <- min(10000L, floor(iterations / 5L))
if (!is.finite(iterations) || iterations <= burn + 1000L) {
  stop("POISTREE_IRJ_VERIFY_ITER must exceed burn-in by at least 1,000.")
}

hyperparameters <- c(
  kappa, 6, a_xi, b_xi, 0.6,
  1, 1, 0.15, 1, 0.1, 0.15,
  gate, 4, 4 / 30, 0.12, 0.5
)

run_kernel <- function(model, seed) {
  set.seed(seed)
  mppstree_multi(
    X, monitor, matrix(numeric(), 0L, 1L), region,
    hyperparameters,
    update_hyper = 0L, update_gate = 0L, model = model,
    al = alpha, eta = eta,
    Dmax = max_depth, nmin = min_leaf_n,
    iters = iterations, burn = burn, thin = 1L,
    nmove = 2L, ncc = 2L, cut_mode = 1L,
    ncand = cut_candidates, chains = 1L, verbose = 0L,
    proposal_mode = 3L
  )
}

cases <- data.frame(
  model = c(0L, 7L),
  seed = c(103L, 203L),
  model_name = c("MPPT", "S-MPPT"),
  kernel = "combined"
)

exact <- list(hard = exact_target(FALSE), soft = exact_target(TRUE))
results <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
  specification <- cases[i, ]
  fit <- run_kernel(
    specification$model,
    specification$seed
  )
  target <- if (specification$model == 0L) exact$hard else exact$soft
  estimated_split <- fit$mean_leaves - 1
  relative_intensity_error <- max(
    abs(as.numeric(fit$mean) - target$mean_intensity) /
      target$mean_intensity
  )
  data.frame(
    model = specification$model_name,
    kernel = specification$kernel,
    seed = specification$seed,
    exact_split = target$split_probability,
    estimated_split = estimated_split,
    split_error = abs(estimated_split - target$split_probability),
    maximum_relative_intensity_error = relative_intensity_error
  )
}))

print(results, row.names = FALSE, digits = 7)

if (any(results$split_error > 0.02) ||
    any(results$maximum_relative_intensity_error > 0.035)) {
  stop("Production iRJ-MCMC failed the depth-one exact-target verification.")
}

cat("Production iRJ-MCMC passed the depth-one exact-target verification.\n")
