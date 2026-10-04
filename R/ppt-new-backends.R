.ppt_retained_per_chain <- function(iter, burn, thin) {
  length(seq.int(as.integer(burn), as.integer(iter) - 1L,
                 by = as.integer(thin)))
}

.ppt_fit_soft_leaf_rjmcmc <- function(
    x, region, predict_at = x, test = NULL,
    a = 0.5, b = NULL,
    gate = 12, a_gate = 36, b_gate = 3,
    sd_gate = 0.07, gate_min = 0,
    gate_family = c("logistic", "compact"),
    gate_structure = c("dimension", "shared"),
    update_gate = TRUE,
    alpha = 0.95, eta = 2,
    max_depth = 8L, min_leaf_n = 1L,
    chains = 4L, iter = 10000L, burn = 2500L, thin = 3L,
    tree_moves = 3L, change_moves = 8L,
    cut_proposal = c("quantile", "uniform", "data"),
    cut_candidates = 50L, seed = 1L, verbose = TRUE) {
  .ppt_fit_soft_leaf_mcmc(
    x, region, predict_at, test, a, b, gate, a_gate, b_gate, sd_gate,
    gate_min, gate_family, gate_structure, update_gate, alpha, eta,
    max_depth, min_leaf_n, chains, iter, burn, thin, tree_moves,
    change_moves, cut_proposal, cut_candidates, seed, verbose,
    informed = FALSE
  )
}

.ppt_fit_soft_leaf_irjmcmc <- function(
    x, region, predict_at = x, test = NULL,
    a = 0.5, b = NULL,
    gate = 12, a_gate = 36, b_gate = 3,
    sd_gate = 0.07, gate_min = 0,
    gate_family = c("logistic", "compact"),
    gate_structure = c("dimension", "shared"),
    update_gate = TRUE,
    alpha = 0.95, eta = 2,
    max_depth = 8L, min_leaf_n = 1L,
    chains = 4L, iter = 10000L, burn = 2500L, thin = 3L,
    tree_moves = 3L, change_moves = 8L,
    cut_proposal = c("quantile", "uniform", "data"),
    cut_candidates = 50L, seed = 1L, verbose = TRUE) {
  .ppt_fit_soft_leaf_mcmc(
    x, region, predict_at, test, a, b, gate, a_gate, b_gate, sd_gate,
    gate_min, gate_family, gate_structure, update_gate, alpha, eta,
    max_depth, min_leaf_n, chains, iter, burn, thin, tree_moves,
    change_moves, cut_proposal, cut_candidates, seed, verbose,
    informed = TRUE
  )
}

.ppt_fit_soft_leaf_pcg <- function(
    x, region, predict_at = x, test = NULL,
    a = 0.5, b = NULL,
    gate = 12, a_gate = 36, b_gate = 3,
    sd_gate = 0.07, gate_min = 0,
    gate_family = c("logistic", "compact"),
    gate_structure = c("dimension", "shared"),
    update_gate = TRUE,
    alpha = 0.95, eta = 2,
    max_depth = 8L, min_leaf_n = 1L,
    chains = 4L, iter = 10000L, burn = 2500L, thin = 3L,
    tree_moves = 3L, change_moves = 8L,
    cut_proposal = c("quantile", "uniform", "data"),
    cut_candidates = 50L, seed = 1L, verbose = TRUE,
    ram_target = 0.234, ram_decay = 0.7, ram_adapt = NULL) {
  .ppt_fit_soft_leaf_mcmc(
    x, region, predict_at, test, a, b, gate, a_gate, b_gate, sd_gate,
    gate_min, gate_family, gate_structure, update_gate, alpha, eta,
    max_depth, min_leaf_n, chains, iter, burn, thin, tree_moves,
    change_moves, cut_proposal, cut_candidates, seed, verbose,
    informed = FALSE, pcg = TRUE, ram_target = ram_target,
    ram_decay = ram_decay, ram_adapt = ram_adapt
  )
}

.ppt_validate_ram_controls <- function(target, decay, adapt, burn) {
  scalar_finite <- function(x) {
    is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x)
  }
  if (!scalar_finite(target) || target <= 0 || target >= 1) {
    stop("`ram_target` must be a finite numeric scalar between 0 and 1.",
         call. = FALSE)
  }
  if (!scalar_finite(decay) || decay <= 0.5 || decay > 1) {
    stop("`ram_decay` must be a finite numeric scalar in (0.5, 1].",
         call. = FALSE)
  }
  if (is.null(adapt)) adapt <- floor(0.8 * burn)
  if (!scalar_finite(adapt) || adapt != floor(adapt) ||
      adapt < 0 || adapt > burn || adapt > .Machine$integer.max) {
    stop("`ram_adapt` must be an integer between 0 and `burn`.",
         call. = FALSE)
  }
  list(target = target, decay = decay, adapt = as.integer(adapt))
}

.ppt_fit_soft_leaf_mcmc <- function(
    x, region, predict_at, test, a, b, gate, a_gate, b_gate, sd_gate,
    gate_min, gate_family, gate_structure, update_gate, alpha, eta,
    max_depth, min_leaf_n, chains, iter, burn, thin, tree_moves,
    change_moves, cut_proposal, cut_candidates, seed, verbose, informed,
    pcg = FALSE, ram_target = 0.234, ram_decay = 0.7, ram_adapt = 0L) {
  gate_family <- match.arg(gate_family, c("logistic", "compact"))
  gate_structure <- match.arg(gate_structure, c("dimension", "shared"))
  cut_proposal <- match.arg(cut_proposal, c("quantile", "uniform", "data"))
  cut_mode <- switch(cut_proposal,
                     data = 0L, quantile = 1L, uniform = 2L)

  x <- .ppt_validate_points(x, name = "x")
  d <- ncol(x)
  region <- .ppt_validate_region(region, d)
  x <- .ppt_validate_points(x, d, region, "x")
  predict_at <- .ppt_validate_points(
    predict_at, d, region, "predict_at", allow_empty = TRUE
  )
  if (is.null(test)) {
    test <- matrix(numeric(), 0L, d)
  } else {
    test <- .ppt_validate_points(
      test, d, region, "test", allow_empty = TRUE
    )
  }

  if (is.null(b)) {
    b <- a * prod(region[, 2L] - region[, 1L]) / nrow(x)
  }
  if (length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    stop("The soft leaf model requires positive scalar `a` and `b`.",
         call. = FALSE)
  }
  if (identical(gate_structure, "shared") &&
      any(lengths(list(gate, a_gate, b_gate, sd_gate, gate_min)) != 1L)) {
    stop("Shared gates require scalar gate-prior parameters.",
         call. = FALSE)
  }
  gate <- .ppt_expand_parameter(gate, d, "gate")
  a_gate <- .ppt_expand_parameter(a_gate, d, "a_gate")
  b_gate <- .ppt_expand_parameter(b_gate, d, "b_gate")
  sd_gate <- .ppt_expand_parameter(sd_gate, d, "sd_gate")
  gate_min <- .ppt_expand_parameter(
    gate_min, d, "gate_min", allow_zero = TRUE
  )
  if (any(gate <= gate_min)) {
    stop("Every initial gate must exceed its `gate_min`.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1 ||
      length(eta) != 1L || !is.finite(eta) || eta < 0) {
    stop("Require `0 < alpha < 1` and `eta >= 0`.", call. = FALSE)
  }
  max_depth <- .ppt_validate_depth(max_depth)
  controls <- c(
    iter, burn, thin, chains, max_depth, min_leaf_n,
    tree_moves, change_moves, cut_candidates
  )
  if (any(!is.finite(controls)) || any(controls != floor(controls)) ||
      iter <= burn || burn < 0 || thin < 1 || chains < 1 ||
      max_depth < 0 || min_leaf_n < 1 || tree_moves < 0 ||
      change_moves < 0 || cut_candidates < 2) {
    stop("Invalid RJ-MCMC or tree controls.", call. = FALSE)
  }
  if (!is.logical(update_gate) || length(update_gate) != 1L ||
      is.na(update_gate)) {
    stop("`update_gate` must be a logical scalar.", call. = FALSE)
  }

  if (pcg && informed) {
    stop("The PCG sampler uses standard collapsed tree proposals.", call. = FALSE)
  }
  ram <- .ppt_validate_ram_controls(ram_target, ram_decay, ram_adapt, burn)

  set.seed(seed)
  raw <- ppstree_multi(
    x, predict_at, test, region,
    a, b, gate, a_gate, b_gate, sd_gate, gate_min,
    as.integer(identical(gate_structure, "shared")),
    alpha, eta, as.integer(max_depth), as.integer(min_leaf_n),
    as.integer(iter), as.integer(burn), as.integer(thin),
    as.integer(tree_moves), as.integer(change_moves),
    cut_mode, as.integer(cut_candidates), as.integer(update_gate),
    as.integer(match(gate_family, c("logistic", "compact")) - 1L),
    as.integer(chains), as.integer(verbose), informed = informed, pcg = pcg,
    ram_target = ram$target, ram_decay = ram$decay, ram_adapt = ram$adapt
  )

  input_names <- colnames(x)
  if (is.null(input_names)) input_names <- paste0("x", seq_len(d))
  gate_mean <- stats::setNames(as.numeric(raw$gate_mean), input_names)
  gate_accept <- stats::setNames(
    as.numeric(raw$gate_accept), input_names
  )
  tree_accept <- stats::setNames(
    as.numeric(raw$accept), c("grow", "prune", "change")
  )
  integrated_intensity <- as.numeric(raw$integrated_intensity)
  prediction_draws <- as.matrix(raw$draws)
  retained_per_chain <- .ppt_retained_per_chain(iter, burn, thin)
  leaf_count_trace <- vapply(raw$state_nodes, function(nodes) {
    sum(nodes[, 2L] < 0)
  }, numeric(1))
  max_depth_trace <- vapply(raw$state_nodes, function(nodes) {
    max(floor(log2(nodes[, 1L])))
  }, numeric(1))

  fit <- structure(
    list(
      call = NULL,
      model = list(
        gating = "soft", gate_family = gate_family,
        scales = "leaf",
        sampler = if (pcg) "pcg" else if (informed) "irjmcmc" else "rjmcmc",
        algorithm = if (pcg) "Partially collapsed Gibbs (RAM)" else
          if (informed) "Informed MH" else "RJ-MCMC",
        label = "S-PPT"
      ),
      data = list(
        x = x, n = nrow(x), dimension = d, region = region,
        test = if (nrow(test)) test else NULL
      ),
      prediction = list(
        locations = predict_at,
        mean = as.numeric(raw$mean),
        median = as.numeric(raw$median),
        lower = as.numeric(raw$lower95),
        upper = as.numeric(raw$upper95),
        level = 0.95,
        draws = prediction_draws
      ),
      posterior = list(
        mean_leaves = as.numeric(raw$mean_leaves),
        mean_max_depth = as.numeric(raw$mean_max_depth),
        mean_gate = if (identical(gate_structure, "shared")) {
          unname(gate_mean[1L])
        } else {
          mean(gate_mean)
        },
        gate_by_dimension = gate_mean,
        mean_log_likelihood = as.numeric(raw$loglik_mean),
        mean_integrated_intensity = mean(integrated_intensity),
        integrated_intensity_draws = integrated_intensity,
        lppd = if (nrow(test)) as.numeric(raw$logpred) else NA_real_,
        log_evidence = NA_real_,
        log_target_normalizer = NA_real_,
        log_relative_normalizer = NA_real_,
        draws = ncol(prediction_draws),
        particle_weights = numeric(),
        tree_draws = list(),
        state = list(
          mode = "heap",
          nodes = raw$state_nodes,
          gate = as.matrix(raw$state_gate),
          gate_mode = if (identical(gate_family, "logistic")) 2L else 1L,
          gate_depth = 0,
          chain = rep(seq_len(chains), each = retained_per_chain)
        )
      ),
      diagnostics = list(
        acceptance = c(
          tree_accept,
          stats::setNames(unname(gate_accept), paste0("gate_", input_names))
        ),
        tree_acceptance = tree_accept,
        gate_acceptance = gate_accept,
        particle_ess = NA_real_, ess_history = numeric(),
        unique_trees = NA_integer_, leaf_count_trace = leaf_count_trace,
        max_depth_trace = max_depth_trace, log_evidence_increment = numeric(),
        log_evidence_running = numeric(),
        expanded_nodes = NA_real_,
        resampling_events = NA_real_
      ),
      prior = list(
        intensity = list(shape = a, rate = b),
        gate = list(
          shape = a_gate, rate = b_gate, lower = gate_min,
          structure = gate_structure, family = gate_family
        ),
        tree = list(alpha = alpha, eta = eta)
      ),
      control = list(
        max_depth = as.integer(max_depth),
        min_leaf_n = as.integer(min_leaf_n),
        chains = as.integer(chains), iter = as.integer(iter),
        burn = as.integer(burn), thin = as.integer(thin),
        retained_per_chain = retained_per_chain,
        tree_moves = as.integer(tree_moves),
        change_moves = as.integer(change_moves),
        cut_proposal = cut_proposal,
        cut_candidates = as.integer(cut_candidates),
        gate_family = gate_family, update_gate = isTRUE(update_gate),
        seed = seed
      ),
      backend = "ppstree_multi"
    ),
    class = "ppt"
  )
  if (pcg) {
    fit$control$sd_gate <- sd_gate
    fit$control$ram_target <- ram$target
    fit$control$ram_decay <- ram$decay
    fit$control$ram_adapt <- ram$adapt
    fit$diagnostics$gate_joint_acceptance <- as.numeric(raw$gate_joint_accept)
    fit$diagnostics$chain_gate_joint_acceptance <-
      as.numeric(raw$chain_gate_joint_accept)
    fit$diagnostics$chain_gate_joint_acceptance_probability <-
      as.numeric(raw$chain_gate_joint_accept_prob)
    fit$diagnostics$ram <- list(
      covariance = raw$ram_covariance, factor = raw$ram_factor,
      updates = as.integer(raw$ram_updates),
      failures = as.integer(raw$ram_failures)
    )
  }
  fit
}
