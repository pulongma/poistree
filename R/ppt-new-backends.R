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
    gate_structure = c("shared", "dimension"),
    update_gate = TRUE,
    alpha = 0.95, eta = 2,
    max_depth = 8L, min_leaf_n = 1L,
    chains = 4L, iter = 10000L, burn = 2500L, thin = 3L,
    tree_moves = 3L, change_moves = 8L,
    cut_proposal = c("quantile", "uniform", "data"),
    cut_candidates = 30L, seed = 1L, verbose = TRUE) {
  gate_family <- match.arg(gate_family)
  gate_structure <- match.arg(gate_structure)
  cut_proposal <- match.arg(cut_proposal)
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
    as.integer(chains), as.integer(verbose)
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

  structure(
    list(
      call = NULL,
      model = list(
        gating = "soft", gate_family = gate_family,
        scales = "leaf", scale_prior = NA_character_,
        sampler = "rjmcmc", label = "S-PPT"
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
        kappa = NA_real_, tau = NA_real_,
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
        unique_trees = NA_integer_, leaf_count_trace = numeric(),
        max_depth_trace = numeric(), log_evidence_increment = numeric(),
        log_evidence_running = numeric()
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
}

.ppt_fit_multiscale_rjmcmc <- function(
    x, region, gating, scale_prior,
    sampler_name = c("rjmcmc", "irjmcmc"),
    predict_at = x, test = NULL,
    kappa = 1, tau = 6,
    a_xi = 2, b_xi = NULL, sd_xi = 0.4,
    update_hyper = TRUE, update_gate = TRUE,
    a_kappa = 1, b_kappa = 1, sd_kappa = 0.15,
    a_tau = 1, b_tau = 0.1, sd_tau = 0.15,
    gate = 30, a_gate = 4, b_gate = 4 / 30,
    sd_gate = 0.12, gate_depth = 0.5,
    gate_family = c("logistic", "compact"),
    gate_structure = c("shared", "dimension"),
    alpha = NULL, eta = NULL,
    max_depth = 12L, min_leaf_n = 1L,
    chains = 8L, iter = 4000L, burn = 1500L, thin = 3L,
    tree_moves = 3L, change_moves = 6L,
    cut_proposal = c("quantile", "uniform", "data"),
    cut_candidates = 30L, seed = 1L, verbose = TRUE) {
  sampler_name <- match.arg(sampler_name)
  gate_family <- match.arg(gate_family)
  gate_structure <- match.arg(gate_structure)
  cut_proposal <- match.arg(cut_proposal)
  cut_mode <- switch(cut_proposal,
                     data = 0L, quantile = 1L, uniform = 2L)

  is_soft <- identical(gating, "soft")
  is_independent <- identical(scale_prior, "independent")
  if (identical(sampler_name, "irjmcmc")) {
    if (!is_independent) {
      stop("iRJ-MCMC currently requires `scale_prior = \"independent\"`.",
           call. = FALSE)
    }
    ## iRJ-MCMC always combines informed cuts with sequential conditional
    ## allocation proposals. Keeping one production kernel avoids an otherwise
    ## unnecessary tuning choice and makes fitted objects directly comparable.
    proposal_mode <- 3L
  } else {
    proposal_mode <- 0L
  }
  model_code <- if (is_soft && is_independent) {
    if (identical(gate_family, "logistic")) 7L else 5L
  } else if (is_soft) {
    if (identical(gate_family, "logistic")) 6L else 4L
  } else if (is_independent) {
    0L
  } else {
    1L
  }
  if (is.null(alpha)) alpha <- if (is_soft) 0.99 else 0.95
  if (is.null(eta)) eta <- if (is_soft) 1 else 2
  ## Retain explicit-NULL compatibility while using the same permissive
  ## minimum child occupancy for every backend.
  if (is.null(min_leaf_n)) min_leaf_n <- 1L

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

  if (is.null(b_xi)) b_xi <- a_xi / nrow(x)
  positive_scalars <- c(
    kappa = kappa, tau = tau, a_xi = a_xi, b_xi = b_xi,
    sd_xi = sd_xi, a_kappa = a_kappa, b_kappa = b_kappa,
    sd_kappa = sd_kappa, a_tau = a_tau, b_tau = b_tau,
    sd_tau = sd_tau
  )
  if (any(!is.finite(positive_scalars)) || any(positive_scalars <= 0)) {
    stop("Scale and hyperprior parameters must be positive finite scalars.",
         call. = FALSE)
  }
  if (length(gate_depth) != 1L || !is.finite(gate_depth) ||
      gate_depth < 0) {
    stop("`gate_depth` must be a nonnegative finite scalar.",
         call. = FALSE)
  }
  if (identical(gate_structure, "shared") &&
      any(lengths(list(gate, a_gate, b_gate, sd_gate)) != 1L)) {
    stop("Shared gates require scalar gate-prior parameters.",
         call. = FALSE)
  }
  if (identical(gate_structure, "dimension")) {
    gate <- .ppt_expand_parameter(gate, d, "gate")
    a_gate <- .ppt_expand_parameter(a_gate, d, "a_gate")
    b_gate <- .ppt_expand_parameter(b_gate, d, "b_gate")
    sd_gate <- .ppt_expand_parameter(sd_gate, d, "sd_gate")
  } else {
    gate <- as.numeric(gate)
    a_gate <- as.numeric(a_gate)
    b_gate <- as.numeric(b_gate)
    sd_gate <- as.numeric(sd_gate)
    if (any(!is.finite(c(gate, a_gate, b_gate, sd_gate))) ||
        any(c(gate, a_gate, b_gate, sd_gate) <= 0)) {
      stop("Gate and gate-prior parameters must be positive.",
           call. = FALSE)
    }
  }
  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1 ||
      length(eta) != 1L || !is.finite(eta) || eta < 0) {
    stop("Require `0 < alpha < 1` and `eta >= 0`.", call. = FALSE)
  }
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
  valid_logical_scalar <- function(z) {
    is.logical(z) && length(z) == 1L && !is.na(z)
  }
  if (!valid_logical_scalar(update_hyper) ||
      !valid_logical_scalar(update_gate)) {
    stop("`update_hyper` and `update_gate` must be logical scalars.",
         call. = FALSE)
  }

  hyperparameters <- c(
    kappa, tau, a_xi, b_xi, sd_xi,
    a_kappa, b_kappa, sd_kappa,
    a_tau, b_tau, sd_tau,
    gate, a_gate, b_gate, sd_gate, gate_depth
  )

  set.seed(seed)
  raw <- mppstree_multi(
    x, predict_at, test, region,
    hyperparameters, as.integer(update_hyper), as.integer(update_gate),
    model_code, alpha, eta,
    as.integer(max_depth), as.integer(min_leaf_n),
    as.integer(iter), as.integer(burn), as.integer(thin),
    as.integer(tree_moves), as.integer(change_moves),
    cut_mode, as.integer(cut_candidates),
    as.integer(chains), as.integer(verbose), proposal_mode
  )

  input_names <- colnames(x)
  if (is.null(input_names)) input_names <- paste0("x", seq_len(d))
  gate_mean <- stats::setNames(as.numeric(raw$gate_mean), input_names)
  gate_accept <- stats::setNames(
    as.numeric(raw$gate_accept), input_names
  )
  if (!is_soft) {
    gate_mean[] <- NA_real_
    gate_accept[] <- NA_real_
  }
  finite_gate_accept <- gate_accept[is.finite(gate_accept)]
  if (length(finite_gate_accept)) {
    names(finite_gate_accept) <- paste0(
      "gate_", names(finite_gate_accept)
    )
  }
  tree_accept <- stats::setNames(
    as.numeric(raw$tree_accept), c("grow", "prune", "change")
  )
  integrated_intensity <- as.numeric(raw$integrated_intensity)
  prediction_draws <- as.matrix(raw$draws)
  retained_per_chain <- .ppt_retained_per_chain(iter, burn, thin)

  structure(
    list(
      call = NULL,
      model = list(
        gating = gating,
        gate_family = if (is_soft) gate_family else NA_character_,
        scales = "multiscale", scale_prior = scale_prior,
        sampler = sampler_name,
        algorithm = if (identical(sampler_name, "irjmcmc")) {
          "iRJ-MCMC"
        } else {
          "RJ-MCMC"
        },
        label = .ppt_model_label(gating, "multiscale", scale_prior)
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
        kappa = as.numeric(raw$kappa_mean),
        tau = as.numeric(raw$tau_mean),
        mean_gate = if (is_soft) mean(gate_mean) else NA_real_,
        gate_by_dimension = gate_mean,
        mean_log_likelihood = as.numeric(raw$loglik_mean),
        mean_integrated_intensity = mean(integrated_intensity),
        integrated_intensity_draws = integrated_intensity,
        lppd = if (nrow(test)) as.numeric(raw$logpred) else NA_real_,
        log_evidence = NA_real_,
        draws = ncol(prediction_draws),
        particle_weights = numeric(),
        tree_draws = list(),
        state = list(
          mode = "heap",
          nodes = raw$state_nodes,
          gate = as.matrix(raw$state_gate),
          gate_mode = if (!is_soft) 0L
                      else if (identical(gate_family, "logistic")) 2L
                      else 1L,
          gate_depth = if (is_soft && identical(gate_family, "compact")) {
            gate_depth
          } else {
            0
          },
          chain = rep(seq_len(chains), each = retained_per_chain)
        )
      ),
      diagnostics = list(
        acceptance = c(tree_accept, finite_gate_accept),
        tree_acceptance = tree_accept,
        gate_acceptance = gate_accept,
        particle_ess = NA_real_, ess_history = numeric(),
        unique_trees = NA_integer_, leaf_count_trace = numeric(),
        max_depth_trace = numeric(), log_evidence_increment = numeric(),
        log_evidence_running = numeric()
      ),
      prior = list(
        increments = list(kappa = kappa),
        scales = list(
          type = scale_prior,
          tau = if (is_independent) NA_real_ else tau,
          shape = a_xi, rate = b_xi
        ),
        gate = if (is_soft) list(
          shape = a_gate, rate = b_gate,
          structure = gate_structure, family = gate_family,
          depth_exponent = if (identical(gate_family, "compact")) {
            gate_depth
          } else {
            NA_real_
          }
        ) else NULL,
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
        informed_cut = identical(sampler_name, "irjmcmc"),
        conditional_labels = identical(sampler_name, "irjmcmc"),
        gate_family = if (is_soft) gate_family else NA_character_,
        update_hyper = isTRUE(update_hyper),
        update_gate = if (is_soft) isTRUE(update_gate) else FALSE,
        seed = seed
      ),
      backend = "mppstree_multi"
    ),
    class = "ppt"
  )
}

.ppt_fit_hard_multiscale_markov_rjmcmc <- function(x, region, ...) {
  .ppt_fit_multiscale_rjmcmc(
    x, region, gating = "hard", scale_prior = "markov", ...
  )
}

.ppt_fit_hard_multiscale_independent_rjmcmc <- function(x, region, ...) {
  .ppt_fit_multiscale_rjmcmc(
    x, region, gating = "hard", scale_prior = "independent", ...
  )
}

.ppt_fit_soft_multiscale_markov_rjmcmc <- function(x, region, ...) {
  .ppt_fit_multiscale_rjmcmc(
    x, region, gating = "soft", scale_prior = "markov", ...
  )
}

.ppt_fit_soft_multiscale_independent_rjmcmc <- function(x, region, ...) {
  .ppt_fit_multiscale_rjmcmc(
    x, region, gating = "soft", scale_prior = "independent", ...
  )
}

.ppt_fit_hard_multiscale_independent_irjmcmc <- function(x, region, ...) {
  .ppt_fit_multiscale_rjmcmc(
    x, region, gating = "hard", scale_prior = "independent",
    sampler_name = "irjmcmc", ...
  )
}

.ppt_fit_soft_multiscale_independent_irjmcmc <- function(x, region, ...) {
  .ppt_fit_multiscale_rjmcmc(
    x, region, gating = "soft", scale_prior = "independent",
    sampler_name = "irjmcmc", ...
  )
}
