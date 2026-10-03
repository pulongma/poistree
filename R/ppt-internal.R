.ppt_backend_key <- function(gating, scales, sampler) {
  paste(gating, scales, sampler, sep = ":")
}

# Backend names are stored as strings so the registry can be created before
# the implementation functions are loaded.
.ppt_backend_registry <- c(
  "hard:leaf:smc" = ".ppt_fit_hard_leaf_smc",
  "hard:leaf:rjmcmc" = ".ppt_fit_hard_leaf_rjmcmc",
  "hard:leaf:pgas" = ".ppt_fit_hard_leaf_pgas",
  "soft:leaf:rjmcmc" = ".ppt_fit_soft_leaf_rjmcmc"
)

.ppt_model_components <- function(model) {
  out <- unlist(model[c("gating", "scales", "sampler")], use.names = TRUE)
  out[!is.na(out)]
}

.ppt_model_label <- function(gating, scales) {
  if (identical(gating, "soft")) "S-PPT" else "PPT"
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.tree_sig <- function(nodes) {
  leaves <- Filter(
    function(node) !is.null(node) && isTRUE(node$is_leaf),
    nodes
  )
  if (!length(leaves)) return("")
  boxes <- vapply(leaves, function(node) {
    paste(round(as.numeric(node$region), 4), collapse = ",")
  }, character(1))
  paste(sort(boxes), collapse = ";")
}

.n_unique_trees <- function(particles) {
  length(unique(vapply(particles, .tree_sig, character(1))))
}

.leaves_per_particle <- function(particles) {
  vapply(particles, function(nodes) {
    sum(vapply(
      nodes,
      function(node) !is.null(node) && isTRUE(node$is_leaf),
      logical(1)
    ))
  }, integer(1))
}

.ppt_validate_region <- function(region, d) {
  region <- as.matrix(region)
  storage.mode(region) <- "double"
  if (!identical(dim(region), c(d, 2L)) ||
      any(!is.finite(region)) ||
      any(region[, 2L] <= region[, 1L])) {
    stop("`region` must be a finite d by 2 matrix with upper > lower.",
         call. = FALSE)
  }
  region
}

.ppt_validate_max_aspect_ratio <- function(max_aspect_ratio) {
  if (!is.numeric(max_aspect_ratio) || length(max_aspect_ratio) != 1L ||
      is.na(max_aspect_ratio) || max_aspect_ratio < 1) {
    stop("`max_aspect_ratio` must be a scalar in [1, Inf].",
         call. = FALSE)
  }
  as.numeric(max_aspect_ratio)
}

.ppt_expand_parameter <- function(x, d, name, allow_zero = FALSE) {
  x <- as.numeric(x)
  if (length(x) == 1L) x <- rep(x, d)
  bad <- length(x) != d || any(!is.finite(x)) ||
    if (allow_zero) any(x < 0) else any(x <= 0)
  if (bad) {
    stop(
      "`", name, "` must have length one or d and be ",
      if (allow_zero) "nonnegative." else "positive.",
      call. = FALSE
    )
  }
  x
}

.ppt_validate_points <- function(x, d = NULL, region = NULL,
                                 name = "x", allow_empty = FALSE) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if ((!allow_empty && !nrow(x)) || !ncol(x) || any(!is.finite(x))) {
    stop("`", name, "` must be a finite numeric matrix",
         if (!allow_empty) " with at least one row." else ".",
         call. = FALSE)
  }
  if (!is.null(d) && ncol(x) != d) {
    stop("`", name, "` must have ", d, " columns.", call. = FALSE)
  }
  if (!is.null(region) && nrow(x)) {
    scale <- max(1, abs(region))
    tolerance <- sqrt(.Machine$double.eps) * scale
    too_low <- sweep(x, 2L, region[, 1L] - tolerance, "<")
    too_high <- sweep(x, 2L, region[, 2L] + tolerance, ">")
    if (any(too_low | too_high)) {
      stop("All observations in `", name,
           "` must lie inside `region`.", call. = FALSE)
    }
    for (j in seq_len(ncol(x))) {
      x[, j] <- pmin(region[j, 2L], pmax(region[j, 1L], x[, j]))
    }
  }
  x
}

.ppt_match_prediction_rows <- function(newdata, locations) {
  key <- function(z) {
    apply(z, 1L, function(row)
      paste(format(row, digits = 17L, scientific = TRUE), collapse = "\r"))
  }
  match(key(newdata), key(locations))
}

.ppt_posterior_weights <- function(object, n_draws) {
  weights <- object$posterior$particle_weights
  if (length(weights) != n_draws || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    return(rep(1 / n_draws, n_draws))
  }
  as.numeric(weights / sum(weights))
}

.ppt_weighted_quantile <- function(x, weights, probability) {
  keep <- is.finite(x) & is.finite(weights) & weights >= 0
  x <- x[keep]
  weights <- weights[keep]
  if (!length(x) || sum(weights) <= 0) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  weights <- weights[ord] / sum(weights)
  x[which(cumsum(weights) >= probability)[1L]]
}

.ppt_tree_marginal <- function(tree, grid, variable, region,
                               average = TRUE) {
  leaves <- Filter(
    function(node) !is.null(node) && isTRUE(node$is_leaf), tree
  )
  if (!length(leaves)) {
    stop("A posterior tree draw contains no terminal nodes.",
         call. = FALSE)
  }
  d <- nrow(region)
  out <- numeric(length(grid))
  upper <- region[variable, 2L]
  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(upper))
  other_axes <- setdiff(seq_len(d), variable)

  for (leaf in leaves) {
    box <- as.matrix(leaf$region)
    lambda <- as.numeric(leaf$lambda)
    if (!identical(dim(box), c(d, 2L)) || length(lambda) != 1L ||
        !is.finite(lambda) || lambda < 0) {
      stop("A posterior tree draw contains an invalid leaf.",
           call. = FALSE)
    }
    is_last <- abs(box[variable, 2L] - upper) <= tolerance
    active <- grid >= box[variable, 1L] &
      (grid < box[variable, 2L] |
       (is_last & grid <= box[variable, 2L] + tolerance))
    other_volume <- if (length(other_axes)) {
      prod(box[other_axes, 2L] - box[other_axes, 1L])
    } else {
      1
    }
    out[active] <- out[active] + lambda * other_volume
  }

  if (isTRUE(average) && length(other_axes)) {
    domain_volume <- prod(
      region[other_axes, 2L] - region[other_axes, 1L]
    )
    out <- out / domain_volume
  }
  out
}

.ppt_fit_hard_leaf_smc <- function(
    x, region, predict_at = x, test = NULL,
    a = 0.5, b = 0, resample_thresh = 0.5,
    max_depth = 8L, min_leaf_n = 1L,
    max_aspect_ratio = Inf,
    particles = 1000L, seed = 1L) {
  x <- .ppt_validate_points(x, name = "x")
  d <- ncol(x)
  region <- .ppt_validate_region(region, d)
  x <- .ppt_validate_points(x, d, region, "x")
  predict_at <- .ppt_validate_points(
    predict_at, d, region, "predict_at", allow_empty = TRUE
  )
  if (is.null(test)) {
    test <- matrix(numeric(0), 0L, d)
  } else {
    test <- .ppt_validate_points(
      test, d, region, "test", allow_empty = TRUE
    )
  }

  controls <- c(max_depth, min_leaf_n, particles)
  if (any(!is.finite(controls)) ||
      max_depth < 0 || min_leaf_n < 1 || particles < 2 ||
      any(controls != floor(controls))) {
    stop("Invalid SMC or tree controls.", call. = FALSE)
  }
  if (length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || a <= 0 || b < 0) {
    stop("Require scalar a > 0 and b >= 0.", call. = FALSE)
  }
  if (length(resample_thresh) != 1L || !is.finite(resample_thresh) ||
      resample_thresh <= 0 || resample_thresh > 1) {
    stop("resample_thresh must lie in (0, 1].", call. = FALSE)
  }
  max_aspect_ratio <- .ppt_validate_max_aspect_ratio(max_aspect_ratio)
  max_depth <- as.integer(max_depth)
  min_leaf_n <- as.integer(min_leaf_n)
  particles <- as.integer(particles)

  evaluation_locations <- rbind(predict_at, test)
  set.seed(seed)
  raw <- PPT_fit_SMC(
    pts = x, grid = evaluation_locations, region = region,
    max_depth = max_depth, P = particles,
    min_leaf_n = min_leaf_n, resample_thresh = resample_thresh,
    a = a, b = b, max_aspect_ratio = max_aspect_ratio
  )

  weights <- as.numeric(raw$weights)
  if (length(weights) != length(raw$particle) ||
      any(!is.finite(weights)) || any(weights < 0) ||
      sum(weights) <= 0) {
    weights <- rep(1 / length(raw$particle), length(raw$particle))
  } else {
    weights <- weights / sum(weights)
  }
  particle_ess <- 1 / sum(weights^2)

  leaf_counts <- .leaves_per_particle(raw$particle)
  max_depths <- vapply(raw$particle, function(particle) {
    leaves <- Filter(
      function(node) !is.null(node) && isTRUE(node$is_leaf),
      particle
    )
    if (!length(leaves)) return(NA_real_)
    max(vapply(leaves, function(node) as.numeric(node$depth), numeric(1)))
  }, numeric(1))
  mean_leaves <- if (length(leaf_counts) == length(weights)) {
    sum(weights * leaf_counts)
  } else {
    NA_real_
  }
  valid_depth <- is.finite(max_depths)
  mean_max_depth <- if (any(valid_depth)) {
    sum(weights[valid_depth] * max_depths[valid_depth]) /
      sum(weights[valid_depth])
  } else {
    NA_real_
  }

  n_prediction <- nrow(predict_at)
  prediction_rows <- if (n_prediction) seq_len(n_prediction) else integer()
  lambda_draws <- as.matrix(raw$lambda$draws)
  integrated_intensity <- vapply(raw$particle, function(particle) {
    leaves <- Filter(
      function(node) !is.null(node) && isTRUE(node$is_leaf),
      particle
    )
    sum(vapply(leaves, function(node) {
      box <- as.matrix(node$region)
      as.numeric(node$lambda) * prod(box[, 2L] - box[, 1L])
    }, numeric(1)))
  }, numeric(1))
  lppd <- NA_real_
  if (nrow(test)) {
    test_rows <- n_prediction + seq_len(nrow(test))
    test_draws <- pmax(
      lambda_draws[test_rows, , drop = FALSE],
      .Machine$double.xmin
    )
    log_predictive_draw <- colSums(log(test_draws)) - integrated_intensity
    log_weighted <- log(weights) + log_predictive_draw
    center <- max(log_weighted)
    lppd <- center + log(sum(exp(log_weighted - center)))
  }

  input_names <- colnames(x)
  if (is.null(input_names)) input_names <- paste0("x", seq_len(d))
  log_evidence <- as.numeric(raw$logZ)
  if (length(log_evidence) != 1L || !is.finite(log_evidence)) {
    log_evidence <- NA_real_
  }

  structure(
    list(
      call = NULL,
      model = list(
        gating = "hard",
        gate_family = NA_character_,
        scales = "leaf",
        sampler = "smc",
        label = "PPT"
      ),
      data = list(
        x = x,
        n = nrow(x),
        dimension = d,
        region = region,
        test = if (nrow(test)) test else NULL
      ),
      prediction = list(
        locations = predict_at,
        mean = as.numeric(raw$lambda$mean)[prediction_rows],
        median = as.numeric(raw$lambda$median)[prediction_rows],
        lower = as.numeric(raw$lambda$lower95)[prediction_rows],
        upper = as.numeric(raw$lambda$upper95)[prediction_rows],
        level = 0.95,
        draws = lambda_draws[prediction_rows, , drop = FALSE]
      ),
      posterior = list(
        mean_leaves = as.numeric(mean_leaves),
        mean_max_depth = as.numeric(mean_max_depth),
        mean_gate = NA_real_,
        gate_by_dimension = stats::setNames(
          rep(NA_real_, d), input_names
        ),
        mean_log_likelihood = as.numeric(raw$loglik),
        mean_integrated_intensity = sum(weights * integrated_intensity),
        integrated_intensity_draws = as.numeric(integrated_intensity),
        lppd = as.numeric(lppd),
        log_evidence = log_evidence,
        draws = length(weights),
        particle_weights = weights,
        tree_draws = raw$particle,
        state = list(mode = "leafbox", chain = NULL)
      ),
      diagnostics = list(
        acceptance = numeric(),
        tree_acceptance = numeric(),
        gate_acceptance = numeric(),
        particle_ess = as.numeric(particle_ess),
        ess_history = as.numeric(raw$ESS),
        unique_trees = .n_unique_trees(raw$particle),
        leaf_count_trace = as.numeric(leaf_counts),
        max_depth_trace = as.numeric(max_depths),
        log_evidence_increment = as.numeric(raw$logZ_inc),
        log_evidence_running = as.numeric(raw$logZ_run)
      ),
      prior = list(
        intensity = list(shape = a, rate = b),
        tree = list(split_probability = 0.5, axis_rate = 1 / d)
      ),
      control = list(
        max_depth = max_depth,
        min_leaf_n = min_leaf_n,
        particles = particles,
        resample_thresh = resample_thresh,
        max_aspect_ratio = max_aspect_ratio,
        seed = seed
      ),
      backend = "PPT_fit_SMC"
    ),
    class = "ppt"
  )
}

.ppt_fit_hard_leaf_rjmcmc <- function(
    x, region, predict_at = x, test = NULL,
    a = 0.5, b = 0, alpha = 0.95, eta = 2,
    max_depth = 8L, min_leaf_n = 1L,
    chains = 4L, iter = 4000L, burn = 1000L,
    cut_candidates = 30L, prediction_draws = 300L,
    seed = 1L, verbose = TRUE) {
  x <- .ppt_validate_points(x, name = "x")
  d <- ncol(x)
  region <- .ppt_validate_region(region, d)
  x <- .ppt_validate_points(x, d, region, "x")
  predict_at <- .ppt_validate_points(
    predict_at, d, region, "predict_at", allow_empty = TRUE
  )
  if (is.null(test)) {
    test <- matrix(numeric(0), 0L, d)
  } else {
    test <- .ppt_validate_points(
      test, d, region, "test", allow_empty = TRUE
    )
  }

  if (length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || a <= 0 || b < 0) {
    stop("Require scalar `a > 0` and `b >= 0`.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1 ||
      length(eta) != 1L || !is.finite(eta) || eta < 0) {
    stop("Require `0 < alpha < 1` and `eta >= 0`.", call. = FALSE)
  }

  controls <- c(
    max_depth, min_leaf_n, chains, iter, burn,
    cut_candidates, prediction_draws
  )
  if (any(!is.finite(controls)) || any(controls != floor(controls)) ||
      max_depth < 0 || min_leaf_n < 1 || chains < 1 ||
      iter <= burn || burn < 0 || cut_candidates < 2 ||
      prediction_draws < 1) {
    stop("Invalid RJ-MCMC or tree controls.", call. = FALSE)
  }
  max_depth <- as.integer(max_depth)
  min_leaf_n <- as.integer(min_leaf_n)
  chains <- as.integer(chains)
  iter <- as.integer(iter)
  burn <- as.integer(burn)
  cut_candidates <- as.integer(cut_candidates)
  prediction_draws <- as.integer(prediction_draws)
  kept <- iter - burn

  evaluation_locations <- rbind(predict_at, test)
  set.seed(seed)
  raw_chains <- vector("list", chains)
  for (chain in seq_len(chains)) {
    raw_chains[[chain]] <- PPT_fit_MCMC(
      x, evaluation_locations, region,
      niter = kept, burnin = burn,
      max_depth = max_depth, min_leaf_n = min_leaf_n,
      cut_grid_n = cut_candidates,
      a = a, b = b, alpha = alpha, eta = eta,
      n_pred = min(prediction_draws, kept)
    )
    if (isTRUE(verbose)) {
      message(
        "PPT [RJ-MCMC] chain ", chain, "/", chains,
        " completed; mean leaves = ",
        format(mean(raw_chains[[chain]]$nleaves), digits = 5L)
      )
    }
  }

  lambda_draws <- do.call(
    cbind, lapply(raw_chains, function(z) as.matrix(z$lambda$draws))
  )
  tree_draws <- unlist(
    lapply(raw_chains, `[[`, "tree_draws"), recursive = FALSE
  )
  integrated_intensity <- unlist(
    lapply(raw_chains, `[[`, "integrated_intensity"), use.names = FALSE
  )
  poisson_loglik <- unlist(
    lapply(raw_chains, `[[`, "poisson_loglik"), use.names = FALSE
  )
  leaf_trace <- unlist(
    lapply(raw_chains, `[[`, "nleaves"), use.names = FALSE
  )
  depth_trace <- unlist(
    lapply(raw_chains, `[[`, "max_depth"), use.names = FALSE
  )
  acceptance <- Reduce(
    "+", lapply(raw_chains, function(z) as.numeric(z$accept))
  ) / chains
  names(acceptance) <- c("grow", "prune", "change")

  n_prediction <- nrow(predict_at)
  prediction_rows <- if (n_prediction) seq_len(n_prediction) else integer()
  prediction_matrix <- lambda_draws[prediction_rows, , drop = FALSE]
  summarize <- function(prob = NULL) {
    if (!n_prediction) return(numeric())
    if (is.null(prob)) return(rowMeans(prediction_matrix))
    apply(
      prediction_matrix, 1L, stats::quantile,
      probs = prob, names = FALSE, type = 8
    )
  }

  lppd <- NA_real_
  if (nrow(test)) {
    test_rows <- n_prediction + seq_len(nrow(test))
    test_draws <- pmax(
      lambda_draws[test_rows, , drop = FALSE],
      .Machine$double.xmin
    )
    log_predictive_draw <-
      colSums(log(test_draws)) - integrated_intensity
    center <- max(log_predictive_draw)
    lppd <- center + log(mean(exp(log_predictive_draw - center)))
  }

  input_names <- colnames(x)
  if (is.null(input_names)) input_names <- paste0("x", seq_len(d))

  structure(
    list(
      call = NULL,
      model = list(
        gating = "hard",
        gate_family = NA_character_,
        scales = "leaf",
        sampler = "rjmcmc",
        label = "PPT"
      ),
      data = list(
        x = x,
        n = nrow(x),
        dimension = d,
        region = region,
        test = if (nrow(test)) test else NULL
      ),
      prediction = list(
        locations = predict_at,
        mean = summarize(),
        median = summarize(0.5),
        lower = summarize(0.025),
        upper = summarize(0.975),
        level = 0.95,
        draws = prediction_matrix
      ),
      posterior = list(
        mean_leaves = mean(leaf_trace),
        mean_max_depth = mean(depth_trace),
        mean_gate = NA_real_,
        gate_by_dimension = stats::setNames(
          rep(NA_real_, d), input_names
        ),
        mean_log_likelihood = mean(poisson_loglik),
        mean_integrated_intensity = mean(integrated_intensity),
        integrated_intensity_draws = integrated_intensity,
        lppd = as.numeric(lppd),
        log_evidence = NA_real_,
        draws = ncol(lambda_draws),
        particle_weights = numeric(),
        tree_draws = tree_draws,
        state = list(
          mode = "leafbox",
          chain = rep(
            seq_along(raw_chains),
            vapply(raw_chains, function(z) length(z$tree_draws), integer(1))
          )
        )
      ),
      diagnostics = list(
        acceptance = acceptance,
        tree_acceptance = acceptance,
        gate_acceptance = numeric(),
        particle_ess = NA_real_,
        ess_history = numeric(),
        unique_trees = .n_unique_trees(tree_draws),
        leaf_count_trace = as.numeric(leaf_trace),
        max_depth_trace = as.numeric(depth_trace),
        log_evidence_increment = numeric(),
        log_evidence_running = numeric()
      ),
      prior = list(
        intensity = list(shape = a, rate = b),
        tree = list(alpha = alpha, eta = eta)
      ),
      control = list(
        max_depth = max_depth,
        min_leaf_n = min_leaf_n,
        chains = chains,
        iter = iter,
        burn = burn,
        cut_candidates = cut_candidates,
        prediction_draws = prediction_draws,
        seed = seed
      ),
      backend = "PPT_fit_MCMC"
    ),
    class = "ppt"
  )
}

.ppt_fit_hard_leaf_pgas <- function(
    x, region, predict_at = x, test = NULL,
    a = 0.5, b = 0,
    max_depth = 8L, min_leaf_n = 1L,
    max_aspect_ratio = Inf,
    particles = 500L, chains = 1L,
    iter = 500L, burn = 100L,
    seed = 1L, verbose = TRUE) {
  x <- .ppt_validate_points(x, name = "x")
  d <- ncol(x)
  region <- .ppt_validate_region(region, d)
  x <- .ppt_validate_points(x, d, region, "x")
  predict_at <- .ppt_validate_points(
    predict_at, d, region, "predict_at", allow_empty = TRUE
  )
  if (is.null(test)) {
    test <- matrix(numeric(0), 0L, d)
  } else {
    test <- .ppt_validate_points(
      test, d, region, "test", allow_empty = TRUE
    )
  }

  if (length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || a <= 0 || b < 0) {
    stop("Require scalar `a > 0` and `b >= 0`.", call. = FALSE)
  }
  controls <- c(
    max_depth, min_leaf_n, particles, chains, iter, burn
  )
  if (any(!is.finite(controls)) || any(controls != floor(controls)) ||
      max_depth < 0 || min_leaf_n < 1 || particles < 2 ||
      chains < 1 || iter <= burn || burn < 0) {
    stop("Invalid Particle-Gibbs or tree controls.", call. = FALSE)
  }
  max_aspect_ratio <- .ppt_validate_max_aspect_ratio(max_aspect_ratio)
  max_depth <- as.integer(max_depth)
  min_leaf_n <- as.integer(min_leaf_n)
  particles <- as.integer(particles)
  chains <- as.integer(chains)
  iter <- as.integer(iter)
  burn <- as.integer(burn)
  retained <- seq.int(burn + 1L, iter)

  evaluation_locations <- rbind(predict_at, test)
  set.seed(seed)
  raw_chains <- vector("list", chains)
  for (chain in seq_len(chains)) {
    raw_chains[[chain]] <- PPT_fit_PG(
      x, evaluation_locations, region,
      max_depth = max_depth, niter = iter, P = particles,
      min_leaf_n = min_leaf_n, resample_thresh = 0.5,
      a = a, b = b, verbose = verbose,
      max_aspect_ratio = max_aspect_ratio
    )
  }

  lambda_draws <- do.call(cbind, lapply(raw_chains, function(z) {
    as.matrix(z$lambda$draws)[, retained, drop = FALSE]
  }))
  tree_draws <- unlist(lapply(raw_chains, function(z) {
    z$particles[retained]
  }), recursive = FALSE)
  loglik_draws <- unlist(lapply(raw_chains, function(z) {
    as.numeric(z$loglik)[retained]
  }), use.names = FALSE)
  conditional_weights <- do.call(cbind, lapply(raw_chains, function(z) {
    as.matrix(z$weights)[, retained, drop = FALSE]
  }))

  leaf_counts <- .leaves_per_particle(tree_draws)
  max_depths <- vapply(tree_draws, function(tree) {
    leaves <- Filter(
      function(node) !is.null(node) && isTRUE(node$is_leaf), tree
    )
    if (!length(leaves)) return(NA_real_)
    max(vapply(leaves, function(node) as.numeric(node$depth), numeric(1)))
  }, numeric(1))
  integrated_intensity <- vapply(tree_draws, function(tree) {
    leaves <- Filter(
      function(node) !is.null(node) && isTRUE(node$is_leaf), tree
    )
    sum(vapply(leaves, function(node) {
      box <- as.matrix(node$region)
      as.numeric(node$lambda) * prod(box[, 2L] - box[, 1L])
    }, numeric(1)))
  }, numeric(1))
  pgas_ess <- apply(conditional_weights, 2L, function(weights) {
    weights <- pmax(as.numeric(weights), 0)
    if (!all(is.finite(weights)) || sum(weights) <= 0) return(NA_real_)
    weights <- weights / sum(weights)
    1 / sum(weights^2)
  })

  n_prediction <- nrow(predict_at)
  prediction_rows <- if (n_prediction) seq_len(n_prediction) else integer()
  prediction_matrix <- lambda_draws[prediction_rows, , drop = FALSE]
  summarize <- function(prob = NULL) {
    if (!n_prediction) return(numeric())
    if (is.null(prob)) return(rowMeans(prediction_matrix))
    apply(
      prediction_matrix, 1L, stats::quantile,
      probs = prob, names = FALSE, type = 8
    )
  }

  lppd <- NA_real_
  if (nrow(test)) {
    test_rows <- n_prediction + seq_len(nrow(test))
    test_draws <- pmax(
      lambda_draws[test_rows, , drop = FALSE],
      .Machine$double.xmin
    )
    log_predictive_draw <-
      colSums(log(test_draws)) - integrated_intensity
    center <- max(log_predictive_draw)
    lppd <- center + log(mean(exp(log_predictive_draw - center)))
  }

  input_names <- colnames(x)
  if (is.null(input_names)) input_names <- paste0("x", seq_len(d))

  structure(
    list(
      call = NULL,
      model = list(
        gating = "hard",
        gate_family = NA_character_,
        scales = "leaf",
        sampler = "pgas",
        label = "PPT"
      ),
      data = list(
        x = x,
        n = nrow(x),
        dimension = d,
        region = region,
        test = if (nrow(test)) test else NULL
      ),
      prediction = list(
        locations = predict_at,
        mean = summarize(),
        median = summarize(0.5),
        lower = summarize(0.025),
        upper = summarize(0.975),
        level = 0.95,
        draws = prediction_matrix
      ),
      posterior = list(
        mean_leaves = mean(leaf_counts),
        mean_max_depth = mean(max_depths),
        mean_gate = NA_real_,
        gate_by_dimension = stats::setNames(
          rep(NA_real_, d), input_names
        ),
        mean_log_likelihood = mean(loglik_draws),
        mean_integrated_intensity = mean(integrated_intensity),
        integrated_intensity_draws = integrated_intensity,
        lppd = as.numeric(lppd),
        log_evidence = NA_real_,
        draws = ncol(lambda_draws),
        particle_weights = numeric(),
        tree_draws = tree_draws,
        state = list(
          mode = "leafbox",
          chain = rep(seq_len(chains), each = length(retained))
        )
      ),
      diagnostics = list(
        acceptance = numeric(),
        tree_acceptance = numeric(),
        gate_acceptance = numeric(),
        particle_ess = mean(pgas_ess, na.rm = TRUE),
        ess_history = as.numeric(pgas_ess),
        unique_trees = .n_unique_trees(tree_draws),
        leaf_count_trace = as.numeric(leaf_counts),
        max_depth_trace = as.numeric(max_depths),
        log_evidence_increment = numeric(),
        log_evidence_running = numeric()
      ),
      prior = list(
        intensity = list(shape = a, rate = b),
        tree = list(split_probability = 0.5, axis_rate = 1 / d)
      ),
      control = list(
        max_depth = max_depth,
        min_leaf_n = min_leaf_n,
        particles = particles,
        chains = chains,
        iter = iter,
        burn = burn,
        conditional_smc = TRUE,
        ancestor_sampling = FALSE,
        resampling_schedule = "tree_level",
        max_aspect_ratio = max_aspect_ratio,
        seed = seed
      ),
      backend = "PPT_fit_PG"
    ),
    class = "ppt"
  )
}
