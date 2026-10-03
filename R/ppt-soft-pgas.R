# Fixed per-axis cut grid for the soft Particle-Gibbs backend: type-1 quantiles
# of each coordinate on [0.05, 0.95], strictly inside the region.  A fixed,
# box-free grid is what makes ancestor sampling valid for trees.
.ppt_soft_cut_grid <- function(x, region, cut_candidates) {
  probs <- seq(0.05, 0.95, length.out = cut_candidates)
  lapply(seq_len(ncol(x)), function(j) {
    buffer <- 1e-3 * (region[j, 2L] - region[j, 1L])
    cuts <- unique(stats::quantile(x[, j], probs = probs, type = 1, names = FALSE))
    cuts <- cuts[cuts > region[j, 1L] + buffer & cuts < region[j, 2L] - buffer]
    if (!length(cuts)) cuts <- mean(region[j, ])
    sort(cuts)
  })
}

.ppt_fit_soft_leaf_pgas <- function(
    x, region, predict_at = x, test = NULL,
    a = 0.5, b = NULL,
    gate = 12, a_gate = 36, b_gate = 3, sd_gate = 0.07, gate_min = 0,
    gate_structure = c("dimension", "shared"),
    update_gate = TRUE,
    alpha = 0.95, eta = 2,
    max_depth = 6L, cut_candidates = 30L, cut_grid = NULL,
    particles = 50L, chains = 1L, iter = 1000L, burn = 200L, thin = 1L,
    label_sweeps = 1L, ancestor_sampling = TRUE,
    exact_max = 150L, defensive = 0,
    resampling = c("node", "level"), ess_threshold = 1,
    allocation = c("sequential", "rates"),
    seed = 1L, verbose = TRUE) {
  gate_structure <- match.arg(gate_structure)
  resampling <- match.arg(resampling)
  allocation <- match.arg(allocation)

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
    test <- .ppt_validate_points(test, d, region, "test", allow_empty = TRUE)
  }

  if (is.null(b)) b <- a * prod(region[, 2L] - region[, 1L]) / nrow(x)
  if (length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    stop("The soft leaf model requires positive scalar `a` and `b`.",
         call. = FALSE)
  }
  if (identical(gate_structure, "shared") &&
      any(lengths(list(gate, a_gate, b_gate, sd_gate, gate_min)) != 1L)) {
    stop("Shared gates require scalar gate-prior parameters.", call. = FALSE)
  }
  gate <- .ppt_expand_parameter(gate, d, "gate")
  a_gate <- .ppt_expand_parameter(a_gate, d, "a_gate")
  b_gate <- .ppt_expand_parameter(b_gate, d, "b_gate")
  sd_gate <- .ppt_expand_parameter(sd_gate, d, "sd_gate")
  gate_min <- .ppt_expand_parameter(gate_min, d, "gate_min", allow_zero = TRUE)
  if (any(gate <= gate_min)) {
    stop("Every initial gate must exceed its `gate_min`.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      length(eta) != 1L || !is.finite(eta) || eta < 0) {
    stop("Require `0 < alpha < 1` and `eta >= 0`.", call. = FALSE)
  }
  controls <- c(max_depth, cut_candidates, particles, chains, iter, burn, thin,
                label_sweeps, exact_max)
  if (any(!is.finite(controls)) || any(controls != floor(controls)) ||
      max_depth < 1 || cut_candidates < 2 || particles < 2 || chains < 1 ||
      iter <= burn || burn < 0 || thin < 1 || label_sweeps < 0 || exact_max < 0) {
    stop("Invalid Particle-Gibbs or tree controls.", call. = FALSE)
  }
  if (length(defensive) != 1L || !is.finite(defensive) ||
      defensive < 0 || defensive >= 1) {
    stop("`defensive` must lie in [0, 1).", call. = FALSE)
  }
  if (length(ess_threshold) != 1L || !is.finite(ess_threshold) ||
      ess_threshold <= 0 || ess_threshold > 1) {
    stop("`ess_threshold` must lie in (0, 1].", call. = FALSE)
  }
  for (flag in list(update_gate, ancestor_sampling)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) {
      stop("`update_gate` and `ancestor_sampling` must be logical scalars.",
           call. = FALSE)
    }
  }
  if (is.null(cut_grid)) {
    cut_grid <- .ppt_soft_cut_grid(x, region, cut_candidates)
  } else {
    if (!is.list(cut_grid) || length(cut_grid) != d) {
      stop("`cut_grid` must be a list with one numeric vector per input.",
           call. = FALSE)
    }
    cut_grid <- lapply(cut_grid, function(g) sort(unique(as.numeric(g))))
  }

  set.seed(seed)
  raw <- lapply(seq_len(chains), function(chain) {
    SPPT_fit_PGAS(
      x, predict_at, test, region, cut_grid,
      a, b, gate, a_gate, b_gate, sd_gate, gate_min,
      identical(gate_structure, "shared"), alpha, eta, as.integer(max_depth),
      as.integer(particles), as.integer(iter), as.integer(burn), as.integer(thin),
      as.integer(label_sweeps), isTRUE(update_gate), isTRUE(ancestor_sampling),
      as.integer(exact_max), defensive,
      identical(resampling, "node"), ess_threshold,
      identical(allocation, "rates"), isTRUE(verbose)
    )
  })
  pull <- function(name) unlist(lapply(raw, function(z) as.numeric(z[[name]])))
  draws <- do.call(cbind, lapply(raw, function(z) as.matrix(z$draws)))
  state_nodes <- unlist(lapply(raw, function(z) z$state_nodes), recursive = FALSE)
  state_gate <- do.call(rbind, lapply(raw, function(z) as.matrix(z$gate)))
  loglik <- pull("loglik")
  integrated_intensity <- pull("integrated_intensity")
  retained_per_chain <- .ppt_retained_per_chain(iter, burn, thin)

  lppd <- NA_real_
  if (nrow(test)) {
    lt <- pull("loglik_test")
    lppd <- max(lt) + log(mean(exp(lt - max(lt))))
  }
  input_names <- colnames(x)
  if (is.null(input_names)) input_names <- paste0("x", seq_len(d))
  gate_mean <- stats::setNames(colMeans(state_gate), input_names)
  gate_accept <- mean(vapply(raw, function(z) as.numeric(z$gate_accept), numeric(1)))
  summarize <- function(prob = NULL) {
    if (!nrow(predict_at)) return(numeric())
    if (is.null(prob)) return(rowMeans(draws))
    apply(draws, 1L, stats::quantile, probs = prob, names = FALSE, type = 8)
  }

  structure(
    list(
      call = NULL,
      model = list(
        gating = "soft", gate_family = "logistic",
        scales = "leaf",
        sampler = "pgas",
        algorithm = if (isTRUE(ancestor_sampling)) "PGAS" else "PG",
        label = "S-PPT"
      ),
      data = list(
        x = x, n = nrow(x), dimension = d, region = region,
        test = if (nrow(test)) test else NULL
      ),
      prediction = list(
        locations = predict_at,
        mean = summarize(), median = summarize(0.5),
        lower = summarize(0.025), upper = summarize(0.975),
        level = 0.95, draws = draws
      ),
      posterior = list(
        mean_leaves = mean(pull("nleaf")),
        mean_max_depth = mean(pull("max_depth")),
        mean_gate = if (identical(gate_structure, "shared")) {
          unname(gate_mean[1L])
        } else {
          mean(gate_mean)
        },
        gate_by_dimension = gate_mean,
        mean_log_likelihood = mean(loglik),
        mean_integrated_intensity = mean(integrated_intensity),
        integrated_intensity_draws = integrated_intensity,
        lppd = lppd,
        log_evidence = NA_real_,
        draws = ncol(draws),
        particle_weights = numeric(),
        tree_draws = list(),
        state = list(
          mode = "heap",
          nodes = state_nodes,
          gate = state_gate,
          gate_mode = 2L,
          gate_depth = 0,
          chain = rep(seq_len(chains), each = retained_per_chain)
        )
      ),
      diagnostics = list(
        acceptance = stats::setNames(gate_accept, "gate"),
        tree_acceptance = numeric(),
        gate_acceptance = stats::setNames(rep(gate_accept, d), input_names),
        particle_ess = mean(pull("ess")),
        ess_history = pull("ess"),
        ancestor_move_rate = mean(pull("as_rate"), na.rm = TRUE),
        expanded_nodes = mean(pull("expanded")),
        resampling_events = mean(pull("resampled")),
        unique_trees = NA_integer_,
        leaf_count_trace = pull("nleaf"),
        max_depth_trace = pull("max_depth"),
        log_evidence_increment = numeric(),
        log_evidence_running = numeric()
      ),
      prior = list(
        intensity = list(shape = a, rate = b),
        gate = list(
          shape = a_gate, rate = b_gate, lower = gate_min,
          structure = gate_structure, family = "logistic"
        ),
        tree = list(alpha = alpha, eta = eta, cut_grid = cut_grid)
      ),
      control = list(
        max_depth = as.integer(max_depth),
        cut_candidates = as.integer(cut_candidates),
        particles = as.integer(particles),
        chains = as.integer(chains), iter = as.integer(iter),
        burn = as.integer(burn), thin = as.integer(thin),
        retained_per_chain = retained_per_chain,
        label_sweeps = as.integer(label_sweeps),
        conditional_smc = TRUE,
        ancestor_sampling = isTRUE(ancestor_sampling),
        resampling_schedule = if (identical(resampling, "node")) "tree_node" else "tree_level",
        ess_threshold = ess_threshold,
        allocation = allocation,
        exact_max = as.integer(exact_max), defensive = defensive,
        update_gate = isTRUE(update_gate), seed = seed
      ),
      backend = "SPPT_fit_PGAS"
    ),
    class = "ppt"
  )
}
