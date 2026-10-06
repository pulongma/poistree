# Small, deterministic posterior states. References below use explicit leaf
# paths and probability formulas, independently of native prediction code.
prediction_region <- matrix(c(-2, 3, 4, 12), ncol = 2, byrow = TRUE)

prediction_tree <- function(multiplier = 1) {
  cbind(heap_id = 1:9,
        axis = c(0, 1, 0, 0, -1, -1, -1, -1, -1),
        cut = c(-0.3, 7, 1.8, -1, rep(NA_real_, 5)),
        lambda = multiplier * c(0, 0, 0, 0, 5, 12, 2, 3, 8))
}

prediction_leaf_path <- function(nodes, leaf, region, gate_depth) {
  ids <- numeric()
  child <- leaf
  while (child > 1) {
    ids <- c(child, ids)
    child <- floor(child / 2)
  }
  box <- region
  path <- list()
  for (child in ids) {
    parent <- nodes[match(floor(child / 2), nodes[, 1]), ]
    axis <- parent[2] + 1L
    width <- box[axis, 2] - box[axis, 1]
    depth <- floor(log2(parent[1]))
    right <- child %% 2 == 1
    path[[length(path) + 1L]] <- list(
      axis = axis, cut = parent[3], right = right,
      width = width, compact_width = width / (1 + depth)^gate_depth
    )
    if (right) box[axis, 1] <- parent[3] else box[axis, 2] <- parent[3]
  }
  # The legacy evaluator sums from the terminal node towards the root.
  rev(path)
}

prediction_reference <- function(nodes, gates, at, region, mode,
                                 gate_depth = 0) {
  answer <- matrix(0, nrow(at), length(nodes))
  for (s in seq_along(nodes)) {
    tree <- nodes[[s]][order(nodes[[s]][, 1]), , drop = FALSE]
    leaves <- which(tree[, 2] < 0)
    for (leaf in leaves) {
      rate <- tree[leaf, 4]
      if (!is.finite(rate) || rate <= 0) next
      path <- prediction_leaf_path(tree, tree[leaf, 1], region, gate_depth)
      log_membership <- rep(0, nrow(at))
      for (node in path) {
        j <- node$axis
        if (mode %in% c(2L, 3L)) {
          width <- if (mode == 3L) node$width else diff(region[j, ])
          z <- gates[s, j] * (at[, j] - node$cut) / width
          log_probability <- stats::plogis(
            if (node$right) z else -z, log.p = TRUE
          )
        } else {
          half_width <- node$compact_width / gates[s, j]
          t <- pmax(0, pmin(1, (at[, j] - (node$cut - half_width)) /
                                  (2 * half_width)))
          probability <- t^2 * (3 - 2 * t)
          if (!node$right) probability <- 1 - probability
          log_probability <- log(probability)
        }
        log_membership <- log_membership + log_probability
      }
      membership <- ifelse(log_membership > -745, exp(log_membership), 0)
      answer[, s] <- answer[, s] + rate * membership
    }
  }
  answer
}

prediction_fixture <- function(mode = 3L, root_only = FALSE) {
  multipliers <- c(0.8, 0.8, 1, 1.9)
  nodes <- lapply(multipliers, prediction_tree)
  if (root_only) {
    nodes <- lapply(c(2, 2, 5, 9), function(rate)
      matrix(c(1, -1, NA_real_, rate), nrow = 1))
  }
  gates <- matrix(rep(c(8, 13), each = 4), nrow = 4)
  locations <- rbind(c(-1.5, 5), c(0, 8), c(2.2, 11), c(0.7, 6))
  colnames(locations) <- c("longitude", "season")
  draws <- prediction_reference(nodes, gates, locations, prediction_region, mode)
  weights <- c(0, 0.2, 0.5, 0.3)
  q <- t(apply(draws, 1, stats::quantile, probs = c(0.5, 0.025, 0.975)))
  # Stored backend summaries deliberately follow their original type-7
  # quantile convention, which ppt_predict on stored rows must preserve.
  structure(list(
    model = list(gating = "soft", scales = "leaf"),
    data = list(x = locations, dimension = 2L, region = prediction_region,
                test = locations),
    prediction = list(locations = locations, draws = draws,
                      mean = as.numeric(draws %*% weights),
                      median = q[, 1], lower = q[, 2], upper = q[, 3]),
    posterior = list(
      tree_draws = list(), particle_weights = weights,
      integrated_intensity_draws = if (root_only) 40 * c(2, 2, 5, 9) else 100 * multipliers,
      mean_integrated_intensity = if (root_only) 40 * sum(weights * c(2, 2, 5, 9)) else 100 * sum(weights * multipliers),
      lppd = NA_real_,
      state = list(mode = "heap", nodes = nodes, gate = gates,
                   gate_mode = mode, gate_depth = 0)
    )
  ), class = "ppt")
}
