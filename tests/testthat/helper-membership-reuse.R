# Independent R gate and tree calculations for the membership reuse tests.
# C++ gate families are 0: root logistic, 1: compact, 2: node logistic.
membership_gate_logs <- function(x, axis, cut, parent_width, region, gate, family) {
  if (family == 1L) {
    h <- parent_width / gate[axis]
    t <- pmax(0, pmin(1, (x[, axis] - (cut - h)) / (2 * h)))
    right <- t * t * (3 - 2 * t)
    return(cbind(log(1 - right), log(right)))
  }
  width <- if (family == 2L) parent_width else diff(region[axis, ])
  z <- gate[axis] * (x[, axis] - cut) / width
  cbind(plogis(z, lower.tail = FALSE, log.p = TRUE),
        plogis(z, log.p = TRUE))
}

membership_tree_log_phi <- function(x, region, splits, gate, family) {
  leaves <- list()
  visit <- function(id, prefix, box) {
    row <- match(id, splits[, 1L])
    if (is.na(row)) {
      leaves[[as.character(id)]] <<- prefix
      return(invisible(NULL))
    }
    axis <- as.integer(splits[row, 2L]) + 1L
    cut <- splits[row, 3L]
    log_gate <- membership_gate_logs(x, axis, cut, diff(box[axis, ]),
                                     region, gate, family)
    left_box <- right_box <- box
    left_box[axis, 2L] <- right_box[axis, 1L] <- cut
    visit(2L * id, prefix + log_gate[, 1L], left_box)
    visit(2L * id + 1L, prefix + log_gate[, 2L], right_box)
  }
  visit(1L, numeric(nrow(x)), region)
  order <- order(as.integer(names(leaves)))
  list(log_phi = do.call(cbind, unname(leaves[order])),
       ids = as.integer(names(leaves)[order]))
}

membership_native_run <- function(args, enabled, diagnostic = FALSE) {
  args$cache_geometry <- enabled
  cut_cache_run(args, TRUE, diagnostic = diagnostic)
}
