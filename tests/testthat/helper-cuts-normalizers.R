cut_cache_native_args <- function() {
  list(
    X = cbind(c(.02, .04, .09, .16, .24, .31, .56, .66, .74, .86, .93, .97),
              c(.08, .15, .67, .22, .75, .29, .55, .89, .65, .42, .87, .94)),
    grid = cbind(c(.03, .29, .63, .98), c(.21, .72, .37, .93)),
    Xtest = cbind(c(.17, .44, .81), c(.32, .67, .54)),
    region = matrix(c(0, 1, 0, 1), 2L, byrow = TRUE),
    a = .5, b = .04, gate = 7, a_gate = 2, b_gate = .3,
    sd_gate = .7, gate_min = 0, gate_shared = 0L,
    alpha = .85, eta = .5, Dmax = 3L, nmin = 1L,
    iters = 72L, burn = 16L, thin = 1L, nmove = 3L, ncc = 2L,
    cut_mode = 1L, ncand = 4L, update_gate = 1L, gate_family = 0L,
    chains = 1L, verbose = 0L, informed = FALSE, pcg = TRUE,
    ram_adapt = 12L
  )
}

cut_cache_run <- function(args, enabled, diagnostic = FALSE, seed = 619L) {
  if (diagnostic) {
    args$mon <- args$grid
    args[c("grid", "Xtest", "chains")] <- NULL
  }
  set.seed(seed)
  value <- do.call(if (diagnostic) poistree:::ppstree_diag else
    poistree:::ppstree_multi, c(args, list(cache_cuts = enabled)))
  list(value = value, rng = .Random.seed)
}

cut_cache_changed_rules <- function(states) {
  seen <- list()
  for (nodes in states) {
    for (i in seq_len(nrow(nodes))) {
      id <- as.character(nodes[i, 1L])
      if (nodes[i, 2L] < 0) next
      rule <- nodes[i, 2:3]
      if (!is.null(seen[[id]]) && !identical(seen[[id]], rule)) return(TRUE)
      seen[[id]] <- rule
    }
  }
  FALSE
}
