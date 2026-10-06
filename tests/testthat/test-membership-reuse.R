test_that("cached parent membership preserves deep split probabilities", {
  region <- rbind(c(-2, 3), c(4, 7))
  x <- cbind(c(-2, -.9, -.4, .1, .2, .3, .4, .5, 2.5, 3),
              c(4, 6.7, 5, 4.3, 5.4, 5.8, 5.9, 6.2, 6.8, 7))
  # Three successive splits on the first coordinate give a deep parent.
  # The compact basis is exactly zero outside some ancestor supports.
  deep_path <- rbind(c(0, .5, 5, -1), c(0, -.4, 2.5, 1),
                     c(1, 5.9, 3, -1), c(0, .1, .9, 1))
  for (family in 0:2) {
    for (gate in list(c(3, 7), c(700, 1400))) {
      for (parent_path in list(matrix(numeric(), 0L, 4L), deep_path)) {
        prefix <- numeric(nrow(x))
        for (j in seq_len(nrow(parent_path))) {
          g <- parent_path[j, ]
          logs <- membership_gate_logs(x, g[1L] + 1L, g[2L], g[3L],
                                      region, gate, family)
          prefix <- prefix + logs[, if (g[4L] < 0) 1L else 2L]
        }
        proposed <- membership_gate_logs(x, 1L, .3, .4, region, gate, family)
        expected <- proposed + prefix
        den <- pmax(expected[, 1L], expected[, 2L])
        left <- exp(expected[, 1L] - den)
        right <- exp(expected[, 2L] - den)
        probability <- left / (left + right)
        args <- list(X = x, region = region, parent_path = parent_path,
                     axis = 0L, cut = .3, parent_width = .4,
                     gate = gate, gate_family = family)
        cached <- do.call(poistree:::ppstree_side_inspect,
                          c(args, list(cache_geometry = TRUE)))
        direct <- do.call(poistree:::ppstree_side_inspect,
                          c(args, list(cache_geometry = FALSE)))
        expect_identical(cached, direct)
        expect_equal(cached$log_children, expected, tolerance = 1e-12)
        expect_equal(as.numeric(cached$p_left), probability, tolerance = 1e-12)
        if (family != 1L) expect_true(all(is.finite(cached$p_left)))
        if (family == 1L && nrow(parent_path))
          expect_true(any(is.infinite(cached$log_children)))
      }
    }
  }
})

test_that("parent reuse preserves complete seeded tree and diagnostic chains", {
  accepted <- rejected <- rep(FALSE, 3L)
  gate_both <- changed_rule <- FALSE
  for (sampler in c("pcg", "rjmcmc")) {
    for (family in 0:2) {
      for (shared in 0:1) {
        args <- cut_cache_native_args()
        args$gate_family <- family
        args$gate_shared <- shared
        args$pcg <- sampler == "pcg"
        for (diagnostic in c(FALSE, TRUE)) {
          cached <- membership_native_run(args, TRUE, diagnostic)
          direct <- membership_native_run(args, FALSE, diagnostic)
          expect_identical(cached, direct)
          accept <- as.numeric(cached$value$accept)
          accepted <- accepted | accept > 0
          rejected <- rejected | accept < 1
          gate_both <- gate_both || any(cached$value$gate_accept > 0 &
                                        cached$value$gate_accept < 1)
          if (!diagnostic) changed_rule <- changed_rule ||
            cut_cache_changed_rules(cached$value$state_nodes)
        }
      }
    }
  }
  expect_true(all(accepted))
  expect_true(all(rejected))
  expect_true(gate_both)
  expect_true(changed_rule)
})
