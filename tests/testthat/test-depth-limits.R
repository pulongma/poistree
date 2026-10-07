depth_test_backends <- function() {
  list(
    list(gating = "hard", sampler = "smc", engine = "shared", particles = 3L),
    list(gating = "hard", sampler = "smc", engine = "dense", particles = 3L),
    list(gating = "hard", sampler = "pgas", particles = 3L,
         chains = 1L, iter = 4L, burn = 1L, verbose = FALSE),
    list(gating = "hard", sampler = "rjmcmc", chains = 1L,
         iter = 4L, burn = 1L, verbose = FALSE),
    list(gating = "hard", sampler = "irjmcmc", chains = 1L,
         iter = 4L, burn = 1L, verbose = FALSE),
    list(gating = "soft", sampler = "rjmcmc", chains = 1L,
         iter = 4L, burn = 1L, update_gate = FALSE, verbose = FALSE),
    list(gating = "soft", sampler = "irjmcmc", chains = 1L,
         iter = 4L, burn = 1L, update_gate = FALSE, verbose = FALSE),
    list(gating = "soft", sampler = "pgas", particles = 3L,
         chains = 1L, iter = 4L, burn = 1L, update_gate = FALSE,
         verbose = FALSE),
    list(gating = "soft", sampler = "pcg", chains = 1L,
         iter = 4L, burn = 1L, update_gate = FALSE, verbose = FALSE)
  )
}

test_that("all fitting backends reject malformed and unsafe maximum depths", {
  x <- matrix(c(0.2, 0.8), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  invalid <- list(NULL, numeric(), c(1, 2), NA_real_, NaN, Inf, -Inf,
                  -1, 1.5, 21, 31, 32, 2^31, 1e100, 2 + 0i, 2 + 1i,
                  TRUE, "2", factor("2"))
  for (backend in depth_test_backends()) {
    for (depth in invalid) {
      args <- c(list(x = x, region = region, max_depth = depth), backend)
      expect_error(do.call(ppt_fit, args), "max_depth")
    }
  }
})

test_that("depth zero remains a valid root-only model where supported", {
  x <- matrix(seq(0.05, 0.95, length.out = 20L), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  for (backend in depth_test_backends()) {
    args <- c(list(x = x, region = region, max_depth = 0L, a = 0.5, b = 1), backend)
    if (backend$gating == "soft" && backend$sampler == "pgas") {
      expect_error(do.call(ppt_fit, args), "max_depth")
    } else {
      fit <- do.call(ppt_fit, args)
      expect_identical(fit$control$max_depth, 0L)
      expect_true(all(fit$diagnostics$max_depth_trace == 0L))
      expect_true(all(fit$diagnostics$leaf_count_trace == 1L))
      expect_true(is.finite(as.numeric(poistree:::ppt_logLik(fit))))
      if (backend$sampler == "smc") {
        expected <- lgamma(0.5 + nrow(x)) - lgamma(0.5) -
          (0.5 + nrow(x)) * log(2)
        expect_equal(fit$posterior$log_target_normalizer, expected,
                     tolerance = 1e-12)
        expect_true(is.na(fit$posterior$log_evidence))
        expect_equal(ppt_lambda(fit, at = x, type = "draws")$draws,
                     fit$prediction$draws, tolerance = 1e-12)
      }
    }
  }
})

depth_native_calls <- function() {
  x <- matrix(c(0.2, 0.8), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  empty <- matrix(numeric(), ncol = 1L)
  smc <- list(pts = x, grid = x, region = region, max_depth = 1,
              P = 3L, min_leaf_n = 1L, resample_thresh = 0.5,
              a = 0.5, b = 0.1, max_aspect_ratio = Inf, cut_grid_n = 3L)
  pg <- c(smc, list(niter = 4L, verbose = FALSE))
  rj <- list(pts = x, grid = x, region = region, max_depth = 1,
             niter = 3L, burnin = 1L)
  soft <- list(X = x, grid = x, Xtest = empty, region = region,
               a = 0.5, b = 0.1, gate = 10, a_gate = 2, b_gate = 0.1,
               sd_gate = 0.1, gate_min = 0, gate_shared = 0L,
               alpha = 0.5, eta = 2, Dmax = 1, nmin = 1L,
               iters = 4L, burn = 1L, thin = 1L, nmove = 1L, ncc = 1L,
               cut_mode = 0L, ncand = 3L, update_gate = 0L,
               gate_family = 0L, chains = 1L, verbose = 0L)
  diag <- soft[!names(soft) %in% c("grid", "Xtest", "chains", "verbose")]
  diag$mon <- x
  soft_pg <- list(X = x, grid = x, Xtest = empty, region = region,
                  cut_grid = list(0.5), a = 0.5, b = 0.1, gate = 10,
                  a_gate = 2, b_gate = 0.1, sd_gate = 0.1, gate_min = 0,
                  gate_shared = FALSE, rho = 0.5, eta = 2, max_depth = 1,
                  P = 3L, niter = 4L, burn = 1L, thin = 1L,
                  label_sweeps = 0L, update_gate = FALSE,
                  ancestor_sampling = TRUE, exact_max = 10L,
                  defensive = 0.01, resample_node = TRUE,
                  ess_threshold = 1, allocation_rates = FALSE, verbose = FALSE)
  hard_inspect <- list(pts = x, region = region,
                       splits = matrix(numeric(), 0L, 3L), max_depth = 1)
  soft_inspect <- list(X = x, region = region,
                       splits = matrix(numeric(), 0L, 3L), labels = c(1L, 1L),
                       gate = 10, a = 0.5, b = 0.1, alpha = 0.5, eta = 2,
                       Dmax = 1, nmin = 1L, cut_mode = 0L, ncand = 3L,
                       gate_family = 0L, kind = 0L)
  list(PPT_fit_MCMC = rj, PPT_fit_IMCMC = rj,
       PPT_IMCMC_transition = hard_inspect,
       PPT_fit_PG = pg, PPT_fit_SMC = smc, PPT_fit_SMC_shared = smc,
       SPPT_fit_PGAS = soft_pg, ppstree_multi = soft,
       ppstree_diag = diag, ppstree_informed_transition = soft_inspect)
}

test_that("native entrypoints validate depth before integer conversion", {
  for (name in names(depth_native_calls())) {
    original <- depth_native_calls()[[name]]
    key <- if ("Dmax" %in% names(original)) "Dmax" else "max_depth"
    fun <- get(name, envir = asNamespace("poistree"))
    for (depth in list(-1, 0.5, NA_real_, NaN, Inf, 21, 31, 32, 2^31, 1e100)) {
      args <- original
      args[[key]] <- depth
      expect_error(do.call(fun, args), "max_depth")
    }
    for (depth in list(numeric(), c(1, 2), "bad")) {
      args <- original
      args[key] <- list(depth)
      expect_error(do.call(fun, args))
    }
  }
})

test_that("eager tree storage is checked before allocating particle arrays", {
  calls <- depth_native_calls()
  for (name in c("PPT_fit_SMC", "PPT_fit_PG", "SPPT_fit_PGAS")) {
    args <- calls[[name]]
    args$max_depth <- 20
    fun <- get(name, envir = asNamespace("poistree"))
    expect_error(do.call(fun, args), "tree storage limit")
  }
  x <- matrix(c(0.2, 0.8), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  for (backend in depth_test_backends()[c(2L, 3L, 8L)]) {
    args <- c(list(x = x, region = region, max_depth = 20L), backend)
    expect_error(do.call(ppt_fit, args), "tree storage limit")
  }

  expect_invisible(poistree:::.ppt_validate_tree_storage(0L, 8388606, "dense"))
  expect_error(poistree:::.ppt_validate_tree_storage(0L, 8388607, "dense"),
               "tree storage limit")
  expect_invisible(poistree:::.ppt_validate_tree_storage(19L, layout = "soft"))
  expect_error(poistree:::.ppt_validate_tree_storage(20L, layout = "soft"),
               "tree storage limit")
})
