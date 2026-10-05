# Usage: Rscript pgas-hard-proposal-benchmark.R LIBRARY NEW_OUTPUT_DIRECTORY
# One installed package; 18 fits with paired seeds and randomized method order.
# ESS estimates from these short chains are rough mixing diagnostics, not
# evidence of convergence, posterior accuracy, or general efficiency gains.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop("Usage: Rscript pgas-hard-proposal-benchmark.R LIBRARY NEW_OUTPUT_DIRECTORY")
  }
  benchmark_library <- normalizePath(args[1], mustWork = TRUE)
  .libPaths(c(benchmark_library, .libPaths()))
  library(poistree, lib.loc = benchmark_library)
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("The benchmark requires the coda package for effectiveSize().")
  }
  output_dir <- args[2]
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                                no.. = TRUE))) {
    stop("Use a new or empty output directory to preserve existing results.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  on.exit(writeLines(capture.output(sessionInfo()),
                     file.path(output_dir, "sessionInfo.txt")), add = TRUE)

  cases <- list(
    n400_d2 = list(n = 400L, d = 2L),
    n400_d10 = list(n = 400L, d = 10L),
    n1500_d10 = list(n = 1500L, d = 10L)
  )
  methods <- data.frame(
    method = c("laplace", "hard"),
    proposal_score = c("laplace", "hard"),
    proposal_temperature = c(NA_real_, 0.5),
    proposal_defensive = c(NA_real_, 0.1),
    stringsAsFactors = FALSE
  )
  fit_seeds <- c(911L, 912L, 913L)
  data_seed <- 241L
  order_seed <- 20261005L
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")

  # Generate all execution orders before fitting, independently of fit RNG use.
  set.seed(order_seed)
  blocks <- list()
  for (case_name in names(cases)) {
    for (repeat_id in seq_along(fit_seeds)) {
      method_order <- sample(seq_len(nrow(methods)))
      blocks[[length(blocks) + 1L]] <- data.frame(
        case = case_name, repeat_id = repeat_id,
        seed = fit_seeds[repeat_id], order_in_pair = seq_along(method_order),
        methods[method_order, ], row.names = NULL
      )
    }
  }
  plan <- do.call(rbind, blocks)
  rownames(plan) <- NULL
  plan$run_id <- seq_len(nrow(plan))
  write.csv(plan, file.path(output_dir, "run-plan.csv"), row.names = FALSE)
  saveRDS(list(cases = cases, methods = methods, fit_seeds = fit_seeds,
               data_seed = data_seed, order_seed = order_seed,
               library = benchmark_library,
               package_path = find.package("poistree"),
               package_version = as.character(packageVersion("poistree"))),
          file.path(output_dir, "benchmark-config.rds"))
  writeLines(c(
    "PGAS Laplace versus hard-routing proposal benchmark, one installed package.",
    "Three fixed data sets; three paired seeds per case; one chain per fit.",
    "Each fit uses 500 sweeps, burn=100, thin=1, and retains 400 states.",
    "Hard proposal: proposal_temperature=0.5 and proposal_defensive=0.1.",
    "The existing defensive control is zero in both methods.",
    "Method execution order is randomized within each case/seed block.",
    "Paired seeds do not imply identical RNG consumption or sampled trajectories.",
    "Timing includes ppt_fit and prediction at 12 fixed points, not output I/O or diagnostics.",
    "ESS is coda::effectiveSize for each individual chain and scalar quantity.",
    "ESS/sec uses total fit time, including discarded sweeps; no chains are pooled.",
    "A constant trace is flagged; its ESS and ESS/sec are unavailable.",
    "expanded_nodes is mean distinct expanded nodes over retained sweeps.",
    "Proposal changes can change visited trees and the amount of computation.",
    "These short-chain ESS estimates are rough mixing diagnostics, not a convergence",
    "or posterior-accuracy proof. Means and SDs describe the retained samples only."
  ), file.path(output_dir, "README.txt"))

  extract_traces <- function(fit, expected_draws) {
    values <- list(
      log_likelihood = fit$diagnostics$log_likelihood_trace,
      integrated_intensity = fit$posterior$integrated_intensity_draws,
      leaves = fit$diagnostics$leaf_count_trace,
      gate_x1 = if (!is.null(fit$posterior$state$gate)) {
        fit$posterior$state$gate[, 1L]
      } else NULL
    )
    available <- vapply(values, function(x) !is.null(x), logical(1))
    for (quantity in names(values)) {
      if (is.null(values[[quantity]])) {
        values[[quantity]] <- rep(NA_real_, expected_draws)
      } else if (length(values[[quantity]]) != expected_draws) {
        stop("Unexpected trace length for ", quantity, ".")
      }
    }
    list(values = values, available = available)
  }

  summarize_trace <- function(x, available, elapsed) {
    status <- "ok"
    ess <- NA_real_
    note <- ""
    if (!available) {
      status <- "unavailable"
    } else if (any(!is.finite(x))) {
      status <- "nonfinite"
    } else if (length(x) < 3L) {
      status <- "too_short"
    } else if (stats::sd(x) == 0) {
      status <- "constant"
    } else {
      warnings <- character()
      result <- tryCatch(withCallingHandlers(
        as.numeric(coda::effectiveSize(coda::mcmc(x))),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }), error = function(e) e)
      if (inherits(result, "error")) {
        status <- "ess_error"
        note <- conditionMessage(result)
      } else if (length(result) != 1L || !is.finite(result) || result <= 0) {
        status <- "ess_unavailable"
      } else {
        ess <- result
        if (length(warnings)) {
          status <- "ess_warning"
          note <- paste(warnings, collapse = " | ")
        }
      }
    }
    finite <- available && all(is.finite(x))
    data.frame(
      retained_draws = length(x),
      mean = if (finite) mean(x) else NA_real_,
      sd = if (finite) stats::sd(x) else NA_real_,
      ESS = ess, ESS_per_second = if (elapsed > 0) ess / elapsed else NA_real_,
      status = status, note = note, stringsAsFactors = FALSE
    )
  }

  timing <- list()
  efficiency <- list()
  for (case_name in names(cases)) {
    config <- cases[[case_name]]
    set.seed(data_seed)
    x <- matrix(runif(config$n * config$d), config$n, config$d)
    x[, 1L] <- qbeta(x[, 1L], 2, 5)
    region <- cbind(rep(0, config$d), rep(1, config$d))
    common <- list(
      x = x, region = region, predict_at = x[seq_len(12L), , drop = FALSE],
      gating = "soft", scales = "leaf", sampler = "pgas",
      a = 0.5, b = 0.5 / config$n,
      gate = 12, a_gate = 36, b_gate = 3, sd_gate = 0.07,
      gate_min = 0, gate_structure = "dimension", update_gate = TRUE,
      particles = 8L, chains = 1L, iter = 500L, burn = 100L, thin = 1L,
      max_depth = 4L, alpha = 0.9, eta = 1,
      exact_max = 150L, cut_candidates = 50L,
      label_sweeps = 1L, ancestor_sampling = TRUE, defensive = 0,
      resampling = "node", ess_threshold = 1, allocation = "sequential",
      verbose = FALSE
    )
    saveRDS(common, file.path(output_dir, paste0(case_name, "-inputs.rds")))
    case_plan <- plan[plan$case == case_name, , drop = FALSE]
    expected_draws <- length(seq.int(common$burn, common$iter - 1L, by = common$thin))
    for (i in seq_len(nrow(case_plan))) {
      run <- case_plan[i, , drop = FALSE]
      proposal <- list(proposal_score = run$proposal_score)
      if (run$proposal_score == "hard") {
        proposal$proposal_temperature <- run$proposal_temperature
        proposal$proposal_defensive <- run$proposal_defensive
      }
      fit_args <- c(common, proposal, list(seed = run$seed))
      elapsed <- system.time(
        fit <- do.call(ppt_fit, fit_args), gcFirst = TRUE
      )
      trace_result <- extract_traces(fit, expected_draws)
      traces <- data.frame(
        chain = 1L, retained_draw = seq_len(expected_draws),
        sweep = seq.int(common$burn, common$iter - 1L, by = common$thin) + 1L,
        trace_result$values, check.names = FALSE
      )
      run_efficiency <- do.call(rbind, lapply(names(trace_result$values), function(quantity) {
        data.frame(run, chain = 1L, quantity = quantity,
                   summarize_trace(trace_result$values[[quantity]],
                                   trace_result$available[[quantity]],
                                   elapsed[["elapsed"]]), row.names = NULL)
      }))
      efficiency[[length(efficiency) + 1L]] <- run_efficiency
      timing[[length(timing) + 1L]] <- data.frame(
        run, n = config$n, d = config$d, iter = common$iter,
        burn = common$burn, thin = common$thin, particles = common$particles,
        cuts = common$cut_candidates, max_depth = common$max_depth,
        update_gate = common$update_gate, exact_max = common$exact_max,
        seconds = unname(elapsed[["elapsed"]]),
        user_seconds = unname(elapsed[["user.self"]]),
        system_seconds = unname(elapsed[["sys.self"]]),
        expanded_nodes = fit$diagnostics$expanded_nodes,
        mean_leaves = fit$posterior$mean_leaves,
        gate_acceptance = unname(fit$diagnostics$acceptance[["gate"]]),
        row.names = NULL
      )
      saveRDS(list(run = run, control = fit$control, prior = fit$prior,
                   prediction = fit$prediction, posterior = fit$posterior,
                   diagnostics = fit$diagnostics, traces = traces,
                   trace_available = trace_result$available,
                   efficiency = run_efficiency, rng = .Random.seed, timing = elapsed),
              file.path(output_dir, sprintf("%s-seed%d-%s.rds",
                                             case_name, run$seed, run$method)))
      write.csv(do.call(rbind, timing), file.path(output_dir, "timing.csv"),
                row.names = FALSE)
      write.csv(do.call(rbind, efficiency), file.path(output_dir, "efficiency.csv"),
                row.names = FALSE)
      cat(sprintf("%d/%d %s seed=%d %s: %.3f seconds\n", run$run_id,
                  nrow(plan), case_name, run$seed, run$method, elapsed[["elapsed"]]))
      flush.console()
      rm(fit)
    }
  }
  invisible(list(timing = do.call(rbind, timing),
                 efficiency = do.call(rbind, efficiency)))
}

main()
