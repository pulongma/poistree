capture_mcmc_progress <- function(args) {
  # RcppProgress writes to the message stream; other fit messages may use stdout.
  output <- capture.output({
    messages <- capture.output({
      fit <- do.call(ppt_fit, args)
      random_seed <- .Random.seed
    }, type = "message")
  })
  list(fit = fit, random_seed = random_seed,
       output = output, messages = messages)
}

for (gating in c("hard", "soft")) {
  for (sampler in c("rjmcmc", "irjmcmc", "pgas",
                    if (gating == "soft") "pcg")) {
    test_that(paste(gating, sampler, "reports progress without changing draws"), {
      args <- list(
        x = matrix(c(0.08, 0.14, 0.31, 0.65, 0.82, 0.93), ncol = 1L),
        region = matrix(c(0, 1), nrow = 1L),
        predict_at = matrix(c(0.2, 0.5, 0.8), ncol = 1L),
        gating = gating, sampler = sampler,
        chains = 2L, iter = 7L, burn = 2L,
        max_depth = 2L, cut_candidates = 3L,
        a = 0.5, b = 0.1, seed = 812L
      )
      if (gating == "soft") {
        # The final iteration is discarded, so thinning must not stall the bar.
        args$thin <- 3L
      } else if (sampler != "pgas") {
        args$prediction_draws <- 2L
      }
      if (sampler == "pgas") args$particles <- 3L

      visible <- capture_mcmc_progress(args)
      quiet <- capture_mcmc_progress(c(args, list(verbose = FALSE)))

      expect_identical(quiet$output, character())
      expect_identical(quiet$messages, character())
      expect_equal(sum(grepl("^0%.*100%$", visible$messages)), 2L)
      expect_equal(sum(grepl("^\\*+\\|$", visible$messages)), 2L)
      progress_text <- paste(c(visible$output, visible$messages), collapse = "\n")
      expect_match(progress_text, "chain 1/2", fixed = TRUE)
      expect_match(progress_text, "chain 2/2", fixed = TRUE)

      expect_identical(visible$fit$prediction, quiet$fit$prediction)
      expect_identical(visible$fit$posterior, quiet$fit$posterior)
      expect_identical(visible$fit$diagnostics, quiet$fit$diagnostics)
      expect_identical(visible$random_seed, quiet$random_seed)
    })
  }
}
