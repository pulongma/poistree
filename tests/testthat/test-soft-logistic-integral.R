logistic_path_value <- function(x, cuts, sides, gate, lower = 0, upper = 1) {
  vapply(x, function(value) {
    right <- stats::plogis(gate * (value - cuts) / (upper - lower))
    prod(ifelse(sides > 0L, right, 1 - right))
  }, numeric(1))
}

logistic_path_integral <- function(cuts, sides, gate,
                                   lower = 0, upper = 1) {
  stats::integrate(
    logistic_path_value,
    lower = lower, upper = upper,
    cuts = cuts, sides = sides, gate = gate,
    rel.tol = 1e-11, abs.tol = 1e-15,
    subdivisions = 5000L
  )$value
}

test_that("deep same-axis logistic exposure agrees with adaptive integration", {
  # This valid depth-14 recursive path is a regression case for catastrophic
  # cancellation in the former simple-pole partial-fraction calculation.  It
  # returned about 0.61 although its exposure is about 1.40e-5.
  cuts <- c(
    0.217256216173992, 0.0717259424743955,
    0.0196324874387363, 0.0322867908044393,
    0.0298023492333524, 0.0212512087672371,
    0.0262508854998823, 0.0271087326184606,
    0.0293653440434951, 0.0289748157170296,
    0.0284457626937788, 0.0280220642046908,
    0.0281667670123337, 0.028122339806845
  )
  sides <- c(
    -1L, -1L, 1L, -1L, -1L, 1L, 1L,
    1L, -1L, -1L, -1L, 1L, -1L, -1L
  )
  gate <- 34.9963777256801
  region <- matrix(c(0, 1), nrow = 1L)
  points <- matrix(c(0.01, 0.03, 0.2), ncol = 1L)
  expected <- logistic_path_integral(cuts, sides, gate)

  mppt <- poistree:::mppstree_logistic_geometry(
    rep(1L, length(cuts)), cuts, sides, points, region, gate
  )
  sppt <- poistree:::ppstree_geometry(
    rep(1L, length(cuts)), cuts, rep(1, length(cuts)), sides,
    points, region, gate
  )

  expect_equal(mppt$H, expected, tolerance = 2e-10)
  expect_equal(sppt$H, expected, tolerance = 2e-10)
  expect_equal(mppt$H, sppt$H, tolerance = 1e-13)
})

test_that("deep logistic child exposures add to their parent", {
  cuts <- c(
    0.20, 0.80, 0.30, 0.70, 0.35, 0.65,
    0.40, 0.60, 0.43, 0.57, 0.45, 0.55
  )
  sides <- rep(c(1L, -1L), 6L)
  split <- 0.50
  region <- matrix(c(0, 1), nrow = 1L)
  points <- matrix(seq(0, 1, length.out = 31L), ncol = 1L)

  for (gate in c(2, 10, 30, 500)) {
    parent <- poistree:::mppstree_logistic_geometry(
      rep(1L, length(cuts)), cuts, sides, points, region, gate
    )
    left <- poistree:::mppstree_logistic_geometry(
      rep(1L, length(cuts) + 1L), c(cuts, split), c(sides, -1L),
      points, region, gate
    )
    right <- poistree:::mppstree_logistic_geometry(
      rep(1L, length(cuts) + 1L), c(cuts, split), c(sides, 1L),
      points, region, gate
    )
    expected <- logistic_path_integral(cuts, sides, gate)

    expect_equal(parent$H, expected, tolerance = 2e-9)
    expect_equal(left$H + right$H, parent$H, tolerance = 2e-9)
    expect_equal(left$phi + right$phi, parent$phi, tolerance = 1e-12)
  }
})
