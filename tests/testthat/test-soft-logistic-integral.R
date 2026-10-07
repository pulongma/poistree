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

  sppt <- poistree:::ppstree_geometry(
    rep(1L, length(cuts)), cuts, rep(1, length(cuts)), sides,
    points, region, gate
  )

  expect_equal(sppt$H, expected, tolerance = 2e-10)
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
  geometry <- function(cuts, sides, gate) {
    poistree:::ppstree_geometry(
      rep(1L, length(cuts)), cuts, rep(1, length(cuts)), sides,
      points, region, gate
    )
  }

  for (gate in c(2, 10, 30, 500)) {
    parent <- geometry(cuts, sides, gate)
    left <- geometry(c(cuts, split), c(sides, -1L), gate)
    right <- geometry(c(cuts, split), c(sides, 1L), gate)
    expected <- logistic_path_integral(cuts, sides, gate)

    expect_equal(parent$H, expected, tolerance = 2e-9)
    expect_equal(left$H + right$H, parent$H, tolerance = 2e-9)
    expect_equal(left$phi + right$phi, parent$phi, tolerance = 1e-12)
  }
})
