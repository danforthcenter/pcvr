library(testthat)
if (!interactive()) pdf(NULL)
wide <- mvSim(
  dists = list(
    rnorm = list(mean = 5, sd = 1),
    runif = list(min = 1, max = 9)
  ),
  n_samples = 10,
  counts = 100,
  min_bin = 0,
  max_bin = 10
)
wide$image <- seq_len(nrow(wide))
long <- as.data.frame(
  data.table::melt(
    data.table::as.data.table(
      wide
    ),
    id.vars = c("group", "image")
  )
)

test_that("wide emd matrix output works", {
  out <- pcv.emd(wide, cols = "sim_", mat = TRUE, plot = TRUE, raiseError = FALSE)
  expect_equal(names(out), c("data", "plot"))
  out <- pcv.euc(wide, cols = "sim_", mat = TRUE, plot = TRUE, raiseError = FALSE)
  expect_equal(names(out), c("data", "plot"))
})

test_that("long emd works", {
  out <- pcv.emd(long,
    cols = "sim_", trait = "variable",
    mat = TRUE, plot = TRUE, raiseError = FALSE
  )
  expect_equal(names(out), c("data", "plot"))
  out <- pcv.euc(long,
    cols = "sim_", trait = "variable",
    mat = TRUE, plot = TRUE, raiseError = FALSE
  )
  expect_equal(names(out), c("data", "plot"))
})

test_that("helpers throw correct errors", {
  expect_error(euc1d(c(1, 1, 1), c(1, 1)))
  expect_error(emd1d(c(1, 1, 1), c(1, 1)))
  small_data <- data.frame(trait = "x")
  expect_message(.emdRaiseError(raiseError = TRUE, df = small_data, parallel = 1, trait = "trait"))
  med_data <- data.frame(trait = rep("x", 1000))
  expect_warning(.emdRaiseError(raiseError = TRUE, df = med_data, parallel = 1, trait = "trait"))
  large_data <- data.frame(trait = rep("x", 10000))
  expect_error(.emdRaiseError(raiseError = TRUE, df = large_data, parallel = 1, trait = "trait"))
})
