if (!interactive()) pdf(NULL)

test_that("conjugate vonmises method is consistent for SV and MV", {
  set.seed(123)
  method <- "vonmises"
  prior <- list(mu = 0, kappa = 1, known_kappa = 1, boundary = c(0, 180), n = 1)
  generating <- list(
    s1 = list(f = "rnorm", n = 20, mean = 20, sd = 5),
    s2 = list(f = "rnorm", n = 20, mean = 170, sd = 5)
  )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  mu1 <- as.numeric(out$posteriors[c(1, 3), "mu"])
  mu2 <- as.numeric(out$posteriors[c(2, 4), "mu"])
  expect_equal(
    abs(diff(mu1) / mean(mu1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.1 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(mu2) / mean(mu2)),
    0,
    tolerance = 0.1
  )
  k1 <- as.numeric(out$posteriors[c(1, 3), "kappa"])
  k2 <- as.numeric(out$posteriors[c(2, 4), "kappa"])
  expect_equal(
    abs(diff(k1) / mean(k1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.5 # but these have more variability
  )
  expect_equal(
    abs(diff(k2) / mean(k2)),
    0,
    tolerance = 0.5
  )
})

test_that("conjugate vonmises2 method is consistent for SV and MV", {
  set.seed(123)
  method <- "vonmises2"
  prior <- list(mu = 0, kappa = 1, known_kappa = 1, boundary = c(0, 180), n = 1)
  generating <- list(
    s1 = list(f = "rnorm", n = 20, mean = 20, sd = 5),
    s2 = list(f = "rnorm", n = 20, mean = 170, sd = 5)
  )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  mu1 <- as.numeric(out$posteriors[c(1, 3), "mu"])
  mu2 <- as.numeric(out$posteriors[c(2, 4), "mu"])
  expect_equal(
    abs(diff(mu1) / mean(mu1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.1 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(mu2) / mean(mu2)),
    0,
    tolerance = 0.1
  )
  k1 <- as.numeric(out$posteriors[c(1, 3), "kappa"])
  k2 <- as.numeric(out$posteriors[c(2, 4), "kappa"])
  expect_equal(
    abs(diff(k1) / mean(k1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.5 # but these have more variability
  )
  expect_equal(
    abs(diff(k2) / mean(k2)),
    0,
    tolerance = 0.5
  )
})

test_that("conjugate T method is consistent for SV and MV", {
  set.seed(324)
  method <- "t"
  prior <- list(mu = c(50, 50), n = c(1, 1), s2 = c(100, 100))
  generating <- list(
    s1 = list(f = "rnorm", n = 10, mean = 50, sd = 10),
    s2 = list(f = "rnorm", n = 10, mean = 60, sd = 15)
  )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  mu1 <- as.numeric(out$posteriors[c(1, 3), "mu"])
  mu2 <- as.numeric(out$posteriors[c(2, 4), "mu"])
  expect_equal(
    abs(diff(mu1) / mean(mu1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.1 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(mu2) / mean(mu2)),
    0,
    tolerance = 0.1
  )
  v1 <- as.numeric(out$posteriors[c(1, 3), "sd"])
  v2 <- as.numeric(out$posteriors[c(2, 4), "sd"])
  expect_equal(
    abs(diff(v1) / mean(v1)),
    0, # diff should be 0, scaled or not
    tolerance = 1 # but these have a lot of variability from the MV traits seeming wider.
    # worth considering how to adapt the MV method if I don't want this to be the case
  )
  expect_equal(
    abs(diff(v2) / mean(v2)),
    0,
    tolerance = 1
  )
})

test_that("conjugate Gaussian method is consistent for SV and MV", {
  set.seed(324)
  method <- "gaussian"
  prior <- list(mu = c(50, 50), n = c(1, 1), s2 = c(100, 100))
  generating <- list(
    s1 = list(f = "rnorm", n = 10, mean = 50, sd = 10),
    s2 = list(f = "rnorm", n = 10, mean = 60, sd = 15)
  )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  mu1 <- as.numeric(out$posteriors[c(1, 3), "mu"])
  mu2 <- as.numeric(out$posteriors[c(2, 4), "mu"])
  expect_equal(
    abs(diff(mu1) / mean(mu1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.1 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(mu2) / mean(mu2)),
    0,
    tolerance = 0.1
  )
  v1 <- as.numeric(out$posteriors[c(1, 3), "sd"])
  v2 <- as.numeric(out$posteriors[c(2, 4), "sd"])
  expect_equal(
    abs(diff(v1) / mean(v1)),
    0, # diff should be 0, scaled or not
    tolerance = 1 # see note for T distribution
  )
  expect_equal(
    abs(diff(v2) / mean(v2)),
    0,
    tolerance = 1
  )
})


test_that("conjugate Beta method is consistent for SV and MV", {
  # MV traits are stronger and this has very loose tolerances. I'd like to tighten it up, but it is
  # not obvious how much the extra information in MV traits should change this vs how much it does
  # change this. These tests are one step towards making that very stable, but it is not totally done.
  set.seed(324)
  method <- "beta"
  prior <- list(a = c(0.5, 0.5), b = c(0.5, 0.5))
  generating <- list(
    s1 = list(f = "rbeta", n = 50, shape1 = 8, shape2 = 4),
    s2 = list(f = "rbeta", n = 50, shape1 = 7, shape2 = 6)
  )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  out$posteriors
  out$summaries
  a1 <- as.numeric(out$posteriors[c(1, 3), "a"])
  a2 <- as.numeric(out$posteriors[c(2, 4), "a"])
  expect_equal(
    abs(diff(a1) / mean(a1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.5 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(a2) / mean(a2)),
    0,
    tolerance = 0.5
  )
  b1 <- as.numeric(out$posteriors[c(1, 3), "b"])
  b2 <- as.numeric(out$posteriors[c(2, 4), "b"])
  expect_equal(
    abs(diff(b1) / mean(b1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.5
  )
  expect_equal(
    abs(diff(b2) / mean(b2)),
    0,
    tolerance = 0.5
  )
})

test_that("conjugate lognormal method is consistent for SV and MV", {
  set.seed(324)
  method <- "lognormal"
  prior <- list(mu = c(3, 3), sd = c(2, 2))
  generating <- list(
    s1 = list(f = "rlnorm", n = 20, meanlog = 3.5, sdlog = log(3)),
    s2 = list(f = "rlnorm", n = 20, meanlog = 3, sdlog = log(4))
  )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  out$posteriors
  out$summaries
  mu1 <- as.numeric(out$posteriors[c(1, 3), "mu"])
  mu2 <- as.numeric(out$posteriors[c(2, 4), "mu"])
  expect_equal(
    abs(diff(mu1) / mean(mu1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.15 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(mu2) / mean(mu2)),
    0,
    tolerance = 0.15
  )
  # really this next part should test sd, but the sd is still very
  # sensitive to MV vs SV traits
  v1 <- as.numeric(out$posteriors[c(1, 3), "lognormal_sigma"])
  v2 <- as.numeric(out$posteriors[c(2, 4), "lognormal_sigma"])
  expect_equal(
    abs(diff(v1) / mean(v1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.25 # but these have more variability
  )
  expect_equal(
    abs(diff(v2) / mean(v2)),
    0,
    tolerance = 0.25
  )
})

#* for local testing on different distributions to get a better
#* evaluation of any bias that exists between the mv and sv methods
#* Sometimes the apparent bias is more related to how mvSim works
#* than the actual updating though.
if (FALSE) {
  library(brms)
  library(ggplot2)
  res_main <- do.call(rbind, lapply(c(10, 25, 50), function(n) {
    n_df <- do.call(rbind, parallel::mclapply(1:1000, function(i) {
      sv_vm <- rvon_mises(1000 * n, 1, 4)
      mv_vm <- mvSim(
        dists = list(
          rvon_mises = list(mu = 1, kappa = 4)
        ),
        n_samples = n,
        min_bin = -3.14,
        max_bin = 3.14,
        binwidth = 0.01
      )
      vm2_sv <- conjugate(
        s1 = sv_vm,
        method = "vonmises2",
        priors = list(mu = 0, kappa = 0.1, boundary = c(-pi, pi), n = 1),
        cred.int.level = 0.95,
        plot = FALSE
      )

      vm2_mv <- conjugate(
        s1 = mv_vm[, -1],
        method = "vonmises2",
        priors = list(mu = 0, kappa = 0.1, boundary = c(-pi, pi), n = 1),
        cred.int.level = 0.95,
        plot = FALSE
      )

      out <- data.frame(
        data = c(rep(c("sv", "mv"), each = 3), "sv", "mv"),
        quantity = c(rep(c("hde", "hdi_low", "hdi_high"), 2), "hde", "hde"),
        param = c(rep("mu", 6), "kappa", "kappa"),
        value = c(
          as.numeric(vm2_sv$summary),
          as.numeric(vm2_mv$summary),
          vm2_sv$posterior[[1]]$kappa,
          vm2_mv$posterior[[1]]$kappa
        ),
        i = i,
        n = n,
        conjugate = "main"
      )
      return(out)
    }, mc.cores = 10))
    return(n_df)
  }))

  ggplot(
    res_main[res_main$param == "kappa", ],
    aes(x = value, group = data, fill = data)
  ) +
    facet_wrap(~n) +
    geom_histogram(alpha = 0.75, position = "identity") +
    pcv_theme() +
    labs(title = "kappa")

  ggplot(
    res_main[res_main$param == "mu" & res_main$quantity == "hde", ],
    aes(x = value, group = data, fill = data)
  ) +
    facet_wrap(~n) +
    geom_histogram(alpha = 0.75, position = "identity") +
    pcv_theme() +
    labs(title = "mu")
}
