devtools::check(
  pkg = "~/pcvr",
  manual = TRUE,
  cran = TRUE,
  remote = TRUE,
  env_vars = c(NOT_CRAN = "false")
)