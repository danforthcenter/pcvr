library(covr)

x <- covr::package_coverage(
  path = "~/pcvr",
  type = c("all"),
  combine_types = TRUE,
  quiet = TRUE,
  clean = TRUE,
  pre_clean = TRUE,
  line_exclusions = list("R/plsr.R", "R/vipPlot.R", "R/pcvsubread.R", "R/readpcv3.R")
)
x

report(x, file = "~/pcvr/tools/pcvrCoverage-report.html")