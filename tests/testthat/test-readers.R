library(pcvr)
library(testthat)

test_that(".readpcvCheckDataState works", {
  df <- data.frame(
    "trait" = "dummy",
    "value" = 1,
    "label" = "pixel"
  )
  out1 <- .readpcvCheckDataState(
    df,
    traitCol = "trait",
    valueCol = "value",
    labelCol = "label",
    mode = NULL
  )
  expect_equal(out1, list("startsLong" = TRUE, "outputMode" = "long"))
  expect_warning(
    out2 <- .readpcvCheckDataState(
      df,
      traitCol = "trait",
      valueCol = "value2",
      labelCol = "label",
      mode = NULL
    )
  )
  expect_equal(out2, list("startsLong" = FALSE, "outputMode" = "wide"))
})

test_that("read pcv raises error", {
  link1 <- "https://gist.githubusercontent.com/seankross/"
  link2 <- "a412dfbd88b3db70b74b/raw/5f23f993cd87c283ce766e7ac6b329ee7cc2e1d1/mtcars.csv"
  file <- paste0(link1, link2)
  expect_error(suppressWarnings(read.pcv(file, awk = list("gear in 4, 3"))))
})
