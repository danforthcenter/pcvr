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
