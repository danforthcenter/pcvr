if (!interactive()) pdf(NULL)
test_that("reading sv github data as long works", {
  sv <- read.pcv(paste0("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/",
                        "main/pcv4-single-value-traits.csv"), mode = "long")
  #* check read in
  expect_equal(dim(sv), c(77058, 20))
  expect_equal(colnames(sv), c(
    "camera", "imgtype", "zoom", "exposure", "gain", "frame", "rotation",
    "lifter", "timestamp", "id", "barcode", "treatment", "velocity",
    "cartag", "measurementlabel", "other", "image", "sample", "trait",
    "value"
  ))

  #* check bw.time
  sv <- bw.time(sv,
    plantingDelay = 7, phenotype = "area_pixels", cutoff = 10, timeCol = "timestamp",
    group = c("barcode", "rotation"), plot = TRUE
  )
  expect_equal(colnames(sv)[21:23], c("DAS", "DAP", "DAE"))
  expect_equal(head(sv$DAS), 4:9)
  expect_equal(head(sv$DAP), 11:16)
  expect_equal(head(sv$DAE), 0:5)

  expect_true(all(sapply(sv, function(c) sum(is.na(c))) == 0))


  #* check bw.outliers
  sv$genotype <- substr(sv$barcode, 3, 5)
  sv$genotype <- ifelse(sv$genotype == "002", "B73",
    ifelse(sv$genotype == "003", "W605S",
      ifelse(sv$genotype == "004", "MM", "Mo17")
    )
  )
  sv$fertilizer <- substr(sv$barcode, 8, 8)
  sv$fertilizer <- ifelse(sv$fertilizer == "A", "100",
    ifelse(sv$fertilizer == "B", "50", "0")
  )

  svNoOutliers <- suppressWarnings(bw.outliers(
    df = sv, phenotype = "area_pixels", group = c("DAS", "genotype", "fertilizer"),
    cutoff = 3, plot = TRUE
  ))
  pct_removed <- nrow(svNoOutliers$data) / nrow(sv)
  expect_equal(pct_removed, 0.997, tolerance = 0.0015)
  expect_s3_class(svNoOutliers$plot, "ggplot")
  #* check cumulativePheno
  csv <- cumulativePheno(sv,
    phenotypes = c("area_pixels", "height_pixels", "width_pixels"),
    group = c("barcode", "rotation")
  )
  #* check relativeTolerance
  rt <- relativeTolerance(
    df = sv, phenotypes = "area_pixels", grouping = c("genotype", "fertilizer")
  )
  expect_equal(dim(rt), c(9L, 9L))
  expect_equal(dim(csv), c(85620, 26))
  expect_equal(sum(csv[csv[["trait"]] == "height_pixels_csum", "value"]), 10646423)
})
