if (!interactive()) pdf(NULL)
test_that("reading mv github data as long works", {
  skip_on_cran()
  library(data.table)
  mv <- read.pcv(paste0(
    "https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/",
    "main/pcv4-multi-value-traits.csv"
  ), mode = "wide", reader = "fread")
  expect_equal(dim(mv), c(2854, 198))

  mv$genotype <- substr(mv$barcode, 3, 5)
  mv$genotype <- ifelse(mv$genotype == "002", "B73",
    ifelse(mv$genotype == "003", "W605S",
      ifelse(mv$genotype == "004", "MM", "Mo17")
    )
  )
  mv$fertilizer <- substr(mv$barcode, 8, 8)
  mv$fertilizer <- ifelse(mv$fertilizer == "A", "100",
    ifelse(mv$fertilizer == "B", "50", "0")
  )
  # test bw.time
  mv <- bw.time(mv, timeCol = "timestamp", group = "barcode", plot = FALSE)
  expect_equal(dim(mv), c(2854, 201))
  expect_equal(colnames(mv)[201], "DAS")

  # test bw.outliers

  phenotypes <- which(grepl("hue_freq", colnames(mv)))

  mvNoOutliers <- suppressWarnings(bw.outliers(
    df = mv, phenotype = phenotypes, naTo0 = FALSE, plot = FALSE,
    group = c("DAS", "genotype", "fertilizer"), cutoff = 3, plotgroup = c("barcode", "rotation")
  ))

  pct_removed <- nrow(mvNoOutliers) / nrow(mv)
  expect_equal(pct_removed, 0.93, tolerance = 0.015)

  mvNoOutliers <- suppressWarnings(bw.outliers(
    df = mv, phenotype = phenotypes, naTo0 = FALSE, plot = TRUE, outlierMethod = "mahalanobis",
    group = c("DAS", "genotype", "fertilizer"), cutoff = 3, plotgroup = c("barcode", "rotation")
  ))

  pct_removed <- nrow(mvNoOutliers$data) / nrow(mv)
  expect_equal(pct_removed, 0.945, tolerance = 0.015)
  expect_s3_class(mvNoOutliers$plot, "ggplot")
  #* test joyplot
  joyplot <- pcv.joyplot(mv[mv$DAS == 18, ],
    index = "hue_frequencies",
    group = c("fertilizer", "genotype")
  )
  expect_s3_class(joyplot, "ggplot")

  #* test mv_ag
  set.seed(123)
  mv$svt <- rnorm(nrow(mv))
  mv_ag1 <- mv_ag(
    df = mv,
    group = c("DAS", "genotype", "fertilizer"),
    n_per_group = 2,
    keep = "svt"
  )
  expect_error(
    mv_ag2 <- mv_ag(
      df = mv,
      group = c("DAS", "genotype", "fertilizer"),
      n_per_group = 2,
      labelCol = "camera" # if only some of the id columns are there then error should be thrown
    )
  )

  expect_equal(dim(mv_ag1), c(460, 184))

  #* test EMD
  images <- unique(mv$image)[1:10]
  emd <- pcv.emd(
    df = mv[mv$image %in% images, ], cols = "hue_frequencies", reorder = c("fertilizer", "genotype"),
    mat = FALSE, plot = TRUE, parallel = 1, raiseError = FALSE
  )
  expect_s3_class(emd$plot, "ggplot")
  expect_equal(dim(emd$data), c(110, 7))
  expect_equal(sum(emd$data$emd), 5684.034, tolerance = 0.01)

  #* test network
  net <- pcv.net(emd$data, filter = 0.05)
  expect_s3_class(net[[3]], "igraph")
  expect_equal(dim(net$nodes), c(10, 12))
  expect_equal(dim(net$edges), c(30, 11))
})
