if (!interactive()) pdf(NULL)
test_that("reading sv github data as wide works", {
  skip_if_offline(host = "r-project.org")
  sv <- read.pcv(paste0(
    "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/",
    "main/pcv4-single-value-traits.csv"
  ), mode = "wide")
  #* check read in
  expect_equal(dim(sv), c(2854, 45))
  expect_equal(colnames(sv), c(
    "camera", "imgtype", "zoom", "exposure", "gain", "frame", "rotation",
    "lifter", "timestamp", "id", "barcode", "treatment", "velocity",
    "cartag", "measurementlabel", "other", "image", "sample", "area_pixels",
    "area_above_reference_pixels", "area_below_reference_pixels",
    "color_chip_size_median", "convex_hull_area_pixels", "convex_hull_vertices_none",
    "ellipse_angle_degrees", "ellipse_eccentricity_none", "ellipse_major_axis_pixels",
    "ellipse_minor_axis_pixels", "height_pixels", "height_above_reference_pixels",
    "height_below_reference_pixels", "horizontal_reference_position_none",
    "hue_circular_mean_degrees", "hue_circular_std_degrees", "hue_median_degrees",
    "in_bounds_none", "longest_path_pixels", "median_color_chip_height_median",
    "median_color_chip_width_median", "object_in_frame_none", "percent_area_above_reference_none",
    "percent_area_below_reference_none", "perimeter_pixels", "solidity_none",
    "width_pixels"
  ))

  #* check pcv.time
  sv <- pcv.time(sv,
    plantingDelay = 7, phenotype = "area_pixels", cutoff = 10, timeCol = "timestamp",
    group = c("barcode", "rotation"), plot = TRUE
  )$data
  expect_equal(colnames(sv)[46:48], c("DAS", "DAP", "DAE"))
  expect_equal(head(sv$DAS), 4:9)
  expect_equal(head(sv$DAP), 11:16)
  expect_equal(head(sv$DAE), 0:5)

  #* check pcv.outliers
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

  #* see notes from 9/6/2023 on why this is done differently, also see test-long_sv_workflow.R
  svNoOutliers <- suppressWarnings(pcv.outliers(
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
  expect_equal(dim(csv), c(2854, 54))
  expect_equal(sum(csv$height_pixels_csum), 10646423)
  #* relative tolerance
  rt <- relativeTolerance(
    df = sv, phenotypes = "area_pixels", grouping = c("genotype", "fertilizer")
  )
  expect_equal(dim(rt), c(9L, 9L))
  #* check growthSS (R CMD might throw a fit about brms and my SUGGESTS vs DEPENDS)
  sv$group <- interaction(sv$fertilizer, sv$genotype)
  sv$area_cm2 <- sv$area_pixels / (42.5^2)
  ss <- growthSS(
    model = "gompertz", form = area_pixels ~ DAS | barcode / group, sigma = "spline", df = sv,
    start = list("A" = 130, "B" = 10, "C" = 0.5)
  )
  expect_type(ss, "list")

  expect_s3_class(ss[["formula"]], "brmsformula")
  expect_s3_class(ss[["prior"]], "brmsprior")
  expect_type(ss[["initfun"]], "closure")
  expect_s3_class(ss[["df"]], "data.frame")
  expect_type(ss[["family"]], "character")
  expect_s3_class(ss[["pcvrForm"]], "formula")

  #* pending
})
