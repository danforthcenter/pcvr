test_that("reading sv github data as long works", {
  
  sv <- read.pcv("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcv4-single-value-traits.csv",
                 mode = "long")
  #* check read in
  expect_equal(dim(sv), c(77058, 20))
  expect_equal(colnames(sv), c("camera", "imgtype", "zoom", "exposure", "gain", "frame", "rotation", 
                               "lifter", "timestamp", "id", "barcode", "treatment", "velocity", 
                               "cartag", "measurementlabel", "other", "image", "sample", "trait", 
                               "value"))

  #* check bw.time
  sv <- bw.time(sv, plantingDelay = 7, phenotype = "area_pixels", cutoff = 10, timeCol = "timestamp",
                group = c("barcode", "rotation"), plot = FALSE)
  expect_equal(colnames(sv)[21:23], c('DAS', 'DAP', 'DAE'))
  expect_equal(head(sv$DAS), 4:9)
  expect_equal(head(sv$DAP), 11:16)
  expect_equal(head(sv$DAE), 0:5)

  expect_true(all(sapply(sv, function(c) sum(is.na(c)))==0))
  
  
  #* check bw.outliers
  sv$genotype = substr(sv$barcode, 3,5)
  sv$genotype = ifelse(sv$genotype == "002", "B73",
                       ifelse(sv$genotype == "003", "W605S",
                              ifelse(sv$genotype == "004", "MM", "Mo17")))
  sv$fertilizer = substr(sv$barcode, 8, 8)
  sv$fertilizer = ifelse(sv$fertilizer == "A", "100",
                         ifelse(sv$fertilizer == "B", "50", "0"))
  
  #* component test
  #* see notes from 9/6/23 for details of what I was working on here and why the bw.outliers output isn't passed on
  #* 
  # subdf <- sv[complete.cases( sv[sv[["trait"]]=="area_pixels", c("value", "trait", "DAS", "genotype", "fertilizer")] ) & sv[["trait"]]=="area_pixels" , ]
  # cooksd <- cooks.distance(glm(data=subdf, as.numeric(value) ~ as.factor(DAS):as.factor(genotype):as.factor(fertilizer)))
  # 
  # expect_equal(quantile(cooksd, seq(0,1, 0.1), na.rm=TRUE), c(`0%` = 5.52344863051265e-34, `10%` = 4.18408115603277e-08, 
  #                                                            `20%` = 2.04530194702589e-07, `30%` = 6.06184338605635e-07, `40%` = 1.49581206806256e-06, 
  #                                                            `50%` = 3.45942904424853e-06, `60%` = 7.15936760830706e-06, `70%` = 1.44450407789928e-05, 
  #                                                            `80%` = 3.12378355203182e-05, `90%` = 8.31906611451925e-05, `100%` = 2.45717262428721
  # ))
  # expect_equal(sum(is.na(cooksd)), 7)
  # 
  # expect_equal(head(cooksd), c(`1` = 1.25816970417831e-07, `2` = 8.69123109207279e-07, `3` = 6.52481789994869e-07, 
  #                              `4` = 1.69313226500173e-07, `5` = 6.1543502355451e-07, `6` = 9.13091057791556e-09))
  
  #* full outliers test
  
  sv_noOutliers<-bw.outliers(df = sv, phenotype="area_pixels", group = c("DAS", "genotype", "fertilizer"),
                  cutoff = 3, plot=FALSE)
  pct_removed <- nrow(sv_noOutliers)/nrow(sv)
  expect_equal( pct_removed , 0.997, tolerance = 0.0015 )

  #* check cumulativePheno
  csv <- cumulativePheno(sv, phenotypes = c("area_pixels", "height_pixels", "width_pixels"),
                         group = c("barcode", "rotation"))
  expect_equal(dim(csv), c(85620, 26))
  expect_equal(sum(csv[csv[["trait"]]=="height_pixels_csum", "value"]), 10646423)

  #* check pcvBox makes a ggplot
  sv_box <- pcvBox(sv[sv$DAS==15 & sv$trait =="area_pixels", ], x="fertilizer", y="value", compare="0", showPoints = T)
  expect_s3_class(sv_box, "ggplot")
})
