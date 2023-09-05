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
  
  #* check bw.outliers
  sv$genotype = substr(sv$barcode, 3,5)
  sv$genotype = ifelse(sv$genotype == "002", "B73",
                       ifelse(sv$genotype == "003", "W605S",
                              ifelse(sv$genotype == "004", "MM", "Mo17")))
  sv$fertilizer = substr(sv$barcode, 8, 8)
  sv$fertilizer = ifelse(sv$fertilizer == "A", "100",
                         ifelse(sv$fertilizer == "B", "50", "0"))
  sv<-bw.outliers(df = sv, phenotype="area_pixels", group = c("DAS", "genotype", "fertilizer"),
                  cutoff = 3, plot=FALSE)
  expect_equal(dim(sv), c(76788, 25))
  
  #* check cumulativePheno
  csv <- cumulativePheno(sv, phenotypes = c("area_pixels", "height_pixels", "width_pixels"),
                         group = c("barcode", "rotation"))
  expect_equal(dim(csv), c(85320, 26))
  expect_equal(sum(csv[csv[["trait"]]=="height_pixels_csum", "value"]), 16182836)
  #* check pcvBox makes a ggplot
  sv_box <- pcvBox(sv[sv$DAS==15 & sv$trait =="area_pixels", ], x="fertilizer", y="value", compare="0", showPoints = T)
  expect_s3_class(sv_box, "ggplot")
 
  #* pending
  
})
