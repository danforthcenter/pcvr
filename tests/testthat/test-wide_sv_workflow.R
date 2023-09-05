test_that("reading sv github data as wide works", {
  
  sv <- read.pcv("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcv4-single-value-traits.csv",
                 mode = "wide")
  #* check read in
  expect_equal(dim(sv), c(2854, 45))
  expect_equal(colnames(sv), c("camera", "imgtype", "zoom", "exposure", "gain", "frame", "rotation", 
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
                               "width_pixels"))
  
  #* check bw.time
   sv <- bw.time(sv, plantingDelay = 7, phenotype = "area_pixels", cutoff = 10, timeCol = "timestamp",
                 group = c("barcode", "rotation"), plot = FALSE)
   expect_equal(colnames(sv)[46:48], c('DAS', 'DAP', 'DAE'))
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
   expect_equal(dim(sv), c(2844, 50))
   
   #* check cumulativePheno
   csv <- cumulativePheno(sv, phenotypes = c("area_pixels", "height_pixels", "width_pixels"),
                          group = c("barcode", "rotation"))
   expect_equal(dim(csv), c(2844, 54))
   expect_equal(sum(csv$height_pixels_csum), 10634085)
   #* check pcvBox makes a ggplot
   sv_box <- pcvBox(sv[sv$DAS==15, ], x="fertilizer", y="area_pixels", compare="0", showPoints = T)
   expect_s3_class(sv_box, "ggplot")
   #* check growthSS (R CMD might throw a fit about brms and my SUGGESTS vs DEPENDS)
   sv$group <- interaction(sv$fertilizer, sv$genotype)
   sv$area_cm2 <- sv$area_pixels / (42.5^2)
   ss <- growthSS(model="gompertz", form =  area_pixels~DAS|barcode/group, sigma="spline", df=sv,
                  priors = list("A" = 130, "B" = 10, "C" = 0.5))
   expect_type(ss, "list")
   
   expect_s3_class(ss[["formula"]], "brmsformula")
   expect_s3_class(ss[["prior"]], "brmsprior")
   expect_type(ss[["initfun"]], "closure")
   expect_s3_class(ss[["df"]], "data.frame")
   expect_type(ss[["family"]], "character")
   expect_s3_class(ss[["pcvrForm"]], "formula")
   
   #* pending
   
  })









