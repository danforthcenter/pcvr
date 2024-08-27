#' Calculate relative tolerance of some phenotype(s) relative to control
#'
#' @description Often in bellwether experiments we are curious about the effect of some
#' treatment vs control. For certain routes in analysing the data this requires considering
#' phenotypes as relative differences compared to a control. Note that the \code{conjugate}
#' function can also be useful in considering the relative tolerance to stress between groups and that
#' growth models are another suggested way to test relative tolerance questions.
#'
#' @param df Dataframe to use, this can be in long or wide format.
#' @param phenotypes A character vector of column names for the phenotypes
#'     that should be compared against control.
#' @param grouping A character vector of column names that identify groups in the data.
#'    These groups will be calibrated separately,
#'    with the exception of the group that identifies a control within the greater hierarchy.
#'    Note that for levels of grouping where the control group does not exist the output will be NA.
#' @param control A column name for the variable to be used to select the control observations.
#'     If left NULL (the default) then this will be taken as the first string in the group argument.
#' @param controlGroup The level of the control variable to compare groups against.
#' @param traitCol Column with phenotype names, defaults to "trait".
#'    This should generally not need to be changed from the default.
#'    If this and valueCol are present in colnames(df) then the data
#'    is assumed to be in long format.
#' @param valueCol Column with phenotype values, defaults to "value".
#'    This should generally not need to be changed from the default.
#' @return A dataframe with relative tolerance columns added.
#' @importFrom stats sd setNames
#' @keywords single-value-trait
#' @examples
#'
#' \donttest{
#' sv <- read.pcv(
#'   "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcv4-single-value-traits.csv",
#'   reader = "fread"
#' )
#' sv$genotype <- substr(sv$barcode, 3, 5)
#' sv$genotype <- ifelse(sv$genotype == "002", "B73",
#'   ifelse(sv$genotype == "003", "W605S",
#'     ifelse(sv$genotype == "004", "MM", "Mo17")
#'   )
#' )
#' sv$fertilizer <- substr(sv$barcode, 8, 8)
#' sv$fertilizer <- ifelse(sv$fertilizer == "A", "100",
#'   ifelse(sv$fertilizer == "B", "50", "0")
#' )
#'
#' sv <- bw.time(sv,
#'   plantingDelay = 0, phenotype = "area_pixels",
#'   cutoff = 10, timeCol = "timestamp", group = c("barcode", "rotation"), plot = FALSE
#' )
#' phenotypes <- colnames(sv)[19:35]
#' phenoForm <- paste0("cbind(", paste0(phenotypes, collapse = ", "), ")")
#' groupForm <- "DAS+DAP+barcode+genotype+fertilizer"
#' form <- as.formula(paste0(phenoForm, "~", groupForm))
#' sv <- aggregate(form, data = sv, mean, na.rm = TRUE)
#' sv <- bw.outliers(sv,
#'   phenotype = "area_pixels",
#'   group = c("DAS", "genotype", "fertilizer"),
#'   plotgroup = c("barcode")
#' )$data
#'
#' pixels_per_cmsq <- 42.5^2 # pixel per cm^2
#' sv$area_cm2 <- sv$area_pixels / pixels_per_cmsq
#' sv$height_cm <- sv$height_pixels / 42.5
#'
#' df <- sv
#' phenotypes <- c("area_cm2", "height_cm")
#' grouping <- c("fertilizer", "genotype", "DAS")
#' controlGroup <- "100"
#' control <- "fertilizer"
#'
#' rt <- relativeTolerance(df, phenotypes, grouping, control, controlGroup)
#' head(rt)
#' sapply(rt, function(c) sum(is.na(c)))
#' }
#'
#' @export
#'
relativeTolerance <- function(df, phenotypes = NULL, grouping = NULL, control = NULL,
                              controlGroup = NULL, traitCol = "trait", valueCol = "value") {
  if (all(c(traitCol, valueCol) %in% colnames(df))) {
    wide <- FALSE
  } else {
    wide <- TRUE
  }
  if (is.null(grouping)) {
    grouping <- control
  }
  if (is.null(control)) {
    control <- grouping[1]
  }
  if (is.null(controlGroup)) {
    controlGroup <- unique(df[[control]])[1]
  }

  if (control %in% grouping) {
    group_no_control <- grouping[grouping != control]
  } else {
    group_no_control <- grouping
  }
  group_no_control_factor <- interaction(df[, group_no_control])
  datsp <- split(x = df, f = group_no_control_factor)

  #* if Z = x/y
  #* x ~ N(mu_1, sd_1)
  #* y ~ N(mu_2, sd_2)
  #* mu_z = mu_1/mu_2
  #* sd_z = sqrt( ((sd_1 / mu_1)^2) + ((sd_2/mu_2)^2) ) * mu_1/mu_2

  #* `Wide`
  if (wide) {
    d2 <- do.call(rbind, lapply(datsp, function(d) {
      d_res <- do.call(rbind, lapply(phenotypes, function(pheno) {
        ctrl_mean <- mean(d[d[[control]] == controlGroup, pheno])
        ctrl_se <- sd(d[d[[control]] == controlGroup, pheno]) /
          length(d[d[[control]] == controlGroup, pheno])
        pheno_res <- do.call(rbind, lapply(setdiff(unique(d[[control]]), controlGroup), function(cg) {
          #* experimental group parameters
          mu_eg <- mean(d[d[[control]] == cg, pheno])
          se_eg <- sd(d[d[[control]] == cg, pheno]) / length(d[d[[control]] == cg, pheno])
          #* relative center and propogated error
          mu_rel <- mu_eg / ctrl_mean
          se_rel <- sqrt(((se_eg / mu_eg)^2) + ((ctrl_se / ctrl_mean)^2)) * mu_rel
          setNames(
            data.frame(
              cg, d[1, group_no_control], pheno, mu_rel, se_rel,
              mu_eg, se_eg, ctrl_mean, ctrl_se
            ),
            c(
              control, group_no_control, "phenotype", "mu_rel", "se_rel",
              "mu_trt", "se_trt", "mu_control", "se_control"
            )
          )
        }))
        pheno_res
      }))
      d_res
    }))
    rownames(d2) <- NULL
  } else { #* `Long`
    d2 <- do.call(rbind, lapply(datsp, function(d) {
      d_res <- do.call(rbind, lapply(phenotypes, function(pheno) {
        ctrl_mean <- mean(d[d[[control]] == controlGroup & d[[traitCol]] == pheno, valueCol])
        ctrl_se <- sd(d[d[[control]] == controlGroup & d[[traitCol]] == pheno, valueCol]) /
          length(d[d[[control]] == controlGroup & d[[traitCol]] == pheno, valueCol])
        pheno_res <- do.call(rbind, lapply(setdiff(unique(d[[control]]), controlGroup), function(cg) {
          #* experimental group parameters
          mu_eg <- mean(d[d[[control]] == cg & d[[traitCol]] == pheno, valueCol])
          se_eg <- sd(d[d[[control]] == cg, pheno]) /
            length(d[d[[control]] == cg & d[[traitCol]] == pheno, valueCol])
          #* relative center and propogated error
          mu_rel <- mu_eg / ctrl_mean
          se_rel <- sqrt(((se_eg / mu_eg)^2) + ((ctrl_se / ctrl_mean)^2)) * mu_rel
          setNames(
            data.frame(
              cg, d[1, group_no_control], pheno, mu_rel, se_rel,
              mu_eg, se_eg, ctrl_mean, ctrl_se
            ),
            c(
              control, group_no_control, "phenotype", "mu_rel", "se_rel",
              "mu_trt", "se_trt", "mu_control", "se_control"
            )
          )
        }))
        pheno_res
      }))
      d_res
    }))
    rownames(d2) <- NULL
  }
  return(d2)
}
