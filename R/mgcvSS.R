#' Ease of use nlrq starter function for 6 growth model parameterizations
#'
#' Internal to growthSS
#'
#' @examples
#'
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#'
#' ss <- .mgcvSS(model = "gam", form = y ~ time | id / group, df = simdf)
#' names(ss) # formula, df, pcvrForm
#'
#' @keywords internal
#' @noRd

.mgcvSS <- function(model = "gam", form, df) {
  #* `parse form argument`
  parsed_form <- .parsePcvrForm(form, df)
  y <- parsed_form$y
  x <- parsed_form$x
  group <- parsed_form$group
  USEGROUP <- parsed_form$USEG
  if (parsed_form$USEID) {
    message(paste0("Individual is not used with type = 'gam'."))
  }
  df <- parsed_form$data
  df[[paste(group, collapse = ".")]] <- interaction(df[, group])
  group <- paste(group, collapse = ".")
  if (USEGROUP) {
    df[, group] <- lapply(group, function(grp) {
      return(factor(df[[grp]]))
    })
    df[, paste0(group, "_numericLabel")] <- lapply(group, function(grp) {
      return(unclass(df[[grp]]))
    })
  }
  #* `assemble gam formula`
  if (USEGROUP) {
    gam_form <- stats::as.formula(paste0(y, "~0+", group, "+s(", x, ", by=",
                                         paste(group, collapse = ":"), ")"))
  } else {
    gam_form <- stats::as.formula(paste0(y, "~0+s(", x, ")"))
  }
  #* `return list`
  out <- list()
  out[["formula"]] <- gam_form
  out[["df"]] <- df
  out[["pcvrForm"]] <- form
  return(out)
}
