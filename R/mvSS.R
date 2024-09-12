#' Ease of use multi-value trait model helper function.
#' 
#' This function provides a simplified interface to modeling multi-value traits using \link{growthSS}.
#' Output from this should be passed to \link{fitGrowth} to fit the specified model.
#'
#' @param form A formula similar to \code{label | value ~ time + id/group} where label is a column
#' of histogram bins, value is the counts within those bins, time is an optional time variable,
#' id identifies an individual, and group contains the treatment groups.
#' @param sigma Distributional models passed to \link{growthSS}.
#' @param df Data passed to \link{growthSS}.
#' @param pars Parameters to vary, passed to \link{growthSS}.
#' @param start Starting values or priors, passed to \link{growthSS}.
#' @param type Backend to use, passed to \link{growthSS}.
#' @param tau Quantile to model, passed to \link{growthSS}.
#' @param hierarchy Formulae describing any hierarchical models, see \link{growthSS}.
#' @param spectral_index Optionally, a spectral index
#' \href{https://plantcv.readthedocs.io/en/stable/spectral_index/}{from those calculated by PlantCV}.
#' If this is given then the appropriate truncation and model family (if applicable)
#' will be included for the index you are using without you having to write it in the formula.
#' @keywords multi-value
#' @return A named list of plots showing prior distributions that \code{growthSS} would use,
#' optionally with a plot of simulated growth curves using draws from those priors.
#'
#' @examples
#' set.seed(123)
#' mv_df <- mvSim(dists = list(rnorm = list(mean = 100, sd = 30)), wide = FALSE)
#' mv_df$group <- rep(c("a", "b"), times = 900)
#' mv_df <- mv_df[mv_df$value > 0, ]
#' mv_df$label <- as.numeric(gsub("sim_", "", mv_df$variable))
#' 
#' ss1 <- mvSS(model = "linear", form = label | value ~ group, df = mv_df,
#'             start = list("A" = 5), type = "brms", spectral_index = "ci_rededge")
#' \donttest{
#' mod1 <- fitGrowth(ss1, backend = "cmdstanr", iter = 1000, chains = 1, cores = 1)
#' growthPlot(mod1, ss1$pcvrForm, df = ss1$df)
#' }
#' 
#' # when the model is longitudinal the same model is possible with growthSS
#'
#' m1 <- mvSim(dists = list(rnorm = list(mean = 100, sd = 30),
#'                          rnorm = list(mean = 110, sd = 25),
#'                          rnorm = list(mean = 120, sd = 20),
#'                          rnorm = list(mean = 135, sd = 15)),
#'                          wide = FALSE, n = 6)
#' m1$time <- rep(1:4, times = 6*180)
#' m2 <- mvSim(dists = list(rnorm = list(mean = 85, sd = 25),
#'                          rnorm = list(mean = 95, sd = 20),
#'                          rnorm = list(mean = 105, sd = 15),
#'                          rnorm = list(mean = 110, sd = 15)),
#'                          wide = FALSE, n = 6)
#' m2$time <- rep(1:4, times = 6*180)
#' mv_df2 <- rbind(m1, m2)
#' mv_df2$group <- rep(c("a", "b"), each = 4320)
#' mv_df2 <- mv_df2[mv_df2$value > 0, ]
#' mv_df2$label <- as.numeric(gsub("sim_", "", mv_df2$variable))
#' 
#' ss_mv1 <- mvSS(model = "linear", form = label | value ~ time | group, df = mv_df2,
#'             start = list("A" = 50), type = "brms", spectral_index = "ci_rededge")
#' ss_mv2 <- growthSS(model = "skew_normal: linear",
#'                    form = label | resp_weights(value) + trunc(lb = -1, ub = Inf) ~ time | group,
#'                    df = mv_df2, start = list("A" = 50))
#' identical(names(ss_mv1), names(ss_mv2))
#' # ignoring environments and other such details these are identical except for the
#' # function call.
#' unlist(lapply(names(ss_mv1), function(nm) {
#' if(!identical(ss_mv1[[nm]], ss_mv2[[nm]], ignore.environment = TRUE,
#'               ignore.srcref = TRUE)) {
#'   if(!identical(as.character(ss_mv1[[nm]]), as.character(ss_mv2[[nm]])  )) {
#'     nm
#'   }
#' }
#' }))
#' 
#' \donttest{
#' m2 <- fitGrowth(ss_mv1, backend = "cmdstanr", iter = 1000, chains = 1, cores = 1)
#' growthPlot(m2, ss_mv1$pcvrForm, df = ss_mv1$df)
#' }
#' 
#' @export

mvSS <- function(model = "linear", form, sigma = NULL, df, start = NULL,
                 pars = NULL, type = "brms", tau = 0.5, hierarchy = NULL,
                 spectral_index = c("none", "ari", "ci_rededge", "cri550", "cri700",
                                    "egi", "evi", "gdvi", "mari", "mcari", "mtci", "ndre",
                                    "ndvi", "pri", "psnd_chlorophyll_a", "psnd_chlorophyll_b",
                                    "psnd_caroteniods", "psri", "pssr_chlorophyll_a",
                                    "pssr_chlorophyll_b", "pssr_caroteniods", "rgri",
                                    "rvsi", "savi", "sipi", "sr", "vari", "vi_green", "wi",
                                    "fvfm", "fqfm")) {
  #* `get spectral index helper function`
  #* spectral index function should just return truncation and family, then there will be a separate
  #* function that applies those changes to the form argument.
  pcvrForm <- form
  spec_helper <- get(paste0(".", spectral_index, "_mvSS_helper"))
  #* `run spectral index helper function`
  spec_helper_res <- spec_helper()
  family <- spec_helper_res$family
  trunc <- spec_helper_res$trunc
  #* `run formula cleaner`
  #* This should return the final, usable formula and model
  form_res <- .mvSS_formula_builder(form, family, trunc, type, model, df)
  form <- form_res$formula
  weights <- form_res$weights
  model <- form_res$model
  df <- form_res$df
  has_x_var <- form_res$has_x_var
  if (has_x_var) {
    #* `if time is a variable, call growthSS with new model formula`
    out <- growthSS(model = model, form = form, sigma = sigma, df = df, start = start,
                    pars = pars, type = type, tau = tau, hierarchy = hierarchy)
  } else {
    #* `if x variable is missing then call mvSS_helpers to make simple model`
    #* this is the tricky part I think. but on the other hand if time is missing then it's
    #* really only one option for what happens next, so maybe it's not bad. Model is basically
    #* ignored in this case.
    form_fun <- get(paste0(".", type, "_mvSS"))
    out <- form_fun(form, df, start, family, model, tau, weights, pcvrForm)
  }
  out$call <- match.call()
  return(out)
}


#' `mvSS_formula_builder`
#' add family, truncation, and weights to formula
#' @keywords internal
#' @noRd

.mvSS_formula_builder <- function(form, family, trunc, type, model, df) {
  form_char <- as.character(form)
  has_x_var <- TRUE
  parsed <- .parsePcvrForm(form, df)
  df <- parsed$data
  if (!is.numeric(df[, parsed$x]) & !parsed$USEG & !parsed$USEID) {
    has_x_var <- FALSE # treating single element of RHS as grouping
  }
  weights <- NULL
  #* for brms backend make weights in formula
  if (type == "brms") {
    form_char[2] <- paste0(gsub("[|]", "| resp_weights(", form_char[2]), ")")
    if (any(!is.infinite(trunc))) {
      form_char[2] <- paste0(form_char[2], " + trunc(lb = ", trunc[1],", ub = ", trunc[2], ")")
    }
    if (!grepl("[:]", model)) {
      model <- paste0(family, ": ", model)
    }
  } else {
    #* for other backends save weights variable
    lhs <- trimws(strsplit(form_char[2], "[|]")[[1]])
    weights <- lhs[2]
    form_char[2] <- lhs[1]
  }
  parsed_form <- as.formula(paste0(form_char[c(2, 1, 3)], collapse = ""))
  return(list("formula" = parsed_form,
              "weights" = weights,
              "model" = model,
              "df" = df,
              "has_x_var" = has_x_var))
}

#' `mvSS specified simple model helper functions`

#' @keywords internal
#' @noRd

.brms_mvSS <- function(form = NULL, df = NULL, start = NULL, family = NULL, model = NULL, tau = NULL,
                       weights = NULL, pcvrform = NULL) {
  out <- list()
  #* `Make bayesian non-linear formula`
  bf1 <- as.formula(paste0(as.character(form)[2], "~ A"))
  bf2 <- as.formula(paste0("A ~ 0 + ", as.character(form)[3]))
  
  out[["formula"]] <- brms::bf(bf1, bf2, nl = TRUE)
  out[["prior"]] <- .makePriors(priors = start,
                                pars = "A", df = df,
                                group = "dummyGroup", # group from parse pcvr form should be dummy
                                USEGROUP = FALSE,
                                sigma = FALSE, family = family,
                                formula = out[["formula"]])
  out[["initfun"]] <- 0 # no fancy initialization here
  out[["df"]] <- df
  out[["family"]] <- family
  out[["pcvrForm"]] <- form
  out[["type"]] <- "brms"
  out[["model"]] <- trimws(gsub(".*:", "", model))
  return(out)
}

#' @keywords internal
#' @noRd

.nls_mvSS <- function(form = NULL, df = NULL, start = NULL, family = NULL, model = NULL, tau = NULL,
                      weights = NULL, pcvrForm = NULL) {
  out <- list()
  out[["formula"]] <- form
  out[["start"]] <- NULL
  out[["df"]] <- df
  out[["pcvrForm"]] <- pcvrForm
  out[["type"]] <- "lm"
  out[["model"]] <- trimws(gsub(".*:", "", model))
  out[["weights"]] <- df[[weights]]
  return(out)
}

#' @keywords internal
#' @noRd

.nlrq_mvSS <- function(form = NULL, df = NULL, start = NULL, family = NULL, model = NULL, tau = NULL,
                       weights = NULL, pcvrform = NULL) {
  out <- list()
  out[["formula"]] <- form
  out[["start"]] <- NULL
  out[["df"]] <- df
  out[["pcvrForm"]] <- pcvrform
  out[["type"]] <- "rq"
  out[["model"]] <- trimws(gsub(".*:", "", model))
  out[["taus"]] <- tau
  out[["weights"]] <- df[[weights]]
  return(out)
}

#' `Spectral Index helpers`

#' @keywords internal
#' @noRd
.none_mvSS_helper <- function() {
  truncation <- c(Inf, Inf) # unbounded
  family <- "student" # if not specified then leave as student T per default
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.ari_mvSS_helper <- function() {
  truncation <- c(Inf, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.ci_rededge_mvSS_helper <- function() {
  truncation <- c(-1, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.cri550_mvSS_helper <- function() {
  truncation <- c(Inf, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.cri700_mvSS_helper <- function() {
  truncation <- c(Inf, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.egi_mvSS_helper <- function() {
  truncation <- c(-1, 2)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.evi_mvSS_helper <- function() {
  truncation <- c(Inf, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.gdvi_mvSS_helper <- function() {
  truncation <- c(-2, 2)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.mari_mvSS_helper <- function() {
  truncation <- c(Inf, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.mcari_mvSS_helper <- function() {
  truncation <- c(Inf, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.mtci_mvSS_helper <- function() {
  truncation <- c(Inf, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.ndre_mvSS_helper <- function() {
  truncation <- c(-1, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.ndvi_mvSS_helper <- function() {
  truncation <- c(-1, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.pri_mvSS_helper <- function() {
  truncation <- c(-1, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.psnd_chlorophyll_a_mvSS_helper <- function() {
  truncation <- c(-1, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.psnd_chlorophyll_b_mvSS_helper <- function() {
  truncation <- c(-1, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.psnd_caroteniods_mvSS_helper <- function() {
  truncation <- c(-1, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.psri_mvSS_helper <- function() {
  truncation <- c(Inf, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.pssr_chlorophyll_a_mvSS_helper <- function() {
  truncation <- c(-1, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.pssr_chlorophyll_b_mvSS_helper <- function() {
  truncation <- c(-1, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.pssr_caroteniods_mvSS_helper <- function() {
  truncation <- c(-1, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.rgri_mvSS_helper <- function() {
  truncation <- c(0, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.rvsi_mvSS_helper <- function() {
  truncation <- c(-1, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.savi_mvSS_helper <- function() {
  truncation <- c(-1.2, 1.2)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.sipi_mvSS_helper <- function() {
  truncation <- c(Inf, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.sr_mvSS_helper <- function() {
  truncation <- c(0, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.vari_mvSS_helper <- function() {
  truncation <- c(Inf, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.vi_green_mvSS_helper <- function() {
  truncation <- c(-1, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.wi_mvSS_helper <- function() {
  truncation <- c(0, Inf)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.fvfm_mvSS_helper <- function() {
  truncation <- c(0, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}
#' @keywords internal
#' @noRd
.fqfm_mvSS_helper <- function() {
  truncation <- c(0, 1)
  family <- "skew_normal"
  return(list("trunc" = truncation, "family" = family))
}


