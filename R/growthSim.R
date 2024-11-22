#' Growth data simulating function
#'
#' @description growthSim can be used to help pick reasonable parameters for common
#'  growth models to use in prior distributions or to simulate data for example models/plots.
#'
#' @param model One of "logistic", "gompertz", "weibull", "frechet", "gumbel", "monomolecular",
#' "exponential", "linear", "power law", "logarithmic", "bragg",
#' "double logistic", or "double gompertz".
#' Alternatively this can be a pseudo formula to generate data from a segmented growth curve by
#' specifying "model1 + model2", see examples and \code{\link{growthSS}}.
#' Decay can be specified by including "decay" as part of the model such as "logistic decay" or
#' "linear + linear decay". Count data can be specified with the "count: " prefix,
#' similar to using "poisson: model" in \link{growthSS}. Similarly intercepts can be added with the
#' "int_" prefix, in which case an "I" parameter should be specified.
#' While "gam" models are supported by \code{growthSS}
#' they are not simulated by this function.
#' @param n Number of individuals to simulate over time per each group in params
#' @param t Max time (assumed to start at 1) to simulate growth to as an integer.
#' @param params A list of numeric parameters. A, B, C notation is used in the order that parameters
#' appear in the formula (see examples). Number of groups is inferred from the length of these vectors
#' of parameters. In the case of the "double" models there are also A2, B2, and C2 terms.
#' Changepoints should be specified as "changePointX" or "fixedChangePointX" as in
#' \code{\link{growthSS}}.
#' @param D If decay is being simulated then this is the starting point for decay. This defaults to 0.
#'
#' @return Returns a dataframe of example growth data following the input parameters.
#'
#' @importFrom stats rnorm setNames
#'
#' @examples
#'
#' library(ggplot2)
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "Logistic")
#'
#' simdf <- growthSim("gompertz",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(0.2, 0.25))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "Gompertz")
#'
#' simdf <- growthSim("weibull",
#'   n = 20, t = 25,
#'   params = list("A" = c(100, 100), "B" = c(1, 0.75), "C" = c(2, 3))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "weibull")
#'
#' simdf <- growthSim("frechet",
#'   n = 20, t = 25,
#'   params = list("A" = c(100, 110), "B" = c(2, 1.5), "C" = c(5, 2))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "frechet")
#'
#' simdf <- growthSim("gumbel",
#'   n = 20, t = 25,
#'   list("A" = c(120, 140), "B" = c(6, 5), "C" = c(4, 3))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "gumbel")
#'
#' simdf <- growthSim("double logistic",
#'   n = 20, t = 70,
#'   params = list(
#'     "A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5),
#'     "A2" = c(400, 300), "B2" = c(35, 40), "C2" = c(3.25, 2.75)
#'   )
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "Double Logistic")
#'
#' simdf <- growthSim("double gompertz",
#'   n = 20, t = 100,
#'   params = list(
#'     "A" = c(180, 140), "B" = c(13, 11), "C" = c(0.2, 0.2),
#'     "A2" = c(400, 300), "B2" = c(50, 50), "C2" = c(0.1, 0.1)
#'   )
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "Double Gompertz")
#'
#' simdf <- growthSim("monomolecular",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(0.08, 0.1))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "Monomolecular")
#'
#' simdf <- growthSim("exponential",
#'   n = 20, t = 25,
#'   params = list("A" = c(15, 20), "B" = c(0.095, 0.095))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "Exponential")
#'
#' simdf <- growthSim("linear",
#'   n = 20, t = 25,
#'   params = list("A" = c(1.1, 0.95))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "Linear")
#'
#' simdf <- growthSim("int_linear",
#'   n = 20, t = 25,
#'   params = list("A" = c(1.1, 0.95), I = c(100, 120))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "Linear with Intercept")
#'
#' simdf <- growthSim("logarithmic",
#'   n = 20, t = 25,
#'   params = list("A" = c(2, 1.7))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "Logarithmic")
#'
#' simdf <- growthSim("power law",
#'   n = 20, t = 25,
#'   params = list("A" = c(16, 11), "B" = c(0.75, 0.7))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "Power Law")
#'
#' simdf <- growthSim("bragg",
#'   n = 20, t = 100,
#'   list("A" = c(10, 15), "B" = c(0.01, 0.02), "C" = c(50, 60))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "bragg")
#'
#' # simulating models from segmented growth models
#'
#' simdf <- growthSim(
#'   model = "linear + linear", n = 20, t = 25,
#'   params = list("linear1A" = c(16, 11), "linear2A" = c(0.75, 0.7), "changePoint1" = c(11, 14))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "linear + linear")
#'
#' simdf <- growthSim(
#'   model = "linear + linear decay", n = 20, t = 25,
#'   params = list("linear1A" = c(16, 11), "linear2A" = c(3, 2), "changePoint1" = c(11, 14))
#' )
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "linear + linear decay")
#'
#' simdf <- growthSim(
#'   model = "linear + linear + logistic", n = 20, t = 50,
#'   params = list(
#'     "linear1A" = c(16, 11), "linear2A" = c(3, 4), # linear slopes, very intuitive
#'     "changePoint1" = c(11, 14), "changePoint2" = c(10, 12),
#'     # changepoint1 is standard, changepoint2 happens relative to changepoint 1
#'     "logistic3A" = c(200, 210), "logistic3B" = c(20, 25), "logistic3C" = c(3, 3)
#'   )
#' )
#' # similar to changepoint2, the asymptote and inflection point are relative to the starting
#' # point of the logistic growth component. This is different than the model output
#' # if you were to fit a curve to this model using `growthSS`.
#' ggplot(simdf, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   labs(title = "linear + linear + logistic")
#'
#' @details
#'     The \code{params} argument requires some understanding of how each growth model is parameterized.
#'     Examples of each are below should help, as will the examples.
#'     \itemize{
#'     \item \bold{Logistic}: `A / (1 + exp( (B-x)/C) )`
#'     Where A is the asymptote, B is the inflection point, C is the growth rate.
#'     \item \bold{Gompertz}: `A * exp(-B * exp(-C*x))`
#'     Where A is the asymptote, B is the inflection point, C is the growth rate.
#'     \item \bold{Weibull}: `A * (1-exp(-(x/C)^B))`
#'     Where A is the asymptote, B is the weibull shape parameter, C is the weibull scale parameter.
#'     \item \bold{Frechet}: `A * exp(-((x-0)/C)^(-B))`
#'     Where A is the asymptote, B is the frechet shape parameter, C is the frechet scale parameter.
#'     Note that the location parameter (conventionally m) is 0 in these models for simplicity but is
#'     still included in the formula.
#'     \item \bold{Gumbel}: `A * exp(-exp(-(x-B)/C))`
#'     Where A is the asymptote, B is the inflection point (location), C is the growth rate (scale).
#'     \item \bold{Double Logistic}: `A / (1+exp((B-x)/C)) + ((A2-A) /(1+exp((B2-x)/C2)))`
#'     Where A is the asymptote, B is the inflection point, C is the growth rate,
#'     A2 is the second asymptote, B2 is the second inflection point, and C2 is the second
#'     growth rate.
#'     \item \bold{Double Gompertz}: `A * exp(-B * exp(-C*x)) + ((A2-A) * exp(-B2 * exp(-C2*(x-B))))`
#'     Where A is the asymptote, B is the inflection point, C is the growth rate,
#'     A2 is the second asymptote, B2 is the second inflection point, and C2 is the second
#'     growth rate.
#'     \item \bold{Monomolecular}: `A-A * exp(-B * x)`
#'     Where A is the asymptote and B is the growth rate.
#'     \item \bold{Exponential}: `A * exp(B * x)`
#'     Where A is the scale parameter and B is the growth rate.
#'     \item \bold{Linear}: `A * x`
#'     Where A is the growth rate.
#'     \item \bold{Logarithmic}: `A * log(x)`
#'     Where A is the growth rate.
#'     \item \bold{Power Law}: `A * x ^ (B)`
#'     Where A is the scale parameter and B is the growth rate.
#'     \item \bold{Bragg}: `A * exp(-B * (x - C) ^ 2)`
#'     This models minima and maxima as a dose-response curve where A is the max response,
#'     B is the "precision" or slope at inflection, and C is the x position of the max response.
#'     \item \bold{Lorentz}: `A / (1 + B * (x - C) ^ 2)`
#'     This models minima and maxima as a dose-response curve where A is the max response,
#'     B is the "precision" or slope at inflection, and C is the x position of the max response.
#'     Generally Bragg is preferred to Lorentz for dose-response curves.
#'     \item \bold{Beta}: `A * (((x - D) / (C - D)) * ((E - x) / (E - C)) ^ ((E - C) / (C - D))) ^ B`
#'     This models minima and maxima as a dose-response curve where A is the Maximum Value,
#'     B is a shape/concavity exponent similar to the sum of alpha and beta in a Beta distribution,
#'     C is the position of maximum value, D is the minimum position where distribution > 0,
#'     E is the maximum position where distribution > 0.
#'     This is a difficult model to fit but can model non-symmetric dose-response relationships which
#'     may sometimes be worth the extra effort.
#'     }
#'     Note that for these distributions parameters generally do not exist in a vacuum.
#'     Changing one will make the others look different in the resulting data.
#'     The examples are a good place to start if you are unsure what parameters to use.
#'
#' @export
#'
#'

growthSim <- function(
    model = c(
      "logistic", "gompertz", "double logistic", "double gompertz",
      "monomolecular", "exponential", "linear", "power law", "frechet",
      "weibull", "gumbel", "logarithmic", "bragg", "lorentz", "beta"
    ),
    n = 20, t = 25, params = list(), D = 0) {
  if (grepl("count:", model)) {
    COUNT <- TRUE
    model <- trimws(gsub("count:", "", model))
  } else {
    COUNT <- FALSE
  }
  int <- FALSE
  if (grepl("^int", model)) {
    int <- TRUE
    model <- trimws(sub("^int_?", "", model))
  }
  if (is.null(names(params))) {
    names(params) <- c(LETTERS[seq_along(params)])
  }
  if (any(names(params) %in% letters)) {
    names(params) <- c(LETTERS[which(letters %in% substr(names(params), 1, 1))])
  }
  params <- as.list(params)
  noise <- lapply(params, function(i) mean(i) / 10)
  names(noise) <- names(params)
  if (any(grepl("fixedChangePoint", names(noise), ignore.case = TRUE))) {
    noise[grepl("fixedChangePoint", names(noise))] <- 0
    nms <- names(noise)
    nms <- sub("fixedC", "c", nms)
    names(noise) <- nms
    nms <- names(params)
    nms <- sub("fixedC", "c", nms)
    names(params) <- nms
  }

  #* check that params are all the same length, if not then rep until they are
  if (!all(unlist(lapply(params, length)) == max(unlist(lapply(params, length))))) {
    message("params are not uniform length, values are being recycled to fit max length")
    diffLengths <- which(!unlist(lapply(params, length)) == max(unlist(lapply(params, length))))
    params[diffLengths] <- lapply(
      diffLengths,
      function(i) {
        rep(params[[i]], length.out = max(unlist(lapply(params, length))))
      }
    )
  }
  #* decide which internal funciton to use
  if (!grepl("\\+", model)) {
    out <- .singleGrowthSim(model, n, t, params, noise, D, int)
  } else {
    out <- .multiGrowthSim(model, n, t, params, noise, D, int)
  }
  if (COUNT) {
    out <- do.call(rbind, lapply(split(out, interaction(out$group, out$id)), function(sub) {
      sub$y <- round(cummax(sub$y))
      sub
    }))
    rownames(out) <- NULL
  }
  return(out)
}

#' Internal helper function to simulate growth from a series of growth models
#' @keywords internal
#' @noRd

.multiGrowthSim <- function(model, n = 20, t = 25, params = list(), noise = NULL, D = 0, int) {
  component_models <- trimws(strsplit(model, "\\+")[[1]])

  firstModel <- component_models[1]
  firstModel <- trimws(firstModel)
  firstModelFindParams <- trimws(gsub("decay", "", firstModel))
  firstParams <- params[grepl(paste0(firstModelFindParams, "1"), names(params))]
  firstChangepoints <- params[["changePoint1"]]
  firstNoise <- noise[grepl(paste0(firstModelFindParams, "1|changePoint1"), names(noise))]
  names(firstNoise) <- sub(paste0(firstModelFindParams, "1|Point."), "", names(firstNoise))

  if (is.null(firstChangepoints)) {
    stop("Simulating segmented data requires 'changePointX' parameters as described in growthSS.")
  }

  df1 <- do.call(rbind, lapply(1:n, function(i) {
    firstChangepointsRand <- lapply(firstChangepoints, function(fc) {
      round(rnorm(1, fc, firstNoise$change))
    })

    n_df <- do.call(rbind, lapply(seq_along(firstChangepointsRand), function(g) {
      .singleGrowthSim(firstModel,
        n = 1, t = firstChangepointsRand[[g]],
        params = stats::setNames(
          lapply(firstParams, function(l) l[[g]]),
          c(sub(paste0(firstModel, "1"), "", names(firstParams)))
        ),
        noise = firstNoise, D
      )
    }))
    n_df$group <- rep(letters[seq_along(firstChangepointsRand)],
      times = unlist(firstChangepointsRand)
    )
    n_df$id <- paste0("id_", i)
    n_df
  }))
  dataList <- list(df1)

  for (u in 2:length(component_models)) {
    iterModel <- component_models[u]
    iterModel <- trimws(iterModel)
    iterModelFindParams <- trimws(gsub("decay", "", iterModel))
    iterParams <- params[grepl(paste0(iterModelFindParams, u), names(params))]

    nextChangepoints <- params[[paste0("changePoint", u)]]
    iterNoise <- noise[grepl(paste0(iterModelFindParams, u, "|changePoint", u), names(noise))]
    names(iterNoise) <- sub(paste0(iterModelFindParams, u, "|Point."), "", names(iterNoise))

    iter_data <- do.call(rbind, lapply(1:n, function(i) {
      if (is.null(nextChangepoints) | u == length(component_models)) {
        iterChangepointsRand <- rep(t, length(iterParams[[1]]))
      } else {
        iterChangepointsRand <- lapply(nextChangepoints, function(fc) {
          round(rnorm(1, fc, iterNoise$change))
        })
      }
      n_df <- do.call(rbind, lapply(seq_along(iterChangepointsRand), function(g) {
        if (u == length(component_models)) {
          grp <- unique(df1$group)[g]
          gt <- t - max(df1[df1$id == paste0("id_", i) & df1$group == grp, "time"])
        } else {
          gt <- iterChangepointsRand[[g]]
        }

        inner_df <- .singleGrowthSim(iterModel,
          n = 1, t = gt,
          params = stats::setNames(
            lapply(iterParams, function(l) l[[g]]),
            c(sub(paste0(iterModelFindParams, u), "", names(iterParams)))
          ),
          noise = iterNoise, D, int
        )
        inner_df$group <- letters[g]
        inner_df
      }))
      n_df$id <- paste0("id_", i)
      n_df
    }))

    prev_data <- dataList[[(u - 1)]]

    new_data <- do.call(rbind, lapply(unique(paste0(iter_data$id, iter_data$group)), function(int) {
      prev_data_sub <- prev_data[paste0(prev_data$id, prev_data$group) == int, ]
      iter_data_sub <- iter_data[paste0(iter_data$id, iter_data$group) == int, ]
      y_end <- prev_data_sub[prev_data_sub$time == max(prev_data_sub$time), "y"]
      iter_data_sub$time <- iter_data_sub$time + max(prev_data_sub$time)
      iter_data_sub$y <- iter_data_sub$y - iter_data_sub$y[1]
      iter_data_sub$y <- y_end + iter_data_sub$y
      iter_data_sub
    }))
    dataList[[(u)]] <- new_data
  }
  out <- do.call(rbind, dataList)
  out <- out[out$time < t, ]
  return(out)
}




#' Internal helper function to simulate growth from a single growth model
#' @keywords internal
#' @noRd

.singleGrowthSim <- function(model, n = 20, t = 25, params = list(), noise = NULL, D, int) {
  models <- c(
    "logistic", "gompertz", "double logistic", "double gompertz",
    "monomolecular", "exponential", "linear", "power law", "frechet", "weibull", "gumbel",
    "logarithmic", "bragg", "lorentz", "beta"
  )

  if (grepl("decay", model)) {
    decay <- TRUE
    model <- trimws(gsub("decay", "", model))
  } else {
    decay <- FALSE
  }

  matched_model <- match.arg(model, models)
  gsi <- match.fun(paste0("gsi_", gsub(" ", "", matched_model)))

  if (decay) {
    gsid <- function(D = 0, ...) {
      D - gsi(...)
    }
  } else {
    gsid <- function(D = 0, ...) {
      0 + gsi(...)
    }
  }

  out <- do.call(rbind, lapply(seq_along(params[[1]]), function(i) {
    pars <- lapply(params, function(p) p[i])
    as.data.frame(rbind(do.call(rbind, lapply(1:n, function(e) {
      iter_data <- data.frame(
        "id" = paste0("id_", e), "group" = letters[i], "time" = 1:t,
        "y" = gsid(D = D, 1:t, pars, noise), stringsAsFactors = FALSE
      )
      if (int) {
        iter_data$y <- iter_data$y + rnorm(1, mean = pars[["I"]], sd = noise[["I"]])
      }
      iter_data
    }))))
  }))

  return(out)
}


#* ************************************************************
#* ***** `gsi functions to simulate individual plants` *****
#* ************************************************************

gsi_logistic <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  c_r <- pars[["C"]] + rnorm(1, mean = 0, sd = noise[["C"]])
  return(a_r / (1 + exp((b_r - x) / c_r)))
}
gsi_gompertz <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  c_r <- pars[["C"]] + rnorm(1, mean = 0, sd = noise[["C"]])
  return(a_r * exp(-b_r * exp(-c_r * x)))
}
gsi_doublelogistic <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  c_r <- pars[["C"]] + rnorm(1, mean = 0, sd = noise[["C"]])
  a2_r <- pars[["A2"]] + rnorm(1, mean = 0, sd = noise[["A2"]])
  b2_r <- pars[["B2"]] + rnorm(1, mean = 0, sd = noise[["B2"]])
  c2_r <- pars[["C2"]] + rnorm(1, mean = 0, sd = noise[["C2"]])
  return(a_r / (1 + exp((b_r - x) / c_r)) + ((a2_r - a_r) / (1 + exp((b2_r - x) / c2_r))))
}
gsi_doublegompertz <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  c_r <- pars[["C"]] + rnorm(1, mean = 0, sd = noise[["C"]])
  a2_r <- pars[["A2"]] + rnorm(1, mean = 0, sd = noise[["A2"]])
  b2_r <- pars[["B2"]] + rnorm(1, mean = 0, sd = noise[["B2"]])
  c2_r <- pars[["C2"]] + rnorm(1, mean = 0, sd = noise[["C2"]])
  return((a_r * exp(-b_r * exp(-c_r * x))) + ((a2_r - a_r) * exp(-b2_r * exp(-c2_r * (x - b_r)))))
}
gsi_monomolecular <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  return(a_r - a_r * exp(-b_r * x))
}
gsi_exponential <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  return(a_r * exp(b_r * x))
}
gsi_linear <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  return(a_r * x)
}
gsi_powerlaw <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  return(a_r * x^(b_r))
}
gsi_logarithmic <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  return(a_r * log(x))
}
gsi_frechet <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- max(c(0, pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])))
  c_r <- max(c(0, pars[["C"]] + rnorm(1, mean = 0, sd = noise[["C"]])))
  # holding location to 0, b is shape parameter, c is scale (growth rate)
  return(a_r * exp(-((x - 0) / c_r)^(-b_r)))
}
gsi_gumbel <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  c_r <- pars[["C"]] + rnorm(1, mean = 0, sd = noise[["C"]])
  # b is location, c is scale (rate)
  return(a_r * exp(-exp(-(x - b_r) / c_r)))
}
gsi_weibull <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  c_r <- pars[["C"]] + rnorm(1, mean = 0, sd = noise[["C"]])
  # c is scale, b is shape
  return(a_r * (1 - exp(-(x / c_r)^b_r)))
}
gsi_bragg <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  c_r <- pars[["C"]] + rnorm(1, mean = 0, sd = noise[["C"]])
  # a is max response, b is precision, c is x position of max response
  return(a_r * exp(-b_r * (x - c_r)^2))
}
gsi_lorentz <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  c_r <- pars[["C"]] + rnorm(1, mean = 0, sd = noise[["C"]])
  # a is max response, b is precision, c is x position of max response
  return(a_r / (1 + b_r * (x - c_r)^2))
}
gsi_beta <- function(x, pars, noise) {
  a_r <- pars[["A"]] + rnorm(1, mean = 0, sd = noise[["A"]])
  b_r <- pars[["B"]] + rnorm(1, mean = 0, sd = noise[["B"]])
  c_r <- pars[["C"]] + rnorm(1, mean = 0, sd = noise[["C"]])
  d_r <- pars[["D"]] + rnorm(1, mean = 0, sd = noise[["D"]])
  e_r <- pars[["E"]] + rnorm(1, mean = 0, sd = noise[["E"]])
  y <- a_r * (((x - d_r) / (c_r - d_r)) * ((e_r - x) / (e_r - c_r))^((e_r - c_r) / (c_r - d_r)))^b_r
  return(y)
}
