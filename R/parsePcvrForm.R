#' Helper function for parsing pcvr formulas used in growthSS and downstream functions
#'
#' @param form The pcvr style formula specifying outcome, predictor, individuals, and groups.
#' @param df The data that will be used to fit the model
#'
#' @keywords internal
#' @noRd

.parsePcvrForm <- function(form, df = NULL) {
  #* `parse form argument`
  y <- as.character(form)[2]
  x <- as.character(form)[3]
  if (grepl("\\|", x) && grepl("\\/", x)) { # Y ~ X per id within group
    x3 <- trimws(strsplit(x, "[|]|[/]")[[1]])
    x <- x3[1]
    individual <- x3[2]
    group <- trimws(strsplit(x3[length(x3)], "[+]")[[1]])
    USEGROUP <- TRUE
    USEINDIVIDUAL <- TRUE
  } else if (grepl("[\\]|[|]", x)) { # Y ~ X by group
    x2 <- trimws(strsplit(x, "[|]")[[1]])
    x <- x2[1]
    individual <- "dummyIndividual"
    if (!is.null(df)) {
      df[[individual]] <- "dummyIndividual"
    }
    group <- trimws(strsplit(x2[length(x2)], "[+]")[[1]])
    USEGROUP <- TRUE
    USEINDIVIDUAL <- FALSE
  } else { # Y ~ X
    x2 <- trimws(strsplit(x, "[|]")[[1]])
    x <- x2[1]
    individual <- "dummyIndividual"
    group <- "dummyGroup"
    if (!is.null(df)) {
      df[[individual]] <- "dummyIndividual"
      df[[group]] <- "dummyGroup"
    }
    USEGROUP <- FALSE
    USEINDIVIDUAL <- FALSE
  }
  if (grepl("\\+", x)) {
    x_components <- lapply(strsplit(x, "\\+"), trimws)
    x <- x_components[[1]][1]
    hierarchical_predictor <- x_components[[1]][2]
  } else {
    hierarchical_predictor <- NULL
  }
  if (!is.null(df)) {
    if (length(unique(interaction(df[, group]))) == 1) {
      USEGROUP <- FALSE
    } else {
      USEGROUP <- TRUE
    } # if there is only one group then ignore grouping
    tryCatch(
      {
        df <- df[complete.cases(df[, c(x, y, individual, group, hierarchical_predictor)]), ]
        df <- df[!is.infinite(df[[x]]), ]
        df <- df[!is.infinite(df[[y]]), ]
        df <- df[!is.infinite(df[[hierarchical_predictor]]), ]
        formatted <- .formatNonIntegerTime(df, x,
                                           format = "%Y-%m-%d %H:%M:%S",
                                           index = NULL, digits = 2)
        df <- formatted$data
        x <- formatted$timeCol
      },
      error = function(err) {}
    )
  }
  return(list(
    "y" = y, "x" = x, "individual" = individual, "group" = group,
    "USEG" = USEGROUP, "USEID" = USEINDIVIDUAL, "data" = df,
    "hierarchical_predictor" = hierarchical_predictor
  ))
}
