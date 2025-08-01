#* styling for pcvr
library(lintr)
library(styler)
devtools::load_all("~/pcvr") # make sure ggplot syntax is recognized

x <- lintr::lint_package(path = "~/pcvr",
                    linters = linters_with_defaults(line_length_linter(length = 105L),
                                                    object_name_linter(styles = c("snake_case", "symbols",
                                                                                  "camelCase", "dotted.case",
                                                                                  "lowercase", "UPPERCASE")),
                                                    return_linter(return_style = "explicit"),
                                                    brace_linter(allow_single_line = TRUE)
                    ))
x
length(x)
#* run styling
if (FALSE) {
  style_pkg("~/pcvr", dry = "off", scope = "line_breaks")
}


x <- lintr::lint(filename = "~/pcvr/vignettes/articles/field_capacity.Rmd",
                         linters = linters_with_defaults(line_length_linter(length = 105L),
                                                         object_name_linter(styles = c("snake_case", "symbols",
                                                                                       "camelCase", "dotted.case",
                                                                                       "lowercase", "UPPERCASE")),
                                                         return_linter(return_style = "explicit"),
                                                         brace_linter(allow_single_line = TRUE)
                         ))
x
styler::style_dir(path = "~/pcvr/tutorials", scope = "tokens")




if(FALSE){
  dput(unique(names(x)))
  c("R/brms_segmentedForm.R", "R/bwoutliers.R", "R/bwtime.R", "R/nlmeSS.R", 
    "tests/testthat/test-brmsModels.R", "tests/testthat/test-growthModels.R", 
    "tests/testthat/test-growthSS_helpers.R")
  file = "~/pcvr/tests/testthat/test-growthSS_helpers.R"
  styler::style_file(file, scope = "line_breaks")
  lintr::lint(file, linters = lintr::linters_with_defaults(lintr::line_length_linter(length = 105L),
                                                           lintr::object_name_linter(styles = c("snake_case", "symbols",
                                                                                                "camelCase", "dotted.case",
                                                                                                "lowercase", "UPPERCASE")),
                                                           lintr::return_linter(return_style = "explicit"),
                                                           lintr::brace_linter(allow_single_line = TRUE)
  ))
  #* or for tokens
  styler::style_file(file, dry = "off", scope = "tokens")
}
