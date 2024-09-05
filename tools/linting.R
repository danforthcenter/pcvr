#* styling for pcvr
library(lintr)
library(styler)
devtools::load_all("~/pcvr") # make sure ggplot syntax is recognized

x <- lintr::lint_package(path = "~/pcvr",
                    linters = linters_with_defaults(line_length_linter(length = 105L),
                                                    object_name_linter(styles = c("snake_case", "symbols",
                                                                                  "camelCase", "dotted.case",
                                                                                  "lowercase", "UPPERCASE")),
                                                    brace_linter(allow_single_line = TRUE)
                    ))
x
length(x)
#* dry run styling
if (FALSE) {
  style_pkg("~/pcvr", dry = "on")
}

if(FALSE){
  file = "~/pcvr/R/emd.R"
  styler::style_file(file, scope = "line_breaks")
  lintr::lint(file, linters = lintr::linters_with_defaults(lintr::line_length_linter(length = 105L),
                                                           lintr::object_name_linter(styles = c("snake_case", "symbols",
                                                                                                "camelCase", "dotted.case",
                                                                                                "lowercase", "UPPERCASE")),
                                                           lintr::brace_linter(allow_single_line = TRUE)
  ))
  #* or for tokens
  styler::style_file(file, dry = "off", scope = "tokens")
}
