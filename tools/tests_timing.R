library(ggplot2)
library(rmarkdown)
library(knitr)
devtools::load_all("~/pcvr")
sv <- NULL
tdf <- do.call(rbind,
               lapply(dir("~/pcvr/tests/testthat", pattern = ".[Rr]$", full.names = TRUE),
                      function(file) {
                 message(paste0("Running ", gsub(".*/", "", file)))
                 x <- system.time({
                   source(file)
                 })
                 data.frame(time = x[["elapsed"]], Test = gsub(".*/", "", file))
               }))
tdf
if (!interactive()) {pdf(NULL)}
ggplot(tdf, aes(x = 1, y = time, fill = Test)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d()+
  theme_light()
