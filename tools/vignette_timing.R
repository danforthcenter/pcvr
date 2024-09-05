library(ggplot2)
library(rmarkdown)
library(knitr)
devtools::load_all("~/pcvr")
sv <- NULL
tdf <- do.call(rbind,
               lapply(dir("~/pcvr/vignettes", pattern = ".[Rr]md$", full.names = TRUE), function(vig) {
                 message(paste0("Rendering ", gsub(".*/", "", vig), " Vignette"))
                 x <- system.time({
                   rmarkdown::render(vig)
                   })
                 data.frame(time = x[["elapsed"]], Vignette = gsub(".*/", "", vig))
                 }))
tdf
if (!interactive()) {pdf(NULL)}
ggplot(tdf, aes(x = 1, y = time, fill = Vignette)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d()+
  theme_light()
