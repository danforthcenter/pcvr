wd <- getwd()
library(ggplot2)
setwd("~/pcvr/R")
devtools::load_all("~/pcvr")
tdf <- do.call(rbind, lapply(dir("~/pcvr/man", pattern = ".[Rr]d$", full.names = TRUE), function(doc) {
  message(paste0("Running ", gsub(".*/", "", doc), " examples"))
  x <- system.time({
    pkgload::run_example(doc, quiet = TRUE)
  })
  data.frame(time = x[["elapsed"]], fun = gsub(".*/", "", doc))
}))
setwd(wd)
tdf
ggplot(tdf, aes(x = 1, y = time, fill = fun)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d()+
  theme_light()
