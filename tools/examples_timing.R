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
setwd("~/pcvr/tools/")
write.csv(tdf, file = "examples_timing.csv", row.names = FALSE)

p1 <- ggplot(tdf, aes(x = 1, y = time, fill = fun)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d()+
  theme_light()
ggsave("examples_timing.png", plot = p1, width = 7, height = 6, dpi = 300, bg = "#ffffff")