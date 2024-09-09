library(ggplot2)
library(quarto)
devtools::load_all("~/pcvr")
tdf <- do.call(rbind,
               lapply(dir("~/pcvr/tutorials", pattern = ".[Qq]md$", recursive = TRUE, 
                          full.names = TRUE), function(tut) {
                            message(paste0("Rendering ", gsub(".*/", "", tut), " Tutorial"))
                            x <- system.time({
                              quarto::quarto_render(tut)
                            })
                            data.frame(time = x[["elapsed"]], Tutorial = gsub(".*/", "", tut))
                          }))
tdf
if (!interactive()) {pdf(NULL)}
ggplot(tdf, aes(x = 1, y = time, fill = Tutorial)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d()+
  theme_light()