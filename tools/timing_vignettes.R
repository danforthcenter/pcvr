devtools::load_all("~/pcvr")

vignette_files <- dir("~/pcvr/vignettes", pattern = ".*.Rmd", full.names = TRUE)
times <- list()
for (file in vignette_files) {
  x <- system.time({
    rmarkdown::render(file)
  })
  nm <- gsub(".*vignettes.", "", file)
  times[[nm]] <- data.frame(test = nm, time = x["elapsed"])
}
out <- do.call(rbind, times)
out
write.csv(out, file = "~/Desktop/timing_vignettes.csv", row.names = FALSE)

out <- read.csv("~/Desktop/timing_vignettes.csv")
library(ggplot2)
ggplot(out, aes(x = 1, y = time, fill = test)) +
  geom_col(position = "stack")+
  scale_fill_viridis_d(option = "plasma")+
  theme_minimal()