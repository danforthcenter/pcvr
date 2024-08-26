devtools::load_all("~/pcvr")
roxygen2::roxygenize("~/pcvr")
pdf(NULL)
docs_files <- sort(dir(path = "~/pcvr/man", pattern = "\\.[Rr]d$", full.names=TRUE))
times <- list()
for (file in docs_files) {
  x <- system.time({
    suppressMessages(pkgload::run_example(file, quiet = TRUE))
  })
  nm <- gsub(".*man/", "", file)
  times[[nm]] <- data.frame(test = nm, time = x["elapsed"])
}
out <- do.call(rbind, times)
out
write.csv(out, file = "~/Desktop/timing_examples.csv", row.names = FALSE)

out <- read.csv("~/Desktop/timing_examples.csv")
library(ggplot2)
ggplot(out, aes(x = 1, y = time, fill = test)) +
  geom_col(position = "stack")+
  scale_fill_viridis_d(option = "plasma")+
  theme_minimal()
dev.off()