devtools::load_all("~/pcvr")
test_files <- dir("~/pcvr/tests/testthat", pattern = "test-.*.R", full.names = TRUE)
times <- list()
for (file in test_files) {
  x <- system.time({
    source(file)
  })
  nm <- gsub(".*test-", "test-", file)
  times[[nm]] <- data.frame(test = nm, time = x["elapsed"])
}
out <- do.call(rbind, times)
out
write.csv(out, file = "~/Desktop/timing_tests.csv", row.names = FALSE)

out <- read.csv("~/Desktop/timing_tests.csv")
library(ggplot2)
ggplot(out, aes(x = 1, y = time, fill = test)) +
  geom_col(position = "stack")+
  scale_fill_viridis_d(option = "plasma")+
  theme_minimal()
