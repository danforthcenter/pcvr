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
write.csv(tdf, file = "~/pcvr/tools/tests_timing.csv", row.names = FALSE)

p1 <- ggplot(tdf, aes(x = 1, y = time, fill = Test)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d()+
  theme_light()
ggsave("~/pcvr/tools/tests_timing.png", plot = p1, width = 7, height = 6, dpi = 300, bg = "#ffffff")


# now with CI and CRAN

Sys.setenv("CI" = "true")
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
write.csv(tdf, file = "~/pcvr/tools/tests_timing_ci.csv", row.names = FALSE)

p2 <- ggplot(tdf, aes(x = 1, y = time, fill = Test)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d()+
  theme_light()
ggsave("~/pcvr/tools/tests_timing_ci.png", plot = p2, width = 7, height = 6, dpi = 300, bg = "#ffffff")

Sys.setenv("NOT_CRAN" = "false")
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
write.csv(tdf, file = "~/pcvr/tools/tests_timing_cran.csv", row.names = FALSE)

p3 <- ggplot(tdf, aes(x = 1, y = time, fill = Test)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d()+
  theme_light()
ggsave("~/pcvr/tools/tests_timing_cran.png", plot = p3, width = 7, height = 6, dpi = 300, bg = "#ffffff")
