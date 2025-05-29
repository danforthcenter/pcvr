#* `FREM figure script`
setwd("~/pcvr/papers/pcvr_v1.x/")
library(devtools)
devtools::load_all("~/pcvr/")

#* `Simulate a dataset`

set.seed(123)
area <- growthSim("logistic",
                  n = 20, t = 25,
                  params = list("A" = c(200, 190, 175, 150),
                                "B" = c(13, 11, 12, 11),
                                "C" = c(3, 3.5, 3, 3.25))
                  )
height <- growthSim("monomolecular",
                  n = 20, t = 25,
                  params = list("A" = c(45, 40, 37, 35),
                                "B" = c(0.1, 0.15))
                  )
width <- growthSim("monomolecular",
                    n = 20, t = 25,
                    params = list("A" = c(45, 45, 44, 43),
                                  "B" = c(0.2, 0.2, 0.2, 0.18))
)
leaf_count <- growthSim("count:logistic",
                   n = 20, t = 25,
                   params = list("A" = c(8, 7, 6, 6),
                                 "B" = c(13, 13, 13, 13),
                                 "C" = c(3, 3, 3, 3))
)
hue <- growthSim("logistic",
                 n = 20, t = 25,
                 params = list("A" = c(16, 15, 12, 11),
                               "B" = c(13, 13, 13, 13),
                               "C" = c(3, 3, 3, 3))
)
hue$y <- hue$y + rep(rnorm(length(unique(interaction(hue[, c("id", "group")]))), 90, 5), each = 25)
hue$y <- round(hue$y)
colnames(area)[4] <- "area"
colnames(height)[4] <- "height"
colnames(width)[4] <- "width"
colnames(leaf_count)[4] <- "leaf_count"
colnames(hue)[4] <- "hue"
df <- Reduce(plyr::join, list(area, height, width, leaf_count, hue))
df$treatment <- ifelse(df$group %in% c("a", "c"),
                    "trt1", "trt2")
df$genotype <- ifelse(df$group %in% c("a", "b"),
                   "geno1", "geno2")
write.csv(df, "~/pcvr/papers/pcvr_v1.x/simdf.csv", row.names = FALSE)

#* `FREM`

p1 <- frem(df, des = c("genotype", "treatment"),
     phenotypes = c("area", "height", "width", "leaf_count", "hue"),
     timeCol = "time", time = NULL) &
  plot_annotation(title = "Cross Sectional Variance Partitioning")

p2 <- frem(df, des = c("genotype", "treatment"),
     phenotypes = c("area", "height", "width", "leaf_count", "hue"),
     timeCol = "time", time = "all", cor = FALSE) +
  ggplot2::labs(title = "Longitudinal Variance Partitioning")

ggsave("~/pcvr/papers/pcvr_v1.x/frem_CS.png", p1,
       width = 7, height = 6, dpi = 300, bg = "#ffffff")
ggsave("~/pcvr/papers/pcvr_v1.x/frem_Longitudinal.png", p2,
       width = 7, height = 6, dpi = 300, bg = "#ffffff")

