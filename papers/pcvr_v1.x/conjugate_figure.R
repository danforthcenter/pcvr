#* `conjugate figure script`
setwd("~/pcvr/papers/pcvr_v1.x/")
library(devtools)
devtools::load_all("~/pcvr/")
df <- read.csv("simdf.csv")
sub <- df[df$time == max(df$time), ]

s2 <- sub[sub$genotype == "geno1", ]

conj <- conjugate(
  s1 = area ~ treatment, s2 = s2,
  method = "t",
  priors = list(mu = 180, sd = 10),
  rope_range = c(-10, 10), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal",
  bayes_factor = c(180, 200)
)

plot(conj)

ggsave("~/pcvr/papers/pcvr_v1.x/conjugate.png",
       width = 8, height = 5, dpi = 300, bg = "#ffffff")
