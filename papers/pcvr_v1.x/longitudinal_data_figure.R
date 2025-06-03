#* `longitudinal figure script`
setwd("~/pcvr/papers/pcvr_v1.x/")
library(devtools)
devtools::load_all("~/pcvr/")
simdf <- read.csv("simdf.csv")
sub <- simdf[simdf$genotype == "geno1" &
               simdf$treatment == "trt1", ]

#* `qualities of longitudinal data`

#* ***** `Autocorrelation Patch`
x1a <- 17
x2a <- 18
y1a <- sub[sub$time == x1a & sub$id == "id_6", "area"]
y2a <- sub[sub$time == x2a & sub$id == "id_6", "area"]

y1b <- sub[sub$time == x1a & sub$id == "id_15", "area"]
y2b <- sub[sub$time == x2a & sub$id == "id_15", "area"]

autocor_plot <- 
  ggplot(sub, aes(time, area, group = interaction(group, id))) +
  geom_line(aes(color = ifelse(id %in% c("id_6", "id_15"), "a", "b")),
            show.legend = FALSE) +
  scale_color_manual(values = c("black", "gray60")) +
  annotate("segment", x = x1a, y = y1a, xend = x1a, yend = y2a, color = "red") +
  annotate("segment", x = x1a, y = y2a, xend = x2a, yend = y2a, color = "red") +
  annotate("segment", x = x1a, y = y1b, xend = x1a, yend = y2b, color = "red") +
  annotate("segment", x = x1a, y = y2b, xend = x2a, yend = y2b, color = "red") +
  labs(title = "Autocorrelation") +
  pcv_theme() +
  theme(axis.title.y.left = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.title.x.bottom = element_blank())

#* ***** `Nonlinearity Patch`

nlPlot <- ggplot(sub, aes(time, area, group = interaction(group, id))) +
  geom_line() +
  geom_abline(slope = 8, intercept = 0, color = "red") +
  labs(title = "Non-Linearity") +
  pcv_theme() +
  theme(axis.title.y.left = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.title.x.bottom = element_blank())

sub$lmResid <- residuals(lm(area ~ time, sub))

residPlot <- ggplot(sub, aes(time, lmResid, group = interaction(group, id))) +
  geom_line(color = "gray70") +
  geom_hline(yintercept = 0, color = "red", linetype = 5) +
  labs(y = "lm Residuals") +
  pcv_theme() +
  theme(axis.text.x.bottom = element_blank(),
        axis.title.x.bottom = element_blank())

#* ***** `Heteroskedasticity Patch`

r1 <- range(sub[sub$time == 1, "area"])
r2 <- range(sub[sub$time == 5, "area"])
r3 <- range(sub[sub$time == 10, "area"])
r4 <- range(sub[sub$time == 20, "area"])

main <- ggplot(sub, aes(time, area, group = interaction(group, id))) +
  geom_line() +
  annotate("segment", x = 1, xend = 1, y = r1[1], yend = r1[2],
           color = "blue", linewidth = 2) +
  annotate("segment", x = 5, xend = 5, y = r2[1], yend = r2[2],
           color = "blue", linewidth = 2) +
  annotate("segment", x = 10, xend = 10, y = r3[1], yend = r3[2],
           color = "blue", linewidth = 2) +
  annotate("segment", x = 20, xend = 20, y = r4[1], yend = r4[2],
           color = "blue", linewidth = 2) +
  labs(title = "Heteroskedasticity") +
  pcv_theme() +
  theme(axis.title.x.bottom = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.title.y.left = element_blank())

sigma_df <- aggregate(area ~ group + time, data = sub, FUN = sd)

sigmaPlot <- ggplot(sigma_df, aes(x = time, y = area, group = group)) +
  geom_line(color = "blue", linetype = 5) +
  pcv_theme() +
  labs(y = "SD of area") +
  theme(plot.title = element_blank())

#* `Patchwork`

final_design <- c(
  patchwork::area(1, 1, 5, 4),
  patchwork::area(6, 1, 9, 4),
  patchwork::area(10, 1, 11, 4),
  patchwork::area(12, 1, 15, 4),
  patchwork::area(16, 1, 17, 4)
)

patch <- autocor_plot + nlPlot + residPlot +
  main + sigmaPlot +
  plot_layout(design = final_design, axes = "collect") +
  plot_annotation(title = "Longitudinal Area")

ggsave("~/pcvr/papers/pcvr_v1.x/longitudinal_data_traits.png",
       width = 7, height = 10, dpi = 300, bg = "#ffffff")

