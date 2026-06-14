## Output: book/Figures/classical/fig-p3c13-loglogistic-hazard.png
## Log-logistic hazard curves with fixed scale 1 and varying shape parameter,
## showing the non-monotonic hazard the log-logistic distribution allows.

suppressPackageStartupMessages({
  library(ggplot2)
})
source("book/experiments/figure-prep.R")

loglogistic_hazard <- function(t, shape, scale = 1) {
  (shape / scale) * (t / scale)^(shape - 1) /
    (1 + (t / scale)^shape)
}

t_grid <- seq(0.001, 5, length.out = 1000)
shapes <- c(0.5, 1, 1.5, 3, 7)

dat <- expand.grid(
  T = t_grid,
  Shape = shapes
)

dat$hazard <- loglogistic_hazard(dat$T, dat$Shape, scale = 1)

g <- ggplot(dat, aes(x = T, y = hazard, color = factor(Shape))) +
  geom_line(linewidth = 0.8) +
  labs(x = "Time", y = "Log-logistic hazard, h(t)", color = "Shape") +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 5)) +
  theme(aspect.ratio = 1)            # square panel, consistent with sibling curves

save_fig(g, "book/Figures/classical/fig-p3c13-loglogistic-hazard.png",
  width = 6, height = 5)
