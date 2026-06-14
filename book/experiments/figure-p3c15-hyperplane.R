## Output: book/Figures/svm/fig-p3c15-hyperplane.png
##
## Two-panel schematic of hyperplanes as prediction surfaces. Left (panel A): a
## single predictor, the fitted model is a line y = alpha + x * beta1, with
## short segments showing prediction error. Right (panel B): two predictors,
## the fitted model is a plane shown as an oblique-projected grid, with points
## and residual segments to the surface. Self-contained; sources figure-prep.R
## for the book-wide export helper (theme_void() is set explicitly on panel B).

source("book/experiments/figure-prep.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

## Fitted-plane coefficients and the oblique 2D projection used throughout.
## The plane is y = alpha + beta1 * x1 + beta2 * x2 and `project()` maps a
## 3D point (x1, x2, y) onto the same oblique 2D axes as the grid in panel B.
alpha <- 1
beta1 <- 0.55
beta2 <- 0.10

project <- function(x1, x2, y) {
  list(
    xp = x1 + 0.45 * x2,
    yp = y + 0.35 * x2 - 0.25 * x2
  )
}

set.seed(1)

pts <- data.frame(
  x1 = runif(10, 0.5, 5.5),
  x2 = runif(10, 0.5, 5.5)
)

pts$y_hat <- alpha + beta1 * pts$x1 + beta2 * pts$x2
pts$resid <- rnorm(10, 0, 0.65)
pts$y_obs <- pts$y_hat + pts$resid

# Project fitted and observed points
fit_proj <- project(pts$x1, pts$x2, pts$y_hat)
obs_proj <- project(pts$x1, pts$x2, pts$y_obs)

pts$xp_fit <- fit_proj$xp
pts$yp_fit <- fit_proj$yp
pts$xp_obs <- obs_proj$xp
pts$yp_obs <- obs_proj$yp
pts$resid_abs <- abs(pts$resid)

# Panel A: line
dat_line <- data.frame(
  x = seq(0, 6, length.out = 30)
)
dat_line$y <- 1 + 0.8 * dat_line$x

pts_line <- data.frame(
  x = seq(0.5, 5.5, length.out = 8)
)
pts_line$y <- 1 + 0.8 * pts_line$x + rnorm(8, 0, 0.45)

p1 <- ggplot() +
  geom_segment(
    data = pts_line,
    aes(
      x = x,
      xend = x,
      y = 1 + 0.8 * x,
      yend = y
    ),
    linewidth = 0.3,
    alpha = 0.7
  ) +
  geom_point(data = pts_line, aes(x, y), size = 2.2, alpha = 0.7) +
  geom_line(data = dat_line, aes(x, y), linewidth = 0.9) +
  annotate("text", x = 4.1, y = 2.0,
           label = "y == alpha + x * beta[1]",
           parse = TRUE, size = 5) +
  labs(
    x = "x", y = "y") +
  coord_cartesian(xlim = c(0, 6), ylim = c(0.5, 6.5))

# Panel B: projected plane / hyperplane
grid <- expand.grid(
  x1 = seq(0, 6, length.out = 13),
  x2 = seq(0, 6, length.out = 13)
)

# Oblique projection into 2D
grid$xp <- grid$x1 + 0.45 * grid$x2
grid$yp <- 1 + 0.55 * grid$x1 + 0.35 * grid$x2 - 0.25 * grid$x2

lines_x1 <- do.call(rbind, lapply(split(grid, grid$x2), function(d) {
  d[order(d$x1), ]
}))

lines_x2 <- do.call(rbind, lapply(split(grid, grid$x1), function(d) {
  d[order(d$x2), ]
}))

p2 <- ggplot() +
  geom_segment(
    data = pts,
    aes(
      x = xp_fit,
      y = yp_fit,
      xend = xp_obs,
      yend = yp_obs
    ),
    linewidth = 0.3,
    alpha = 0.7
  ) +
  geom_path(
    data = lines_x1,
    aes(x = xp, y = yp, group = x2),
    linewidth = 0.35,
    alpha = 0.4
  ) +
  geom_path(
    data = lines_x2,
    aes(x = xp, y = yp, group = x1),
    linewidth = 0.35,
    alpha = 0.4
  ) +
  annotate("text", x = 5.2, y = 1.25,
           label = "y == alpha + x[1] * beta[1] + x[2] * beta[2]",
           parse = TRUE, size = 5) +
  annotate("segment", x = 0, xend = 6, y = 1, yend = 4.3,
           linewidth = 0.6, arrow = arrow(length = unit(0.12, "inches"))) +
  annotate("segment", x = 0, xend = 2.7, y = 1, yend = -0.5,
           linewidth = 0.6, arrow = arrow(length = unit(0.12, "inches"))) +
  annotate("segment", x = 0, xend = 0, y = 1, yend = 5.7,
           linewidth = 0.6, arrow = arrow(length = unit(0.12, "inches"))) +
  annotate("text", x = 5.3, y = 4.5, label = expression(x[1]), size = 5) +
  annotate("text", x = 3.1, y = -0.67, label = expression(x[2]), size = 5) +
  annotate("text", x = -0.2, y = 5.9, label = expression(y), size = 5) +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_equal(xlim = c(-0.5, 8.8), ylim = c(-0.9, 6.2), clip = "off") +
  theme_void() +
  geom_point(
    data = pts,
    aes(x = xp_obs, y = yp_obs, size = resid_abs), alpha = 0.7
  ) +
  scale_size_continuous(
    range = c(0.8, 2.2),
    guide = "none"
  ) +
  theme(legend.position = "none")

g_hyperplane <- p1 + p2

save_fig(g_hyperplane, "book/Figures/svm/fig-p3c15-hyperplane.png",
         width = 8, height = 3.5, trim = FALSE)
