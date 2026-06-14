## Output: book/Figures/classical/fig-p3c13-ph-vs-aft.png
## Increasing a covariate by log(2): left panel shows the PH model doubling the
## hazard h(t) (multiplicative on the y-axis / risk); right panel shows the AFT
## model reaching S(t) in half the time (multiplicative on the x-axis / time).

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)  # if_else
})
source("book/experiments/figure-prep.R")

t <- seq(0.1, 2.5, by = 0.02)

hweib <- function(shape, scale, eta, t, form) {
  if (form == "AFT") {
    mod_scale <- scale * exp(-eta)
  } else if (form == "PH") {
    mod_scale <- scale * exp(-eta / shape)
  }
  (shape / mod_scale) * (t / mod_scale)^(shape - 1)
}

sweib <- function(shape, scale, eta, t, form) {
  if (form == "AFT") {
    mod_scale <- scale * exp(-eta)
  } else if (form == "PH") {
    mod_scale <- scale * exp(-eta / shape)
  }
  exp(-(t / mod_scale)^shape)
}

# Modify shape and scale parameters for better intercepts
plotWeib <- function(type = c("hazard", "survival"), shape = 3, scale = 2) {
  type <- match.arg(type)

  fun <- ifelse(type == "hazard", hweib, sweib)
  baseline <- fun(shape, scale, 0, t, "PH")
  PH <- fun(shape, scale, log(2), t, "PH")
  AFT <- fun(shape, scale, log(2), t, "AFT")
  ylabel <- ifelse(type == "hazard", "hazard", "Survival Probability")

  df <- data.frame(
    y = c(baseline, PH, AFT),
    t = rep(t, 3),
    Model = factor(rep(c("Baseline", "PH", "AFT"),
      each = length(t)), levels = c("Baseline", "PH", "AFT"))
  )

  ggplot(df, aes(x = t, y = y, color = Model)) +
    geom_line(aes(alpha = if_else(
      type == "hazard" & Model == "AFT" |
        type == "survival" & Model == "PH",
      0, 1))) +
    theme(legend.position = "right") +
    guides(alpha = "none") +
    ylab(ylabel) +
    xlab("Time") +
    scale_color_manual(values = c("Baseline" = "black", "AFT" = "red", "PH" = "blue"),
      aesthetics = c("color", "fill"))
}

segment <- function(start, form) {
  if (form == "PH") {
    map <- aes(
      x = start,
      xend = start,
      y = hweib(3, 2, 0, start, "PH"),
      yend = hweib(3, 2, log(2), start, "PH")
    )
  } else if (form == "AFT") {
    map <- aes(
      x = start,
      xend = start * 2,
      y = sweib(3, 2, 0, start * 2, "AFT"),
      yend = sweib(3, 2, log(2), start, "AFT")
    )
  }

  geom_segment(map,
    arrow = arrow(ends = "both", length = unit(0.1, "in")),
    inherit.aes = FALSE, size = 0.3,
    color = "#6f6f6f"
  )
}

p1 <- plotWeib("hazard") +
  ylim(0, 5) +
  segment(1, "PH") + segment(1.5, "PH") + segment(2, "PH") +
  geom_label(aes(x = x, y = y), data.frame(x = 1.45, y = hweib(3, 2, log(2), 2, "PH")),
    label = expression(h[PH](t) == 2 * h[0](t)),
    inherit.aes = FALSE)

p2 <- plotWeib("survival") +
  ylim(0, 1) +
  segment(0.5, "AFT") + segment(0.75, "AFT") +
  segment(1, "AFT") +
  geom_label(aes(x = x, y = y), data.frame(x = 2.08, y = sweib(3, 2, 0, 1.5, "AFT")),
    label = expression(S[AFT](t) == S[0]("2t")),
    inherit.aes = FALSE)

g <- (p1 + p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", aspect.ratio = 1)   # keep both panels square

save_fig(g, "book/Figures/classical/fig-p3c13-ph-vs-aft.png",
  width = 8, height = 4, trim = FALSE)
