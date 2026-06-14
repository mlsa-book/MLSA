## Output: book/Figures/svm/fig-p3c15-regression.png
##
## Schematic support-vector-regression figure: an epsilon-tube around the
## fitted line g(x) = x, with red circles inside the tube (non-support
## observations) and blue diamonds outside (support vectors), annotated with
## the two slack parameters. Self-contained; sources figure-prep.R for the
## book-wide theme/export helpers.

source("book/experiments/figure-prep.R")
suppressPackageStartupMessages(library(ggplot2))

df <- data.frame(
  x = c(1.8, 2.4, 2.7, 3.0, 3.3, 4.4, 4.2, 4.5, 5.7, 6.2, 6.5, 5),
  y = c(1.4, 2.2, 2.6, 3.4, 1.2, 2.0, 2.8, 7, 3.8, 5.2, 6.2, 4.2)
)

eps <- 1
line <- function(x) x
upper <- function(x) x + eps
lower <- function(x) x - eps

df$support <- with(
  df,
  y >= upper(x) | y <= lower(x)
)

tube <- data.frame(
  x = seq(0.5, 8.5, length.out = 200)
)

g_svm <- ggplot(df, aes(x, y)) +
  geom_ribbon(
    data = tube,
    aes(
      x = x,
      ymin = lower(x),
      ymax = upper(x)
    ),
    inherit.aes = FALSE,
    fill = "grey90",
    alpha = 0.5
  ) +

  geom_abline(intercept = 0, slope = 1, linewidth = 0.8) +
  geom_abline(intercept = eps, slope = 1, linetype = "dashed", color = "#0cb702", linewidth = 0.75) +
  geom_abline(intercept = -eps, slope = 1, linetype = "dashed", color = "#ed68ed", linewidth = 0.75) +

  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18)) +
  scale_color_manual(values = c(`FALSE` = "#f8766d", `TRUE` = "#619cff")) +

  # xi*: point above upper epsilon tube
  geom_segment(
    aes(x = 4.5, xend = 4.5, y = upper(4.5), yend = 7),
    arrow = arrow(ends = "both", length = unit(0.12, "inches")),
  ) +
  annotate("text", x = 4.65, y = 6.5, label = expression(zeta[i]^"*"), color = "black", size = 6) +

  # xi': point below lower epsilon tube
  geom_segment(
    aes(x = 5.7, xend = 5.7, y = 3.8, yend = lower(5.7)),
    arrow = arrow(ends = "both", length = unit(0.12, "inches")),
  ) +
  annotate("text", x = 5.9, y = 4.3, label = expression(zeta[i]^minute), color = "black", size = 6) +

  annotate("text", x = 6.5, y = 8.2, label = expression(y + epsilon), hjust = 0, color = "#0cb702", size = 6) +
  annotate("text", x = 7.1, y = line(7.55),  label = "y", hjust = 0, size = 6) +
  annotate("text", x = 6.9, y = lower(7.55), label = expression(y - epsilon), hjust = 0, color = "#ed68ed", size = 6) +

  geom_point(
    aes(shape = support, color = support),
    size = 2.7,
    show.legend = FALSE
  ) +

  coord_cartesian(xlim = c(0.5, 8), ylim = c(0.5, 8.7), expand = FALSE) +
  labs(x = "x", y = "g(x)") +
  theme(aspect.ratio = 1)            # square panel -> diagonal reads at ~45 deg

save_fig(g_svm, "book/Figures/svm/fig-p3c15-regression.png",
         width = 6, height = 5)
