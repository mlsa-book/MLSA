## Brier and log loss scoring rules for @fig-p2c10-brier-logloss (P2C10_rules.qmd)
##
## Output: book/Figures/evaluation/fig-p2c10-brier-logloss.png
##
## Two-panel figure: Brier score (left) and log loss (right) over a probabilistic
## prediction in [0, 1], for true outcomes y_i = 1 (blue) and y_i = 0 (red).
## Both losses are minimized when the prediction matches the observed outcome.
##
## Self-contained extract of the relevant block in book/experiments/code.R.

source("book/experiments/figure-prep.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

##------------------
## Logloss and Brier
##------------------
p <- seq(0.001, 0.999, length.out = 500)

df <- data.frame(
  Prediction = rep(p, 4),
  Value = c(
    (1 - p)^2, p^2,
    -log(p), -log(1 - p)
  ),
  Class = rep(c("y1", "y0", "y1", "y0"), each = length(p)),
  Loss = rep(c("Brier score", "Log loss"), each = 2 * length(p))
)


brier <- ggplot(subset(df, Loss == "Brier score"),
                aes(Prediction, Value, colour = Class)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  annotate(
    "label", x = 0.10, y = 1.00,
    label = "y[i] == 1",
    parse = TRUE
  ) +
  annotate(
    "label", x = 0.90, y = 1.00,
    label = "y[i] == 0",
    parse = TRUE,
  ) +
  labs(
    x = expression(Prediction ~ hat(p)[i]),
    y = "Brier score"
  ) +
  theme(legend.position = "none")

logloss <- ggplot(subset(df, Loss == "Log loss"),
                  aes(Prediction, Value, colour = Class)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    limits = c(0, 3.15),
    breaks = c(0, log(2), 1, 2, 3),
    labels = c("0", "0.69", "1", "2", "3")
  ) +
  annotate(
    "label", x = 0.15, y = 3.00,
    label = "y[i] == 1",
    parse = TRUE
  ) +
  annotate(
    "label", x = 0.85, y = 3.00,
    label = "y[i] == 0",
    parse = TRUE,
  ) +
  labs(
    x = expression(Prediction ~ hat(p)[i]),
    y = "Log loss"
  ) +
  theme(legend.position = "none")


save_fig((brier + logloss) & theme(text = element_text(size = 13)),
         "book/Figures/evaluation/fig-p2c10-brier-logloss.png",
         width = 8, height = 4)
