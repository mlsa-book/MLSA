## Output: book/Figures/forests/fig-p3c14-bootstrap.png
##
## Four-panel illustration of bootstrapping Kaplan-Meier estimators across three
## decision trees (red, blue, green):
##   a) the individual per-tree step-function estimates,
##   b) the union of event times to aggregate over (dashed verticals),
##   c) each tree's survival probability at every event time,
##   d) the averaged survival probabilities connected by a step function.
## Self-contained; extracted from book/experiments/code.R.

library(ggplot2)
library(patchwork)
source("book/experiments/figure-prep.R")

x <- 0:13
y1 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.9, 0.8, 0.1, 0.1)
y2 <- c(1, 1, 0.6, 0.6, 0.6, 0.6, 0.3, 0.3, 0.2, 0.2, 0, 0, 0, 0)
y3 <- c(1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2)
group <- rep(c("blue", "red", "green"), each = 14)
df <- data.frame(x = x, y = c(y1, y2, y3), group = group)

p0 <- ggplot(df, aes(x = x, y = y, color = group))
p1 <- p0 + geom_step() + labs(x = "Time", y = "Survival Probability")
p2 <- p0 + geom_vline(xintercept = x, lty = 2, color = "gray") +
  geom_step() + labs(x = "Time", y = "Survival Probability")
p3 <- p0 + geom_vline(xintercept = x, lty = 2, color = "gray") +
  geom_point() +
  geom_label(aes(x = x, y = y),
    data.frame(x = 5, y = c(0.9, 0.2, 0.4),
               group = c("blue", "red", "green")),
    label = sprintf("S(6) = %s", c(1, 0.3, 0.5))) +
  labs(x = "Time", y = "Survival Probability")

df2 <- data.frame(x = x, y = apply(data.frame(y1, y2, y3), 1, mean))
p4 <- ggplot(df2, aes(x = x, y = y)) + geom_step() + geom_point() +
  geom_label(aes(x = x, y = y),
    data.frame(x = 6, y = 0.5), label = "S(6) = 0.6") +
  labs(x = "Time", y = "Survival Probability")

ybreaks <- seq.int(0, 1, 0.25)
xbreaks <- seq.int(0, 12, 3)
g <- p1 + p2 + p3 + p4 &
  guides(color = "none") &
  scale_y_continuous(limits = c(0, 1), breaks = ybreaks, labels = ybreaks) &
  scale_x_continuous(limits = c(0, 14), breaks = xbreaks, labels = xbreaks,
                     expand = c(0, 0)) &
  plot_annotation(tag_levels = "a", tag_suffix = ")")

## Multi-panel patchwork: keep raw ggsave margins (trim = FALSE).
save_fig(g, "book/Figures/forests/fig-p3c14-bootstrap.png",
         width = 7, height = 6, trim = FALSE)
