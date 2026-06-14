## figure-p3c15-comparable.R
## P3C15 (SVM) · fig-p3c15-comparable
## Survival-SVM as a ranking reduction: comparable (i, j(i)) observation pairs and
## the constraints between them. Schematic (no external data).
## The per-pair arrow colours are MEANINGFUL (one colour per origin observation i),
## so they are kept faithful to the original rather than remapped to a book palette.

source("book/experiments/figure-prep.R")

set.seed(2)
dat <- data.frame(
  i = 1:6, t = 1:6,
  status = c("censored", "uncensored", "censored", "uncensored", "censored", "censored")
)
pairs <- data.frame(
  i = c(2, 3, 4, 5, 6), j = c(1, 2, 2, 4, 4),
  x = c(2, 3, 4, 5, 6), xj = c(1, 2, 2, 4, 4),
  y_drop = c(0.55, 0.9, 1.4, 1.6, 2.3)
)
pairs$xj_arrow <- pairs$xj + runif(nrow(pairs), -0.07, 0.07)
pairs$lab_y <- pairs$y_drop - 0.17

g_cp <- ggplot(dat, aes(t, i)) +
  geom_segment(data = pairs, aes(x = x, xend = x, y = i, yend = y_drop, colour = factor(i)),
               inherit.aes = FALSE, linewidth = 0.8) +
  geom_point(aes(shape = status, fill = status, colour = status), size = 6, stroke = 1.1) +
  geom_segment(data = pairs, aes(x = x, xend = xj_arrow, y = y_drop, yend = y_drop, colour = factor(i)),
               inherit.aes = FALSE, linewidth = 0.8) +
  geom_segment(data = pairs, aes(x = xj_arrow, xend = xj_arrow, y = y_drop, yend = j - 0.2, colour = factor(i)),
               inherit.aes = FALSE, linewidth = 0.8, arrow = arrow(length = unit(0.12, "inches"))) +
  geom_text(data = pairs, aes(x = (x + xj) / 2 + 0.05, y = lab_y,
                              label = paste0("(i=", i, ", j(i)=", j, ")"), colour = factor(i)),
            inherit.aes = FALSE, size = 5, hjust = 0.5) +
  scale_shape_manual(values = c(censored = 21, uncensored = 22)) +
  scale_fill_manual(values = c(censored = "#d9e8ff", uncensored = "#ffd6d6")) +
  scale_colour_manual(values = c(
    censored = "#7aa6d9", uncensored = "#e06666",
    # one distinct colour per origin observation i (was: i=3 and i=4 both green)
    `2` = "#D55E00", `3` = "#009E73", `4` = "#E69F00", `5` = "#CC79A7", `6` = "#0072B2"
  )) +
  labs(x = "Outcome time", y = "Observation") +     # y title sits alongside the 1..6 ticks
  theme(legend.position = "none", aspect.ratio = 1,
        text = element_text(size = 13))             # square panel + axis text/title +2

save_fig(g_cp, "book/Figures/svm/fig-p3c15-comparable.png", width = 6, height = 5)
