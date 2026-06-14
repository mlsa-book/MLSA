## Output: book/Figures/svm/fig-p3c15-survival.png
##
## Schematic survival-time SSVM figure: fitted line g(x) = x + 0.2 with red
## circles (non-influential observations) and blue diamonds (support vectors)
## above and below the line, annotated with the upper- and lower-bound slack
## parameters. Self-contained; sources figure-prep.R for the book-wide
## theme/export helpers.

source("book/experiments/figure-prep.R")
suppressPackageStartupMessages(library(ggplot2))

df <- data.frame(
  x = c(1.7, 2.2, 2.4, 2.6, 3.0, 3.3, 3.9, 4.2, 5.2, 5.7, 6.1, 6.5),
  y = c(1.4, 2.2, 2.8, 5.2, 3.5, 2.0, 6.8, 2.4, 7.0, 1.2, 4.7, 7.1),
  support = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
              TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
)

line <- function(x) x + 0.2

g_svm_surv <- ggplot(df, aes(x, y)) +
  geom_abline(intercept = 0.2, slope = 1, linewidth = 0.8) +

  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18)) +
  scale_color_manual(values = c(`FALSE` = "#f8766d", `TRUE` = "#619cff")) +

  # xi double dagger: finite upper bound, above fitted line
  geom_segment(
    aes(x = 3.9, xend = 3.9, y = line(3.9), yend = 6.8),
    color = "#0cb702",
    arrow = arrow(ends = "both", length = unit(0.12, "inches"))
  ) +
  annotate(
    "text",
    x = 4.4, y = 6.45,
    label = expression(xi[i]^"*" * "," ~ i %in% U),
    size = 5
  ) +

  # xi dagger: finite lower bound, below fitted line
  geom_segment(
    aes(x = 5.7, xend = 5.7, y = 1.2, yend = line(5.7)),
    color = "#ed68ed",
    arrow = arrow(ends = "both", length = unit(0.12, "inches"))
  ) +
  annotate(
    "text",
    x = 6.2, y = 2.15,
    label = expression(xi[i]^minute * "," ~ i %in% L),
    size = 5
  ) +

  annotate("text", x = 7.2, y = 7.7, label = "t", size = 6) +

  geom_point(
    aes(shape = support, color = support),
    size = 2.8,
    show.legend = FALSE
  ) +

  coord_cartesian(xlim = c(0.6, 8), ylim = c(0.5, 8.4), expand = FALSE) +
  labs(x = "x", y = "g(x)") +
  theme(aspect.ratio = 1)            # square panel -> diagonal reads at ~45 deg

save_fig(g_svm_surv, "book/Figures/svm/fig-p3c15-survival.png",
         width = 6, height = 5)
