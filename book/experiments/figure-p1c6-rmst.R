## figure-p1c6-rmst.R
## P1C6 (Survival task) · fig-survtsk-rmst
## Output: book/Figures/survtsk/fig-p1c6-rmst.png
##
## Illustrates the RMST at tau = 6 as the area under the survival curve for two
## groups (i, j), using a left-endpoint Riemann-sum approximation (shaded
## rectangles). The in-plot RMST annotation is enlarged for legibility.
## Self-contained: sources figure-prep.R for shared theme/export conventions.

suppressPackageStartupMessages({
  library(dplyr)
  library(patchwork)
})

source("book/experiments/figure-prep.R")

## --- Data ----------------------------------------------------------------
## Two stylised survival step-functions on the time grid 0..9.
yi <- c(1, 0.8, 0.75, 0.75, 0.7, rep(0.6, 5))
yj <- c(1, 0.9, 0.85, 0.6, 0.5, rep(0.1, 5))
df <- data.frame(
  x     = rep(0:9, 2),
  y     = c(yi, yj),
  Group = rep(c("i", "j"), each = 10)
)

## --- Per-group panel -----------------------------------------------------
plot_rmst <- function(df, group, tau = 6) {
  rect_df <- df %>%
    filter(Group == group, x < tau) %>%
    arrange(x) %>%
    mutate(xmin = x, xmax = x + 1, ymin = 0, ymax = y)

  rmst_hat <- sum(rect_df$y)

  ggplot(df, aes(x = x, y = y, group = Group)) +
    geom_rect(
      data = rect_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      alpha = 0.3,
      inherit.aes = FALSE
    ) +
    geom_step(aes(linetype = Group), linewidth = 1) +
    geom_vline(xintercept = tau, linetype = "dotted", linewidth = 0.5) +
    annotate(
      "text",
      x = 0.25,
      y = 0.4,
      label = as.expression(
        bquote(RMST[.(group)] * "(" * .(tau) * ")" == .(format(round(rmst_hat, 2), nsmall = 2)))
      ),
      parse = TRUE,
      ## Left-anchored (hjust = 0) just inside the y-axis. Kept moderate so the
      ## label does not run into the dotted tau guide line at x = tau.
      size = 4,
      hjust = 0
    ) +
    labs(
      x = "Time",
      y = "Survival Probability",
      title = paste0("RMST(", tau, ") for group ", group)
    ) +
    theme(legend.position = "right")
}

p1 <- plot_rmst(df, "i")
p2 <- plot_rmst(df, "j")

p_rmst_survtsk <- (p1 + p2) + plot_layout(guides = "collect") &
  theme(axis.title   = element_text(size = 13),
        axis.text    = element_text(size = 11),
        legend.title = element_text(size = 13),
        legend.text  = element_text(size = 11))

save_fig(p_rmst_survtsk, "book/Figures/survtsk/fig-p1c6-rmst.png",
         width = 7.5, height = 3.5, dpi = 600)
