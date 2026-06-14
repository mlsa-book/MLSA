## Output: book/Figures/forests/fig-p3c14-logrank.png
##
## Three-panel Kaplan-Meier illustration of the log-rank splitting rule on the
## `lung` dataset: (a) overall KM, (b) split at age > 50, (c) split at age > 75.
## Higher chi-squared / lower p-value (panel c) => greater dissimilarity between
## branches => better split. Self-contained; extracted from book/experiments/code.R.

library(survival)
library(ggplot2)
library(patchwork)
source("book/experiments/figure-prep.R")

set.seed(20241125)

sf0 <- survfit(Surv(time, status) ~ 1, lung)
p0 <- ggplot(
  rbind(data.frame(x = sf0$time, y = sf0$surv), data.frame(x = 0, y = 1)),
  aes(x = x, y = y)
) +
  geom_step() +
  labs(x = "Time", y = "Survival Probability")

df_1 <- cbind(
  lung,
  split = factor(lung$age > 50,
                 levels = c(FALSE, TRUE),
                 labels = c("Age ≤ 50", "Age > 50"))
)

df_2 <- cbind(
  lung,
  split = factor(lung$age > 75,
                 levels = c(FALSE, TRUE),
                 labels = c("Age ≤ 75", "Age > 75"))
)

logrank_1 <- survdiff(Surv(time, status) ~ split, data = df_1)
logrank_2 <- survdiff(Surv(time, status) ~ split, data = df_2)
sf_1 <- survfit(Surv(time, status) ~ split, df_1)
sf_2 <- survfit(Surv(time, status) ~ split, df_2)

make_plot <- function(sf, logrank, split_lab) {
  strata <- rep(names(sf$strata), sf$strata)
  ggplot(rbind(
    data.frame(x = sf$time, y = sf$surv, group = strata),
    data.frame(x = 0, y = 1, group = names(sf$strata))),
    aes(x = x, y = y, color = group)
  ) +
    geom_step(linewidth = 0.8) +
    geom_label(data = data.frame(x = 700, y = 0.80), aes(x = x, y = y),
      label = sprintf("%s\nχ² = %.2f\np = %.2f",
                      split_lab, logrank$chisq, logrank$pvalue),
      inherit.aes = FALSE, size = 3.9
    ) +
    scale_colour_manual(values = c("#0072B2", "#D55E00")) +   # usual 2-factor scheme
    labs(x = "Time", y = "Survival Probability", color = NULL)
}

p1 <- make_plot(sf_1, logrank_1, "Split: age > 50")
p2 <- make_plot(sf_2, logrank_2, "Split: age > 75")

g <- (p0 / (p1 + p2)) &
  guides(color = "none") &
  scale_x_continuous(limits = c(0, 1000), expand = c(0, 0)) &
  plot_annotation(tag_levels = "a", tag_suffix = ")")

## Multi-panel patchwork: keep raw ggsave margins (trim = FALSE).
save_fig(g, "book/Figures/forests/fig-p3c14-logrank.png",
         width = 8, height = 6, trim = FALSE)
