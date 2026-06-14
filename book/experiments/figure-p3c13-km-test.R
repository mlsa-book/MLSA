## Output: book/Figures/classical/fig-p3c13-km-test.png
## Kaplan-Meier estimator as a predictive tool (rats dataset): top panel shows
## the step function with event/censoring times, bottom panel demonstrates
## reading off a survival probability at a given time.

suppressPackageStartupMessages({
  library(ggplot2)
  library(survival)
  library(patchwork)
})
source("book/experiments/figure-prep.R")

fit <- survfit(Surv(rats$time, rats$status) ~ 1)

g <- ggplot(data.frame(x = fit$time, y = fit$surv), aes(x = x, y = y)) +
  geom_step() + labs(x = "Time", y = "Survival Probability") +
  scale_x_continuous(expand = c(0, 0)) +
  theme(aspect.ratio = 1, text = element_text(size = 13))

g1 <- g +
  geom_vline(xintercept = fit$time[5:7], lty = 2, alpha = 1, color = 3, lwd = 1) +
  geom_vline(xintercept = fit$time[9:10], lty = 3, alpha = 1, color = 4, lwd = 1)

## filled-triangle arrowheads (type = "closed") for the read-off demonstration
g2 <- g +
  geom_segment(x = 60, y = 0, yend = fit$surv[12], color = 2, lwd = 1,
    arrow = arrow(type = "closed", length = unit(0.12, "inches"))) +
  geom_segment(x = 23, xend = fit$time[12], y = fit$surv[12], color = 2, lwd = 1,
    arrow = arrow(ends = "first", type = "closed", length = unit(0.12, "inches")))

g3 <- g1 | g2          # side by side (was stacked top/bottom)

save_fig(g3, "book/Figures/classical/fig-p3c13-km-test.png",
  width = 9, height = 4.5)
