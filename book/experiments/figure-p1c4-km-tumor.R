## figure-p1c4-km-tumor.R
## P1C4 (Survival analysis) · fig-km-tumor
## Kaplan-Meier estimate of overall survival on the tumor data, with the
## median survival time read off via dotted guide lines. Single panel -> square.
## Data: pammtools::tumor.

source("book/experiments/figure-prep.R")   # theme_bw + WHITE strips + legend right + save_fig
library(survival)

data(tumor, package = "pammtools")

km  <- survfit(Surv(days, status) ~ 1, data = tumor)
bkm <- data.frame(time = km$time, estimate = km$surv)
med <- as.numeric(quantile(km, probs = 0.5)$quantile)
df_med <- data.frame(x = c(0, med), y = c(0.5, 0),
                     xend = c(med, med), yend = c(0.5, 0.5))

p_km <- ggplot(bkm, aes(time, estimate)) +
  geom_step() +
  geom_segment(data = df_med, aes(x = x, xend = xend, y = y, yend = yend), lty = 3) +
  ylim(0, 1) +
  labs(x = "Time", y = "Survival Probability") +
  theme_square_panel   # single-panel survival curve -> square coordinate box

save_fig(p_km, "book/Figures/survival/fig-p1c4-km-tumor.png", width = 5, height = 5)
