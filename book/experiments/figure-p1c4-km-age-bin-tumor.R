## figure-p1c4-km-tumor-age-bin.R
## P1C4 (Survival analysis) · fig-km-tumor-age-bin
## Kaplan-Meier estimates of survival on the tumor data, stratified into two
## age groups (< 50 vs. >= 50). Median guide lines for the older group.
## Single panel, two coloured curves -> square coordinate box.
## Data: pammtools::tumor.

source("book/experiments/figure-prep.R")   # theme_bw + WHITE strips + legend right + save_fig
library(survival)
library(dplyr)

data(tumor, package = "pammtools")

tumor <- tumor |>
  mutate(age_bin = factor(age < 50, levels = c(TRUE, FALSE),
                          labels = c("< 50", ">= 50")))

km_age  <- survfit(Surv(days, status) ~ age_bin, data = tumor)
bkm_age <- data.frame(time = km_age$time, estimate = km_age$surv,
                      strata = rep(names(km_age$strata), km_age$strata)) |>
  mutate(age = sub("age_bin=", "", strata))

med_age <- as.numeric(quantile(km_age, probs = 0.5)$quantile)[2]
df_age  <- data.frame(x = c(0, med_age), y = c(0.5, 0),
                      xend = c(med_age, med_age), yend = c(0.5, 0.5))

## Age groups -> colour from the Okabe-Ito palette.
age_cols <- c("< 50" = okabe_ito[5], ">= 50" = okabe_ito[6])  # blue / vermillion

p_km_age <- ggplot(bkm_age, aes(time, estimate)) +
  geom_step(aes(colour = age)) +
  geom_segment(data = df_age, aes(x = x, xend = xend, y = y, yend = yend), lty = 3) +
  geom_hline(yintercept = 0.5, lty = 3) +
  ylim(0, 1) +
  scale_colour_manual(values = age_cols, name = "Age") +
  labs(x = "Time", y = "Survival Probability") +
  theme_square_panel   # single-panel survival curve -> square coordinate box

save_fig(p_km_age, "book/Figures/survival/fig-p1c4-km-age-bin-tumor.png",
         width = 6, height = 5)
