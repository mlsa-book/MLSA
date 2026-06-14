## figure-p4c21-pseudo-lm.R
## Pseudo-value regression for survival probabilities: stratified Kaplan-Meier
## curves (by complications) with linear-model pseudo-value predictions overlaid
## as points at tau in {1000, 2000, 3000} days (tumor data).
## Output: book/Figures/survival/fig-p4c21-pseudo-lm.png

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(survival)
  library(pammtools)
  library(pseudo)
  library(broom)
})

source("book/experiments/figure-prep.R")

data("tumor", package = "pammtools")
tumor_comp <- tumor |> select(days, status, complications)
tau <- c(1000, 2000, 3000)

## Pseudo-values for survival probabilities at each tau
pseudo_dfs <- purrr::map_dfr(tau, function(.x) {
  pseudo_mat <- pseudosurv(time = tumor$days, event = tumor$status, tmax = .x)
  cbind(tumor_comp, data.frame(pseudo = pseudo_mat$pseudo, time = pseudo_mat$time,
                               age = tumor$age))
}) |> mutate(time = factor(time))

ndf <- pammtools::make_newdata(pseudo_dfs, complications = unique(complications),
                               time = unique(time))
lm_fit <- lm(formula = pseudo ~ complications * time, data = pseudo_dfs)

ndf$estimate <- predict(lm_fit, newdata = ndf)
ndf <- ndf |>
  mutate(model = "pseudo-values (LM)") |>
  mutate(time = as.numeric(as.character(time)))

## Stratified KM wrt complications
km_complications <- survfit(Surv(days, status) ~ complications, data = tumor)
bkm_complications <- broom::tidy(km_complications) |>
  mutate(complications = gsub("complications=", "", strata)) |>
  mutate(model = "KM")

p_km_complications <- ggplot(bkm_complications, aes(x = time, y = estimate)) +
  geom_step(aes(col = complications)) +
  geom_vline(xintercept = tau, lty = 3) +
  geom_point(data = ndf, pch = 19, size = 2) +
  scale_colour_manual(values = col_complications, name = "Complications") +
  ylim(c(0, 1)) +
  ylab("Survival Probability") +
  xlab("Time") +
  theme_square_panel

save_fig(p_km_complications,
         "book/Figures/survival/fig-p4c21-pseudo-lm.png",
         width = 5, height = 4)
