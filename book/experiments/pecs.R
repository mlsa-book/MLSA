## Prediction error curves (PECs) for @fig-eval-pecs (P2C10_rules.qmd)
##
## Recreates book/Figures/evaluation/pecs.png:
##   x = time, y = survival Brier / Graf score (the "loss") at each time point
##   Cox PH (red) vs survival SVM (blue); Cox PH is the better (lower) curve.
##
## Theme: repo standard theme_bw(); y-axis labelled "Loss" (not "Score").
##
## Models / score route
## --------------------
##   * tgen("simsurv") simulated right-censored single-event task (mlr3proba).
##   * Cox PH: lrn("surv.coxph") -- predicts a survival distribution natively.
##   * Survival SVM: lrn("surv.svm", type = "regression") from {survivalsvm}.
##     surv.svm only returns a relative risk ranking (crank); it has no `lp`,
##     so mlr3proba's distrcompositor leaves the prediction without a distr
##     (its .predict early-returns when `pred$lp` is NULL). We therefore
##     compose the SVM survival distribution by hand using the standard
##     proportional-hazards form  S_i(t) = S0(t)^exp(eta_i)  with a Kaplan-Meier
##     baseline S0 and eta = standardised SVM crank. This yields a valid `distr`
##     and a model that performs consistently WORSE than the Cox PH, as desired.
##   * Time-resolved Brier via msr("surv.graf", integrated = FALSE, times = t)
##     evaluated on a held-out test set over a grid of time points.

rm(list = ls())
suppressMessages({
  library(mlr3)
  library(mlr3proba)
  library(mlr3extralearners)
  library(survivalsvm)
  library(simsurv)
  library(distr6)
  library(ggplot2)
  library(dplyr)
})

theme_set(theme_bw())
set.seed(20231211)

## ---- Simulate data + train/test split -------------------------------------
task  <- tgen("simsurv")$generate(500)
split <- partition(task)

## ---- Cox PH (native distr prediction) --------------------------------------
pcox <- lrn("surv.coxph")$train(task, split$train)$predict(task, split$test)

## ---- Survival SVM + manual PH distribution composition ---------------------
psvm <- lrn("surv.svm", type = "regression", gamma = 0.1)$
  train(task, split$train)$predict(task, split$test)

## Kaplan-Meier baseline S0(t) (population-average survival on the test set).
pkm       <- lrn("surv.kaplan")$train(task, split$train)$predict(task, split$test)
base_surv <- colMeans(pkm$data$distr)
times     <- as.numeric(names(base_surv))

## Use the SVM risk ranking (crank) as the PH linear predictor.
eta <- psvm$crank
eta <- (eta - mean(eta)) / sd(eta) * 0.7   # centre + modest spread
nr  <- length(eta); nc <- length(times)
survmat  <- matrix(base_surv, nr, nc, byrow = TRUE)
surv_svm <- survmat^exp(matrix(eta, nr, nc))   # S_i(t) = S0(t)^exp(eta_i)

psvm <- PredictionSurv$new(
  row_ids  = psvm$row_ids,
  truth    = psvm$truth,
  crank    = psvm$crank,
  response = psvm$response,
  distr    = surv_return(times = times, surv = surv_svm)$distr
)

## ---- Time-resolved Brier (Graf) score over a grid --------------------------
test_times <- task$truth(split$test)[, 1L]
## keep the grid strictly inside the observed test-time range so the curve
## does not show a spurious tail artefact at the largest evaluation time.
grid <- seq(0, quantile(test_times, 0.95), length.out = 60)[-1]
grid <- grid[grid < max(test_times)]

brier_at <- function(pred, ts) {
  vapply(ts, function(tt) {
    as.numeric(pred$score(msr("surv.graf", integrated = FALSE, times = tt)))
  }, numeric(1))
}

df <- bind_rows(
  data.frame(time = grid, Loss = brier_at(pcox, grid), Model = "Cox PH"),
  data.frame(time = grid, Loss = brier_at(psvm, grid), Model = "SVM")
)

## ---- Plot ------------------------------------------------------------------
g <- ggplot(df, aes(x = time, y = Loss, color = Model)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Cox PH" = "#D7191C", "SVM" = "#2C7BB6")) +
  labs(x = "Time", y = "Loss", color = "Model") +
  ylim(0, NA)

ggsave("book/Figures/evaluation/pecs.png", g,
       width = 7, height = 4.5, units = "in", dpi = 300)

## keep the rendered _book copy in sync
if (dir.exists("book/_book/Figures/evaluation")) {
  ggsave("book/_book/Figures/evaluation/pecs.png", g,
         width = 7, height = 4.5, units = "in", dpi = 300)
}

cat("Cox PH integrated Brier:", as.numeric(pcox$score(msr("surv.graf"))), "\n")
cat("SVM    integrated Brier:", as.numeric(psvm$score(msr("surv.graf"))), "\n")
