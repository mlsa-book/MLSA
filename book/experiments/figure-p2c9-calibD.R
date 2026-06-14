## figure-p2c9-calibD.R
## P2C9 (Calibration) · fig-p2c9-calibD
## D-calibration plot: predicted vs. true quantiles for three candidate models
## (CPH, RF, RRT), with the diagonal as the perfectly-calibrated reference and a
## label reporting D-calibration chi-squared statistics / p-values.
## Data: mlr3proba tgen("simsurv") simulated single-event task.

source("book/experiments/figure-prep.R")
suppressMessages({
  library(mlr3)
  library(mlr3proba)
  library(mlr3extralearners)
  library(mlr3pipelines)
  library(dplyr)
})

set.seed(20231211)
t <- tgen("simsurv")$generate(400)
s <- partition(t)

prpart <- lrn("surv.rpart")$train(t, s$train)$predict(t, s$test)
pkm    <- lrn("surv.kaplan")$train(t, s$train)$predict(t, s$test)
prp <- PredictionSurv$new(row_ids = prpart$row_ids, truth = prpart$truth,
                          crank = prpart$crank, lp = prpart$crank)
prrt <- po("distrcompose",
           param_vals = list(form = "aft", overwrite = FALSE, scale_lp = FALSE))$
  predict(list(base = pkm, pred = prp))[[1]]

pcox <- lrn("surv.coxph")$train(t, s$train)$predict(t, s$test)
pran <- lrn("surv.ranger")$train(t, s$train)$predict(t, s$test)

drrt <- autoplot(prrt, "dcalib")$data |> mutate(Group = "RRT")
dcox <- autoplot(pcox, "dcalib")$data |> mutate(Group = "CPH")
dran <- autoplot(pran, "dcalib")$data |> mutate(Group = "RF")

dcal_coxp <- as.numeric(pcox$score(msr("surv.dcalib", truncate = Inf, chisq = TRUE)))
dcal_cox  <- as.numeric(pcox$score(msr("surv.dcalib", truncate = Inf, chisq = FALSE)))
dcal_ranp <- as.numeric(pran$score(msr("surv.dcalib", truncate = Inf, chisq = TRUE)))
dcal_ran  <- as.numeric(pran$score(msr("surv.dcalib", truncate = Inf, chisq = FALSE)))
dcal_rrtp <- as.numeric(prrt$score(msr("surv.dcalib", truncate = Inf, chisq = TRUE)))
dcal_rrt  <- as.numeric(prrt$score(msr("surv.dcalib", truncate = Inf, chisq = FALSE)))

# format p-values element-wise: very small ones print as "<0.001" rather than
# e-notation (format() on the whole vector would force e-notation on all of them)
fmt_p <- function(p) vapply(p, function(x)
  if (x < 0.001) "<0.001" else format(signif(x, 2)), character(1))
scores <- paste0(
  sprintf("   %s = %s (%s)", c("CPH", "RF", "RRT"),
          signif(c(dcal_cox, dcal_ran, dcal_rrt), 2),
          fmt_p(c(dcal_coxp, dcal_ranp, dcal_rrtp))),
  collapse = "\n")
scores <- paste0("DCal (p-values):\n", scores)

dat <- rbind(dcox, dran, drrt)
dat$Group <- factor(dat$Group, levels = c("CPH", "RF", "RRT"))

g <- ggplot(dat, aes(x = p, y = q, colour = Group)) +
  geom_abline(slope = 1, intercept = 0, colour = "lightgray", linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = okabe_ito, name = "Model") +
  xlim(0, 1) + ylim(0, 1) +
  labs(x = "True (p)", y = "Predicted") +
  geom_label(aes(x = x, y = y), data.frame(x = 0.55, y = 0.12), label = scores,
             inherit.aes = FALSE, hjust = "left", size = 4.5) +
  theme_square_panel +
  theme(text = element_text(size = 13))

save_fig(g, "book/Figures/evaluation/fig-p2c9-calibD.png", width = 6, height = 5)
