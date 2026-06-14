## figure-p2c9-calibKM.R
## P2C9 (Calibration) · fig-p2c9-calibKM
## Kaplan-Meier-style calibration: predicted average survival S(T) over time for
## three candidate models (CPH, RF, RRT) plus the KM reference curve.
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

# RRT: surv.rpart crank reused as lp, then AFT distr-composed on a KM baseline.
prpart <- lrn("surv.rpart")$train(t, s$train)$predict(t, s$test)
pkm    <- lrn("surv.kaplan")$train(t, s$train)$predict(t, s$test)
prp <- PredictionSurv$new(row_ids = prpart$row_ids, truth = prpart$truth,
                          crank = prpart$crank, lp = prpart$crank)
prrt <- po("distrcompose",
           param_vals = list(form = "aft", overwrite = FALSE, scale_lp = FALSE))$
  predict(list(base = pkm, pred = prp))[[1]]

pcox <- lrn("surv.coxph")$train(t, s$train)$predict(t, s$test)
pran <- lrn("surv.ranger")$train(t, s$train)$predict(t, s$test)

drrt <- autoplot(prrt, "calib")$data |> filter(Group == "Pred") |> mutate(Group = "RRT")
dcox <- autoplot(pcox, "calib")$data |> mutate(Group = if_else(Group == "KM", "KM", "CPH"))
dran <- autoplot(pran, "calib")$data |> filter(Group == "Pred") |> mutate(Group = "RF")

dat <- rbind(dcox, dran, drrt)
dat$Group <- factor(dat$Group, levels = c("KM", "CPH", "RF", "RRT"))

g <- ggplot(dat, aes(x = x, y = y, colour = Group)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = okabe_ito, name = "Model") +
  ylim(0, 1) +
  labs(x = "Time", y = "Survival Probability") +
  theme_square_panel +
  theme(text = element_text(size = 13))

save_fig(g, "book/Figures/evaluation/fig-p2c9-calibKM.png", width = 6, height = 5)
