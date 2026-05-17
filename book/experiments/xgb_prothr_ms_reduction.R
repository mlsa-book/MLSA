## Standalone reproducer for the multi-state reduction example in P4C23:
##   AJ (split by treatment) vs XGBoost-Poisson trained on the multi-state
##   PED via the two-stage reduction
##     (1) multi-state -> transition-specific single-event with left-truncation
##     (2) single-event LT -> Poisson regression on PED (partition-based)
##   producing book/Figures/reductions/tp-prothr-cmp.png.
##
## The integrated version is in code.R (the prothr MS comparison block);
## this script duplicates the dependencies inline so it can be run on its
## own as a quick check.

suppressPackageStartupMessages({
  library(mstate); library(survival)
  library(dplyr); library(tidyr); library(ggplot2)
  library(pammtools); library(xgboost)
})

set.seed(20260516)

## ----- prep ------------------------------------------------------------------
data("prothr", package = "mstate")
my.prothr <- prothr |>
  mutate(transition = as.factor(paste0(from, "->", to))) |>
  filter(Tstart != Tstop)
tmat <- matrix(NA, 3, 3)
tmat[1, 2] <- 1; tmat[1, 3] <- 2; tmat[2, 1] <- 3; tmat[2, 3] <- 4

## ----- AJ baseline split by treatment ---------------------------------------
fit_aj_arm <- function(arm) {
  d <- my.prothr |> filter(treat == arm)
  probtrans(
    msfit(coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = d,
                method = "breslow"), trans = tmat),
    predt = 0, direction = "forward")
}
tp_aj <- list(Placebo    = fit_aj_arm("Placebo"),
              Prednisone = fit_aj_arm("Prednisone"))

## ----- XGBoost-Poisson on the multi-state PED -------------------------------
prothr_ped <- prothr |>
  filter(Tstart != Tstop) |>
  mutate(transition = as.factor(paste0(from, "->", to)),
         treat      = as.factor(treat)) |>
  rename(tstart = Tstart, tstop = Tstop) |>
  select(-trans)
ped <- as_ped(
  data       = prothr_ped,
  formula    = Surv(tstart, tstop, status) ~ .,
  transition = "transition",
  id         = "id",
  timescale  = "calendar")

## Features: tend + treat + transition. Offset = log(intlen) goes in as
## base_margin during TRAINING only; left unset at prediction time so the
## prediction returns the hazard (rate per unit time) directly.
make_X <- function(df) {
  model.matrix(~ tend + treat + transition, data = df)[, -1, drop = FALSE]
}
ped_df <- as.data.frame(ped)
dtrain <- xgb.DMatrix(data = make_X(ped_df),
                      label = as.numeric(ped_df$ped_status))
setinfo(dtrain, "base_margin", ped_df$offset)
xgb_fit <- xgb.train(
  params = list(
    objective        = "count:poisson",
    base_score       = 1,
    eta              = 0.01,
    max_depth        = 3,
    lambda           = 0,
    alpha            = 0,
    gamma            = 0,
    min_child_weight = 1
  ),
  data                  = dtrain,
  nrounds               = 5000,
  watchlist             = list(train = dtrain),
  early_stopping_rounds = 10,
  verbose               = 0
)

## Prediction grid + cumulative hazard + product integral (via pammtools).
ndf <- make_newdata(ped, tend = unique(tend), treat = unique(treat),
                    transition = unique(transition))
brks <- attr(ped, "breaks")
ndf$intlen <- setNames(diff(c(0, brks)), brks)[as.character(ndf$tend)]
ndf$hazard <- predict(xgb_fit, xgb.DMatrix(data = make_X(as.data.frame(ndf))))
ndf <- ndf |>
  group_by(treat, transition) |>
  arrange(treat, transition, tend) |>
  mutate(cumu_hazard = cumsum(hazard * intlen)) |>
  ungroup()
tp_xgb_df <- ndf |>
  group_by(treat) |>
  add_trans_prob(object = NULL, ci = FALSE) |>
  ungroup() |>
  transmute(treat,
            transition = as.character(transition),
            time       = tend,
            trans_prob)

## ----- combine + plot -------------------------------------------------------
build_aj_long <- function(tp, arm) {
  ts <- tp[[1]]$time
  data.frame(
    time       = rep(ts, 4),
    transition = rep(c("0->1", "0->2", "1->0", "1->2"), each = length(ts)),
    Treatment  = factor(arm, levels = c("Placebo", "Prednisone")),
    method     = factor("AJ", levels = c("AJ", "XGBoost")),
    trans_prob = c(tp[[1]]$pstate2, tp[[1]]$pstate3,
                   tp[[2]]$pstate1, tp[[2]]$pstate3))
}
relabel_mstate <- c("1->2" = "0->1", "1->3" = "0->2",
                    "2->1" = "1->0", "2->3" = "1->2")
df_tp_cmp <- bind_rows(
  build_aj_long(tp_aj$Placebo,    "Placebo"),
  build_aj_long(tp_aj$Prednisone, "Prednisone"),
  tp_xgb_df |>
    mutate(transition = unname(relabel_mstate[transition]),
           Treatment  = factor(treat, levels = c("Placebo", "Prednisone")),
           method     = factor("XGBoost", levels = c("AJ", "XGBoost"))) |>
    select(time, transition, Treatment, method, trans_prob)
)

p_tp_cmp <- ggplot(df_tp_cmp,
                   aes(x = time, y = trans_prob,
                       colour = Treatment, linetype = method)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ transition) +
  scale_colour_manual(values = c("steelblue", "firebrick4")) +
  scale_linetype_manual(values = c("AJ" = "solid", "XGBoost" = "dotted")) +
  coord_cartesian(xlim = c(0, 4000), ylim = c(0, 1)) +
  labs(x = "Time", y = "Transition probability", linetype = "Method") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("Figures/reductions/tp-prothr-cmp.png", p_tp_cmp,
       width = 8, height = 6.5, dpi = 300)
cat("Saved book/Figures/reductions/tp-prothr-cmp.png\n")
