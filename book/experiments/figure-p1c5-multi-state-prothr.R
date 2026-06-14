## figure-p1c5-multi-state-prothr.R
## P1C5 (Event history analysis) · fig-multi-state-prothr
## Transition probabilities by treatment (Placebo vs Prednisone) for the
## liver-cirrhosis illness-death model. Data: mstate::prothr; mstate::probtrans.

source("book/experiments/figure-prep.R")   # theme_bw + white strips + legend right
library(survival)
library(mstate)
library(dplyr)

data(prothr, package = "mstate")
my.prothr <- prothr |>
  mutate(transition = as.factor(paste0(from, "->", to))) |>
  filter(Tstart != Tstop)
tmat <- matrix(NA, 3, 3); tmat[1, 2] <- 1; tmat[1, 3] <- 2; tmat[2, 1] <- 3; tmat[2, 3] <- 4

fit_tp <- function(trt) {
  cox <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans),
               data = my.prothr |> filter(treat == trt), method = "breslow")
  tp <- probtrans(msfit(cox, trans = tmat), predt = 0, direction = "forward")
  tt <- tp[[1]]$time
  data.frame(time = rep(tt, 4),
             transition = rep(c("0->1", "0->2", "1->0", "1->2"), each = length(tt)),
             Treatment = factor(trt, levels = c("Placebo", "Prednisone")),
             trans_prob = c(tp[[1]]$pstate2, tp[[1]]$pstate3, tp[[2]]$pstate1, tp[[2]]$pstate3))
}
overall_df <- rbind(fit_tp("Placebo"), fit_tp("Prednisone"))

## math strip labels 0 -> 1 etc. via plotmath (label_parsed)
overall_df$transition <- factor(overall_df$transition,
  levels = c("0->1", "0->2", "1->0", "1->2"),
  labels = c("0 %->% 1", "0 %->% 2", "1 %->% 0", "1 %->% 2"))

p_ms <- ggplot(overall_df, aes(time, trans_prob)) +
  facet_wrap(~ transition, labeller = label_parsed) +
  geom_step(aes(colour = Treatment), linewidth = 1) +
  scale_colour_manual(values = col_treatment, name = "Treatment") +   # Placebo blue, Prednisone vermillion
  coord_cartesian(xlim = c(0, 4000), ylim = c(0, 1)) +
  labs(x = "Time", y = "Transition probability") +
  theme_square_tile() +
  theme(text = element_text(size = 13))

save_fig(p_ms, "book/Figures/survival/fig-p1c5-multi-state-prothr.png", width = 7, height = 5.6)
