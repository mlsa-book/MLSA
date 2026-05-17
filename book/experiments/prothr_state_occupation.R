## ---------------------------------------------------------------------------
## State-occupation probabilities for the `prothr` data (mstate)
##
## Computes the Aalen-Johansen plug-in estimator of the marginal
## state-occupation probabilities  P_q(tau) = P(E(tau) = q)  for the
## three-state cirrhosis multi-state process. The state labels follow the
## convention used in the book chapter (lines ~253ff of P1C5_eha.qmd):
##
##   state 0 = normal prothrombin   (= mstate code 1)
##   state 1 = abnormal prothrombin (= mstate code 2)
##   state 2 = death (absorbing)    (= mstate code 3)
##
## with possible transitions 0->1, 0->2, 1->0, 1->2 (cf. Fig. multi-state-
## prothr-states in book/Figures/survival/). Marginalisation uses the
## empirical initial-state distribution observed at study entry, separately
## for each treatment arm (prednisone vs. placebo). Most subjects enter the
## trial in state 1 (abnormal prothrombin).
##
## Output: stacked-area figure at
##   book/Figures/survival/prothr-state-occupation.png
## plus a small numeric summary table printed to stdout.
##
## This script mirrors the setup of the existing transition-probability
## script in book/experiments/code.R (lines ~643-704).
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mstate)
  library(survival)
})

## resolve the figure output path relative to this script's location so the
## script can be sourced from anywhere (./, ./book, ./book/experiments, ...)
fig_dir <- (function() {
  for (cand in c("Figures/survival",
                 "book/Figures/survival",
                 "../Figures/survival")) {
    if (dir.exists(cand)) return(cand)
  }
  "Figures/survival"
})()

## --- 1. Data and transition structure --------------------------------------
data(prothr, package = "mstate")

my.prothr <- prothr |>
  mutate(transition = as.factor(paste0(from, "->", to))) |>
  filter(Tstart != Tstop)  # drop instantaneous transitions

## transition matrix matching the existing chapter setup:
##   states are 1, 2, 3 in mstate (= 0, 1, 2 in the book chapter)
##   transitions: 1->2 (0->1), 1->3 (0->2), 2->1 (1->0), 2->3 (1->2)
tmat <- matrix(NA, nrow = 3, ncol = 3)
tmat[1, 2] <- 1
tmat[1, 3] <- 2
tmat[2, 1] <- 3
tmat[2, 3] <- 4
dimnames(tmat) <- list(
  from = c("0_normal", "1_abnormal", "2_death"),
  to   = c("0_normal", "1_abnormal", "2_death")
)

## --- 2. Aalen-Johansen transition probabilities per treatment arm ----------
fit_aj <- function(arm) {
  d <- my.prothr |> dplyr::filter(treat == arm)
  cox <- coxph(
    Surv(Tstart, Tstop, status) ~ strata(trans),
    data   = d,
    method = "breslow"
  )
  haz <- msfit(cox, trans = tmat)
  probtrans(haz, predt = 0, direction = "forward")
}

tp.placebo    <- fit_aj("Placebo")
tp.prednisone <- fit_aj("Prednisone")

## --- 3. Empirical initial-state distribution -------------------------------
## For each subject, the entry state is the value of `from` in the row with
## minimum Tstart. We take the unique (id, entry_state) per treatment arm.
init_state <- function(arm) {
  prothr |>
    dplyr::filter(treat == arm) |>
    group_by(id) |>
    slice_min(Tstart, with_ties = FALSE) |>
    ungroup() |>
    transmute(id = id, entry = from)
}

init_placebo    <- init_state("Placebo")
init_prednisone <- init_state("Prednisone")

pi_hat <- function(init_tbl) {
  ## returns a length-2 vector (pi_1, pi_2) = (P(start in state 0),
  ## P(start in state 1)) using mstate's 1/2 coding.
  tab <- table(factor(init_tbl$entry, levels = c(1, 2)))
  as.numeric(tab) / sum(tab)
}

pi_placebo    <- pi_hat(init_placebo)
pi_prednisone <- pi_hat(init_prednisone)

cat("Empirical initial-state distribution (entry state in mstate coding):\n")
cat(sprintf("  Placebo:    pi(state 0) = %.3f, pi(state 1) = %.3f, n = %d\n",
            pi_placebo[1], pi_placebo[2], nrow(init_placebo)))
cat(sprintf("  Prednisone: pi(state 0) = %.3f, pi(state 1) = %.3f, n = %d\n\n",
            pi_prednisone[1], pi_prednisone[2], nrow(init_prednisone)))

## --- 4. State-occupation probabilities -------------------------------------
## tp[[s]] holds P_{s, .}(0, tau) on the grid tp[[s]]$time.
## Columns are pstate1, pstate2, pstate3 = P(E(tau)=0), =1, =2.
## We marginalise: P_q(tau) = sum_l pi_l * P_{l q}(0, tau).
soccup <- function(tp, pi) {
  stopifnot(all.equal(tp[[1]]$time, tp[[2]]$time))
  data.frame(
    time = tp[[1]]$time,
    P0   = pi[1] * tp[[1]]$pstate1 + pi[2] * tp[[2]]$pstate1,
    P1   = pi[1] * tp[[1]]$pstate2 + pi[2] * tp[[2]]$pstate2,
    P2   = pi[1] * tp[[1]]$pstate3 + pi[2] * tp[[2]]$pstate3
  )
}

so_placebo    <- soccup(tp.placebo,    pi_placebo)
so_prednisone <- soccup(tp.prednisone, pi_prednisone)

so_placebo$Treatment    <- factor("Placebo",
                                  levels = c("Placebo", "Prednisone"))
so_prednisone$Treatment <- factor("Prednisone",
                                  levels = c("Placebo", "Prednisone"))

so_long <- bind_rows(so_placebo, so_prednisone) |>
  pivot_longer(cols = c(P0, P1, P2),
               names_to  = "state",
               values_to = "prob") |>
  mutate(
    state = factor(state,
                   levels = c("P2", "P1", "P0"),  # death on top of stack
                   labels = c("State 2 (death)",
                              "State 1 (abnormal)",
                              "State 0 (normal)"))
  )

## --- 5. Plot ---------------------------------------------------------------
## colour palette: blue-greens for transient states, red for absorbing
pal <- c(
  "State 0 (normal)"   = "#8FB996",  # muted green
  "State 1 (abnormal)" = "#7E9CC1",  # muted blue
  "State 2 (death)"    = "#B0413E"   # muted red
)

p_soccup <- ggplot(so_long, aes(x = time, y = prob, fill = state)) +
  geom_area(position = "stack", alpha = 0.92) +
  facet_wrap(~Treatment) +
  scale_fill_manual(values = pal,
                    breaks = c("State 0 (normal)",
                               "State 1 (abnormal)",
                               "State 2 (death)")) +
  coord_cartesian(xlim = c(0, 4000), ylim = c(0, 1), expand = FALSE) +
  labs(
    x    = expression(tau),
    y    = expression(hat(P)[q](tau)),
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

fig_path <- file.path(fig_dir, "prothr-state-occupation.png")
ggsave(
  filename = fig_path,
  plot     = p_soccup,
  width    = 6.5,
  height   = 4.2,
  dpi      = 150
)

## --- 6. Numeric summary at selected horizons -------------------------------
summarise_at_time <- function(so_df, taus) {
  out <- lapply(taus, function(tt) {
    idx <- max(which(so_df$time <= tt))
    data.frame(
      tau = tt,
      P0  = so_df$P0[idx],
      P1  = so_df$P1[idx],
      P2  = so_df$P2[idx]
    )
  })
  do.call(rbind, out)
}

taus <- c(365, 1825)

cat("State-occupation probabilities at selected horizons:\n\n")
cat("Placebo arm:\n")
print(round(summarise_at_time(so_placebo, taus), 3))
cat("\nPrednisone arm:\n")
print(round(summarise_at_time(so_prednisone, taus), 3))

cat(sprintf("\nFigure written to %s\n", normalizePath(fig_path)))
