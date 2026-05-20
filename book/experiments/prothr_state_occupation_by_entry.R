## ---------------------------------------------------------------------------
## State-occupation probabilities for the `prothr` cohort, CONDITIONAL on
## starting state.
##
## P(E(tau) = q | E(0) = ell) is, by definition, the pairwise transition
## probability P_{ell, q}(0, tau). This script just re-visualises the rows
## of the AJ transition-probability matrix as stacked-area plots, so that
## each ell-row sums to 1 across q. (No new estimand: same numbers as in
## the line-plot Fig @fig-multi-state-prothr, different presentation.)
##
## Layout: 2 x 2 facet
##   rows: entry state (state 0 = normal, state 1 = abnormal)
##   cols: treatment arm (Placebo, Prednisone)
##
## Output:
##   book/Figures/survival/prothr-state-occupation-by-entry.png
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mstate)
  library(survival)
})

fig_dir <- (function() {
  for (cand in c("Figures/survival",
                  "book/Figures/survival",
                  "../Figures/survival")) {
    if (dir.exists(cand)) return(cand)
  }
  "Figures/survival"
})()

## --- 1. Data + transition matrix (mirrors prothr_state_occupation.R) -------
data(prothr, package = "mstate")
my.prothr <- prothr |>
  mutate(transition = as.factor(paste0(from, "->", to))) |>
  filter(Tstart != Tstop)

tmat <- matrix(NA, nrow = 3, ncol = 3)
tmat[1, 2] <- 1
tmat[1, 3] <- 2
tmat[2, 1] <- 3
tmat[2, 3] <- 4

## --- 2. Aalen-Johansen TPs per arm -----------------------------------------
fit_aj <- function(arm) {
  d <- my.prothr |> dplyr::filter(treat == arm)
  cox <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans),
                data = d, method = "breslow")
  probtrans(msfit(cox, trans = tmat), predt = 0, direction = "forward")
}

tp.placebo    <- fit_aj("Placebo")
tp.prednisone <- fit_aj("Prednisone")

## --- 3. Reshape rows of the TP matrix into long format ---------------------
## tp[[ell]]  : data frame with $time, pstate1, pstate2, pstate3
##              = P_{ell, 0}(0, tau), P_{ell, 1}(0, tau), P_{ell, 2}(0, tau)
##              (mstate states 1/2/3 correspond to book states 0/1/2).
extract_row <- function(tp, entry_state_book, arm) {
  ## entry_state_book is 0 or 1; mstate row is entry_state_book + 1.
  df <- tp[[entry_state_book + 1L]]
  data.frame(
    time      = df$time,
    P0        = df$pstate1,
    P1        = df$pstate2,
    P2        = df$pstate3,
    entry     = sprintf("Entry state %d (%s)",
                        entry_state_book,
                        c("normal", "abnormal")[entry_state_book + 1L]),
    Treatment = arm
  )
}

cond_long <- bind_rows(
  extract_row(tp.placebo,    0, "Placebo"),
  extract_row(tp.placebo,    1, "Placebo"),
  extract_row(tp.prednisone, 0, "Prednisone"),
  extract_row(tp.prednisone, 1, "Prednisone")
) |>
  pivot_longer(cols = c(P0, P1, P2),
                names_to  = "state",
                values_to = "prob") |>
  mutate(
    state     = factor(state,
                        levels = c("P2", "P1", "P0"),
                        labels = c("State 2 (death)",
                                    "State 1 (abnormal)",
                                    "State 0 (normal)")),
    Treatment = factor(Treatment, levels = c("Placebo", "Prednisone")),
    entry     = factor(entry,
                        levels = c("Entry state 0 (normal)",
                                    "Entry state 1 (abnormal)"))
  )

## --- 4. Plot ---------------------------------------------------------------
pal <- c(
  "State 0 (normal)"   = "#8FB996",
  "State 1 (abnormal)" = "#7E9CC1",
  "State 2 (death)"    = "#B0413E"
)

p <- ggplot(cond_long, aes(x = time, y = prob, fill = state)) +
  geom_area(position = "stack", alpha = 0.92) +
  facet_grid(entry ~ Treatment) +
  scale_fill_manual(values = pal,
                     breaks = c("State 0 (normal)",
                                "State 1 (abnormal)",
                                "State 2 (death)")) +
  coord_cartesian(xlim = c(0, 4000), ylim = c(0, 1), expand = FALSE) +
  labs(
    x    = expression(tau),
    y    = expression(hat(P)[paste(ell, q)](0, tau)),
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

fig_path <- file.path(fig_dir, "prothr-state-occupation-by-entry.png")
ggsave(filename = fig_path, plot = p,
        width = 7.5, height = 5.5, dpi = 150)

## --- 5. Numeric summary at tau = 1825 (5 years), entry state 1 -------------
cat("Conditional state-occupation probabilities at tau = 1825 days,\n")
cat("given entry state 1 (abnormal prothrombin):\n\n")

snapshot <- function(tp, ell_book, tau) {
  df <- tp[[ell_book + 1L]]
  idx <- max(which(df$time <= tau))
  data.frame(
    P0_normal   = df$pstate1[idx],
    P1_abnormal = df$pstate2[idx],
    P2_death    = df$pstate3[idx]
  )
}

cat("Placebo:\n")
print(round(snapshot(tp.placebo,    1, 1825), 3))
cat("\nPrednisone:\n")
print(round(snapshot(tp.prednisone, 1, 1825), 3))

cat(sprintf("\nFigure written to %s\n", normalizePath(fig_path)))
