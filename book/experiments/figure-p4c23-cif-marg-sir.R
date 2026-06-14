## Output: book/Figures/reductions/fig-p4c23-cif-marg-sir.png
##
## Marginal CIF verification for the cause-specific reduction (P4C23,
## sir.adm example). Three estimators, covariate = pneu only:
##   - non-parametric Aalen-Johansen (cmprsk::cuminc)
##   - cause-specific Cox PH (one model per cause, combined via the reduction)
##   - reduction-based RSF (one single-event RSF per cause, combined)
## Faceted: discharge (top) / death (bottom) x pneu = no / yes (columns).
##
## Reduction procedure for Cox PH and RSF:
##   1. one model per cause-specific dataset,
##   2. extract cause-specific cumulative hazards H_q(t | x),
##   3. all-cause survival S(t | x) = exp(-(H_1 + H_2)),
##   4. CIF via discrete sum F_q(t_k) = sum_{j<=k} S(t_{j-1}) (H_q(t_j) - H_q(t_{j-1})).

source("book/experiments/figure-prep.R")

suppressPackageStartupMessages({
  library(randomForestSRC)
  library(survival)
  library(cmprsk)         # cuminc() for non-parametric AJ
  library(mvna)           # sir.adm
  library(dplyr)
  library(ggplot2)
})

data(sir.adm)

dat <- sir.adm |>
  mutate(
    status = as.integer(status),
    pneu   = factor(pneu, levels = c(0, 1), labels = c("no", "yes")),
    sex    = factor(sex,  levels = c("F", "M"))
  )

step_interp <- function(t_eval, t_known, h_known) {
  idx <- findInterval(t_eval, t_known)
  out <- numeric(length(t_eval))
  out[idx > 0] <- h_known[idx[idx > 0]]
  out
}

combine_cif <- function(H1, H2) {
  S_all  <- exp(-(H1 + H2))
  S_left <- cbind(1, S_all[, -ncol(S_all), drop = FALSE])
  dH1    <- cbind(H1[, 1, drop = FALSE], H1[, -1, drop = FALSE] - H1[, -ncol(H1), drop = FALSE])
  dH2    <- cbind(H2[, 1, drop = FALSE], H2[, -1, drop = FALSE] - H2[, -ncol(H2), drop = FALSE])
  list(
    cif_c1 = t(apply(S_left * dH1, 1, cumsum)),
    cif_c2 = t(apply(S_left * dH2, 1, cumsum))
  )
}

## Cox PH: cause-specific Breslow cumulative hazard for each newdata row.
cox_chf <- function(fit, newdata, t_grid) {
  sf <- survfit(fit, newdata = newdata)
  if (is.null(dim(sf$surv))) sf$surv <- matrix(sf$surv, ncol = 1)
  H <- t(-log(pmax(sf$surv, 1e-12)))           # [n_newdata x n_time]
  t(apply(H, 1, function(r) step_interp(t_grid, sf$time, r)))
}

## RSF (single-event): per-row CHF on a common grid.
rsf_chf <- function(fit, newdata, t_grid) {
  pr <- predict(fit, newdata = newdata, importance = "none")
  t(apply(pr$chf, 1, function(r) step_interp(t_grid, pr$time.interest, r)))
}

## Long-format helper for plotting.
mk_long <- function(mat, profile_labels, method, cause, t_grid) {
  data.frame(
    time    = rep(t_grid, each = length(profile_labels)),
    cif     = as.vector(mat),
    profile = rep(profile_labels, times = length(t_grid)),
    method  = method,
    cause   = cause
  )
}

set.seed(20260515)

t_grid <- sort(unique(c(0, dat$time)))
t_max  <- 200                       # cap x-axis at 200 days for readability

newdata_marg  <- data.frame(pneu = factor(c("no", "yes"), levels = c("no", "yes")))
profiles_marg <- c("no", "yes")

## AJ via cmprsk
aj <- with(dat, cuminc(ftime = time, fstatus = status, group = pneu, cencode = 0))
aj_cif <- function(group, cause, t_grid) {
  key <- paste(group, cause)
  step_interp(t_grid, aj[[key]]$time, aj[[key]]$est)
}
aj_c1 <- rbind(aj_cif("no", 1, t_grid), aj_cif("yes", 1, t_grid))
aj_c2 <- rbind(aj_cif("no", 2, t_grid), aj_cif("yes", 2, t_grid))

## Cox PH cause-specific
cox_c1_marg <- coxph(Surv(time, as.integer(status == 1L)) ~ pneu, data = dat)
cox_c2_marg <- coxph(Surv(time, as.integer(status == 2L)) ~ pneu, data = dat)
H1_cox <- cox_chf(cox_c1_marg, newdata_marg, t_grid)
H2_cox <- cox_chf(cox_c2_marg, newdata_marg, t_grid)
cif_cox <- combine_cif(H1_cox, H2_cox)

## Reduction-based RSF
dat_c1 <- mutate(dat, status_c = as.integer(status == 1L))
dat_c2 <- mutate(dat, status_c = as.integer(status == 2L))
rsf_c1_marg <- rfsrc(Surv(time, status_c) ~ pneu, data = dat_c1, ntree = 500, splitrule = "logrank")
rsf_c2_marg <- rfsrc(Surv(time, status_c) ~ pneu, data = dat_c2, ntree = 500, splitrule = "logrank")
H1_rsf_red <- rsf_chf(rsf_c1_marg, newdata_marg, t_grid)
H2_rsf_red <- rsf_chf(rsf_c2_marg, newdata_marg, t_grid)
cif_rsf_red <- combine_cif(H1_rsf_red, H2_rsf_red)

## Build long format (three estimators, matching the caption)
df_marg <- bind_rows(
  mk_long(aj_c1,              profiles_marg, "Aalen-Johansen",      "Cause 1: discharge", t_grid),
  mk_long(aj_c2,              profiles_marg, "Aalen-Johansen",      "Cause 2: death",     t_grid),
  mk_long(cif_cox$cif_c1,     profiles_marg, "Cox PH",              "Cause 1: discharge", t_grid),
  mk_long(cif_cox$cif_c2,     profiles_marg, "Cox PH",              "Cause 2: death",     t_grid),
  mk_long(cif_rsf_red$cif_c1, profiles_marg, "Reduction-based RSF", "Cause 1: discharge", t_grid),
  mk_long(cif_rsf_red$cif_c2, profiles_marg, "Reduction-based RSF", "Cause 2: death",     t_grid)
) |>
  filter(time <= t_max) |>
  mutate(
    method  = factor(method, levels = c("Aalen-Johansen", "Cox PH", "Reduction-based RSF")),
    profile = factor(profile, levels = c("no", "yes"))
  )

## AJ black dotted; Cox red solid; reduction-based RSF blue solid.
method_cols <- c("Aalen-Johansen"      = "#000000",
                 "Cox PH"              = "#D55E00",
                 "Reduction-based RSF" = "#0072B2")
method_lty  <- c("Aalen-Johansen"      = "dotted",
                 "Cox PH"              = "solid",
                 "Reduction-based RSF" = "solid")

p_marg <- ggplot(df_marg, aes(x = time, y = cif, colour = method, linetype = method)) +
  geom_step(linewidth = 1) +
  facet_grid(cause ~ profile,
             labeller = labeller(profile = function(x) paste0("pneu = ", x))) +
  scale_colour_manual(values = method_cols, name = "Method") +
  scale_linetype_manual(values = method_lty,  name = "Method") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time", y = "cumulative incidence") +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        text = element_text(size = 15))

save_fig(p_marg, "book/Figures/reductions/fig-p4c23-cif-marg-sir.png",
         width = 10, height = 7, trim = FALSE)
cat("Saved fig-p4c23-cif-marg-sir.png\n")
