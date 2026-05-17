## Generate the two CIF comparison figures for P4C23 §"Application to ICU
## mortality (sir.adm)":
##
##   book/Figures/reductions/cif-marg-sir.png  — marginal verification:
##     AJ, cause-specific Cox PH, reduction-based RSF, native CR RSF;
##     covariate = pneu only.
##
##   book/Figures/reductions/cif-adj-sir.png   — covariate-adjusted CIFs:
##     cause-specific Cox PH, reduction-based RSF, native CR RSF;
##     covariates = pneu + age + sex; 4 profiles (pneu x age = 40/75, sex = M).
##
## The reduction-based RSF and Cox PH BOTH follow the explicit cause-specific
## reduction procedure (@sec-cr-reduction-procedure):
##   1. one model per cause-specific dataset,
##   2. extract cause-specific cumulative hazards H_e(t | x),
##   3. compute all-cause survival   S(t | x) = exp(-(H_1 + H_2)),
##   4. assemble CIFs via discrete Aalen-Johansen-style sum
##        F_e(t_k | x) = sum_{j<=k} S(t_{j-1} | x) * (H_e(t_j) - H_e(t_{j-1})).
##
## The native CR RSF gets CIFs directly from randomForestSRC.
suppressPackageStartupMessages({
  library(randomForestSRC)
  library(survival)
  library(cmprsk)         # cuminc() for non-parametric AJ
  library(mvna)           # sir.adm
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

data(sir.adm)

dat <- sir.adm |>
  mutate(
    status = as.integer(status),
    pneu   = factor(pneu, levels = c(0, 1), labels = c("no", "yes")),
    sex    = factor(sex,  levels = c("F", "M"))
  )

step_interp <- function(t_eval, t_known, h_known) {
  ## right-continuous step function; values are 0 (or h_known[1]?) for t before t_known[1].
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
  ## sf$surv may be a vector (1 newdata row) or matrix (>1 row); normalise to matrix
  if (is.null(dim(sf$surv))) sf$surv <- matrix(sf$surv, ncol = 1)
  H <- t(-log(pmax(sf$surv, 1e-12)))           # [n_newdata x n_time]
  t(apply(H, 1, function(r) step_interp(t_grid, sf$time, r)))
}

## RSF (single-event): per-row CHF on a common grid.
rsf_chf <- function(fit, newdata, t_grid) {
  pr <- predict(fit, newdata = newdata, importance = "none")
  t(apply(pr$chf, 1, function(r) step_interp(t_grid, pr$time.interest, r)))
}

## Long-format helper for plotting:
##   mat is [n_profiles x n_time]; profile_labels has length n_profiles.
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

## ============================================================================
## MARGINAL: only pneu as covariate, two profiles (no / yes)
## ============================================================================
newdata_marg <- data.frame(
  pneu = factor(c("no", "yes"), levels = c("no", "yes"))
)
profiles_marg <- c("no", "yes")

## AJ via cmprsk
aj <- with(dat, cuminc(ftime = time, fstatus = status, group = pneu, cencode = 0))
aj_cif <- function(group, cause, t_grid) {
  key <- paste(group, cause)
  step_interp(t_grid, aj[[key]]$time, aj[[key]]$est)
}
aj_c1 <- rbind(aj_cif("no", 1, t_grid),
                aj_cif("yes", 1, t_grid))
aj_c2 <- rbind(aj_cif("no", 2, t_grid),
                aj_cif("yes", 2, t_grid))

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

## Native CR RSF
rsf_nat_marg <- rfsrc(Surv(time, status) ~ pneu, data = dat, ntree = 500, splitrule = "logrankCR")
pr_nat <- predict(rsf_nat_marg, newdata = newdata_marg, importance = "none")
cif_nat_c1 <- t(apply(pr_nat$cif[, , 1], 1, function(r) step_interp(t_grid, pr_nat$time.interest, r)))
cif_nat_c2 <- t(apply(pr_nat$cif[, , 2], 1, function(r) step_interp(t_grid, pr_nat$time.interest, r)))

## Build long format
df_marg <- bind_rows(
  mk_long(aj_c1,              profiles_marg, "Aalen-Johansen",      "Cause 1: discharge", t_grid),
  mk_long(aj_c2,              profiles_marg, "Aalen-Johansen",      "Cause 2: death",     t_grid),
  mk_long(cif_cox$cif_c1,     profiles_marg, "Cox PH",              "Cause 1: discharge", t_grid),
  mk_long(cif_cox$cif_c2,     profiles_marg, "Cox PH",              "Cause 2: death",     t_grid),
  mk_long(cif_rsf_red$cif_c1, profiles_marg, "Reduction-based RSF", "Cause 1: discharge", t_grid),
  mk_long(cif_rsf_red$cif_c2, profiles_marg, "Reduction-based RSF", "Cause 2: death",     t_grid),
  mk_long(cif_nat_c1,         profiles_marg, "Native CR RSF",       "Cause 1: discharge", t_grid),
  mk_long(cif_nat_c2,         profiles_marg, "Native CR RSF",       "Cause 2: death",     t_grid)
) |>
  filter(time <= t_max) |>
  mutate(
    method  = factor(method, levels = c("Aalen-Johansen", "Cox PH",
                                          "Reduction-based RSF", "Native CR RSF")),
    profile = factor(profile, levels = c("no", "yes"))
  )

method_cols <- c("Aalen-Johansen"      = "#2E3440",
                  "Cox PH"              = "#BF616A",
                  "Reduction-based RSF" = "#5E81AC",
                  "Native CR RSF"       = "#88C0D0")
method_lty  <- c("Aalen-Johansen"      = "solid",
                  "Cox PH"              = "dashed",
                  "Reduction-based RSF" = "dotted",
                  "Native CR RSF"       = "dotdash")

p_marg <- ggplot(df_marg, aes(x = time, y = cif, colour = method, linetype = method)) +
  geom_step(linewidth = 1) +
  facet_grid(cause ~ profile,
              labeller = labeller(profile = function(x) paste0("pneu = ", x))) +
  scale_colour_manual(values = method_cols, name = NULL) +
  scale_linetype_manual(values = method_lty,  name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "time (days)", y = "cumulative incidence") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("book/Figures/reductions/cif-marg-sir.png", p_marg,
        width = 10, height = 7, units = "in", dpi = 300)
cat("Saved cif-marg-sir.png\n")


## ============================================================================
## ADJUSTED: pneu + age + sex; 4 profiles (pneu x age 40/75, sex = M)
## ============================================================================
newdata_adj <- expand.grid(
  pneu = factor(c("no", "yes"), levels = c("no", "yes")),
  age  = c(40, 75),
  sex  = factor("M", levels = c("F", "M"))
)
newdata_adj$profile <- with(newdata_adj,
                            paste0("pneu=", pneu, ", age=", age))
profiles_adj <- newdata_adj$profile

## Cox PH
cox_c1_adj <- coxph(Surv(time, as.integer(status == 1L)) ~ pneu + age + sex, data = dat)
cox_c2_adj <- coxph(Surv(time, as.integer(status == 2L)) ~ pneu + age + sex, data = dat)
H1_cox_a   <- cox_chf(cox_c1_adj, newdata_adj, t_grid)
H2_cox_a   <- cox_chf(cox_c2_adj, newdata_adj, t_grid)
cif_cox_a  <- combine_cif(H1_cox_a, H2_cox_a)

## Reduction-based RSF
rsf_c1_adj <- rfsrc(Surv(time, status_c) ~ pneu + age + sex,
                     data = dat_c1, ntree = 500, splitrule = "logrank")
rsf_c2_adj <- rfsrc(Surv(time, status_c) ~ pneu + age + sex,
                     data = dat_c2, ntree = 500, splitrule = "logrank")
H1_rsf_a   <- rsf_chf(rsf_c1_adj, newdata_adj, t_grid)
H2_rsf_a   <- rsf_chf(rsf_c2_adj, newdata_adj, t_grid)
cif_rsf_a  <- combine_cif(H1_rsf_a, H2_rsf_a)

## Native CR RSF
rsf_nat_adj <- rfsrc(Surv(time, status) ~ pneu + age + sex,
                      data = dat, ntree = 500, splitrule = "logrankCR")
pr_nat_a    <- predict(rsf_nat_adj, newdata = newdata_adj, importance = "none")
cif_nat_c1_a <- t(apply(pr_nat_a$cif[, , 1], 1, function(r) step_interp(t_grid, pr_nat_a$time.interest, r)))
cif_nat_c2_a <- t(apply(pr_nat_a$cif[, , 2], 1, function(r) step_interp(t_grid, pr_nat_a$time.interest, r)))

df_adj <- bind_rows(
  mk_long(cif_cox_a$cif_c1,  profiles_adj, "Cox PH",              "Cause 1: discharge", t_grid),
  mk_long(cif_cox_a$cif_c2,  profiles_adj, "Cox PH",              "Cause 2: death",     t_grid),
  mk_long(cif_rsf_a$cif_c1,  profiles_adj, "Reduction-based RSF", "Cause 1: discharge", t_grid),
  mk_long(cif_rsf_a$cif_c2,  profiles_adj, "Reduction-based RSF", "Cause 2: death",     t_grid),
  mk_long(cif_nat_c1_a,      profiles_adj, "Native CR RSF",       "Cause 1: discharge", t_grid),
  mk_long(cif_nat_c2_a,      profiles_adj, "Native CR RSF",       "Cause 2: death",     t_grid)
) |>
  filter(time <= t_max) |>
  mutate(
    method  = factor(method, levels = c("Cox PH",
                                          "Reduction-based RSF",
                                          "Native CR RSF")),
    profile = factor(profile, levels = profiles_adj)
  )

p_adj <- ggplot(df_adj, aes(x = time, y = cif, colour = profile, linetype = method)) +
  geom_step(linewidth = 0.9) +
  facet_wrap(~ cause, ncol = 2) +
  scale_colour_brewer(palette = "Dark2", name = "profile") +
  scale_linetype_manual(values = method_lty[c("Cox PH", "Reduction-based RSF", "Native CR RSF")],
                         name = "method") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "time (days)", y = "cumulative incidence") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("book/Figures/reductions/cif-adj-sir.png", p_adj,
        width = 12, height = 6.5, units = "in", dpi = 300)
cat("Saved cif-adj-sir.png\n")
