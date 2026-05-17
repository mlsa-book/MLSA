## Competing Risks Reductions: PAM and RSF on sir.adm data
## Generates CIF comparison figures for Chapter P4C24

library(survival)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(tidyr)
library(pammtools)
library(randomForestSRC)
library(mvna)
library(patchwork)

set.seed(2025)

# ---- Load and prepare data ----
data(sir.adm, package = "mvna")
sir.adm$pneu = factor(sir.adm$pneu, labels = c("No pneumonia", "Pneumonia"))
sir.adm$sex = factor(sir.adm$sex)

# Create age group for covariate profiles
sir.adm$age_group = ifelse(sir.adm$age >= 60, "Age >= 60", "Age < 60")
sir.adm$age_group = factor(sir.adm$age_group, levels = c("Age < 60", "Age >= 60"))

# ---- Define covariate profiles for prediction ----
newdata_profiles = expand.grid(
  pneu = factor(c("No pneumonia", "Pneumonia"),
                levels = c("No pneumonia", "Pneumonia")),
  age = c(40, 75),
  sex = factor("M", levels = levels(sir.adm$sex))
)
newdata_profiles$age_group = ifelse(newdata_profiles$age >= 60,
                                    "Age >= 60", "Age < 60")
newdata_profiles$profile = paste0(newdata_profiles$pneu, ", ",
                                  newdata_profiles$age_group)

# ---- PAM approach: Cause-specific models ----
# Cause 1: Discharge (status == 1)
sir_discharge = sir.adm |>
  mutate(status_discharge = ifelse(status == 1, 1, 0))
# Cause 2: Death (status == 2)
sir_death = sir.adm |>
  mutate(status_death = ifelse(status == 2, 1, 0))

# Transform to PED format for each cause
ped_discharge = as_ped(
  data = sir_discharge,
  formula = Surv(time, status_discharge) ~ pneu + age + sex,
  id = "id"
)

ped_death = as_ped(
  data = sir_death,
  formula = Surv(time, status_death) ~ pneu + age + sex,
  id = "id"
)

# Fit PAMs using mgcv::gam
pam_discharge = mgcv::gam(
  ped_status ~ s(tend) + pneu + s(age) + sex,
  family = poisson(),
  offset = offset,
  data = ped_discharge
)

pam_death = mgcv::gam(
  ped_status ~ s(tend) + pneu + s(age) + sex,
  family = poisson(),
  offset = offset,
  data = ped_death
)

# Predict hazards on a time grid for each profile
time_grid = seq(0.5, max(sir.adm$time), length.out = 100)

# Create prediction data for each cause and profile
predict_cif_pam = function(pam_discharge, pam_death, ped_discharge,
                           newdata_profiles, time_grid) {
  results = list()
  for (i in seq_len(nrow(newdata_profiles))) {
    prof = newdata_profiles[i, ]

    # Create prediction PED for discharge model
    pred_ped_dis = ped_discharge |>
      make_newdata(tend = unique(tend), pneu = prof$pneu,
                   age = prof$age, sex = prof$sex)
    pred_ped_dis$hazard_discharge = predict(pam_discharge,
      newdata = pred_ped_dis, type = "response")

    # Create prediction PED for death model
    pred_ped_death = ped_discharge |>
      make_newdata(tend = unique(tend), pneu = prof$pneu,
                   age = prof$age, sex = prof$sex)
    pred_ped_death$hazard_death = predict(pam_death,
      newdata = pred_ped_death, type = "response")

    # Combine: all-cause cumulative hazard
    times = pred_ped_dis$tend
    h_dis = pred_ped_dis$hazard_discharge
    h_death = pred_ped_death$hazard_death
    intlen = pred_ped_dis$intlen

    # Cumulative hazards
    H_dis = cumsum(h_dis * intlen)
    H_death = cumsum(h_death * intlen)
    H_all = H_dis + H_death

    # All-cause survival: S(t) = exp(-H(t))
    S_all = exp(-H_all)
    S_all_prev = c(1, S_all[-length(S_all)])

    # CIF: F_e(t) = sum S(t_{k-1}) * h_e(t_k) * intlen_k
    cif_dis = cumsum(S_all_prev * h_dis * intlen)
    cif_death = cumsum(S_all_prev * h_death * intlen)

    results[[i]] = data.frame(
      time = times,
      CIF_discharge = cif_dis,
      CIF_death = cif_death,
      profile = prof$profile
    )
  }
  bind_rows(results)
}

pam_cifs = predict_cif_pam(pam_discharge, pam_death, ped_discharge,
                           newdata_profiles, time_grid)

# Reshape for plotting
pam_cifs_long = pam_cifs |>
  pivot_longer(
    cols = c(CIF_discharge, CIF_death),
    names_to = "event",
    values_to = "cif"
  ) |>
  mutate(
    event = ifelse(event == "CIF_discharge", "Discharge", "Death"),
    method = "PAM"
  )

# ---- RSF approach ----
# rfsrc handles competing risks natively when status has values 0, 1, 2
sir_rfsrc = sir.adm
sir_rfsrc$pneu_num = as.numeric(sir_rfsrc$pneu) - 1
sir_rfsrc$sex_num = as.numeric(sir_rfsrc$sex) - 1

rfsrc_fit = rfsrc(
  Surv(time, status) ~ pneu + age + sex,
  data = sir.adm,
  ntree = 1000,
  cause = c(1, 2),
  seed = 2025
)

# Predict for the covariate profiles
rfsrc_pred = predict(rfsrc_fit, newdata = newdata_profiles)

# Extract CIFs: rfsrc stores CIF for each cause
rfsrc_times = rfsrc_pred$time.interest

rfsrc_cifs = list()
for (i in seq_len(nrow(newdata_profiles))) {
  rfsrc_cifs[[i]] = data.frame(
    time = rfsrc_times,
    CIF_discharge = rfsrc_pred$cif[i, , 1],
    CIF_death = rfsrc_pred$cif[i, , 2],
    profile = newdata_profiles$profile[i]
  )
}
rfsrc_cifs = bind_rows(rfsrc_cifs)

rfsrc_cifs_long = rfsrc_cifs |>
  pivot_longer(
    cols = c(CIF_discharge, CIF_death),
    names_to = "event",
    values_to = "cif"
  ) |>
  mutate(
    event = ifelse(event == "CIF_discharge", "Discharge", "Death"),
    method = "RSF"
  )

# ---- Combined data for comparison ----
all_cifs = bind_rows(pam_cifs_long, rfsrc_cifs_long)

# ---- PAM-only CIF figure ----
p_pam = ggplot(pam_cifs_long, aes(x = time, y = cif, color = profile)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~event) +
  labs(
    x = "Time (days)",
    y = "Cumulative Incidence",
    color = "Profile"
  ) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))

dir.create("book/Figures/reductions", showWarnings = FALSE, recursive = TRUE)

ggsave("book/Figures/reductions/cif-pam-sir.png", p_pam,
       height = 4.5, width = 9, units = "in", dpi = 600)
cat("Saved: book/Figures/reductions/cif-pam-sir.png\n")

# ---- RSF-only CIF figure ----
p_rfsrc = ggplot(rfsrc_cifs_long, aes(x = time, y = cif, color = profile)) +
  geom_step(linewidth = 0.8) +
  facet_wrap(~event) +
  labs(
    x = "Time (days)",
    y = "Cumulative Incidence",
    color = "Profile"
  ) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))

ggsave("book/Figures/reductions/cif-rfsrc-sir.png", p_rfsrc,
       height = 4.5, width = 9, units = "in", dpi = 600)
cat("Saved: book/Figures/reductions/cif-rfsrc-sir.png\n")

# ---- Comparison: PAM vs RSF (side-by-side) ----
p_comparison = ggplot(all_cifs, aes(x = time, y = cif, color = profile,
                                     linetype = method)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~event) +
  labs(
    x = "Time (days)",
    y = "Cumulative Incidence",
    color = "Profile",
    linetype = "Method"
  ) +
  theme(legend.position = "bottom") +
  guides(
    color = guide_legend(nrow = 2, order = 1),
    linetype = guide_legend(order = 2)
  )

ggsave("book/Figures/reductions/cif-comparison-sir.png", p_comparison,
       height = 5, width = 9, units = "in", dpi = 600)
cat("Saved: book/Figures/reductions/cif-comparison-sir.png\n")
