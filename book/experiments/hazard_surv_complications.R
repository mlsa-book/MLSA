## Two-panel figure: hazard vs. survival probability for complications groups
## Illustrates that hazard curves can cross while survival probabilities may not
library(pammtools)
library(mgcv)
library(ggplot2)
library(dplyr)
library(patchwork)

theme_set(theme_bw())

# Load tumor data
data("tumor", package = "pammtools")

# Transform to PED format, keeping complications as covariate
ped <- as_ped(
  data    = tumor,
  formula = Surv(days, status) ~ complications,
  id      = "id"
) |>
  mutate(complications_fac = factor(complications, levels = c("no", "yes")))

# Fit PAMM: separate smooth baseline hazard per complications group
pamm_fit <- gam(
  ped_status ~ s(tend, by = complications_fac, k = 15) + complications_fac + offset(offset),
  data   = ped,
  family = poisson()
)

# Prediction grid: all unique time points x both groups
ndf <- ped |>
  make_newdata(tend = unique(tend), complications_fac = unique(complications_fac)) |>
  add_hazard(pamm_fit) |>
  arrange(complications_fac, tend) |>
  group_by(complications_fac) |>
  mutate(
    dt         = tend - lag(tend, default = 0),
    cum_hazard = cumsum(hazard * dt),
    surv       = exp(-cum_hazard)
  ) |>
  ungroup() |>
  mutate(
    group = factor(
      complications_fac,
      levels = c("no", "yes"),
      labels = c("No complications", "Complications")
    )
  )

p_haz <- ggplot(ndf, aes(x = tend, y = hazard, color = group)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Days", y = "h(t)", color = NULL) +
  theme(legend.position = "bottom")

p_surv <- ggplot(ndf, aes(x = tend, y = surv, color = group)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Days", y = "S(t)", color = NULL) +
  theme(legend.position = "bottom")

p_combined <- (p_haz | p_surv) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dir.create("book/Figures/evaluation", showWarnings = FALSE, recursive = TRUE)

ggsave(
  "book/Figures/evaluation/hazard-surv-complications.png",
  p_combined,
  height = 3.5, width = 8, units = "in", dpi = 600
)

cat("Figure saved to book/Figures/evaluation/hazard-surv-complications.png\n")
