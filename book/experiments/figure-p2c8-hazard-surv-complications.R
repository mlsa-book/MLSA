## figure-p2c8-hazard-surv-complications.R
## P2C8 (Discrimination) · fig-p2c8-hazard-surv-complications
## Two-panel figure for the tumor data, stratified by surgical complications:
##   left  = estimated hazard h(t)         (curves cross)
##   right = estimated survival S(t)       (curves stay separated)
## Illustrates that crossing hazards need not imply crossing survival curves.
## Data: pammtools::tumor; baseline hazard estimated per group with a PAMM.

source("book/experiments/figure-prep.R")
suppressMessages({
  library(pammtools)
  library(mgcv)
  library(dplyr)
  library(patchwork)
})

data("tumor", package = "pammtools")

ped <- as_ped(
  data    = tumor,
  formula = Surv(days, status) ~ complications,
  id      = "id"
) |>
  mutate(complications_fac = factor(complications, levels = c("no", "yes")))

# PAMM: separate smooth baseline hazard per complications group.
pamm_fit <- gam(
  ped_status ~ s(tend, by = complications_fac, k = 15) + complications_fac + offset(offset),
  data   = ped,
  family = poisson()
)

ndf <- ped |>
  make_newdata(tend = unique(tend), complications_fac = unique(complications_fac)) |>
  mutate(offset = 0) |>
  add_hazard(pamm_fit) |>
  arrange(complications_fac, tend) |>
  group_by(complications_fac) |>
  mutate(
    dt         = tend - lag(tend, default = 0),
    cum_hazard = cumsum(hazard * dt),
    surv       = exp(-cum_hazard)
  ) |>
  ungroup() |>
  mutate(complications = as.character(complications_fac))

# Fixed book-wide complications mapping: no = blue, yes = vermillion.
lab_comp <- c(no = "no", yes = "yes")

p_haz <- ggplot(ndf, aes(x = tend, y = hazard, colour = complications)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = col_complications, labels = lab_comp,
                      name = "Complications") +
  labs(x = "Time", y = "Hazard") +
  theme_square_panel +
  theme(text = element_text(size = 15))

p_surv <- ggplot(ndf, aes(x = tend, y = surv, colour = complications)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = col_complications, labels = lab_comp,
                      name = "Complications") +
  labs(x = "Time", y = "Survival Probability") +
  theme_square_panel +
  theme(text = element_text(size = 15))

p_combined <- (p_haz | p_surv) + plot_layout(guides = "collect")

save_fig(p_combined, "book/Figures/evaluation/fig-p2c8-hazard-surv-complications.png",
         width = 9, height = 3.9)
