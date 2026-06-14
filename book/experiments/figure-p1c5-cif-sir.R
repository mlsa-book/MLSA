## figure-p1c5-cif-sir.R
## P1C5 (Event history analysis) · fig-cif-sir
## Cause-specific cumulative incidence (discharge vs. death) by pneumonia status,
## ICU admission data. Data: mvna::sir.adm; CIF via cmprsk::cuminc.

source("book/experiments/figure-prep.R")   # theme_bw + WHITE facet strips + legend right
library(dplyr)
library(stringr)
library(purrr)
library(cmprsk)

data(sir.adm, package = "mvna")
sir_cif <- sir.adm |>
  mutate(to = dplyr::case_when(status == 0 ~ "cens", .default = as.character(status)))
cif <- cuminc(sir_cif$time, sir_cif$to, group = sir_cif$pneu, cencode = "cens")
cif_b <- purrr::imap_dfr(cif[1:4], function(.x, .y) {
  as.data.frame(.x) |>
    cbind(pneumonia = str_sub(.y, 1, 1), transition = str_sub(.y, 3, 3))
}) |>
  rename(cif = est) |>
  mutate(pneumonia  = factor(pneumonia,  labels = c("no", "yes")),
         transition = factor(transition, labels = c("discharge", "death")))

p_cif <- ggplot(cif_b, aes(time, cif)) +
  geom_step(aes(colour = pneumonia), linewidth = 1) +
  facet_wrap(~ transition) +
  geom_vline(xintercept = 120, lty = 3) +
  scale_colour_manual(values = col_complications, name = "Pneumonia") +   # no=blue, yes=vermillion
  labs(x = expression("Time, " * tau), y = "Cumulative incidence") +
  coord_cartesian(xlim = c(0, 125), ylim = c(0, 1)) +
  theme_square_tile() +   # faceted -> square the whole tile (coord + strip)
  theme(text = element_text(size = 13))

save_fig(p_cif, "book/Figures/survival/fig-p1c5-cif-sir.png", width = 8, height = 4)
