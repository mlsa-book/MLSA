## figure-p1c5-cens-vs-cr.R
## P1C5 (Event history analysis) · fig-cens-vs-cr
## Death CIF on the ICU admission data under two assumptions:
##   - competing risks (Aalen-Johansen via cmprsk::cuminc)   -> solid
##   - independent right-censoring (1 - exp(-H) from a Cox)  -> dotted
## Pneumonia status -> colour (no = blue, yes = vermillion).
## Data: mvna::sir.adm.

source("book/experiments/figure-prep.R")   # theme_bw + WHITE strips + legend right + save_fig
library(survival)
library(cmprsk)
library(dplyr)
library(stringr)
library(purrr)

data(sir.adm, package = "mvna")
sir <- sir.adm

## --- competing-risks CIF (Aalen-Johansen) ---------------------------------
sir_cif <- sir |>
  mutate(to = dplyr::case_when(status == 0 ~ "cens", .default = as.character(status)))
cif <- cuminc(sir_cif$time, sir_cif$to, group = sir_cif$pneu, cencode = "cens")
cif_b <- purrr::imap_dfr(cif[1:4], function(.x, .y) {
  as.data.frame(.x) |>
    cbind(pneumonia = str_sub(.y, 1, 1), transition = str_sub(.y, 3, 3))
}) |>
  rename(cif = est) |>
  mutate(pneumonia  = factor(pneumonia,  labels = c("no", "yes")),
         transition = factor(transition, labels = c("discharge", "death")))

## --- independent-censoring CIF (1 - exp(-H) from cause-specific Cox) -------
cox_sir <- purrr::map_dfr(1:2, function(.x) {
  tmp <- sir; tmp$status <- 1L * (tmp$status == .x)
  m <- coxph(Surv(time, status) ~ strata(pneu), data = tmp)
  basehaz(m) |>
    rename(pneumonia = strata) |>
    mutate(cif = 1 - exp(-hazard),
           pneumonia = dplyr::case_when(pneumonia == "pneu=0" ~ "no",
                                        pneumonia == "pneu=1" ~ "yes"),
           transition = ifelse(.x == 1, "discharge", "death"))
})

p_cvc <- ggplot(filter(cox_sir, transition == "death"), aes(time, cif)) +
  geom_step(aes(colour = pneumonia, lty = "independent censoring"), linewidth = 1) +
  geom_step(data = filter(cif_b, transition == "death"),
            aes(colour = pneumonia, lty = "competing risks"), linewidth = 1) +
  scale_colour_manual(values = col_complications, name = "Pneumonia") +  # no=blue, yes=vermillion
  scale_linetype_manual(values = c("competing risks"      = "solid",
                                   "independent censoring" = "dotted"),
                        name = "Assumption") +
  geom_vline(xintercept = 120, lty = 3) +
  labs(x = expression("Time, " * tau), y = "Cumulative incidence") +
  coord_cartesian(xlim = c(0, 125), ylim = c(0, 1)) +
  theme_square_panel +   # single-panel survival-curve plot -> square coordinate box
  theme(text = element_text(size = 14))

save_fig(p_cvc, "book/Figures/survival/fig-p1c5-cens-vs-cr.png", width = 6.5, height = 5)
