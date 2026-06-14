## figure-p2c8-cindex.R
## P2C8 (Discrimination) · fig-p2c8-cindex
## Time-dependent C-index vs. time cutoff under three time-weightings,
## Cox model on lung data. Data: survival::lung.

source("book/experiments/figure-prep.R")
library(survival)
library(dplyr)
library(tidyr)
library(purrr)

set.seed(260607)
train <- sample(nrow(lung), nrow(lung) * 2 / 3)
test  <- setdiff(seq(nrow(lung)), train)
fit1  <- coxph(Surv(time, status) ~ ., lung[train, ])
tmax  <- as.numeric(quantile(lung$time, probs = seq.int(0.1, 1, 0.1), na.rm = TRUE))
t60 <- quantile(lung$time, 0.6, na.rm = TRUE); t70 <- quantile(lung$time, 0.7, na.rm = TRUE)
t80 <- quantile(lung$time, 0.8, na.rm = TRUE); t90 <- quantile(lung$time, 0.9, na.rm = TRUE)
timewts <- c("n", "S", "n/G2")

cindex_dat <- crossing(tmax = tmax, timewt = timewts) |>
  mutate(cindex = map2_dbl(tmax, timewt, \(tm, wt)
    concordance(fit1, ymax = tm, timewt = wt, newdata = lung[test, ])$concordance))

g <- ggplot(cindex_dat, aes(x = tmax, y = cindex, colour = timewt)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  geom_vline(xintercept = c(t60, t70, t80, t90), linetype = 2, colour = "gray40") +
  scale_colour_manual(values = okabe_ito, name = "Weighting") +
  labs(x = "Time cutoff", y = "C-index") +
  coord_cartesian(ylim = c(0.5, 1)) +
  theme_square_panel +
  theme(text = element_text(size = 15))

save_fig(g, "book/Figures/evaluation/fig-p2c8-cindex.png", width = 5.5, height = 4.5)
