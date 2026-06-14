## figure-p1c4-km-infants.R
## P1C4 (Survival) · fig-infants-lt
## Kaplan-Meier for infant survival by mother's vital status: naive KM (left)
## vs. left-truncation-adjusted KM (right). Data: eha::infants.

source("book/experiments/figure-prep.R")
library(survival)
library(broom)
library(patchwork)

data("infants", package = "eha")

## 2-level group -> colour (Okabe-Ito blue / vermillion)
mother_cols <- c(alive = "#0072B2", dead = "#D55E00")

prep <- function(fit) {
  d <- broom::tidy(fit)
  d$Mother <- factor(sub("mother=", "", d$strata), levels = c("alive", "dead"))
  d
}

km_naive <- prep(survfit(Surv(exit, event) ~ mother, data = infants))
km_lt    <- prep(survfit(Surv(enter, exit, event) ~ mother, data = infants))

km_panel <- function(d, title) {
  ggplot(d, aes(x = time, y = estimate, colour = Mother)) +
    geom_step() +
    scale_colour_manual(name = "Mother", values = mother_cols) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title, x = "Time", y = "Survival Probability") +
    theme_square_panel
}

p <- km_panel(km_naive, "Naive (ignoring truncation)") +
     km_panel(km_lt,    "Left-truncation adjusted") +
     plot_layout(guides = "collect")

save_fig(p, "book/Figures/survival/fig-p1c4-km-infants.png", width = 8.5, height = 4.2)
