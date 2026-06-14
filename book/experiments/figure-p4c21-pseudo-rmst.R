## figure-p4c21-pseudo-rmst.R
## RMST visualization at tau = 1000 days, faceted by complications status
## (tumor data). Each panel shows the group Kaplan-Meier survival curve with the
## area under the curve up to tau shaded; the numerical RMST value is annotated.
## Output: book/Figures/reductions/fig-p4c21-pseudo-rmst.png

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
  library(pammtools)
  library(broom)
})

source("book/experiments/figure-prep.R")

data("tumor", package = "pammtools")

## KM fits per group for plotting
km_complications <- survfit(Surv(days, status) ~ complications, data = tumor)
bkm_complications <- broom::tidy(km_complications) |>
  mutate(complications = gsub("complications=", "", strata))

tau_rmst <- 1000

## RMST per group via trapezoidal integration of the KM curve up to tau
calculate_rmst_at_tau <- function(km_fit, tau) {
  times <- km_fit$time[km_fit$time <= tau]
  surv  <- km_fit$surv[km_fit$time <= tau]
  if (length(times) == 0 || times[1] > 0) {
    times <- c(0, times)
    surv  <- c(1, surv)
  }
  if (max(times) < tau) {
    times <- c(times, tau)
    surv  <- c(surv, tail(surv, 1))
  }
  rmst <- 0
  for (i in 2:length(times)) {
    width      <- times[i] - times[i - 1]
    avg_height <- (surv[i] + surv[i - 1]) / 2
    rmst <- rmst + width * avg_height
  }
  rmst
}

km_no  <- survfit(Surv(days, status) ~ 1, data = tumor |> filter(complications == "no"))
km_yes <- survfit(Surv(days, status) ~ 1, data = tumor |> filter(complications == "yes"))

rmst_no  <- calculate_rmst_at_tau(km_no,  tau_rmst)
rmst_yes <- calculate_rmst_at_tau(km_yes, tau_rmst)

## Ribbon data: shaded area under each curve up to tau
bkm_for_ribbon <- bkm_complications |>
  filter(time <= tau_rmst) |>
  group_by(complications) |>
  arrange(time) |>
  group_modify(~ {
    df <- .x
    if (min(df$time) > 0) {
      df <- bind_rows(data.frame(time = 0, estimate = 1), df)
    }
    if (max(df$time) < tau_rmst) {
      df <- bind_rows(df, data.frame(time = tau_rmst, estimate = tail(df$estimate, 1)))
    }
    df
  }) |>
  ungroup()

rmst_annotations <- data.frame(
  complications = c("no", "yes"),
  rmst_value    = c(rmst_no, rmst_yes),
  x = tau_rmst / 2,
  y = 0.25
) |>
  mutate(label = paste0("RMST [", round(rmst_value, 1), "d]"))

p_rmst_complications <- ggplot(bkm_complications, aes(x = time, y = estimate)) +
  geom_ribbon(data = bkm_for_ribbon,
              aes(x = time, ymin = 0, ymax = estimate, fill = complications),
              alpha = 0.3, inherit.aes = FALSE) +
  geom_step(aes(color = complications), linewidth = 1.3) +
  geom_vline(xintercept = tau_rmst, lty = 2, alpha = 0.7, linewidth = 0.8) +
  geom_text(data = rmst_annotations,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 4.5, fontface = "bold") +
  facet_wrap(~ complications,
             labeller = labeller(complications = c("no" = "No Complications",
                                                   "yes" = "Complications"))) +
  scale_colour_manual(values = col_complications) +
  scale_fill_manual(values = col_complications) +
  ylab("Survival Probability") +
  xlab("Time") +
  labs(color = "Complications", fill = "Complications") +
  ylim(c(0, 1)) +
  xlim(c(0, max(bkm_complications$time))) +
  theme(legend.position = "none",
        text = element_text(size = 13),
        strip.text = element_text(size = 13, face = "bold"))

save_fig(p_rmst_complications,
         "book/Figures/reductions/fig-p4c21-pseudo-rmst.png",
         width = 10, height = 5, trim = FALSE)
