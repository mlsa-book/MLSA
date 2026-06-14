## Output: book/Figures/reductions/fig-p4c22-pem-interaction.png
## Piecewise exponential model (Poisson regression with interval x complications
## interaction) vs Kaplan-Meier on the `tumor` data.

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(pammtools)
  library(broom)
})

source("book/experiments/figure-prep.R")

data("tumor", package = "pammtools")
tumor_comp <- tumor |> select(days, status, complications)

max_time <- max(tumor_comp$days)
cut_points <- seq(0, max_time, length.out = 101)

ped_data <- as_ped(
  data = tumor_comp,
  formula = Surv(days, status) ~ complications,
  cut = cut_points[-1],
  id = "id"
) |>
  mutate(interval = factor(interval))

## Poisson regression with offset (offset added automatically by pammtools)
pem_fit <- glm(
  ped_status ~ interval * complications + offset(offset),
  family = poisson(), data = ped_data
)

pred_data <- ped_data |>
  distinct(interval, complications, .keep_all = TRUE) |>
  select(interval, complications, tstart, tend, offset)
pred_data$expected_events <- predict(pem_fit, newdata = pred_data, type = "response")
pred_data$hazard <- pred_data$expected_events / exp(pred_data$offset)

surv_data <- pred_data |>
  arrange(complications, interval) |>
  group_by(complications) |>
  mutate(cumhaz = cumsum(hazard * (tend - tstart)),
         survival = exp(-cumhaz),
         time = tend) |>
  ungroup() |>
  select(complications, time, survival) |>
  mutate(model = "Piecewise Exponential (PEM)")

surv_data <- bind_rows(
  data.frame(complications = unique(tumor_comp$complications),
             time = 0, survival = 1, model = "Piecewise Exponential (PEM)"),
  surv_data
) |>
  arrange(complications, time)

km <- survfit(Surv(days, status) ~ complications, data = tumor_comp)
bkm <- broom::tidy(km) |>
  mutate(complications = gsub("complications=", "", strata),
         model = "Kaplan-Meier") |>
  select(complications, time, estimate, model) |>
  rename(survival = estimate)

plot_data <- bind_rows(surv_data, bkm)

p <- ggplot(plot_data, aes(x = time, y = survival,
                           color = complications, linetype = model)) +
  geom_step(data = filter(plot_data, model == "Kaplan-Meier"), linewidth = 1.2) +
  geom_line(data = filter(plot_data, model == "Piecewise Exponential (PEM)"),
            linewidth = 1.2) +
  ylab("Survival Probability") +
  xlab("Time") +
  ylim(c(0, 1)) +
  scale_color_manual(
    values = col_complications,
    labels = c(no = "no", yes = "yes")
  ) +
  scale_linetype_manual(
    values = c("Kaplan-Meier" = "dotted", "Piecewise Exponential (PEM)" = "solid")
  ) +
  guides(color = guide_legend(title = "Complications"),
         linetype = guide_legend(title = "Model")) +
  theme_square_panel

save_fig(p, "book/Figures/reductions/fig-p4c22-pem-interaction.png",
         width = 6, height = 6)
