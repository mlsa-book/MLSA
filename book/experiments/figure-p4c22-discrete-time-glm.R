## Output: book/Figures/reductions/fig-p4c22-discrete-time-glm.png
## Discrete-time survival reduction (logistic regression, proportional odds)
## vs Kaplan-Meier on the `tumor` data, split by complications status.

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(pammtools)
  library(broom)
})

source("book/experiments/figure-prep.R")

data("tumor", package = "pammtools")
tumor_comp <- tumor |> select(days, status, complications)

## 100 equidistant cut points -> 100 intervals
max_time <- max(tumor_comp$days)
cut_points <- seq(0, max_time, length.out = 101)

ped_data <- as_ped(
  data = tumor_comp,
  formula = Surv(days, status) ~ complications,
  cut = cut_points[-1],
  id = "id"
) |>
  mutate(interval = factor(interval))

transformed_data <- ped_data |> rename(delta_ij = ped_status)

## Proportional-odds logistic regression (common complications effect)
glm_fit <- glm(delta_ij ~ interval + complications,
               family = binomial(), data = transformed_data)

pred_data <- ped_data |>
  distinct(interval, complications, .keep_all = TRUE) |>
  select(interval, complications, tend)
pred_data$hazard <- predict(glm_fit, newdata = pred_data, type = "response")

surv_data <- pred_data |>
  arrange(complications, interval) |>
  group_by(complications) |>
  mutate(survival = cumprod(1 - hazard), time = tend) |>
  ungroup() |>
  select(complications, time, survival) |>
  mutate(model = "Discrete Time (GLM)")

surv_data <- bind_rows(
  data.frame(complications = unique(tumor_comp$complications),
             time = 0, survival = 1, model = "Discrete Time (GLM)"),
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
  geom_line(data = filter(plot_data, model == "Discrete Time (GLM)"), linewidth = 1.2) +
  ylab("Survival Probability") +
  xlab("Time") +
  ylim(c(0, 1)) +
  scale_color_manual(
    values = col_complications,
    labels = c(no = "no", yes = "yes")
  ) +
  scale_linetype_manual(
    values = c("Kaplan-Meier" = "dotted", "Discrete Time (GLM)" = "solid")
  ) +
  guides(color = guide_legend(title = "Complications"),
         linetype = guide_legend(title = "Model")) +
  theme_square_panel +
  theme(text = element_text(size = 13))

save_fig(p, "book/Figures/reductions/fig-p4c22-discrete-time-glm.png",
         width = 6, height = 6)
