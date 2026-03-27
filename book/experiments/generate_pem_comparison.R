## PEM interval comparison figure
library(survival)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(pammtools)
library(patchwork)

# Set seed for reproducibility
set.seed(20250115)

# Generate true hazard function (Weibull hazard with increasing shape)
# h(t) = (shape/scale) * (t/scale)^(shape-1)
shape = 2.5
scale = 1000
max_time = 3000
t_true = seq(0, max_time, length.out = 1000)
h_true = (shape/scale) * (t_true/scale)^(shape-1)
h_true[1] = 0  # Set h(0) = 0

# Simulate survival data from Weibull distribution
n = 500
# Generate event times from Weibull
event_times = rweibull(n, shape = shape, scale = scale)
# Generate censoring times (uniform)
censor_times = runif(n, min = max_time * 0.5, max = max_time * 1.5)
# Observed times and status
obs_times = pmin(event_times, censor_times)
status = as.numeric(event_times <= censor_times)

# Create data frame
sim_data = data.frame(
  id = 1:n,
  time = obs_times,
  status = status
)

# Configuration 1: Few intervals with equidistant boundaries
J_few = 10
cut_points_few_eq = seq(0, max_time, length.out = J_few + 1)

# Configuration 2: Many intervals with equidistant boundaries
J_many = 50
cut_points_many_eq = seq(0, max_time, length.out = J_many + 1)

# Configuration 3: Few intervals with data-driven boundaries (quantiles of event times)
J_few = 10
event_times_only = obs_times[status == 1]
if (length(event_times_only) > 0) {
  quantiles = quantile(event_times_only, probs = seq(0, 1, length.out = J_few + 1))
  quantiles[1] = 0
  quantiles[length(quantiles)] = max_time
  cut_points_few_data = unique(sort(quantiles))
} else {
  cut_points_few_data = seq(0, max_time, length.out = J_few + 1)
}

# Configuration 4: Many intervals with data-driven boundaries
J_many = 50
if (length(event_times_only) > 0) {
  quantiles = quantile(event_times_only, probs = seq(0, 1, length.out = J_many + 1))
  quantiles[1] = 0
  quantiles[length(quantiles)] = max_time
  cut_points_many_data = unique(sort(quantiles))
} else {
  cut_points_many_data = seq(0, max_time, length.out = J_many + 1)
}

# Function to fit PEM and extract hazard estimates
fit_pem_hazard = function(data, cut_points, label) {
  # Transform data
  ped_data = as_ped(
    data = data,
    formula = Surv(time, status) ~ 1,
    cut = cut_points[-1],
    id = "id"
  ) |>
    mutate(interval = factor(interval))
  
  # Fit Poisson regression
  pem_fit = glm(
    ped_status ~ interval - 1 + offset(offset),
    family = poisson(),
    data = ped_data
  )
  
  # Get predicted hazards
  pred_data = ped_data |>
    distinct(interval, .keep_all = TRUE) |>
    select(interval, tstart, tend, offset)
  
  # Predict expected number of events
  pred_data$expected_events = predict(pem_fit, newdata = pred_data, type = "response")
  # Convert to hazard rate
  pred_data$hazard = pred_data$expected_events / exp(pred_data$offset)
  
  # Create data frame with time and hazard
  result = data.frame(
    time = pred_data$tend,
    hazard = pred_data$hazard,
    config = label
  )
  
  # Add time 0
  result = rbind(
    data.frame(time = 0, hazard = result$hazard[1], config = label),
    result
  )
  
  return(result)
}

# Fit all configurations
haz_few_eq = fit_pem_hazard(sim_data, cut_points_few_eq, "Few intervals, equidistant")
haz_many_eq = fit_pem_hazard(sim_data, cut_points_many_eq, "Many intervals, equidistant")
haz_few_data = fit_pem_hazard(sim_data, cut_points_few_data, "Few intervals, data-driven")
haz_many_data = fit_pem_hazard(sim_data, cut_points_many_data, "Many intervals, data-driven")

# Create true hazard data frame (to overlay on all panels)
true_haz_df = data.frame(
  time = t_true,
  hazard = h_true
)

# Combine all PEM estimates
pem_data = bind_rows(
  haz_few_eq,
  haz_many_eq,
  haz_few_data,
  haz_many_data
)

# Create plot with facets
p_pem_comparison = ggplot(pem_data, aes(x = time, y = hazard)) +
  # True hazard (overlay on all panels)
  geom_line(data = true_haz_df, 
            linewidth = 1.2, color = "black", alpha = 0.7) +
  # PEM estimates
  geom_step(aes(color = config), linewidth = 0.8, alpha = 0.8) +
  facet_wrap(~ config, nrow = 2, ncol = 2,
             labeller = labeller(config = c(
               "Few intervals, equidistant" = "(a) Few intervals, equidistant",
               "Many intervals, equidistant" = "(b) Many intervals, equidistant",
               "Few intervals, data-driven" = "(c) Few intervals, data-driven",
               "Many intervals, data-driven" = "(d) Many intervals, data-driven"
             ))) +
  ylab("h(t)") +
  xlab("time") +
  scale_color_manual(
    values = c(
      "Few intervals, equidistant" = "steelblue",
      "Many intervals, equidistant" = "darkgreen",
      "Few intervals, data-driven" = "orange",
      "Many intervals, data-driven" = "purple"
    ),
    name = "Configuration"
  ) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10)) +
  xlim(c(0, max_time))

# Create directory if it doesn't exist
dir.create("book/Figures/reductions", showWarnings = FALSE, recursive = TRUE)

ggsave("book/Figures/reductions/pem-interval-comparison.png", p_pem_comparison,
       height = 6, width = 12, units = "in", dpi = 600)

cat("\nFigure saved to book/Figures/reductions/pem-interval-comparison.png\n")
