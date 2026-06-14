## figure-p2c8-intervals.R
## P2C8 (Discrimination / ranking) · fig-p2c8-intervals
## Output: book/Figures/evaluation/fig-p2c8-intervals.png
##
## Six combinations of overlapping interval-censoring intervals (l, r] for a
## pair of observations (i, j), used to motivate the C-index under interval
## censoring. Self-contained; decision = keep new (figure must not change).

suppressPackageStartupMessages({
  library(dplyr)
})

source("book/experiments/figure-prep.R")

## --- Data: the six interval configurations -------------------------------
## Strip labels are plotmath expressions (rendered via label_parsed) so the
## subscripts in l_i, r_i, l_j, r_j typeset properly.
cases <- tibble::tribble(
  ~case, ~obs, ~l, ~r,
  "'1.'~~l[i]*' < '*r[i]*' < '*l[j]*' < '*r[j]", "i", 1, 2,
  "'1.'~~l[i]*' < '*r[i]*' < '*l[j]*' < '*r[j]", "j", 3, 4,

  "'2.'~~l[j]*' < '*r[j]*' < '*l[i]*' < '*r[i]", "i", 3, 4,
  "'2.'~~l[j]*' < '*r[j]*' < '*l[i]*' < '*r[i]", "j", 1, 2,

  "'3.'~~l[j]*' < '*l[i]*' < '*r[i]*' < '*r[j]", "i", 2, 3,
  "'3.'~~l[j]*' < '*l[i]*' < '*r[i]*' < '*r[j]", "j", 1, 4,

  "'4.'~~l[i]*' < '*l[j]*' < '*r[j]*' < '*r[i]", "i", 1, 4,
  "'4.'~~l[i]*' < '*l[j]*' < '*r[j]*' < '*r[i]", "j", 2, 3,

  "'5.'~~l[j]*' < '*l[i]*' < '*r[j]*' < '*r[i]", "i", 2, 4,
  "'5.'~~l[j]*' < '*l[i]*' < '*r[j]*' < '*r[i]", "j", 1, 3,

  "'6.'~~l[i]*' < '*l[j]*' < '*r[i]*' < '*r[j]", "i", 1, 3,
  "'6.'~~l[i]*' < '*l[j]*' < '*r[i]*' < '*r[j]", "j", 2, 4
) %>%
  mutate(
    y    = ifelse(obs == "i", 2, 1),
    case = factor(case, levels = unique(case))
  )

## --- Plot ----------------------------------------------------------------
g_intervals <- ggplot(cases) +
  geom_segment(aes(x = l, xend = r, y = y, yend = y, lty = obs),
               linewidth = 1) +
  facet_wrap(~case, ncol = 2, labeller = label_parsed) +
  scale_y_continuous(breaks = c(1, 2), labels = c("j", "i")) +
  labs(x = "Time", y = "") +
  theme(text = element_text(size = 18), legend.position = "n",
        strip.text = element_text(size = 19),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 16))

save_fig(g_intervals, "book/Figures/evaluation/fig-p2c8-intervals.png",
         width = 8, height = 6, dpi = 600)
