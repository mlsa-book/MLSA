# figure-p1c6-heavisides.R  (fig-survtsk-heaviside)
# Faithful reproduction of book/experiments/figure-review/redo-heavisides.R:
# four-panel Heaviside / step-function illustration of survival curves built
# from indicator (Heaviside) steps, plus fitted Weibull survivor curves
# (unconditional, and conditional on sex). The conditional panel uses the
# Okabe-Ito categorical palette (blue / vermillion) rather than the ggplot
# default hcl palette.
# Output: book/Figures/survtsk/fig-p1c6-heavisides.png

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(flexsurv)
  library(survival)
})
theme_set(theme_bw())

okabe_ito_2 <- c("#0072B2", "#D55E00")  # blue, vermillion

set.seed(20250822)

sex <- rbinom(10, 1, 0.5)
x_base <- ceiling(runif(10, ifelse(sex == 1, 20, 1), ifelse(sex == 1, 50, 30)))
x <- unlist(lapply(x_base, function(val) c(0, val, 50)))

df <- data.frame(x = x, y = c(1, 0, 0),
                 group = as.factor(rep(1:10, each = 3)),
                 sex = rep(as.factor(sex), each = 3)) %>%
  mutate(alpha = if_else(group == 10, 1, 0.1))

events_df <- data.frame(time = x_base, event = 1L,
                        sex = factor(sex, levels = levels(df$sex)))

fit_uncond <- flexsurv::flexsurvreg(Surv(time, event) ~ 1,
                                    data = events_df, dist = "weibull")
fit_cond   <- flexsurv::flexsurvreg(Surv(time, event) ~ sex,
                                    data = events_df, dist = "weibull")

tgrid <- seq(0, 50, length.out = 200)

pred_uncond <- summary(fit_uncond, t = tgrid, type = "survival", tidy = TRUE)
pred_uncond_df <- data.frame(x = pred_uncond$time, y = pred_uncond$est)

pred_cond_df <- do.call(rbind, lapply(levels(events_df$sex), function(s) {
  nd  <- data.frame(sex = factor(s, levels = levels(events_df$sex)))
  out <- summary(fit_cond, newdata = nd, t = tgrid, type = "survival", tidy = TRUE)
  data.frame(x = out$time, y = out$est,
             sex = factor(s, levels = levels(events_df$sex)))
}))

g <- ggplot(df, aes(x = x, y = y, group = group))
g1 <- g + geom_step(linewidth = 1.3, color = "gray") +
  labs(x = "Time", y = "Survival Probability")
g2 <- g + geom_step(aes(alpha = alpha), linewidth = 1.3) +
  scale_alpha_identity() + labs(x = "Time", y = "Survival Probability")
g3 <- g1 + geom_line(aes(x = x, y = y), data = pred_uncond_df,
                     inherit.aes = FALSE, color = "black", linewidth = 1) +
  labs(x = "Time", y = "Survival Probability")
g4 <- g +
  geom_step(aes(color = sex), linewidth = 1.3, alpha = 0.5) +
  geom_line(aes(x = x, y = y, group = sex, color = sex),
            data = pred_cond_df, inherit.aes = FALSE, linewidth = 1.5) +
  scale_color_manual(values = okabe_ito_2) +
  labs(x = "Time", y = "Survival Probability")

g_final <- (g1 + g2 + g3 + g4) & ylim(0, 1) & xlim(0, 50) &
  guides(color = "none") & theme(text = element_text(size = 13))

ggsave("book/Figures/survtsk/fig-p1c6-heavisides.png",
       g_final, height = 5, width = 7, units = "in", dpi = 600)
## crop surrounding white margin, keep a small uniform border
system2("convert", c("book/Figures/survtsk/fig-p1c6-heavisides.png",
                     "-trim", "+repage", "-bordercolor", "white", "-border", "12x12",
                     "book/Figures/survtsk/fig-p1c6-heavisides.png"))
