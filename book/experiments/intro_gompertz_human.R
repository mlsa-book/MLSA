# Generate the four-panel Gompertz figure (PDF, hazard, CDF, survival) with
# parameters fit to real Swedish age-specific mortality data from 2019
# (pre-COVID), both sexes combined, ages 30-100.
#
# Data: eha::swedeaths and eha::swepop
# Source: https://www.scb.se/ (Statistics Sweden), bundled in the R package
# `eha` (covers 1969-2020).
#
# Output: book/Figures/introduction/gompertz.png

suppressPackageStartupMessages({
  library(eha)
  library(ggplot2)
})

year_use <- 2019
age_lo   <- 30L
age_hi   <- 100L

d <- subset(swedeaths, year == year_use & age >= age_lo & age <= age_hi)
p <- subset(swepop,    year == year_use & age >= age_lo & age <= age_hi)

d_agg <- aggregate(deaths ~ age, data = d, sum)
p_agg <- aggregate(pop    ~ age, data = p, sum)
lt <- merge(d_agg, p_agg)
lt$mx <- lt$deaths / lt$pop
lt$x  <- lt$age - age_lo

## Fit Gompertz hazard h(x) = a * exp(b * x) by population-weighted
## log-linear regression on the age-specific mortality rate.
fit <- lm(log(mx) ~ x, data = lt, weights = pop)
a <- unname(exp(coef(fit)[1]))
b <- unname(coef(fit)[2])

cat(sprintf("Gompertz fit: a = %.4g, b = %.4f (age origin = %d)\n",
            a, b, age_lo))

ages <- seq(age_lo, 110, length.out = 400)
x  <- ages - age_lo
hx <- a * exp(b * x)
Hx <- (a / b) * (exp(b * x) - 1)
Sx <- exp(-Hx)
Fx <- 1 - Sx
fx <- hx * Sx

fun_levels <- c("f(t): density",
                "h(t): hazard",
                "F(t): CDF",
                "S(t): survival")

df <- rbind(
  data.frame(age = ages, y = fx, func = "f(t): density"),
  data.frame(age = ages, y = hx, func = "h(t): hazard"),
  data.frame(age = ages, y = Fx, func = "F(t): CDF"),
  data.frame(age = ages, y = Sx, func = "S(t): survival")
)
df$func <- factor(df$func, levels = fun_levels)

g <- ggplot(df, aes(x = age, y = y)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ func, scales = "free_y", nrow = 2) +
  theme_bw(base_size = 11) +
  labs(x = "age (years)", y = NULL)

ggsave("book/Figures/introduction/gompertz.png", g,
       height = 3.5, width = 7, units = "in", dpi = 600)
