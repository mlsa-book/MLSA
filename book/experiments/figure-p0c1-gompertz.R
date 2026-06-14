## figure-p0c1-gompertz.R
## P0C1 (Introduction) · #fig-distrfunctions
## Density / hazard / CDF / survival of a Gompertz distribution fit to Swedish
## age-specific mortality (2019, ages 30-100).
## Data: eha::swedeaths + eha::swepop (Statistics Sweden, bundled in eha).

source("book/experiments/figure-prep.R")
library(eha)

year_use <- 2019
age_lo   <- 30L
age_hi   <- 100L

d_int <- subset(swedeaths, year == year_use & age >= age_lo & age <= age_hi)
p_int <- subset(swepop,    year == year_use & age >= age_lo & age <= age_hi)

d_agg <- aggregate(deaths ~ age, data = d_int, sum)
p_agg <- aggregate(pop    ~ age, data = p_int, sum)
lt    <- merge(d_agg, p_agg)
lt$mx <- lt$deaths / lt$pop
lt$x  <- lt$age - age_lo

## Gompertz hazard h(x) = a * exp(b * x), population-weighted log-linear fit.
fit_gp <- lm(log(mx) ~ x, data = lt, weights = pop)
a_gp   <- unname(exp(coef(fit_gp)[1]))
b_gp   <- unname(coef(fit_gp)[2])
cat(sprintf("Gompertz fit: a = %.4g, b = %.4f (age origin = %d)\n",
            a_gp, b_gp, age_lo))

ages_gp <- seq(age_lo, 110, length.out = 400)
x_gp    <- ages_gp - age_lo
hx_gp   <- a_gp * exp(b_gp * x_gp)
Hx_gp   <- (a_gp / b_gp) * (exp(b_gp * x_gp) - 1)
Sx_gp   <- exp(-Hx_gp)
Fx_gp   <- 1 - Sx_gp
fx_gp   <- hx_gp * Sx_gp

## Keep the symbol-based strip headers (author preference for this figure).
fun_levels_gp <- c("f(t): density", "h(t): hazard", "F(t): CDF", "S(t): survival")
df_gp <- rbind(
  data.frame(age = ages_gp, y = fx_gp, func = "f(t): density"),
  data.frame(age = ages_gp, y = hx_gp, func = "h(t): hazard"),
  data.frame(age = ages_gp, y = Fx_gp, func = "F(t): CDF"),
  data.frame(age = ages_gp, y = Sx_gp, func = "S(t): survival")
)
df_gp$func <- factor(df_gp$func, levels = fun_levels_gp)

g_gp <- ggplot(df_gp, aes(x = age, y = y)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ func, scales = "free_y", nrow = 2) +
  labs(x = "Age (years)", y = NULL) +
  theme_square_tile() +              # square facet tiles, including strip headers
  theme(text = element_text(size = 13))

save_fig(g_gp, "book/Figures/introduction/fig-p0c1-gompertz.png", width = 6, height = 6)

