## Output: book/Figures/classical/fig-p3c13-dogs.png
## Human vs dog lifespans: human (black, Gompertz), dog under AFT (red, 5x time
## acceleration), dog under PH with HR=5 (solid blue) and the much larger HR
## (~395, dashed blue) needed to match the AFT dog's median survival.

suppressPackageStartupMessages({
  library(ggplot2)
  library(extraDistr)  # pgompertz
})
source("book/experiments/figure-prep.R")

age <- seq.int(1, 100, 1)
surv <- pgompertz(age, 0.00005, 0.09, FALSE)
ph_surv <- surv^5
aft_surv <- round(pgompertz(age * 5, 0.00005, 0.09, FALSE), 2)

# A 5x acceleration of time (AFT) and a 5x hazard ratio (PH) act on different
# scales, so the same factor 5 is not comparable. The PH hazard ratio that
# matches the AFT dog's *median* survival is much larger (~395), because the
# Gompertz baseline hazard is near zero at young ages.
med_aft <- uniroot(function(a) pgompertz(a * 5, 0.00005, 0.09, FALSE) - 0.5, c(1, 100))$root
cstar <- log(0.5) / log(pgompertz(med_aft, 0.00005, 0.09, FALSE))
ph_match_surv <- surv^cstar
lab_match <- sprintf("Dog (PH, HR=%.0f)", cstar)

df <- data.frame(
  age = rep(age, 4),
  survival = c(surv, ph_surv, aft_surv, ph_match_surv),
  Species = factor(rep(c("Human", "Dog (PH, HR=5)", "Dog (AFT)", lab_match), each = 100),
    levels = c("Human", "Dog (AFT)", "Dog (PH, HR=5)", lab_match)))

g <- ggplot(df, aes(x = age, y = survival, color = Species, linetype = Species)) +
  geom_line() + xlim(0, 80) + labs(x = "Time", y = "Survival Probability") +
  scale_color_manual(values = setNames(c("black", "red", "blue", "blue"),
    c("Human", "Dog (AFT)", "Dog (PH, HR=5)", lab_match))) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "solid", "dashed"),
    c("Human", "Dog (AFT)", "Dog (PH, HR=5)", lab_match))) +
  theme(aspect.ratio = 1, text = element_text(size = 13))  # square + text +2

save_fig(g, "book/Figures/classical/fig-p3c13-dogs.png",
  width = 6, height = 5)
