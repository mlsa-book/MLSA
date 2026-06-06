remotes::install_github("mlr-org/mlr3proba", ref = 'v0.5.7', upgrade = "never")
remotes::install_github("mlr-org/mlr3", ref = 'v0.16.1', upgrade = "never")
remotes::install_github("mlr-org/paradox", ref = 'v0.11.1', upgrade = "never")

rm(list = ls())
library(tidyverse)
theme_set(theme_bw())
library(distr6)
library(mlr3)
library(mlr3proba)
library(mlr3pipelines)
library(mlr3extralearners)
library(survival)
library(patchwork)
library(survminer)
library(party)
library(pammtools)
library(extraDistr)
library(mvna)
library(etm)
library(cmprsk)
library(mstate) #prothr dataset
library(ggpubr)
library(pseudo)
library(mgcv)
library(broom)
library(eha) # swedeaths / swepop

## Intro - Gompertz fit to real Swedish age-specific mortality (2019, ages 30-100).
## Data source: eha::swedeaths and eha::swepop (Statistics Sweden, bundled in eha).
## Output: book/Figures/introduction/gompertz.png
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

# Gompertz hazard h(x) = a * exp(b * x), fit by population-weighted log-linear
# regression on the age-specific mortality rate.
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

fun_levels_gp <- c("f(t): density",
                   "h(t): hazard",
                   "F(t): CDF",
                   "S(t): survival")
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
  theme_bw(base_size = 11) +
  labs(x = "age (years)", y = NULL)

ggsave("book/Figures/introduction/gompertz.png", g_gp,
       height = 3.5, width = 7, units = "in", dpi = 600)

## Ranking
s_t = tsk("whas")
time = s_t$unique_times()
c_t = s_t$data() %>% mutate(status = 1 - status) %>% as_task_surv(event = "status")

s_d = data.frame(t = time, surv = s_t$kaplan()$surv, W = "KMS")
c_d = data.frame(t = time, surv = c_t$kaplan()$surv, W = "KMG")
w_d = data.frame(t = time, surv = 1 / (c_t$kaplan()$surv^2), W = "KMG^-2")

cutoff = time[which(s_t$kaplan()$surv < 0.6)[1]]

d = rbind(s_d, c_d, w_d) %>% as.data.frame()
g = ggplot(d, aes(x = t, y = surv, color = W)) +
  geom_line() +
  ylim(0, 5) +
  theme_classic() +
  labs(
    y = "W(t)",
    title = "Kaplan-Meier estimates and weighting on 'whas' data",
    color = "Weight function") +
  geom_vline(xintercept = cutoff, lty = 2, color = "gray") +
  scale_color_discrete(labels = expression(hat(G)[KM], hat(G)[KM]^{-2}, hat(S)[KM]))
ggsave("book/Figures/evaluation/weights.png", g, height = 3, units = "in",
  dpi = 600)

ids = c("W=1", "W=G^-1", "W=G^-2")
m_inf = c(msr("surv.cindex"),
msr("surv.cindex", weight_meth = "G"),
msr("surv.cindex", weight_meth = "G2"))

m_80 = c(msr("surv.cindex", id = "W=1", cutoff = cutoff),
msr("surv.cindex", weight_meth = "G", id = "W=G^-1", cutoff = cutoff),
msr("surv.cindex", weight_meth = "G2", id = "W=G^-2", cutoff = cutoff))

m = c(m_inf, m_80)

set.seed(20231207)
round(resample(s_t, lrn("surv.coxph"), rsmp("cv", folds = 3))$aggregate(m), 2)

## Calibration

set.seed(20231211)
t = tgen("simsurv")$generate(400)
s = partition(t)

prrt = as_learner(ppl("distrcompositor", lrn("surv.rpart")))$
  train(t, s$train)$predict(t, s$test)
pcox = lrn("surv.coxph")$train(t, s$train)$predict(t, s$test)
pran = lrn("surv.ranger")$train(t, s$train)$predict(t, s$test)

drrt = autoplot(prrt, "calib", t, s$test)$data
drrt = drrt %>% filter(Group == "Pred") %>% mutate(Group = "RRT")
dcox = autoplot(pcox, "calib", t, s$test)$data
dcox = dcox %>% mutate(Group = if_else(Group == "KM", "KM", "CPH"))
dran = autoplot(pran, "calib", t, s$test)$data
dran = dran %>% filter(Group == "Pred") %>% mutate(Group = "RF")

g = ggplot(rbind(dcox, dran, drrt), aes(x = x, y = y, color = Group)) +
  geom_line() + theme_classic() + ylim(0, 1) +
  labs(x = "T", y = "S(T)", color = "Model")
ggsave("book/Figures/evaluation/calibKM.png", g, height = 3, units = "in", dpi = 600)

drrt = autoplot(prrt, "dcalib", t, s$test, extend_quantile = TRUE)$data
drrt = drrt %>% mutate(Group = "RRT")
dcox = autoplot(pcox, "dcalib", t, s$test, extend_quantile = TRUE)$data
dcox = dcox %>% mutate(Group = "CPH")
dran = autoplot(pran, "dcalib", t, s$test, extend_quantile = TRUE)$data
dran = dran %>% mutate(Group = "RF")

dcal_coxp = as.numeric(pcox$score(msr("surv.dcalib", truncate = Inf, chisq = TRUE)))
dcal_cox = as.numeric(pcox$score(msr("surv.dcalib", truncate = Inf, chisq = FALSE)))

dcal_ranp = as.numeric(pran$score(msr("surv.dcalib", truncate = Inf, chisq = TRUE)))
dcal_ran = as.numeric(pran$score(msr("surv.dcalib", truncate = Inf, chisq = FALSE)))

dcal_rrtp = as.numeric(prrt$score(msr("surv.dcalib", truncate = Inf, chisq = TRUE)))
dcal_rrt = as.numeric(prrt$score(msr("surv.dcalib", truncate = Inf, chisq = FALSE)))

scores = paste0(sprintf("   %s = %s (%s)", c("CPH", "RF", "RRT"), signif(c(dcal_cox, dcal_ran, dcal_rrt), 2), signif(c(dcal_coxp, dcal_ranp, dcal_rrtp), 2)), collapse = "\n")
scores = paste0("DCal (p-values):\n", scores)

g = ggplot(rbind(dcox, dran, drrt), aes(x = p, y = q, color = Group)) +
  geom_line() + theme_classic() + ylim(0, 1) + xlim(0, 1) +
  geom_abline(slope = 1, intercept = 0, color = "lightgray", lty = "dashed") +
  labs(x = "True (p)", y = "Predicted", color = "Model") +
  geom_label(aes(x = x, y = y), data.frame(x = 0.75, y = 0.1), label = scores,
  inherit.aes = FALSE, hjust = "left", size = 2.5)
ggsave("book/Figures/evaluation/calibD.png", g, height = 3, units = "in",
  dpi = 600)

## Decision trees

set.seed(20241104)

data <- read.csv("book/experiments/cars.csv")[, -1]
colnames(data)[c(2, 3, 5)] <- c("price", "km", "seller")
# convert INR to $000
data$price <- data$price/73000
data$km <- data$km/1000
data <- data[data$fuel %in% c("Diesel", "Petrol"), ]

train <- sample(nrow(data), nrow(data) * 2/3)
train_set <- data[train, ]
test_set <- data[setdiff(seq(nrow(data)), train), ]

fit <- rpart::rpart(price ~ ., train_set, maxdepth = 2)

# test model is 'good' before plotting
pred <- predict(fit, test_set)
base <- mean(train_set$price)
truth <- test_set$price
rmse <- function(yhat) sqrt(mean((truth - yhat)^2))
mae <- function(yhat) mean(abs(truth - yhat))
matrix(
  c(rmse(base), rmse(pred), mae(base), mae(pred)),
  2, 2, FALSE,
  list(c("Baseline", "Prediction"), c("RMSE", "MAE"))
)

png("book/Figures/forests/cars.png", height = 400, width = 600)
rattle::fancyRpartPlot(fit, caption = "")
dev.off()

# Log-rank test
set.seed(20241125)
sf0 <- survfit(Surv(time, status) ~ 1, lung)
p0 <- ggplot(rbind(data.frame(x=sf0$time, y=sf0$surv), data.frame(x=0,y=1)), aes(x=x,y=y))+geom_step()

opt_1 <- lung$age > 50
opt_2 <- lung$age > 75
df_1 <- cbind(lung, split = opt_1)
df_2 <- cbind(lung, split = opt_2)
logrank_1 <- survdiff(Surv(time, status) ~ split, data = df_1)
logrank_2 <- survdiff(Surv(time, status) ~ split, data = df_2)
sf_1 <- survfit(Surv(time, status) ~ split, df_1)
sf_2 <- survfit(Surv(time, status) ~ split, df_2)
p1 <- ggplot(
  rbind(data.frame(x=sf_1$time, y=sf_1$surv,group=c(rep(FALSE, sf_1$strata[1]), rep(TRUE, sf_1$strata[2]))), data.frame(x = 0, y = 1, group=c(TRUE, FALSE))),
  aes(x=x,y=y,color=group)) +
  geom_step() +
  geom_label(aes(x=x,y=y), data.frame(x = 750, y = 0.9), label =sprintf("χ²=%.2f (p=%.2f)", logrank_1$chisq, logrank_1$pvalue), inherit.aes = FALSE)
p2 <- ggplot(
  rbind(data.frame(x=sf_2$time, y=sf_2$surv,group=c(rep(FALSE, sf_2$strata[1]), rep(TRUE, sf_2$strata[2]))), data.frame(x = 0, y = 1, group=c(TRUE, FALSE))),
  aes(x=x,y=y,color=group)) +
  geom_step() +
  geom_label(aes(x=x,y=y), data.frame(x = 750, y = 0.9), label = sprintf("χ²=%.2f (p=%.2f)", logrank_2$chisq, logrank_2$pvalue), inherit.aes = FALSE)

g <- p0 / (p1 + p2) &
  labs(x = "t", y = "S(t)") &
  theme_classic() &
  guides(color = "none") &
  scale_x_continuous(limits = c(0, 1000),  expand = c(0, 0)) &
  plot_annotation(tag_levels = 'a', tag_suffix = ")")

ggsave("book/Figures/forests/logrank.png", g, height = 6, units = "in",
  dpi = 600)


# Random forest plot

set.seed(20241109)
fit <- party::ctree(Surv(time, status) ~ ., lung)
png("book/Figures/forests/lung.png", height = 400, width = 600)
plot(fit)
dev.off()

## bootstrapped rsfs



x = 0:13
y1 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.9, 0.8, 0.1, 0.1)
y2 = c(1, 1, 0.6, 0.6, 0.6, 0.6, 0.3, 0.3, 0.2, 0.2, 0, 0, 0, 0)
y3 = c(1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2)
group = rep(c("blue", "red", "green"), each = 14)
df = data.frame(x = x, y = c(y1, y2, y3), group = group)

p0 <- ggplot(df, aes(x = x, y = y, color = group))
p1 <-  p0 + geom_step()
p2 <- p0 + geom_vline(xintercept = x, lty = 2, color = "gray") + geom_step()
p3 <- p0 + geom_vline(xintercept = x, lty = 2, color = "gray") +  geom_point() +
  geom_label(aes(x=x,y=y),
    data.frame(x=5,y=c(0.9, 0.2, 0.4),group=c("blue", "red", "green")),
    label = sprintf("S(6) = %s", c(1, 0.3, 0.5)))

df2 = data.frame(x = x, y = apply(data.frame(y1, y2, y3), 1, mean))
p4 <- ggplot(df2, aes(x = x, y = y)) + geom_step() + geom_point() +
  geom_label(aes(x=x,y=y),
    data.frame(x=6,y=0.5), label = "S(6) = 0.6")

ybreaks = seq.int(0, 1, 0.25)
xbreaks = seq.int(0, 12, 3)
g = p1 + p2 + p3 + p4 &
  labs(x = "t", y = "S(t)") &
  theme_classic() &
  guides(color = "none") &
  scale_y_continuous(limits = c(0, 1), breaks = ybreaks, labels = ybreaks) &
  scale_x_continuous(limits = c(0, 14), breaks = xbreaks, labels = xbreaks, expand = c(0, 0)) &
  plot_annotation(tag_levels = 'a', tag_suffix = ")")

ggsave("book/Figures/forests/bootstrap.png", g, height = 6, units = "in",
  dpi = 600)

## Kaplan Meier

data("tumor", package = "pammtools")
tumor <- cbind(id = seq_len(nrow(tumor)), tumor)
tumor_duplicated = tumor |>
  filter(days %in% days[duplicated(days)]) |>
  arrange(days)

## Table for illustration of right-censored data
tab_surv_tumor = tumor_duplicated |>
  filter(id %in% c(13, 62, 185, 230, 431, 719)) |>
  select(id, age, sex, complications, days, status) |>
  arrange(id)
knitr::kable(tab_surv_tumor)

km = survfit(Surv(days, status)~1, data = tumor)
bkm = broom::tidy(km)

df_med = data.frame(
  x = c(0, median(km)), # Starting x-coordinates
  y = c(0.5, 0),        # Starting y-coordinates
  xend = c(median(km), median(km)),    # Ending x-coordinates
  yend = c(0.5, 0.5)       # Ending y-coordinates
)

p_km = ggplot(bkm, aes(x = time, y = estimate)) +
  geom_step() +
  geom_segment(data=df_med, aes(x=x, xend=xend, y=y, yend=yend), lty = 3)+
  ylim(c(0, 1)) +
  ylab("S(t)") +
  xlab("days") +
  theme_bw()
p_km
ggsave("book/Figures/survival/km-tumor.png", p_km, height=3, units="in", dpi=600)

# stratified KM wrt age group
tumor = tumor |>
  mutate(age_bin = factor(age < 50, levels = c(TRUE, FALSE), labels = c("< 50", "≥ 50")))
km_age_bin = survfit(Surv(days, status)~age_bin, data = tumor)
bkm_age_bin = broom::tidy(km_age_bin) |>
  mutate(age = sub("age_bin=", "", strata))
med_km_age_bin = as.numeric(median(km_age_bin))
df_age_bin = data.frame(
  x = c(0, med_km_age_bin[2]), # Starting x-coordinates
  y = c(0.5, 0),        # Starting y-coordinates
  xend = c(med_km_age_bin[2], med_km_age_bin[2]),    # Ending x-coordinates
  yend = c(0.5, 0.5)       # Ending y-coordinates
)

p_km_age_bin = ggplot(bkm_age_bin, aes(x = time, y = estimate)) +
  geom_step(aes(col = age)) +
  geom_segment(data=df_age_bin, aes(x=x, xend=xend, y=y, yend=yend), lty = 3)+
  geom_hline(yintercept =  .5, lty = 3) +
  ylim(c(0, 1)) +
  ylab("S(t)") +
  xlab("days") +
  scale_color_discrete(name = "age") +
  theme_bw()
p_km_age_bin
ggsave("book/Figures/survival/km-age-bin-tumor.png", p_km_age_bin, height=3, units="in", dpi=600)


## Left-truncation
data("infants", package = "eha")

# KM for infants with dead/alive mothers
km_infants = survfit(Surv(exit, event)~mother, data = infants)
bkm_infants = broom::tidy(km_infants)

p_km_infants = ggplot(bkm_infants, aes(x = time, y = estimate)) +
  geom_step(aes(col = strata)) +
  ylim(c(0, 1)) +
  ylab("S(t)") +
  xlab("time")
# adjusted for left-truncation
km_infants_lt = survfit(Surv(enter, exit, event)~mother, data = infants)
bkm_infants_lt = broom::tidy(km_infants_lt)

p_km_infants_lt = ggplot(bkm_infants_lt, aes(x = time, y = estimate)) +
  geom_step(aes(col = strata)) +
  ylim(c(0, 1)) +
  ylab("S(t)") +
  xlab("time")


p_km_infants_joined = p_km_infants + p_km_infants_lt + plot_layout(guides =  "collect")
ggsave("book/Figures/survival/km-infants.png", p_km_infants_joined, height=3, width=7, units="in", dpi=600)


## table infant data
inf_sub = infants |>
  filter(stratum %in% c(1, 2, 4)) |>
  select(stratum, enter, exit, event, mother)

inf_sub |> knitr::kable()

# Chapter 13 - Traditional models
set.seed(2029)


t <- tsk("rats")$filter(sample(tsk("rats")$nrow, 5))
t$kaplan()$surv
knitr::kable(t$data()[, c(3,4,5,1,2)],align = "l")

# PH vs AFT
set.seed(290125)
t = seq(0.1, 2.5, by = 0.02)

hweib = function(shape, scale, eta, t, form) {
  if (form == "AFT") {
    mod_scale = scale * exp(-eta)
  } else if (form == "PH") {
    mod_scale = scale * exp(-eta/shape)
  }
  (shape/mod_scale) * (t/mod_scale)^(shape-1)
}

sweib = function(shape, scale, eta, t, form) {
  if (form == "AFT") {
    mod_scale = scale * exp(-eta)
  } else if (form == "PH") {
    mod_scale = scale * exp(-eta/shape)
  }
  exp(-(t/mod_scale)^shape)
}



# Modify shape and scale parameters for better intercepts
plotWeib = function(type = c("hazard", "survival"), shape = 3, scale = 2) {
  type = match.arg(type)

  fun = ifelse (type == "hazard", hweib, sweib)
  baseline = fun(shape, scale, 0, t, "PH")
  PH = fun(shape, scale, log(2), t, "PH")
  AFT = fun(shape, scale, log(2), t, "AFT")
  ylabel = ifelse(type == "hazard", "h(t)", "S(t)")
  
  df = data.frame(
    y = c(baseline, PH, AFT),
    t = rep(t, 3),
    Model = factor(rep(c("Baseline", "PH", "AFT"), 
                      each = length(t)), levels = c("Baseline", "PH", "AFT"))
  )

  g <- ggplot(df, aes(x = t, y = y, color = Model)) +
    geom_line(aes(alpha = if_else(
      type == "hazard" & Model == "AFT" |
      type == "survival" & Model == "PH",
      0, 1))) +
    theme(legend.position = "bottom") +
    guides(alpha = "none") +
    ylab(ylabel) +
    scale_color_manual(values = c("Baseline" = "black", "AFT" = "red", "PH" = "blue"), aesthetics = c("color","fill"))
}


segment = function(start, form) {
  if (form == "PH") {
    map <- aes(
      x = start,
      xend = start,
      y = hweib(3, 2, 0, start, "PH"),
      yend = hweib(3, 2, log(2), start, "PH")
    )
  } else if (form == "AFT") {
    map <- aes(
      x = start,
      xend = start * 2,
      y = sweib(3, 2, 0, start * 2, "AFT"), # sense check
      yend = sweib(3, 2, log(2), start, "AFT")
    )
  }

  geom_segment(map,
    arrow = arrow(ends = "both", length = unit(0.1, "in")),
    inherit.aes = FALSE, size = 0.3,
    color = "#6f6f6f"
  )
}

p1 = plotWeib("hazard") +
  ylim(0, 5) +
  segment(1, "PH") + segment(1.5, "PH") + segment(2, "PH") +
  geom_label(aes(x = x, y = y), data.frame(x = 1.45, y = hweib(3, 2, log(2), 2, "PH")), 
    label = expression(h[PH](t)==2*h[0](t)),
    inherit.aes = FALSE)

p2 = plotWeib("survival") +
  ylim(0, 1) +
  segment(0.5, "AFT") + segment(0.75, "AFT") +
  segment(1, "AFT") +
  geom_label(aes(x = x, y = y), data.frame(x = 2.08, y = sweib(3, 2, 0, 1.5, "AFT")), 
    label = expression(S[AFT](t)==S[0]("2t")),
    inherit.aes = FALSE)

g <- (p1 + p2) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave("book/Figures/classical/compare.png", g, height = 4, units = "in",
  dpi = 600)

## Humans vs dogs
age = seq.int(1, 100, 1)
surv = pgompertz(age, 0.00005, 0.09, FALSE)
ph_surv = surv^5
aft_surv = round(pgompertz(age*5, 0.00005, 0.09, FALSE), 2)
df = data.frame(age, survival = c(surv, ph_surv, aft_surv), Species = rep(c("Human", "Dog (PH)", "Dog (AFT)"), each = 100))

g <- ggplot(df, aes(x = age, y = survival, group = Species, color = Species)) + geom_line() + xlim(0, 80) + labs(x = "T", y = "S(T)") +
  scale_color_manual(values = c("Human" = "black", "Dog (AFT)" = "red", "Dog (PH)" = "blue"), aesthetics = c("color","fill"))

ggsave("book/Figures/classical/dogs.png", g, height = 4, units = "in",
  dpi = 600)


## KM for testing


fit = survfit(Surv(rats$time, rats$status) ~ 1)
fit$time[1:4]
fit$surv
fit$time
g = ggplot(data.frame(x = fit$time,y = fit$surv), aes(x = x, y =y)) +
  geom_step() + labs(x = "t", y = "S(t)") +
  scale_x_continuous(expand = c(0, 0))

g1 = g +
  geom_vline(xintercept = fit$time[5:7], lty = 2, alpha = 1, color = 3, lwd = 1) +
  geom_vline(xintercept = fit$time[9:10], lty = 3, alpha = 1,color = 4,lwd=1)

g2 = g +
  geom_segment(x = 60, y = 0, yend = fit$surv[12], color = 2, lwd = 1, arrow = arrow()) +
  geom_segment(x = 23, xend = fit$time[12], y = fit$surv[12], color = 2, lwd = 1, arrow = arrow(ends = "first")) 

g3 = g1 / g2  

ggsave("book/Figures/classical/km_test.png", g3, height = 6.5, units = "in",
  dpi = 600)


## competing risks
set.seed(241206)
### table
data(sir.adm, package = "mvna")
tab_sir = sir.adm |> group_by(status) |> sample_n(2) |>
  ungroup() |>
  mutate(id = row_number()) |>
  select(id, time, status, pneu)
tab_sir |> knitr::kable()


### AJ CIF estimates

sir_data_cif = sir.adm |>
  mutate(
    from = 0,
    to = dplyr::case_when(
      status == 0 ~ "cens",
      .default = as.character(status)))

sir_data_cif |> dplyr::pull(status) |> table()

cif = cuminc(
  sir_data_cif$time,
  sir_data_cif$to,
  group = sir_data_cif$pneu,
  cencode="cens")

cif_sir_b = purrr::imap_dfr(cif[1:4], function(.x, .y) {
  as.data.frame(.x) |>
    cbind(
      pneumonia = stringr::str_sub(.y, 1, 1),
      transition = stringr::str_sub(.y, 3, 3)
    )
}) |>
  mutate(
    method = "Aalen-Johansen",
    assumption = "competing risks censoring"
  ) |>
  rename(cif = est) |>
  mutate(
    pneumonia = factor(pneumonia, labels = c("no", "yes")),
    transition = factor(transition, labels = c("discharge", "death")))


km_sir_b = broom::tidy(survfit(Surv(time, status!=0)~pneu, data = sir.adm)) |>
  rename(pneumonia = strata) |>
  mutate(pneumonia= dplyr::case_when(
    pneumonia == "pneu=0" ~ "no",
    pneumonia == "pneu=1" ~ "yes"
  )) |>
  mutate(
    method = "KM",
    assumption="independent right-censoring",
    transition = factor("death", levels=c("discharge", "death")),
    cif = 1-estimate
  )


p_sir_cifs = ggplot(cif_sir_b, aes(x = time, y = cif)) +
  geom_step(aes(col = pneumonia)) +
  facet_wrap(~transition) +
  geom_vline(xintercept = 120, lty = 3) +
  # geom_step(data=km_sir_b, aes(col = pneumonia), lty = 2) +
  labs(
    y = expression(P(Y <= tau~ "," ~ E(Y) == e))
  ) +
  coord_cartesian(xlim = c(0, 125), ylim=c(0, 1))


ggsave("book/Figures/survival/cif-sir.png", p_sir_cifs,
  height=3, width=6, dpi=300)


### Independent censoring vs. Competing Risks

cox_sir = purrr::map_dfr(
  .x = 1:2,
  .f = function(.x) {

    tmp = sir.adm
    tmp$status = 1L*(tmp$status == .x)
    m = coxph(Surv(time, status)~strata(pneu), data = tmp)
    bm = basehaz(m) |>
      rename(pneumonia = strata) |>
      mutate(
        cif = 1 - exp(-hazard),
        pneumonia = dplyr::case_when(
          pneumonia == "pneu=0" ~ "no",
          pneumonia == "pneu=1" ~ "yes"
        ),
        transition = ifelse(.x==1, "discharge", "death")
      )
  }
)

p_cens_vs_cr = ggplot(
  filter(cox_sir, transition == "death"),
  aes(x = time, y = cif)) +
    geom_step(aes(col = pneumonia, lty="independent censoring")) +
    geom_step(
      data = filter(cif_sir_b, transition == "death"),
      aes(col=pneumonia, lty="competing risks")
    ) +
    coord_cartesian(xlim = c(0, 125), ylim=c(0, 1)) +
    geom_vline(xintercept = 120, lty = 3) +
    labs(
      y = expression(P(Y <= tau~ "," ~ E(Y) == 2)),
      linetype = "assumption"
    )


ggsave("book/Figures/survival/cens-vs-cr.png", p_cens_vs_cr, height=3, width=6, dpi=300)



### multi-state example
data(prothr, package = "mstate")
prothr |>
  filter(id %in% c(1, 8, 46)) |>
  mutate(from = from -1, to = to -1) |>
  knitr::kable()

my.prothr <- prothr |>
  mutate(transition = as.factor(paste0(from, "->", to))) |>
  filter(Tstart != Tstop) # filter instantaneous transitions

#' We use cox estimation without covariates because this is essentially `Nelson-Aalen`!!!!!)
tmat = matrix(NA, nrow = 3, ncol = 3)
tmat[1,2] <- 1
tmat[1,3] <- 2
tmat[2,1] <- 3
tmat[2,3] <- 4
tmat
# placebo
cox.bm.placebo = coxph(
  Surv(Tstart,Tstop,status) ~ strata(trans),
  data = my.prothr |> dplyr::filter(treat == "Placebo"),
  method = "breslow")
#' `msfit` has `newdata` arg => patient-specific transition hazards
haz.bm.placebo = msfit(cox.bm.placebo, trans = tmat) # cumhaz
tp.placebo = probtrans(haz.bm.placebo, predt=0, direction = "forward")
# prednisone
cox.bm.prednisone = coxph(
  Surv(Tstart,Tstop,status) ~ strata(trans),
  data = my.prothr |> dplyr::filter(treat == "Prednisone"),
  method = "breslow")
#' `msfit` has `newdata` arg => patient-specific transition hazards
haz.bm.prednisone = msfit(cox.bm.prednisone, trans = tmat) # cumhaz
tp.prednisone = probtrans(haz.bm.prednisone, predt=0, direction = "forward")

time.placebo = tp.placebo[[1]]$time
time.prednisone = tp.prednisone[[1]]$time
placebo.df = data.frame(
  time = rep(c(time.placebo), 4),
  transition = rep(c("0->1", "0->2", "1->0", "1->2"), each = length(time.placebo)),
  Treatment = factor("Placebo", levels = c("Placebo", "Prednisone"))
)
placebo.df$trans_prob = c(tp.placebo[[1]]$pstate2, tp.placebo[[1]]$pstate3, tp.placebo[[2]]$pstate1, tp.placebo[[2]]$pstate3)
prednisone.df = data.frame(
  time = rep(c(time.prednisone), 4),
  transition = rep(c("0->1", "0->2", "1->0", "1->2"), each = length(time.prednisone)),
  Treatment = factor("Prednisone", levels = c("Placebo", "Prednisone"))
)
prednisone.df$trans_prob = c(tp.prednisone[[1]]$pstate2, tp.prednisone[[1]]$pstate3, tp.prednisone[[2]]$pstate1, tp.prednisone[[2]]$pstate3)

overall_df = rbind(placebo.df, prednisone.df)

p_trans_prob_prothr = ggplot(overall_df, aes(x = time, y = trans_prob)) +
  facet_wrap(~transition) +
  geom_step(aes(col = Treatment)) +
  scale_color_manual(values = c("steelblue", "firebrick4")) +
  coord_cartesian(xlim = c(0, 4000), ylim = c(0, 1)) +
  theme(legend.position = "bottom") +
  labs(x = "Time", y = "Transition probability")

ggsave("Figures/survival/multi-state-prothr.png",
  p_trans_prob_prothr, width=5, height=5.2)


## P4C23 — multi-state reduction comparison on prothr (with treatment).
## Two transition-probability estimators, both stratified by treatment arm:
##   AJ      : non-parametric Aalen-Johansen per arm via
##             coxph(~ strata(trans)) + mstate::msfit + mstate::probtrans
##   XGBoost : single XGBoost-Poisson learner trained on the stacked
##             multi-state PED with features tend + treat + transition.
##             First reduction (multi-state -> transition-specific
##             single-event LT data, see @sec-ms-transition-datasets) is
##             folded into the same as_ped() call as the second reduction
##             (single-event LT -> Poisson regression on PED via the
##             partition-based reduction, @sec-partition-based-reductions).
##             Offset = log(intlen) is set as base_margin during TRAINING
##             ONLY; at prediction time base_margin is unset so the model
##             returns the hazard directly. No regularisation; small eta +
##             many rounds + early stopping (defensive — runs to nrounds on
##             training-set watchlist).
## Output: book/Figures/reductions/tp-prothr-cmp.png

## Reuse my.prothr / tmat from the AJ-baseline block above.

# --- AJ split by treatment ---------------------------------------------------
fit_aj_arm <- function(arm) {
  d <- my.prothr |> dplyr::filter(treat == arm)
  probtrans(
    msfit(coxph(Surv(Tstart, Tstop, status) ~ strata(trans),
                data = d, method = "breslow"), trans = tmat),
    predt = 0, direction = "forward")
}
tp_aj <- list(Placebo    = fit_aj_arm("Placebo"),
              Prednisone = fit_aj_arm("Prednisone"))

# --- XGBoost-Poisson via the partition-based reduction (PED) ----------------
library(pammtools); library(xgboost)
set.seed(20260516)

prothr_ped <- prothr |>
  dplyr::filter(Tstart != Tstop) |>
  mutate(transition = as.factor(paste0(from, "->", to)),
         treat      = as.factor(treat)) |>
  rename(tstart = Tstart, tstop = Tstop) |>
  select(-trans)
ped <- as_ped(
  data       = prothr_ped,
  formula    = Surv(tstart, tstop, status) ~ .,
  transition = "transition",
  id         = "id",
  timescale  = "calendar")

## Features: tend + treat + transition.
make_X <- function(df) {
  model.matrix(~ tend + treat + transition, data = df)[, -1, drop = FALSE]
}
ped_df <- as.data.frame(ped)
dtrain <- xgb.DMatrix(data = make_X(ped_df),
                      label = as.numeric(ped_df$ped_status))
setinfo(dtrain, "base_margin", ped_df$offset)
xgb_fit <- xgb.train(
  params = list(
    objective        = "count:poisson",
    base_score       = 1,
    eta              = 0.01,
    max_depth        = 3,
    lambda           = 0,
    alpha            = 0,
    gamma            = 0,
    min_child_weight = 1
  ),
  data                  = dtrain,
  nrounds               = 5000,
  watchlist             = list(train = dtrain),
  early_stopping_rounds = 10,
  verbose               = 0
)

ndf <- make_newdata(ped,
                    tend       = unique(tend),
                    treat      = unique(treat),
                    transition = unique(transition))
brks <- attr(ped, "breaks")
ndf$intlen <- setNames(diff(c(0, brks)), brks)[as.character(ndf$tend)]
ndf$hazard <- predict(xgb_fit, xgb.DMatrix(data = make_X(as.data.frame(ndf))))
ndf <- ndf |>
  group_by(treat, transition) |>
  arrange(treat, transition, tend) |>
  mutate(cumu_hazard = cumsum(hazard * intlen)) |>
  ungroup()
tp_xgb_df <- ndf |>
  group_by(treat) |>
  add_trans_prob(object = NULL, ci = FALSE) |>
  ungroup() |>
  transmute(treat,
            transition = as.character(transition),
            time       = tend,
            trans_prob)

# --- Combine into a long-format frame and plot -------------------------------
build_aj_long <- function(tp, arm) {
  ts <- tp[[1]]$time
  data.frame(
    time       = rep(ts, 4),
    transition = rep(c("0->1", "0->2", "1->0", "1->2"), each = length(ts)),
    Treatment  = factor(arm, levels = c("Placebo", "Prednisone")),
    method     = factor("AJ", levels = c("AJ", "XGBoost")),
    trans_prob = c(tp[[1]]$pstate2, tp[[1]]$pstate3,
                   tp[[2]]$pstate1, tp[[2]]$pstate3))
}
relabel_mstate <- c("1->2" = "0->1", "1->3" = "0->2",
                    "2->1" = "1->0", "2->3" = "1->2")

df_tp_cmp <- bind_rows(
  build_aj_long(tp_aj$Placebo,    "Placebo"),
  build_aj_long(tp_aj$Prednisone, "Prednisone"),
  tp_xgb_df |>
    mutate(transition = unname(relabel_mstate[transition]),
           Treatment  = factor(treat, levels = c("Placebo", "Prednisone")),
           method     = factor("XGBoost", levels = c("AJ", "XGBoost"))) |>
    select(time, transition, Treatment, method, trans_prob)
)

p_tp_cmp <- ggplot(df_tp_cmp,
                   aes(x = time, y = trans_prob,
                       colour = Treatment, linetype = method)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ transition) +
  scale_colour_manual(values = c("steelblue", "firebrick4")) +
  scale_linetype_manual(values = c("AJ" = "solid", "XGBoost" = "dotted")) +
  coord_cartesian(xlim = c(0, 4000), ylim = c(0, 1)) +
  labs(x = "Time", y = "Transition probability", linetype = "Method") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("Figures/reductions/tp-prothr-cmp.png", p_tp_cmp,
       width = 8, height = 6.5, dpi = 300)


## Survival task viz


set.seed(20250822)

sex <- rbinom(10, 1, 0.5)
x_base <- ceiling(runif(10, ifelse(sex == 1, 20, 1), ifelse(sex == 1, 50, 30)))
x <- unlist(lapply(x_base, function(val) c(0, val, 50)))

df <- data.frame(x = x, y = c(1, 0, 0), group = as.factor(rep(1:10, each = 3)), sex = rep(as.factor(sex), each  = 3)) %>%
  mutate(alpha = if_else(group == 10, 1, 0.1))

# Fit Weibull models to the simulated event times (treating them as fully
# observed events) so the smooth aggregate curves in g3/g4 are derived from
# the same data as the step functions, not hard-coded.
events_df <- data.frame(time  = x_base,
                        event = 1L,
                        sex   = factor(sex, levels = levels(df$sex)))

fit_uncond <- flexsurv::flexsurvreg(Surv(time, event) ~ 1,
                                     data = events_df, dist = "weibull")
fit_cond   <- flexsurv::flexsurvreg(Surv(time, event) ~ sex,
                                     data = events_df, dist = "weibull")

tgrid <- seq(0, 50, length.out = 200)

pred_uncond <- summary(fit_uncond, t = tgrid, type = "survival", tidy = TRUE)
pred_uncond_df <- data.frame(x = pred_uncond$time, y = pred_uncond$est)

pred_cond_df <- do.call(rbind, lapply(levels(events_df$sex), function(s) {
  nd  <- data.frame(sex = factor(s, levels = levels(events_df$sex)))
  out <- summary(fit_cond, newdata = nd, t = tgrid,
                 type = "survival", tidy = TRUE)
  data.frame(x = out$time, y = out$est,
             sex = factor(s, levels = levels(events_df$sex)))
}))

g <- ggplot(df, aes(x = x, y = y, group = group))
g1 <- g + geom_step(linewidth = 1.3, color = "gray")
g2 <- g + geom_step(aes(alpha = alpha), linewidth = 1.3) + scale_alpha_identity()
g3 <- g1 + geom_line(aes(x = x, y = y), data = pred_uncond_df,
                     inherit.aes = FALSE,
                     color = "black", linewidth = 1)
g4 <- g +
  geom_step(aes(color = sex), linewidth = 1.3, alpha = 0.5) +
  geom_line(aes(x = x, y = y, group = sex, color = sex),
            data = pred_cond_df,
            inherit.aes = FALSE, linewidth = 1.5)

g_final <- (g1 + g2 + g3 + g4) + ylim(0, 1) + xlim(0, 50) & labs(x = "t", y = "S(t)") & guides(color  = "none")

ggsave("book/Figures/survtsk/heavisides.png",
  g_final, height=5, width=7, units="in", dpi=600)

## Survival prediction types
set.seed(1)

n <- 5
dat <- data.frame(
  id    = 1:n,
  age   = round(rnorm(n, mean = 62, sd = 10)),
  sex   = sample(0:1, n, replace = TRUE),
  trt   = sample(0:1, n, replace = TRUE),
  x    = round(rnorm(n), 2),
  time  = round(runif(n, pmax(0, 3 + rnorm(n, sd = 3)), 15 + rnorm(n, sd = 3))),
  event = rbinom(n, 1, prob = 0.65)
)
# data plot
row_cols <- scales::hue_pal()(n)
p_table <- ggtexttable(
  dat,
  rows = NULL,
  theme = ttheme(
    base_size = 9,
    base_colour = "black",
    padding = unit(c(2, 2), "mm"),
    tbody.style = tbody_style(
      fill = rep(row_cols, ncol(dat)),  # repeat each row color across columns
      col = "black"
    )
  )
)

# survival time 'predictions'
times = dat$time + rnorm(n, 0, 2)
p_bars = ggplot(data.frame(y=factor(dat$id), x=times)) +
  geom_col(aes(x,y,fill=y), width = 0.7, alpha = 0.9) +
  labs(
    title = "Predicted survival times",
    x = "Time",
    y = "Subject"
  ) + theme(legend.position = "none")

# relative risk 'predictions'
risks = as.numeric(scale(-times + rnorm(n, 0, 2), TRUE, FALSE))
p_risks <- ggplot(data.frame(x = factor(dat$id), y = risks)) +
  geom_bar(aes(x=x, y=y, fill = x), stat='identity', position=position_dodge(width = 0)) +
  labs(
    title = "Predicted relative risk scores",
    x = "Subject",
    y = "Relative risk"
  ) + theme(legend.position = "none")

# -----------------------------
# Plot D (bottom-right): survival curves by treatment
# -----------------------------
p_curves <- data.frame(Subject = factor(1:n), x = rep(0:20, each = n), y = punif(rep(0:20, each=n), 3 + rnorm(n, 3), 15 + rnorm(n, 3), FALSE)) %>%
  ggplot(aes(x = x, y = y, group = Subject, color = Subject)) +
  ylim(0, 1) + xlim(5, 20) +
  geom_step(lwd = 1) +
  labs(
    title = "Predicted distributions",
    x = "Time",
    y = "Survival probability"
  ) + theme(legend.position = "none")


g <- ((p_table | p_bars) /
  (p_risks | p_curves)) 

ggsave("book/Figures/survtsk/predict_types.png",
  g, height=5, width=7, units="in", dpi=600)

######################## pseudo-values ########################################
data("tumor", package = "pammtools")
tumor_comp = tumor |> select(days, status, complications)
tau = c(1000, 2000, 3000)

## Pseudo-value calculation table for illustration
# Calculate pseudo-values for time point 1000 for first 4 subjects
tau_table = 1000
n = nrow(tumor)

# Calculate overall KM at time tau
km_all = survfit(Surv(days, status) ~ 1, data = tumor)
km_all_tau = summary(km_all, times = tau_table, extend = TRUE)$surv[1]
if (is.na(km_all_tau)) {
  # If tau is beyond last time, use last survival value
  km_all_tau = tail(km_all$surv, 1)
}

# Select 4 subjects for illustration
subjects = 1:4
results_table = data.frame(
  Subject = integer(),
  Days = integer(),
  Status = integer(),
  KM_Overall = numeric(),
  KM_Minus_i = numeric(),
  Pseudo_Value = numeric()
)

for (i in subjects) {
  # Calculate KM without subject i
  tumor_minus_i = tumor[-i, ]
  km_minus_i = survfit(Surv(days, status) ~ 1, data = tumor_minus_i)
  km_minus_i_tau = summary(km_minus_i, times = tau_table, extend = TRUE)$surv[1]
  if (is.na(km_minus_i_tau)) {
    km_minus_i_tau = tail(km_minus_i$surv, 1)
  }
  
  # Calculate pseudo-value: θ_i(τ) = n*S_KM(τ) - (n-1)*S_KM^{-i}(τ)
  pseudo = n * km_all_tau - (n - 1) * km_minus_i_tau
  
  results_table = rbind(results_table, data.frame(
    Subject = i,
    Days = tumor$days[i],
    Status = tumor$status[i],
    KM_Overall = round(km_all_tau, 4),
    KM_Minus_i = round(km_minus_i_tau, 4),
    Pseudo_Value = round(pseudo, 4)
  ))
}

results_table = results_table |> cbind(tumor[1:4, c("age", "complications")])

knitr::kable(results_table, format = "markdown")



### survival probability


pseudo_dfs = purrr::map_dfr(tau, function(.x) {
  pseudo_mat = pseudosurv(time = tumor$days, event = tumor$status, tmax = .x)
  cbind(tumor_comp, data.frame(pseudo = pseudo_mat$pseudo, time = pseudo_mat$time, age = tumor$age))
}) |> mutate(time = factor(time))

ndf = pammtools::make_newdata(pseudo_dfs, complications = unique(complications), time = unique(time))
lm <- lm(formula = pseudo~complications*time, data = pseudo_dfs)
summary(lm)

## Additive GAM: pseudo ~ complications + s(age) at each time point; store models in a list
gam_models = purrr::map(tau, function(.x) {
  d = pseudo_dfs |> filter(as.character(time) == as.character(.x))
  gam(pseudo ~ complications + age, data = d)
})
names(gam_models) = paste0("tau_", tau)

## Table of age and complications coefficients from GAMs at each time point
coef_age_table = purrr::map_dfr(seq_along(tau), function(i) {
  pt = summary(gam_models[[i]])$p.table
  data.frame(
    tau = tau[i],
    age = pt["age", "Estimate"],
    complications = pt["complicationsyes", "Estimate"]
  )
})
knitr::kable(coef_age_table, format = "markdown", digits = 4,
             col.names = c("τ (days)", "age", "complications"))

ndf$estimate = predict(lm, newdata = ndf)
ndf <- ndf |> 
  mutate(model = "pseudo-values (LM)") |> 
  mutate(time = as.numeric(as.character(time)))

# stratified KM wrt complications
km_complications = survfit(Surv(days, status)~complications, data = tumor)
bkm_complications = broom::tidy(km_complications) |> 
  mutate(model = "KM")


p_km_complications = ggplot(bkm_complications, aes(x = time, y = estimate)) +
  geom_step(aes(col = strata)) +
  geom_vline(xintercept = tau, lty = 3) +
  ylim(c(0, 1)) +
  ylab("S(t)") +
  xlab("time") + 
  geom_point(data = ndf, pch=19, size = 2)
p_km_complications
ggsave("book/Figures/survival/pseudo-complications-lm.png", p_km_complications, height=3, units="in", dpi=600)

### restricted mean survival time
data("tumor", package = "pammtools")
tumor_comp = tumor |> select(days, status, complications)
tau = c(1000, 2000, 3000)


rmst_dfs = purrr::map_dfr(tau, function(.x) {
  rmst = pseudomean(time = tumor$days, event = tumor$status, tmax = .x)
  cbind(tumor_comp, data.frame(pseudo = rmst, time = .x))
}) |> mutate(time = factor(time))

ndf = pammtools::make_newdata(rmst_dfs, complications = unique(complications), time = unique(time))
lm <- lm(formula = pseudo~complications*time, data = rmst_dfs)
summary(lm)

ndf$estimate = predict(lm, newdata = ndf)
ndf <- ndf |> 
  mutate(model = "pseudo-values (LM)") |> 
  mutate(time = as.numeric(as.character(time)))

tumor$pseudo = dplyr::filter(rmst_dfs, time == 1000)$pseudo
lm1000 = lm(pseudo ~ complications + age, data = tumor)
summary(lm1000)

# Get KM fits for each group and prepare data for plotting
km_complications = survfit(Surv(days, status)~complications, data = tumor)
bkm_complications = broom::tidy(km_complications) |>
  mutate(complications = gsub("complications=", "", strata))

# Use single tau for RMST visualization
tau_rmst = 1000

# Calculate RMST for each group at tau_rmst
calculate_rmst_at_tau = function(km_fit, tau) {
  times = km_fit$time[km_fit$time <= tau]
  surv = km_fit$surv[km_fit$time <= tau]
  
  # Add time 0 if needed
  if (length(times) == 0 || times[1] > 0) {
    times = c(0, times)
    surv = c(1, surv)
  }
  
  # Add tau if beyond last time point
  if (max(times) < tau) {
    times = c(times, tau)
    surv = c(surv, tail(surv, 1))
  }
  
  # Calculate RMST using trapezoidal rule
  rmst = 0
  for (i in 2:length(times)) {
    width = times[i] - times[i-1]
    avg_height = (surv[i] + surv[i-1]) / 2
    rmst = rmst + width * avg_height
  }
  
  return(rmst)
}

# Get KM fits for each group separately
km_no = survfit(Surv(days, status)~1, data = tumor |> filter(complications == "no"))
km_yes = survfit(Surv(days, status)~1, data = tumor |> filter(complications == "yes"))

rmst_no = calculate_rmst_at_tau(km_no, tau_rmst)
rmst_yes = calculate_rmst_at_tau(km_yes, tau_rmst)

# Prepare data for ribbon (shaded area up to tau_rmst)
bkm_for_ribbon = bkm_complications |>
  filter(time <= tau_rmst) |>
  group_by(complications) |>
  arrange(time) |>
  group_modify(~ {
    df = .x
    # Add time 0 if not present
    if (min(df$time) > 0) {
      df = bind_rows(data.frame(time = 0, estimate = 1), df)
    }
    # Add tau_rmst if needed
    if (max(df$time) < tau_rmst) {
      df = bind_rows(df, data.frame(time = tau_rmst, estimate = tail(df$estimate, 1)))
    }
    df
  }) |>
  ungroup()

# Create annotation data for RMST values
rmst_annotations = data.frame(
  complications = c("no", "yes"),
  rmst_value = c(rmst_no, rmst_yes),
  x = tau_rmst / 2,
  y = 0.25
) |>
  mutate(
    label = paste0("RMST [", round(rmst_value, 1), "d]")
  )

# Create plot with facets
p_rmst_complications = ggplot(bkm_complications, aes(x = time, y = estimate)) +
  # Shade area under curves (RMST) - up to tau_rmst
  geom_ribbon(data = bkm_for_ribbon, 
              aes(x = time, ymin = 0, ymax = estimate, fill = complications),
              alpha = 0.3, inherit.aes = FALSE) +
  # Survival curves
  geom_step(aes(color = complications), linewidth = 1.3) +
  # Vertical line at tau_rmst
  geom_vline(xintercept = tau_rmst, lty = 2, alpha = 0.7, linewidth = 0.8) +
  # RMST annotation
  geom_text(data = rmst_annotations, 
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  facet_wrap(~ complications, 
             labeller = labeller(complications = c("no" = "No Complications", "yes" = "Complications"))) +
  ylab("S(t)") +
  xlab("time (days)") +
  labs(color = "Complications", fill = "Complications") +
  ylim(c(0, 1)) +
  xlim(c(0, max(bkm_complications$time))) +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, face = "bold"))

# p_rmst_complications
ggsave("book/Figures/reductions/pseudo-rmst-complications-lm.png", p_rmst_complications, height=6, width=12, units="in", dpi=600)


## Discrete time model example




data("tumor", package = "pammtools")
tumor_comp = tumor |> select(days, status, complications)

# Create 100 equidistant cut points
max_time = max(tumor_comp$days)
cut_points = seq(0, max_time, length.out = 101)  # 101 points create 100 intervals
J = length(cut_points) - 1

# Transform data to long format using pammtools
ped_data = as_ped(
  data = tumor_comp,
  formula = Surv(days, status) ~ complications,
  cut = cut_points[-1],  # pammtools expects cut points without 0
  id = "id"
) |>
  mutate(interval = factor(interval))  # Convert interval to factor for separate baseline hazards

# For discrete time, we use ped_status as the outcome (event indicator in each interval)
# Rename for clarity in discrete time context
transformed_data = ped_data |>
  rename(delta_ij = ped_status)

# Fit logistic regression model
glm_fit = glm(delta_ij ~ interval + complications, 
              family = binomial(), 
              data = transformed_data)

# Create new data for predictions using distinct combinations from ped_data
# This is equivalent to make_newdata for creating prediction grids
pred_data = ped_data |>
  distinct(interval, complications, .keep_all = TRUE) |>
  select(interval, complications, tend)

# Add hazard predictions on response scale manually
pred_data$hazard = predict(glm_fit, newdata = pred_data, type = "response")

# Calculate survival probabilities from discrete hazards
# S(j) = prod_{k=1}^j (1 - h(k))
surv_data = pred_data |>
  arrange(complications, interval) |>
  group_by(complications) |>
  mutate(
    survival = cumprod(1 - hazard),
    time = tend  # Use tend from ped_data structure
  ) |>
  ungroup() |>
  select(complications, time, survival) |>
  mutate(model = "Discrete Time (GLM)")

# Add time 0 with survival = 1 and ensure proper ordering
surv_data = bind_rows(
  data.frame(
    complications = rep(unique(tumor_comp$complications), each = 1),
    time = 0,
    survival = 1,
    model = "Discrete Time (GLM)"
  ),
  surv_data
) |>
  arrange(complications, time)

# Get KM estimates for comparison
km_complications = survfit(Surv(days, status) ~ complications, data = tumor_comp)
bkm_complications = broom::tidy(km_complications) |>
  mutate(
    complications = gsub("complications=", "", strata),
    model = "Kaplan-Meier"
  ) |>
  select(complications, time, estimate, model) |>
  rename(survival = estimate)

# Combine data for plotting
plot_data = bind_rows(
  surv_data,
  bkm_complications
)

# Create plot - use color for model, linetype for complications
p_discrete_time = ggplot(plot_data, aes(x = time, y = survival, color = model, linetype = complications)) +
  geom_step(data = filter(plot_data, model == "Kaplan-Meier"), linewidth = 1.2) +
  geom_line(data = filter(plot_data, model == "Discrete Time (GLM)"), linewidth = 1.2) +
  ylab("Survival probability") +
  xlab("time (days)") +
  ylim(c(0, 1)) +
  scale_color_manual(
    values = c("Kaplan-Meier" = "black", "Discrete Time (GLM)" = "steelblue"),
    labels = c("Kaplan-Meier" = "Kaplan-Meier", "Discrete Time (GLM)" = "Discrete Time (GLM)")
  ) +
  scale_linetype_manual(
    values = c("no" = "solid", "yes" = "dashed"),
    labels = c("no" = "No Complications", "yes" = "Complications")
  ) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Model"), linetype = guide_legend(title = "Complications"))

ggsave("book/Figures/reductions/discrete-time-complications-glm.png", p_discrete_time, 
       height=4, width=8, units="in", dpi=600)


# Fit logistic regression model with interaction to allow separate baseline hazards for each group
glm_fit_interaction = glm(delta_ij ~ interval * complications, 
                          family = binomial(), 
                          data = transformed_data)

# Create new data for predictions using distinct combinations from ped_data
# This is equivalent to make_newdata for creating prediction grids
pred_data_interaction = ped_data |>
  distinct(interval, complications, .keep_all = TRUE) |>
  select(interval, complications, tend)

# Add hazard predictions on response scale manually
pred_data_interaction$hazard = predict(glm_fit_interaction, newdata = pred_data_interaction, type = "response")

# Calculate survival probabilities from discrete hazards
surv_data_interaction = pred_data_interaction |>
  arrange(complications, interval) |>
  group_by(complications) |>
  mutate(
    survival = cumprod(1 - hazard),
    time = tend  # Use tend from ped_data structure
  ) |>
  ungroup() |>
  select(complications, time, survival) |>
  mutate(model = "Discrete Time (GLM, Interaction)")

# Add time 0 with survival = 1 and ensure proper ordering
surv_data_interaction = bind_rows(
  data.frame(
    complications = rep(unique(tumor_comp$complications), each = 1),
    time = 0,
    survival = 1,
    model = "Discrete Time (GLM, Interaction)"
  ),
  surv_data_interaction
) |>
  arrange(complications, time)

# Combine data for plotting
plot_data_interaction = bind_rows(
  surv_data_interaction,
  bkm_complications
)

# Create plot - use color for model, linetype for complications
p_discrete_time_interaction = ggplot(plot_data_interaction, aes(x = time, y = survival, color = model, linetype = complications)) +
  geom_step(data = filter(plot_data_interaction, model == "Kaplan-Meier"), linewidth = 1.2) +
  geom_line(data = filter(plot_data_interaction, model == "Discrete Time (GLM, Interaction)"), linewidth = 1.2) +
  ylab("Survival probability") +
  xlab("time (days)") +
  ylim(c(0, 1)) +
  scale_color_manual(
    values = c("Kaplan-Meier" = "black", "Discrete Time (GLM, Interaction)" = "steelblue"),
    labels = c("Kaplan-Meier" = "Kaplan-Meier", "Discrete Time (GLM, Interaction)" = "Discrete Time (GLM)")
  ) +
  scale_linetype_manual(
    values = c("no" = "solid", "yes" = "dashed"),
    labels = c("no" = "No Complications", "yes" = "Complications")
  ) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Model"), linetype = guide_legend(title = "Complications"))

ggsave("book/Figures/reductions/discrete-time-complications-glm-interaction.png", p_discrete_time_interaction, 
       height=4, width=8, units="in", dpi=600)


## Piecewise Exponential Model (PEM) example
data("tumor", package = "pammtools")
tumor_comp = tumor |> select(days, status, complications)

# Create 100 equidistant cut points
max_time = max(tumor_comp$days)
cut_points = seq(0, max_time, length.out = 101)  # 101 points create 100 intervals
J = length(cut_points) - 1

# Transform data to piecewise exponential data format using pammtools
ped_data = as_ped(
  data = tumor_comp,
  formula = Surv(days, status) ~ complications,
  cut = cut_points[-1],  # pammtools expects cut points without 0
  id = "id"
) |>
  mutate(interval = factor(interval))  # Convert interval to factor

# Fit Poisson regression model with interaction
# The offset is automatically included by pammtools
pem_fit = glm(
  ped_status ~ interval * complications + offset(offset),
  family = poisson(),
  data = ped_data
)

# Get predicted hazards for each interval and complication group
pred_data_pem = ped_data |>
  distinct(interval, complications, .keep_all = TRUE) |>
  select(interval, complications, tstart, tend, offset)

# Predict expected number of events (which is hazard * time at risk)
pred_data_pem$expected_events = predict(pem_fit, newdata = pred_data_pem, type = "response")
# Convert to hazard rate: hazard = expected_events / time_at_risk
pred_data_pem$hazard = pred_data_pem$expected_events / exp(pred_data_pem$offset)

# Calculate survival probabilities from piecewise constant hazards
# S(t) = exp(-integral of hazard from 0 to t)
surv_data_pem = pred_data_pem |>
  arrange(complications, interval) |>
  group_by(complications) |>
  mutate(
    # Cumulative hazard up to end of each interval
    cumhaz = cumsum(hazard * (tend - tstart)),
    survival = exp(-cumhaz),
    time = tend
  ) |>
  ungroup() |>
  select(complications, time, survival) |>
  mutate(model = "Piecewise Exponential (PEM)")

# Add time 0 with survival = 1
surv_data_pem = bind_rows(
  data.frame(
    complications = rep(unique(tumor_comp$complications), each = 1),
    time = 0,
    survival = 1,
    model = "Piecewise Exponential (PEM)"
  ),
  surv_data_pem
) |>
  arrange(complications, time)

# Get KM estimates for comparison
km_complications = survfit(Surv(days, status) ~ complications, data = tumor_comp)
bkm_complications = broom::tidy(km_complications) |>
  mutate(
    complications = gsub("complications=", "", strata),
    model = "Kaplan-Meier"
  ) |>
  select(complications, time, estimate, model) |>
  rename(survival = estimate)

# Combine data for plotting
plot_data_pem = bind_rows(
  surv_data_pem,
  bkm_complications
)

# Create plot - use color for model, linetype for complications
p_pem = ggplot(plot_data_pem, aes(x = time, y = survival, color = model, linetype = complications)) +
  geom_step(data = filter(plot_data_pem, model == "Kaplan-Meier"), linewidth = 1.2) +
  geom_line(data = filter(plot_data_pem, model == "Piecewise Exponential (PEM)"), linewidth = 1.2) +
  ylab("Survival probability") +
  xlab("time (days)") +
  ylim(c(0, 1)) +
  scale_color_manual(
    values = c("Kaplan-Meier" = "black", "Piecewise Exponential (PEM)" = "steelblue"),
    labels = c("Kaplan-Meier" = "Kaplan-Meier", "Piecewise Exponential (PEM)" = "Piecewise Exponential (PEM)")
  ) +
  scale_linetype_manual(
    values = c("no" = "solid", "yes" = "dashed"),
    labels = c("no" = "No Complications", "yes" = "Complications")
  ) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Model"), linetype = guide_legend(title = "Complications"))

ggsave("book/Figures/reductions/pem-complications-interaction.png", p_pem, 
       height=4, width=8, units="in", dpi=600)


## Comparison: Pooled Logistic Regression vs Nelson-Aalen increments




data("tumor", package = "pammtools")
tumor_comp = tumor |> select(days, status, complications)

# Get unique event times (only where events occurred)
unique_event_times = sort(unique(tumor_comp$days[tumor_comp$status == 1]))

# Create stacked data using intervals at unique event times
# This creates one row per subject per interval up to their event/censoring time
ped_data_na = as_ped(
  data = tumor_comp,
  formula = Surv(days, status) ~ 1,  # No covariates
  cut = unique_event_times,  # Use unique event times as cut points
  id = "id"
) |>
  mutate(interval = factor(interval)) |> 
  group_by(id) |>
  mutate(time = cumsum(exp(offset))) |> 
  filter(time >= tend) |> # <- this needed in order to achieve stacked survival data set
  ungroup()

# Fit pooled logistic regression without covariates
# This estimates the discrete hazard at each event time
glm_na = glm(
  ped_status ~ interval - 1,  # -1 removes intercept, estimates separate hazard for each interval
  family = binomial(),
  data = ped_data_na
)

# Get predicted hazards at each event time
pred_na = ped_data_na |>
  distinct(interval, .keep_all = TRUE) |>
  select(interval, tend)

pred_na$hazard_glm = predict(glm_na, newdata = pred_na, type = "response")

# Calculate Nelson-Aalen increments manually
# h^d(t_k) = d_{t_k} / n_{t_k}
na_increments = data.frame(
  time = unique_event_times,
  hazard_na = NA_real_
)

for (i in seq_along(unique_event_times)) {
  t_k = unique_event_times[i]
  # Number at risk at time t_k
  n_k = sum(tumor_comp$days >= t_k)
  # Number of events at time t_k
  d_k = sum(tumor_comp$days == t_k & tumor_comp$status == 1)
  # Nelson-Aalen increment
  na_increments$hazard_na[i] = d_k / n_k
}

# Match times and compare
comparison = pred_na |>
  mutate(time = tend) |>
  select(time, hazard_glm) |>
  inner_join(na_increments, by = "time") |>
  arrange(time)

# Print comparison
cat("Comparison of Pooled Logistic Regression vs Nelson-Aalen increments:\n")
cat("========================================\n")
print(comparison)

# Calculate differences
comparison = comparison |>
  mutate(difference = hazard_glm - hazard_na,
         relative_diff = (hazard_glm - hazard_na) / hazard_na * 100)

cat("\nSummary statistics:\n")
cat("Mean absolute difference:", mean(abs(comparison$difference)), "\n")
cat("Max absolute difference:", max(abs(comparison$difference)), "\n")
cat("Mean relative difference (%):", mean(abs(comparison$relative_diff)), "\n")

# Create comparison plot
p_comparison = ggplot(comparison, aes(x = time, y = hazard_glm)) +
  geom_point(aes(color = "Pooled Logistic"), size = 2, alpha = 0.7) +
  geom_point(aes(y = hazard_na, color = "Nelson-Aalen"), size = 2, alpha = 0.7) +
  geom_line(aes(y = hazard_glm, color = "Pooled Logistic"), alpha = 0.5) +
  geom_line(aes(y = hazard_na, color = "Nelson-Aalen"), alpha = 0.5) +
  ylab("Discrete hazard h^d(t)") +
  xlab("time (days)") +
  scale_color_manual(
    values = c("Pooled Logistic" = "steelblue", "Nelson-Aalen" = "black"),
    name = "Method"
  ) +
  theme(legend.position = "bottom")

ggsave("book/Figures/reductions/pooled-logistic-vs-na.png", p_comparison, 
       height=4, width=8, units="in", dpi=600)

cat("\nFigure saved to book/Figures/reductions/pooled-logistic-vs-na.png\n")


## PEM interval comparison figure

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

ggsave("book/Figures/reductions/pem-interval-comparison.png", p_pem_comparison,
       height = 6, width = 12, units = "in", dpi = 600)

cat("\nFigure saved to book/Figures/reductions/pem-interval-comparison.png\n")

### AUC
x = c(0, seq.int(0, 1, 0.1))
y_perfect = c(0,  rep(1, 11))   # perfect prediction
y_guess   = c(0, seq.int(0, 1, 0.1))   # random guess (diagonal)
y_good     = x^0.6
y_better    = x^0.3

df = data.frame(
  x = rep(x, 4),
  y = c(y_guess, y_perfect, y_good, y_better),
  group = factor(rep(c("Random guess", "Perfect classifier", "Good model", "Better model"),
                     each = length(x)))
)

g_auc = ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line(lwd = 1) +
  geom_point() +
  labs(
    x = "False positive rate",
    y = "True positive rate",
    color = "Model"
  )

ggsave("book/Figures/evaluation/rocs.png", g_auc,
       height = 5, width = 8, units = "in", dpi = 600)

## Survtsk chapter RMST comparison
yi = c(1,0.8,0.75,0.75,0.7,rep(0.6, 5))
yj = c(1,0.9,0.85,0.6,0.5,rep(0, 5))
df <- data.frame(x = rep(0:9,2), y = c(yi, yj), Patient=rep(c("i", "j"), each = 10))

plot_rmst <- function(df, patient) {
  filtered_df <- df %>%
    filter(Patient == patient, x <= 5)
  rect_df <-  filtered_df %>%
    arrange(x) %>%
    mutate(xmin = x, xmax = lead(x), ymin = 0, ymax = y) %>%
    filter(!is.na(xmax))

  ggplot(df, aes(x = x, y = y, group = Patient)) +
    geom_step(aes(linetype = Patient), lwd = 1) +
    geom_rect(
      data = rect_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      alpha = 0.3,
      inherit.aes = FALSE) +
    annotate("text", x = 2, y = 0.4,
      label = as.expression(bquote(RMST[.(patient)] * "(" * 5 * ")" == .(round(sum(filtered_df$y), 2)))),
      parse = TRUE) +
    labs(x = "Time", y = "Survival probability", title = paste0("RMST(5) for patient ", patient))    
}

p1 <- plot_rmst(df, "i")
p2 <- plot_rmst(df, "j")
p_rmst_survtsk <- (p1 + p2) + plot_layout(guides = "collect")

ggsave("book/Figures/survtsk/rmst.png", p_rmst_survtsk,
       height = 4, width = 9, units = "in", dpi = 600)

## C-index interval censoring
cases <- tibble::tribble(
  ~case, ~obs, ~l, ~r,
  "1. li < ri < lj < rj", "i", 1, 2,
  "1. li < ri < lj < rj", "j", 3, 4,

  "2. lj < rj < li < ri", "i", 3, 4,
  "2. lj < rj < li < ri", "j", 1, 2,

  "3. lj < li < ri < rj", "i", 2, 3,
  "3. lj < li < ri < rj", "j", 1, 4,

  "4. li < lj < rj < ri", "i", 1, 4,
  "4. li < lj < rj < ri", "j", 2, 3,

  "5. lj < li < rj < ri", "i", 2, 4,
  "5. lj < li < rj < ri", "j", 1, 3,

  "6. li < lj < ri < rj", "i", 1, 3,
  "6. li < lj < ri < rj", "j", 2, 4
) %>%
  mutate(y = ifelse(obs == "i", 2, 1))

g_intervals = ggplot(cases) +
  geom_segment(aes(x = l, xend = r, y = y, yend = y, lty = obs),
               linewidth = 1) +
  facet_wrap(~case, ncol = 2) +
  scale_y_continuous(breaks = c(1,2), labels = c("j","i")) +
  labs(x = "Time", y = "") +
  theme(text = element_text(size = 18), legend.position = "n", axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 16))

ggsave("book/Figures/evaluation/intervals.png", g_intervals,
      height = 6, width = 8, units = "in", dpi = 600)

## =============================================================================
## EHA Reductions: cause-specific hazard models on sir.adm
## Generates figures for P4C23 (Reductions for Event-History Analysis):
##   - book/Figures/reductions/cif-marg-sir.png  (marginal: AJ vs PAM/Cox/RSF)
##   - book/Figures/reductions/cif-adj-sir.png   (adjusted: PAM/Cox/RSF for 4 profiles)
## =============================================================================

suppressPackageStartupMessages({
  library(survival); library(mvna); library(cmprsk); library(pammtools)
  library(mgcv); library(randomForestSRC); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(purrr)
})
theme_set(theme_bw())
set.seed(20260505)

## ---- Data prep ----
data(sir.adm, package = "mvna")
sir <- sir.adm |>
  mutate(
    pneu = factor(pneu, labels = c("No pneumonia", "Pneumonia")),
    sex  = factor(sex, levels = c("F","M"))
  )

tau <- max(sir$time)
common_grid <- sort(unique(c(0, sir$time[sir$status > 0], tau)))

## ---- Cause-specific datasets (status==1: discharge, status==2: death) ----
sir_disch <- sir |> mutate(stat_e = as.integer(status == 1))
sir_death <- sir |> mutate(stat_e = as.integer(status == 2))

## ---- (1) Aalen-Johansen stratified by pneumonia (non-parametric baseline) ----
sir_aj <- sir |>
  mutate(to = case_when(status == 0 ~ "cens", .default = as.character(status)))
aj_fit <- cuminc(sir_aj$time, sir_aj$to, group = sir_aj$pneu, cencode = "cens")

aj_df <- imap_dfr(aj_fit[1:4], function(.x, .y) {
  cause_chr <- sub(".* ", "", .y); pneu_chr <- sub(" [0-9]+$", "", .y)
  data.frame(time = .x$time, cif = .x$est, pneu = pneu_chr, cause = cause_chr)
}) |>
  mutate(
    cause  = factor(cause, levels = c("1","2"), labels = c("Discharge","Death")),
    pneu   = factor(pneu,  levels = c("No pneumonia","Pneumonia")),
    method = "Aalen-Johansen")

## ---- (2) Cause-specific PAM (PED format + mgcv::gam, Poisson w/ log offset) ----
ped_disch <- as_ped(sir_disch, Surv(time, stat_e) ~ pneu + age + sex,
                    cut = common_grid, id = "id")
ped_death <- as_ped(sir_death, Surv(time, stat_e) ~ pneu + age + sex,
                    cut = common_grid, id = "id")

pam_disch_marg <- gam(ped_status ~ s(tend) + pneu, family = poisson(),
                      offset = offset, data = ped_disch)
pam_death_marg <- gam(ped_status ~ s(tend) + pneu, family = poisson(),
                      offset = offset, data = ped_death)
pam_disch_adj  <- gam(ped_status ~ s(tend) + pneu + s(age) + sex, family = poisson(),
                      offset = offset, data = ped_disch)
pam_death_adj  <- gam(ped_status ~ s(tend) + pneu + s(age) + sex, family = poisson(),
                      offset = offset, data = ped_death)

# Helper: build a per-profile prediction grid on the common_grid
make_pred_grid <- function(prof_df, ped) {
  tends <- sort(unique(ped$tend))
  out   <- prof_df[rep(1, length(tends)), , drop = FALSE]
  out$tend   <- tends
  out$intlen <- diff(c(0, tends))
  out$offset <- log(out$intlen)
  out
}
# CIF via cause-specific reduction: F_e(t|x) = sum S(t_{k-1}|x) * h_e(t_k|x) * intlen
pred_pam_cif <- function(prof_df, pam_d, pam_e, ped_d) {
  g <- make_pred_grid(prof_df, ped_d)
  g$h_d <- predict(pam_d, newdata = g, type = "response")
  g$h_e <- predict(pam_e, newdata = g, type = "response")
  H_all <- cumsum((g$h_d + g$h_e) * g$intlen)
  S_prev <- c(1, head(exp(-H_all), -1))
  cif_d <- cumsum(S_prev * g$h_d * g$intlen)
  cif_e <- cumsum(S_prev * g$h_e * g$intlen)
  bind_rows(
    data.frame(time = g$tend, cif = cif_d, cause = "Discharge"),
    data.frame(time = g$tend, cif = cif_e, cause = "Death")
  ) |>
    bind_cols(g[rep(seq_len(nrow(g)), 2),
                setdiff(colnames(g), c("tend","intlen","offset","h_d","h_e")),
                drop = FALSE])
}

## ---- (3) Cause-specific Cox PH ----
cox_disch_marg <- coxph(Surv(time, stat_e) ~ pneu, data = sir_disch, x = TRUE)
cox_death_marg <- coxph(Surv(time, stat_e) ~ pneu, data = sir_death, x = TRUE)
cox_disch_adj  <- coxph(Surv(time, stat_e) ~ pneu + age + sex, data = sir_disch, x = TRUE)
cox_death_adj  <- coxph(Surv(time, stat_e) ~ pneu + age + sex, data = sir_death, x = TRUE)

# CIF for cause-specific Cox: combine Breslow baseline with PH risk score
cox_cif <- function(prof_df, cox_d, cox_e, grid) {
  bh_d <- basehaz(cox_d, centered = FALSE)
  bh_e <- basehaz(cox_e, centered = FALSE)
  eta_d <- predict(cox_d, newdata = prof_df, type = "lp", reference = "zero")
  eta_e <- predict(cox_e, newdata = prof_df, type = "lp", reference = "zero")
  step_eval <- function(bh, t_eval)
    approxfun(c(0, bh$time), c(0, bh$hazard), method = "constant", rule = 2, f = 0)(t_eval)
  H_d <- step_eval(bh_d, grid) * exp(eta_d)
  H_e <- step_eval(bh_e, grid) * exp(eta_e)
  intlen <- diff(c(0, grid))
  h_d <- ifelse(intlen > 0, diff(c(0, H_d)) / intlen, 0)
  h_e <- ifelse(intlen > 0, diff(c(0, H_e)) / intlen, 0)
  H_all <- cumsum((h_d + h_e) * intlen)
  S_prev <- c(1, head(exp(-H_all), -1))
  cif_d <- cumsum(S_prev * h_d * intlen)
  cif_e <- cumsum(S_prev * h_e * intlen)
  bind_rows(
    data.frame(time = grid, cif = cif_d, cause = "Discharge"),
    data.frame(time = grid, cif = cif_e, cause = "Death")
  ) |>
    bind_cols(prof_df[rep(1, length(grid) * 2), , drop = FALSE])
}

## ---- (4) Random Survival Forest (native competing risks) ----
rsf_marg <- rfsrc(Surv(time, status) ~ pneu, data = sir, ntree = 500,
                  cause = c(1, 2), seed = 42)
rsf_adj  <- rfsrc(Surv(time, status) ~ pneu + age + sex, data = sir, ntree = 500,
                  cause = c(1, 2), seed = 42)

rsf_cif <- function(prof_df, rsf_fit) {
  pred <- predict(rsf_fit, newdata = prof_df)
  bind_rows(lapply(seq_len(nrow(prof_df)), function(i) {
    bind_rows(
      data.frame(time = pred$time.interest, cif = pred$cif[i,,1], cause = "Discharge"),
      data.frame(time = pred$time.interest, cif = pred$cif[i,,2], cause = "Death")
    ) |>
      bind_cols(prof_df[rep(i, length(pred$time.interest) * 2), , drop = FALSE])
  }))
}

## ---- Marginal verification figure: PAM/Cox/RSF recover AJ ----
profs_marg <- data.frame(
  pneu = factor(c("No pneumonia","Pneumonia"),
                levels = c("No pneumonia","Pneumonia"))
)
profs_marg_pam <- profs_marg
profs_marg_pam$age <- median(sir$age)
profs_marg_pam$sex <- factor("M", levels = c("F","M"))

pam_marg_df <- bind_rows(lapply(seq_len(nrow(profs_marg_pam)), function(i)
  pred_pam_cif(profs_marg_pam[i,], pam_disch_marg, pam_death_marg, ped_disch))) |>
  select(time, cif, cause, pneu) |>
  mutate(method = "PAM", cause = factor(cause, levels = c("Discharge","Death")))

cox_marg_df <- bind_rows(lapply(seq_len(nrow(profs_marg)), function(i)
  cox_cif(profs_marg[i,, drop = FALSE], cox_disch_marg, cox_death_marg, common_grid))) |>
  mutate(method = "Cox PH", cause = factor(cause, levels = c("Discharge","Death")))

rsf_marg_df <- rsf_cif(profs_marg, rsf_marg) |>
  mutate(method = "RSF", cause = factor(cause, levels = c("Discharge","Death")))

marg_df <- bind_rows(aj_df, pam_marg_df, cox_marg_df, rsf_marg_df) |>
  mutate(method = factor(method,
                         levels = c("Aalen-Johansen","PAM","Cox PH","RSF")))

p_marg <- ggplot(marg_df, aes(time, cif, color = method, linetype = method)) +
  geom_step(linewidth = 0.7) +
  facet_grid(cause ~ pneu) +
  scale_color_manual(values = c("Aalen-Johansen"="black","PAM"="#377eb8",
                                "Cox PH"="#e41a1c","RSF"="#4daf4a")) +
  scale_linetype_manual(values = c("Aalen-Johansen"="solid","PAM"="dashed",
                                   "Cox PH"="dotdash","RSF"="dotted")) +
  labs(x = "Time (days)", y = "Cumulative incidence",
       color = NULL, linetype = NULL) +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(0, 120), ylim = c(0, 1))

dir.create("book/Figures/reductions", showWarnings = FALSE, recursive = TRUE)
ggsave("book/Figures/reductions/cif-marg-sir.png", p_marg,
       height = 5.5, width = 8.5, units = "in", dpi = 600)

## ---- Covariate-adjusted figure: PAM/Cox/RSF on a 4-profile grid ----
profs_adj <- expand.grid(
  pneu = factor(c("No pneumonia","Pneumonia"),
                levels = c("No pneumonia","Pneumonia")),
  age  = c(40, 75),
  sex  = factor("M", levels = c("F","M"))
)
profs_adj$profile <- with(profs_adj, paste0(pneu, ", age ", age))
profs_adj$profile <- factor(profs_adj$profile, levels = unique(profs_adj$profile))

pam_adj_df <- bind_rows(lapply(seq_len(nrow(profs_adj)), function(i)
  pred_pam_cif(profs_adj[i,], pam_disch_adj, pam_death_adj, ped_disch))) |>
  select(time, cif, cause, pneu, age, sex, profile) |>
  mutate(method = "PAM", cause = factor(cause, levels = c("Discharge","Death")))

cox_adj_df <- bind_rows(lapply(seq_len(nrow(profs_adj)), function(i)
  cox_cif(profs_adj[i,, drop = FALSE], cox_disch_adj, cox_death_adj, common_grid))) |>
  mutate(method = "Cox PH", cause = factor(cause, levels = c("Discharge","Death")))

rsf_adj_df <- rsf_cif(profs_adj, rsf_adj) |>
  mutate(method = "RSF", cause = factor(cause, levels = c("Discharge","Death")))

adj_df <- bind_rows(pam_adj_df, cox_adj_df, rsf_adj_df) |>
  mutate(method = factor(method, levels = c("PAM","Cox PH","RSF")))

p_adj <- ggplot(adj_df, aes(time, cif, color = profile, linetype = method)) +
  geom_step(linewidth = 0.6) +
  facet_wrap(~cause) +
  labs(x = "Time (days)", y = "Cumulative incidence",
       color = "Profile", linetype = "Method") +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 1)) +
  coord_cartesian(xlim = c(0, 120), ylim = c(0, 1))

ggsave("book/Figures/reductions/cif-adj-sir.png", p_adj,
       height = 5.5, width = 9, units = "in", dpi = 600)

## =========================================================================
## Censoring + truncation schematic figures (Ch. 3)
##  Outputs (both .svg and .png):
##   - book/Figures/survival/censoring.{svg,png}        (right-censoring)
##   - book/Figures/survival/left-truncation.{svg,png}  (left-truncation)
##   - book/Figures/survival/right-truncation.{svg,png} (right-truncation)
## =========================================================================

shape_vals   <- c("observed event" = 21, "unobserved event" = 23, "censoring" = 21)
fill_vals    <- c("observed event" = "black", "unobserved event" = "black", "censoring" = "white")
shape_breaks <- c("observed event", "unobserved event", "censoring")
line_vals    <- c("observed" = "solid", "unobserved" = "dotted")
line_breaks  <- c("observed", "unobserved")

shape_legend <- guide_legend(override.aes = list(
  fill   = c("black", "black", "white"),
  shape  = c(21, 23, 21),
  colour = "black"))

ensure_levels <- function(df, level_col, all_levels) {
  missing <- setdiff(all_levels, as.character(df[[level_col]]))
  if (length(missing) == 0) return(df)
  ph <- df[rep(1L, length(missing)), , drop = FALSE]
  ph[[level_col]] <- factor(missing, levels = all_levels)
  for (nm in setdiff(names(ph), level_col)) ph[[nm]] <- NA
  rbind(df, ph)
}

save_pair <- function(p, base, h = 3.5, w = 8) {
  ggsave(paste0(base, ".svg"), p, height = h, width = w, units = "in")
  ggsave(paste0(base, ".png"), p, height = h, width = w, units = "in", dpi = 300)
}

dy <- 0.10  # half-offset between stacked bars

bigger_text <- theme(
  axis.title  = element_text(size = 17),
  axis.text   = element_text(size = 15),
  legend.text = element_text(size = 15)
)

## =========================================================================
## RIGHT-CENSORING (4 subjects)
## =========================================================================
study_end <- 8
cens <- tibble::tribble(
  ~subject, ~obs_end, ~event_end, ~status,
  1,        7,        7.0,        "event_observed",
  2,        4,        6.0,        "event_unobserved_during",
  3,        1,        9.0,        "event_unobserved_after",
  4,        8,        9.5,        "event_unobserved_after"
)
seg_solid  <- cens |> transmute(subject, x = 0, xend = obs_end, line = "observed")
seg_dashed <- cens |> filter(status != "event_observed") |>
  transmute(subject, x = obs_end, xend = event_end, line = "unobserved")
segs <- bind_rows(seg_solid, seg_dashed) |>
  mutate(line = factor(line, levels = line_breaks))
pts <- bind_rows(
  cens |> filter(status != "event_observed") |> transmute(subject, x = obs_end,   type = "censoring"),
  cens |> filter(status == "event_observed")  |> transmute(subject, x = event_end, type = "observed event"),
  cens |> filter(status != "event_observed") |> transmute(subject, x = event_end, type = "unobserved event")
) |> mutate(type = factor(type, levels = shape_breaks))
pts <- ensure_levels(pts, "type", shape_breaks)

p_cens <- ggplot() +
  geom_vline(xintercept = study_end, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_segment(data = segs,
               aes(x = x, xend = xend, y = subject, yend = subject, linetype = line),
               linewidth = 1.1) +
  geom_point(data = pts,
             aes(x = x, y = subject, shape = type, fill = type),
             colour = "black", size = 3.5, stroke = 0.7) +
  scale_linetype_manual(values = line_vals, breaks = line_breaks, drop = FALSE, name = NULL) +
  scale_shape_manual(values = shape_vals, breaks = shape_breaks, drop = FALSE, name = NULL) +
  scale_fill_manual(values = fill_vals, breaks = shape_breaks, drop = FALSE, name = NULL) +
  scale_y_continuous(breaks = 1:4) +
  scale_x_continuous(breaks = c(0, study_end), labels = c("study start", "study end"),
                     limits = c(0, 10), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "Subject") +
  guides(shape = shape_legend, linetype = guide_legend(order = 1)) +
  theme_bw(base_size = 15) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  bigger_text
save_pair(p_cens, "book/Figures/survival/censoring", h = 3.5, w = 8)

## =========================================================================
## LEFT-TRUNCATION (3 infants on age 0-365 days)
##  steelblue = truncation time; black = time-to-event
##  solid = subject in data; dashed = invisible (t_i < t_i^L)
##  Filled steelblue circle/diamond at end of truncation bar marks
##  whether the subject is in/out of sample.
## =========================================================================
lt <- tibble::tribble(
  ~subject, ~t_L,  ~t_event, ~status,            ~visible,
  1,        100,   50,       "unobserved event", FALSE,
  2,        80,    200,      "observed event",   TRUE,
  3,        150,   365,      "censoring",        TRUE
)
period_breaks_lt <- c("truncation time", "time-to-event")
period_vals_lt   <- c("truncation time" = "steelblue", "time-to-event" = "black")

lt_segs <- bind_rows(
  lt |> transmute(subject, y = subject + dy, x = 0, xend = t_L,
                  period = "truncation time",
                  line   = if_else(visible, "observed", "unobserved")),
  lt |> transmute(subject, y = subject - dy, x = 0, xend = t_event,
                  period = "time-to-event",
                  line   = if_else(visible, "observed", "unobserved"))
) |> mutate(period = factor(period, levels = period_breaks_lt),
            line   = factor(line, levels = line_breaks))
lt_pts <- lt |> transmute(subject, y = subject - dy, x = t_event,
                          type = factor(status, levels = shape_breaks))
lt_pts <- ensure_levels(lt_pts, "type", shape_breaks)

lt_labels <- bind_rows(
  lt |> transmute(x = t_L     - 7, y = subject + dy + 0.25,
                  label = sprintf("t[%d]^'ℓ'", subject)),
  lt |> transmute(x = t_event - 7, y = subject - dy - 0.25,
                  label = sprintf("t[%d]",          subject))
)

p_lt <- ggplot() +
  geom_vline(xintercept = 365, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_segment(data = lt_segs,
               aes(x = x, xend = xend, y = y, yend = y,
                   colour = period, linetype = line),
               linewidth = 1.5) +
  geom_point(data = filter(lt, visible),
             aes(x = t_L, y = subject + dy),
             shape = 21, fill = "steelblue", colour = "steelblue",
             size = 3.5, stroke = 0.7) +
  geom_point(data = filter(lt, !visible),
             aes(x = t_L, y = subject + dy),
             shape = 23, fill = "steelblue", colour = "steelblue",
             size = 3.5, stroke = 0.7) +
  geom_point(data = lt_pts,
             aes(x = x, y = y, shape = type, fill = type),
             colour = "black", size = 3.5, stroke = 0.7) +
  geom_text(data = lt_labels, aes(x = x, y = y, label = label),
            parse = TRUE, size = 5, hjust = 1, vjust = 0.5) +
  scale_colour_manual(values = period_vals_lt, breaks = period_breaks_lt, name = NULL) +
  scale_linetype_manual(values = line_vals, breaks = line_breaks, drop = FALSE, name = NULL) +
  scale_shape_manual(values = shape_vals, breaks = shape_breaks, drop = FALSE, name = NULL) +
  scale_fill_manual(values = fill_vals, breaks = shape_breaks, drop = FALSE, name = NULL) +
  scale_y_continuous(breaks = 1:3, limits = c(0.55, 3.45)) +
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365),
                     labels = c("0", "90", "180", "270", "365\n(study end)"),
                     limits = c(0, 400), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Age (days)", y = "Subject") +
  guides(colour   = guide_legend(order = 1, override.aes = list(linewidth = 1.2)),
         linetype = guide_legend(order = 2),
         shape    = shape_legend) +
  theme_bw(base_size = 15) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank()) +
  bigger_text
save_pair(p_lt, "book/Figures/survival/left-truncation", h = 3.5, w = 8.5)

## =========================================================================
## RIGHT-TRUNCATION (3 subjects, calendar time 0-30 days)
##  steelblue = truncation time (recruit -> db query)
##  black     = time-to-event (recruit -> event)
##  solid = in registry (event before query); dashed = absent from registry
##  Vertical dashed line marks the database-query date.
##  Filled steelblue circle/diamond at db query marks in/out of registry.
## =========================================================================
db_query <- 20
rt <- tibble::tribble(
  ~subject, ~recruit, ~t_event, ~status,            ~visible,
  1,        0,        8,        "observed event",   TRUE,
  2,        10,       28,       "unobserved event", FALSE,
  3,        5,        10,       "observed event",   TRUE
)
period_breaks_rt <- c("truncation time", "time-to-event")
period_vals_rt   <- c("truncation time" = "steelblue", "time-to-event" = "black")

rt_segs <- bind_rows(
  rt |> transmute(subject, y = subject + dy, x = recruit, xend = db_query,
                  period = "truncation time",
                  line   = if_else(visible, "observed", "unobserved")),
  rt |> transmute(subject, y = subject - dy, x = recruit, xend = t_event,
                  period = "time-to-event",
                  line   = if_else(visible, "observed", "unobserved"))
) |> mutate(period = factor(period, levels = period_breaks_rt),
            line   = factor(line, levels = line_breaks))
rt_pts <- rt |> transmute(subject, y = subject - dy, x = t_event,
                          type = factor(status, levels = shape_breaks))
rt_pts <- ensure_levels(rt_pts, "type", shape_breaks)

rt_labels <- bind_rows(
  rt |> transmute(x = db_query - 0.6, y = subject + dy + 0.25, label = sprintf("t[%d]^r", subject)),
  rt |> transmute(x = t_event  - 0.6, y = subject - dy - 0.25, label = sprintf("t[%d]",   subject))
)

p_rt <- ggplot() +
  geom_vline(xintercept = db_query, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_segment(data = rt_segs,
               aes(x = x, xend = xend, y = y, yend = y,
                   colour = period, linetype = line),
               linewidth = 1.5) +
  geom_point(data = filter(rt, visible),
             aes(x = db_query, y = subject + dy),
             shape = 21, fill = "steelblue", colour = "steelblue",
             size = 3.5, stroke = 0.7) +
  geom_point(data = filter(rt, !visible),
             aes(x = db_query, y = subject + dy),
             shape = 23, fill = "steelblue", colour = "steelblue",
             size = 3.5, stroke = 0.7) +
  geom_point(data = rt_pts,
             aes(x = x, y = y, shape = type, fill = type),
             colour = "black", size = 3.5, stroke = 0.7) +
  geom_text(data = rt_labels, aes(x = x, y = y, label = label),
            parse = TRUE, size = 5, hjust = 1, vjust = 0.5) +
  scale_colour_manual(values = period_vals_rt, breaks = period_breaks_rt, name = NULL) +
  scale_linetype_manual(values = line_vals, breaks = line_breaks, drop = FALSE, name = NULL) +
  scale_shape_manual(values = shape_vals, breaks = shape_breaks, drop = FALSE, name = NULL) +
  scale_fill_manual(values = fill_vals, breaks = shape_breaks, drop = FALSE, name = NULL) +
  scale_y_continuous(breaks = 1:3, limits = c(0.55, 3.45)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30),
                     labels = c("0\n(start data\ncollection)", "10", "20\n(database queried)", "30"),
                     limits = c(0, 32), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Calendar time (days)", y = "Subject") +
  guides(colour   = guide_legend(order = 1, override.aes = list(linewidth = 1.2)),
         linetype = guide_legend(order = 2),
         shape    = shape_legend) +
  theme_bw(base_size = 15) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank()) +
  bigger_text
save_pair(p_rt, "book/Figures/survival/right-truncation", h = 3.5, w = 8.5)

cat("Saved.\n")

## =========================================================================
## Prediction-types overview figure (P1C6 survtsk, @fig-survtsk-overview)
##  Output: book/Figures/survtsk/predict_types.{svg,png}
##  Data:   six patients from the Ch.3 tumor table (P1C4)
##  Model:  Weibull AFT fitted on the full `tumor` data via flexsurvreg,
##          with shape ~ complications and scale ~ age + sex + complications.
## =========================================================================

suppressPackageStartupMessages({
  library(flexsurv); library(survival); library(pammtools)
  library(patchwork); library(svglite)
})

pastel_pt <- c("#A1C9F4", "#FFB482", "#8DE5A1", "#FF9F9B", "#D0BBFF", "#FAB0E4")

data(tumor, package = "pammtools")
fit_pt <- flexsurvreg(
  Surv(days, status) ~ age + sex + complications,
  anc  = list(shape = ~ complications),
  data = tumor,
  dist = "weibull"
)

patients_pt <- tibble::tibble(
  id            = 1:6,
  age           = c(71, 70, 67, 58, 39, 59),
  sex           = factor(c("female","male","female","male","female","female"),
                          levels = c("male","female")),
  complications = factor(c("no","no","yes","no","yes","no"),
                          levels = c("no","yes")),
  days          = c(1217, 519, 2414, 397, 1217, 268),
  status        = c(1, 0, 0, 1, 0, 1)
)

mean_pred_pt    <- summary(fit_pt, newdata = patients_pt, type = "mean", tidy = TRUE)
patients_pt$E_t  <- mean_pred_pt$est
patients_pt$lcl  <- mean_pred_pt$lcl
patients_pt$ucl  <- mean_pred_pt$ucl
patients_pt$risk <- - log(patients_pt$E_t)
patients_pt$risk <- patients_pt$risk - mean(patients_pt$risk)
patients_pt$id_f <- factor(patients_pt$id, levels = patients_pt$id)

tgrid_pt     <- seq(0, 3000, length.out = 250)
surv_pred_pt <- summary(fit_pt, newdata = patients_pt, type = "survival",
                         t = tgrid_pt, tidy = TRUE)
surv_pred_pt$id <- factor(rep(patients_pt$id, each = length(tgrid_pt)),
                           levels = patients_pt$id)

## Panel A: tabular data drawn as a ggplot for visual consistency
col_names_pt  <- c("id","age","sex","complications","days","status")
col_widths_pt <- c(0.5, 0.5, 0.85, 1.4, 0.7, 0.85)
col_starts_pt <- c(0, cumsum(col_widths_pt)[-length(col_widths_pt)])
col_ends_pt   <- cumsum(col_widths_pt)

build_cell_pt <- function(row_y, col_idx, value, fill, fontface, border) {
  tibble::tibble(row_y, col_idx,
         xmin = col_starts_pt[col_idx], xmax = col_ends_pt[col_idx],
         ymin = -row_y - 0.5, ymax = -row_y + 0.5,
         xmid = (col_starts_pt[col_idx] + col_ends_pt[col_idx]) / 2,
         ymid = -row_y,
         value, fill, fontface, border)
}

header_cells_pt <- do.call(rbind, lapply(seq_along(col_names_pt), function(j)
  build_cell_pt(0, j, col_names_pt[j], "grey92", "bold", "grey75")))

tbl_long_pt <- patients_pt |>
  dplyr::transmute(id, age, sex = as.character(sex),
                   complications = as.character(complications), days, status)

data_cells_pt <- do.call(rbind, lapply(seq_len(nrow(tbl_long_pt)), function(i)
  do.call(rbind, lapply(seq_along(col_names_pt), function(j) {
    val  <- as.character(tbl_long_pt[[col_names_pt[j]]][i])
    fill <- if (j == 1) pastel_pt[i] else "white"
    build_cell_pt(i, j, val, fill, "plain", "grey85")
  }))))

cells_pt <- dplyr::bind_rows(header_cells_pt, data_cells_pt)

p_table_pt <- ggplot(cells_pt) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = fill, colour = border), linewidth = 0.4) +
  geom_text(aes(x = xmid, y = ymid, label = value, fontface = fontface),
            size = 4.2) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_x_continuous(limits = c(min(col_starts_pt), max(col_ends_pt)),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(-nrow(tbl_long_pt) - 0.6, 0.6),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = NULL, title = "Tabular data") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        axis.title       = element_blank(),
        panel.grid       = element_blank())

## Panel B: predicted survival times with 95% CIs
p_time_pt <- ggplot(patients_pt, aes(x = E_t, y = id_f, colour = id_f)) +
  geom_linerange(aes(xmin = lcl, xmax = ucl), linewidth = 1.2) +
  geom_point(size = 4) +
  scale_colour_manual(values = pastel_pt) +
  scale_x_continuous(limits = c(0, max(patients_pt$ucl) * 1.05),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Predicted survival time (days)", y = "Subject",
       title = "Predicted survival times") +
  guides(colour = "none") +
  theme_bw(base_size = 13) +
  theme(plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())

## Panel C: predicted relative-risk scores (lollipop)
risk_lim_pt <- max(abs(patients_pt$risk)) * 1.30
p_risk_pt <- ggplot(patients_pt, aes(x = id_f, y = risk, colour = id_f)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_segment(aes(xend = id_f, y = 0, yend = risk), linewidth = 1.4) +
  geom_point(size = 4.5) +
  scale_colour_manual(values = pastel_pt) +
  scale_y_continuous(limits = c(-risk_lim_pt, risk_lim_pt),
                     expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = "Subject", y = "Relative risk (centred)",
       title = "Predicted relative risk scores") +
  guides(colour = "none") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
        panel.grid.minor = element_blank())

## Panel D: predicted survival distributions
p_surv_pt <- ggplot(surv_pred_pt, aes(x = time, y = est,
                                       colour = id, group = id)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = pastel_pt) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Time (days)", y = "Survival probability",
       title = "Predicted survival distributions") +
  guides(colour = "none") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
        panel.grid.minor = element_blank())

final_pt <- (p_table_pt | p_time_pt) / (p_risk_pt | p_surv_pt) +
  plot_layout(widths = c(1, 1), heights = c(1, 1))

ggsave("book/Figures/survtsk/predict_types.svg", final_pt,
       width = 9.5, height = 6.5, units = "in")
ggsave("book/Figures/survtsk/predict_types.png", final_pt,
       width = 9.5, height = 6.5, units = "in", dpi = 300)


##------------------
## Logloss and Brier
##------------------
p <- seq(0.001, 0.999, length.out = 500)

df <- data.frame(
  Prediction = rep(p, 4),
  Value = c(
    (1 - p)^2, p^2,
    -log(p), -log(1 - p)
  ),
  Class = rep(c("y1", "y0", "y1", "y0"), each = length(p)),
  Loss = rep(c("Brier score", "Log loss"), each = 2 * length(p))
)


brier <- ggplot(subset(df, Loss == "Brier score"),
                aes(Prediction, Value, colour = Class)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  annotate(
    "label", x = 0.10, y = 1.00,
    label = "y[i] == 1",
    parse = TRUE
  ) +
  annotate(
    "label", x = 0.90, y = 1.00,
    label = "y[i] == 0",
    parse = TRUE,
  ) +
  labs(
    x = expression(Prediction ~ hat(y)[i]),
    y = "Brier score"
  ) +
  theme(legend.position = "none")

logloss <- ggplot(subset(df, Loss == "Log loss"),
                  aes(Prediction, Value, colour = Class)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    limits = c(0, 3.15),
    breaks = c(0, log(2), 1, 2, 3),
    labels = c("0", "0.69", "1", "2", "3")
  ) +
  annotate(
    "label", x = 0.15, y = 3.00,
    label = "y[i] == 1",
    parse = TRUE
  ) +
  annotate(
    "label", x = 0.85, y = 3.00,
    label = "y[i] == 0",
    parse = TRUE,
  ) +
  labs(
    x = expression(Prediction ~ hat(y)[i]),
    y = "Log loss"
  ) +
  theme(legend.position = "none")



ggsave("book/Figures/evaluation/brier_logloss.png", brier + logloss,
      width = 5.5, height = 3, units = "in")
