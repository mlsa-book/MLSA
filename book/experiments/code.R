remotes::install_github("mlr-org/mlr3proba", ref = 'v0.5.7', upgrade = "never")
remotes::install_github("mlr-org/mlr3", ref = 'v0.16.1', upgrade = "never")
remotes::install_github("mlr-org/paradox", ref = 'v0.11.1', upgrade = "never")

library(ggplot2)
theme_set(theme_bw())

library(distr6)
library(ggplot2)
g = dstr("Gompertz", shape = 2, decorators = "ExoticStatistics")
t = seq.int(0, 1.5, length.out = 100)
d = data.frame(t = t, fun = factor(rep(c("Density", "Hazard", "Cumulative Density", "Survival"), each = 100), levels = c("Density", "Hazard", "Cumulative Density", "Survival")), y = c(g$pdf(t), g$hazard(t), g$cdf(t), g$survival(t)))
g = ggplot(d, aes(x = t, y = y, color = fun)) +
  geom_line() +
  facet_wrap(~fun, scales = "free", nrow = 2) +
  theme_bw() +
  theme(legend.position = "n")
ggsave("book/Figures/introduction/gompertz.png", g, height = 3, units = "in",
  dpi = 600)

## Ranking
rm(list = ls())
library(dplyr)
library(mlr3)
library(mlr3proba)

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
rm(list = ls())
library(mlr3proba)
library(mlr3pipelines)
library(mlr3extralearners)
library(ggplot2)
library(dplyr)
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
rm(list = ls())
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
rm(list = ls())
set.seed(20241125)
library(survival)
library(ggplot2)
library(patchwork)
library(survminer)
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
rm(list = ls())
set.seed(20241109)
library(survival)
library(party)
fit <- party::ctree(Surv(time, status) ~ ., lung)
png("book/Figures/forests/lung.png", height = 400, width = 600)
plot(fit)
dev.off()

## bootstrapped rsfs
rm(list = ls())
library(ggplot2)
library(patchwork)

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
library(survival)
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
  xlab("time")
p_km
ggsave("book/Figures/survival/km-tumor.png", p_km, height=3, units="in", dpi=600)

# stratified KM wrt complications
tumor = tumor |>
  mutate(age_bin = factor(age < 50, levels = c(TRUE, FALSE), labels = c("age < 50", "age >= 50")))
km_age_bin = survfit(Surv(days, status)~age_bin, data = tumor)
bkm_age_bin = broom::tidy(km_age_bin)
med_km_age_bin = as.numeric(median(km_age_bin))
df_age_bin = data.frame(
  x = c(0, med_km_age_bin[2]), # Starting x-coordinates
  y = c(0.5, 0),        # Starting y-coordinates
  xend = c(med_km_age_bin[2], med_km_age_bin[2]),    # Ending x-coordinates
  yend = c(0.5, 0.5)       # Ending y-coordinates
)

p_km_age_bin = ggplot(bkm_age_bin, aes(x = time, y = estimate)) +
  geom_step(aes(col = strata)) +
  geom_segment(data=df_age_bin, aes(x=x, xend=xend, y=y, yend=yend), lty = 3)+
  geom_hline(yintercept =  .5, lty = 3) +
  ylim(c(0, 1)) +
  ylab("S(t)") +
  xlab("time")
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
library(patchwork)

p_km_infants_joined = p_km_infants + p_km_infants_lt + plot_layout(guides =  "collect")
ggsave("book/Figures/survival/km-infants.png", p_km_infants_joined, height=3, width=7, units="in", dpi=600)


## table infant data
inf_sub = infants |>
  filter(stratum %in% c(1, 2, 4)) |>
  select(stratum, enter, exit, event, mother)

inf_sub |> knitr::kable()

# Chapter 13 - Traditional models
set.seed(2029)
library(survival)
library(mlr3proba)
t <- tsk("rats")$filter(sample(tsk("rats")$nrow, 5))
t$kaplan()$surv
knitr::kable(t$data()[, c(3,4,5,1,2)],align = "l")

# PH vs AFT
set.seed(290125)
library(survival)
library(distr6)
library(ggplot2)
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
library(extraDistr)
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
library(survival)
library(patchwork)
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
library(mvna)
library(etm)
library(cmprsk)
set.seed(241206)
### table
data(sir.adm, package = "mvna")
tab_sir = sir.adm |> group_by(status) |> sample_n(2) |>
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

library(mstate) #prothr dataset
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


## Survival task viz
library(dplyr)
library(patchwork)
set.seed(20250822)

sex <- rbinom(10, 1, 0.5)
x_base <- ceiling(runif(10, ifelse(sex == 1, 20, 1), ifelse(sex == 1, 50, 30)))
x <- unlist(lapply(x_base, function(val) c(0, val, 50)))

df <- data.frame(x = x, y = c(1, 0, 0), group = as.factor(rep(1:10, each = 3)), sex = rep(as.factor(sex), each  = 3)) %>%
  mutate(alpha = if_else(group == 10, 1, 0.1))

(df %>% group_by(sex) %>% summarise(x = mean(x)))

g <- ggplot(df, aes(x = x, y = y, group = group))
g1 <- g + geom_step(linewidth = 1.3, color = "gray")
g2 <- g + geom_step(aes(alpha = alpha), linewidth = 1.3) + scale_alpha_identity()
g3 <- g1 + geom_smooth(aes(x = x, y = y), data.frame(x = c(0, mean(df$x), 50), y = c(1, 0, 0)),
    inherit.aes = FALSE, method = "loess", se = FALSE, color = "black", linewidth = 1)
g4 <- g +
  geom_step(aes(color = sex), linewidth = 1.3, alpha = 0.5) +
  geom_smooth(aes(x = x, y = y, group = sex, color = sex),
  data.frame(x = c(0, 23.3, 50, 0, 40, 50), y = c(1, 0, 0, 1, 0, 0), sex = as.factor(c(1, 1, 1, 0, 0, 0))),
    inherit.aes = FALSE, method = "loess", se = FALSE, linewidth = 1.5)

g_final <- (g1 + g2 + g3 + g4) + ylim(0, 1) + xlim(0, 50) & labs(x = "t", y = "S(t)") & guides(color  = "none")

ggsave("book/Figures/survtsk/heavisides.png",
  g_final, height=5, width=7, units="in", dpi=600)

## Survival prediction types
library(tidyverse)
library(survival)
library(survminer)
library(ggpubr)
library(patchwork)

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
