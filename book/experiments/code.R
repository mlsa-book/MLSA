remotes::install_github("mlr-org/mlr3proba", ref = 'v0.5.7', upgrade = "never")
remotes::install_github("mlr-org/mlr3", ref = 'v0.16.1', upgrade = "never")
remotes::install_github("mlr-org/paradox", ref = 'v0.11.1', upgrade = "never")

library(ggplot2)
theme_set(theme_bw())

library(distr6)
g = dstr("Gompertz", shape = 2, decorators = "ExoticStatistics")
t = seq.int(0, 1.5, length.out = 100)
functions <- c("Probability density function, y=f(t)", "Cumulative distribution function, y=F(t)", "Hazard function, y=h(t)", "Survival function, y=S(t)")
d = data.frame(t = t, fun = factor(rep(functions, each = 100), levels = functions), y = c(g$pdf(t), g$hazard(t), g$cdf(t), g$survival(t)))
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

######################## pseudo-values ########################################
library(survival)
library(pseudo)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(patchwork)
library(mgcv)
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
library(survival)
library(pseudo)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(patchwork)
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
rm(list = ls())
library(survival)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(pammtools)
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
  ylab("S(t)") +
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
  ylab("S(t)") +
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
rm(list = ls())
library(survival)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(broom)
library(pammtools)
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
  ylab("S(t)") +
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
rm(list = ls())
library(survival)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(pammtools)
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
rm(list = ls())
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

ggsave("book/Figures/reductions/pem-interval-comparison.png", p_pem_comparison,
       height = 6, width = 12, units = "in", dpi = 600)

cat("\nFigure saved to book/Figures/reductions/pem-interval-comparison.png\n")
