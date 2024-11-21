## Ranking
library(dplyr)
library(ggplot2)
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
  theme_bw() +
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
  geom_line() + theme_bw() + ylim(0, 1) +
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
  geom_line() + theme_bw() + ylim(0, 1) + xlim(0, 1) +
  geom_abline(slope = 1, intercept = 0, color = "lightgray", lty = "dashed") +
  labs(x = "True (p)", y = "Predicted", color = "Model") +
  geom_label(aes(x = x, y = y), data.frame(x = 0.75, y = 0.1), label = scores,
  inherit.aes = FALSE, hjust = "left", size = 2.5)
ggsave("book/Figures/evaluation/calibD.png", g, height = 3, units = "in",
  dpi = 600)

# Calibration point processes experiment
set.seed(42)
library(survival)

event = rbinom(100, 1, 0.7)
times = runif(100)
H = survfit(Surv(times, event) ~ 1)$cumhaz
c("Num deaths" = sum(event), "Sum H" = sum(H))
