# Calibration Chapter Exercise 1
set.seed(1)
event = rbinom(100, 1, 0.7)
times = runif(100)
H = survfit(Surv(times, event) ~ 1)$cumhaz
c("Num deaths" = sum(event), "Sum H" = sum(H))
