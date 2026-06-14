## figure-p2c8-rocs.R
## P2C8 (Discrimination) · fig-p2c8-rocs
## Illustrative ROC curves: random guess, perfect classifier, and two models.
## Self-contained (no external data; curves are constructed analytically).

source("book/experiments/figure-prep.R")

x        <- c(0, seq.int(0, 1, 0.1))
y_perfect <- c(0, rep(1, 11))        # perfect prediction
y_guess   <- c(0, seq.int(0, 1, 0.1)) # random guess (diagonal)
y_good    <- x^0.6
y_better  <- x^0.3

df <- data.frame(
  x = rep(x, 4),
  y = c(y_guess, y_perfect, y_good, y_better),
  group = factor(
    rep(c("Random guess", "Perfect classifier", "Good model", "Better model"),
        each = length(x)),
    levels = c("Random guess", "Good model", "Better model", "Perfect classifier")
  )
)

g_auc <- ggplot(df, aes(x = x, y = y, colour = group)) +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_colour_manual(values = okabe_ito, name = "Model") +
  labs(x = "False positive rate", y = "True positive rate") +
  theme_square_panel +
  theme(text = element_text(size = 15))

save_fig(g_auc, "book/Figures/evaluation/fig-p2c8-rocs.png", width = 6.5, height = 5)
