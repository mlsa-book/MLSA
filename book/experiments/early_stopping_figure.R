# Early-stopping figure for the "Early stopping" section of the ML chapter (P1C3).
#
# Reuses the gradient-descent example dataset (teengamb.csv): a small subset is
# used for training and the remainder for validation. Gradient descent is run on
# the training set and the per-observation training/validation losses (MSE) are
# tracked over iterations. Because the fully converged (OLS) fit overfits the
# small training set, the validation loss is U-shaped, illustrating early stopping.
#
# Output: book/Figures/introduction/early-stopping.png

library(ggplot2)

exp_dir <- "book/experiments"
fig_dir <- "book/Figures/introduction"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Data: reuse the gradient-descent example dataset --------------------
df <- read.csv(file.path(exp_dir, "teengamb.csv"))
preds <- c("sex", "status", "income", "verbal")
X0 <- as.matrix(df[, preds])
y <- df$gamble

## Degree-2 polynomial expansion (linear terms + squares + pairwise
## interactions). The extra capacity lets the model overfit on a *balanced*
## train/validation split, so both sets have a comparable outcome variance
## (the loss curves start at a similar level) while the validation loss is
## still U-shaped.
poly2 <- function(M) {
  cols <- list(M)
  for (j in seq_len(ncol(M))) cols[[length(cols) + 1]] <- M[, j]^2
  for (a in 1:(ncol(M) - 1)) for (b in (a + 1):ncol(M)) {
    cols[[length(cols) + 1]] <- M[, a] * M[, b]
  }
  do.call(cbind, cols)
}
X <- poly2(X0)

## Balanced split into training and validation sets.
set.seed(5)
n <- nrow(df)
n_train <- 32
idx <- sample(n)
tr <- idx[seq_len(n_train)]
va <- idx[(n_train + 1):n]

## Standardise predictors and centre the response using training statistics.
mu <- colMeans(X[tr, ])
sdv <- apply(X[tr, ], 2, sd)
sdv[sdv == 0] <- 1
Xtr <- scale(X[tr, ], center = mu, scale = sdv)
Xva <- scale(X[va, ], center = mu, scale = sdv)
ybar <- mean(y[tr])
ytr <- y[tr] - ybar
yva <- y[va] - ybar
Atr <- cbind(1, Xtr)
Ava <- cbind(1, Xva)

## ---- Gradient descent from zero (early iterates are implicitly shrunk) ---
theta <- rep(0, ncol(Atr))
lr <- 0.02
n_steps <- 8000
train_loss <- numeric(n_steps)
val_loss <- numeric(n_steps)
for (t in seq_len(n_steps)) {
  grad <- 2 * crossprod(Atr, Atr %*% theta - ytr) / length(ytr)
  theta <- theta - lr * as.numeric(grad)
  train_loss[t] <- mean((Atr %*% theta - ytr)^2)
  val_loss[t] <- mean((Ava %*% theta - yva)^2)
}
stop_it <- which.min(val_loss)

dat <- rbind(
  data.frame(iter = seq_len(n_steps), loss = train_loss, set = "Training loss"),
  data.frame(iter = seq_len(n_steps), loss = val_loss, set = "Validation loss")
)
dat$set <- factor(dat$set, levels = c("Training loss", "Validation loss"))

## ---- Plot ----------------------------------------------------------------
cols <- c("Training loss" = "#2C7FB8", "Validation loss" = "#D7301F")

g <- ggplot(dat, aes(iter, loss, colour = set)) +
  annotate("rect", xmin = stop_it, xmax = n_steps, ymin = -Inf, ymax = Inf,
           fill = "#D7301F", alpha = 0.05) +
  geom_vline(xintercept = stop_it, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_line(linewidth = 1) +
  geom_point(data = data.frame(iter = stop_it, loss = val_loss[stop_it]),
             aes(iter, loss), inherit.aes = FALSE,
             colour = "#D7301F", size = 2.6) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = "Iteration", y = "Loss (MSE)") +
  theme_bw() +
  theme(legend.position = c(0.015, 0.985), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = "grey80"))

ggsave(file.path(fig_dir, "early-stopping.png"), g,
       width = 6, height = 3, units = "in", dpi = 600)

cat(sprintf("Early-stopping iteration: %d (of %d)\n", stop_it, n_steps))
cat(sprintf("Val MSE at stop = %.1f; final val MSE = %.1f\n",
            val_loss[stop_it], val_loss[n_steps]))
cat(sprintf("Train MSE at stop = %.1f; final train MSE = %.1f\n",
            train_loss[stop_it], train_loss[n_steps]))
