## figure-p1c6-predict-types.R
## P1C6 (Survival task) · fig-survtsk-overview
## Output: book/Figures/survtsk/predict_types.{svg,png}
##
## Prediction-types overview: six tumor patients and a Weibull AFT model,
## showing tabular data (A), predicted survival times with CIs (B), relative
## risk scores (C), and predicted survival distributions (D).
## Self-contained reproduction; decision = KEEP (figure must not change).

suppressPackageStartupMessages({
  library(flexsurv)
  library(survival)
  library(pammtools)
  library(patchwork)
  library(svglite)
  library(dplyr)
})

source("book/experiments/figure-prep.R")

## Pastel per-subject palette (one colour per patient, used across all panels).
pastel_pt <- c("#A1C9F4", "#FFB482", "#8DE5A1", "#FF9F9B", "#D0BBFF", "#FAB0E4")

## --- Model ---------------------------------------------------------------
data(tumor, package = "pammtools")
fit_pt <- flexsurvreg(
  Surv(days, status) ~ age + sex + complications,
  anc  = list(shape = ~ complications),
  data = tumor,
  dist = "weibull"
)

## --- Six example patients ------------------------------------------------
patients_pt <- tibble::tibble(
  id            = 1:6,
  age           = c(71, 70, 67, 58, 39, 59),
  sex           = factor(c("female", "male", "female", "male", "female", "female"),
                         levels = c("male", "female")),
  complications = factor(c("no", "no", "yes", "no", "yes", "no"),
                         levels = c("no", "yes")),
  days          = c(1217, 519, 2414, 397, 1217, 268),
  status        = c(1, 0, 0, 1, 0, 1)
)

mean_pred_pt    <- summary(fit_pt, newdata = patients_pt, type = "mean", tidy = TRUE)
patients_pt$E_t  <- mean_pred_pt$est
patients_pt$lcl  <- mean_pred_pt$lcl
patients_pt$ucl  <- mean_pred_pt$ucl
patients_pt$risk <- -log(patients_pt$E_t)
patients_pt$risk <- patients_pt$risk - mean(patients_pt$risk)
patients_pt$id_f <- factor(patients_pt$id, levels = patients_pt$id)

tgrid_pt     <- seq(0, 3000, length.out = 250)
surv_pred_pt <- summary(fit_pt, newdata = patients_pt, type = "survival",
                        t = tgrid_pt, tidy = TRUE)
surv_pred_pt$id <- factor(rep(patients_pt$id, each = length(tgrid_pt)),
                          levels = patients_pt$id)

## --- Panel A: tabular data drawn as a ggplot -----------------------------
col_names_pt  <- c("id", "age", "sex", "complications", "days", "status")
col_widths_pt <- c(0.5, 0.5, 0.85, 1.75, 0.7, 0.85)
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
            size = 5) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_x_continuous(limits = c(min(col_starts_pt), max(col_ends_pt)),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(-nrow(tbl_long_pt) - 0.6, 0.6),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = NULL, title = "Input data") +
  theme_bw(base_size = 15) +
  theme(plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        axis.title       = element_blank(),
        panel.grid       = element_blank())

## --- Panel B: predicted survival times with 95% CIs ----------------------
p_time_pt <- ggplot(patients_pt, aes(x = E_t, y = id_f, colour = id_f)) +
  geom_linerange(aes(xmin = lcl, xmax = ucl), linewidth = 1.2) +
  geom_point(size = 4) +
  scale_colour_manual(values = pastel_pt) +
  scale_x_continuous(limits = c(0, max(patients_pt$ucl) * 1.05),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Predicted survival time (days)", y = "Subject",
       title = "Predicted survival times") +
  guides(colour = "none") +
  theme_bw(base_size = 15) +
  theme(plot.title         = element_text(face = "bold", size = 16, hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())

## --- Panel C: predicted relative-risk scores (lollipop) ------------------
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
  theme_bw(base_size = 15) +
  theme(plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
        panel.grid.minor = element_blank())

## --- Panel D: predicted survival distributions ---------------------------
p_surv_pt <- ggplot(surv_pred_pt, aes(x = time, y = est,
                                      colour = id, group = id)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = pastel_pt) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Time", y = "Survival Probability",
       title = "Predicted survival distributions") +
  guides(colour = "none") +
  theme_bw(base_size = 15) +
  theme(plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
        panel.grid.minor = element_blank())

final_pt <- (p_table_pt | p_time_pt) / (p_risk_pt | p_surv_pt) +
  plot_layout(widths = c(1, 1), heights = c(1, 1))

ggsave("book/Figures/survtsk/fig-p1c6-predict-types.svg", final_pt,
       width = 9.5, height = 6.5, units = "in")
ggsave("book/Figures/survtsk/fig-p1c6-predict-types.png", final_pt,
       width = 9.5, height = 6.5, units = "in", dpi = 300)
trim_png("book/Figures/survtsk/fig-p1c6-predict-types.png")
