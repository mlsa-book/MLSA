# figure-p1c4-censoring.R
# Right-censoring schematic (4 subjects). Faithful reproduction of the
# censoring block in code.R (~L2152-2236). Reproduced EXACTLY (theme_bw,
# base_size 15, dpi 300) -- predates the book ggplot conventions on purpose.
# Outputs: book/Figures/survival/fig-p1c4-censoring.{svg,png}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

## --- shared schematic styling (censoring / left- / right-truncation) -------
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
  ## crop surrounding white margin, keep a small uniform border
  system2("convert", c(paste0(base, ".png"), "-trim", "+repage",
                       "-bordercolor", "white", "-border", "12x12",
                       paste0(base, ".png")))
}

dy <- 0.10  # half-offset between stacked bars

bigger_text <- theme(
  axis.title  = element_text(size = 17),
  axis.text   = element_text(size = 15),
  legend.text = element_text(size = 15)
)

## --- RIGHT-CENSORING (4 subjects) ------------------------------------------
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
save_pair(p_cens, "book/Figures/survival/fig-p1c4-censoring", h = 3.5, w = 8)
