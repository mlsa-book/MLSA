# figure-p1c4-left-truncation.R
# Left-truncation schematic (3 infants on age 0-365 days). Faithful
# reproduction of the left-truncation block in code.R (~L2238-2309).
# Reproduced EXACTLY (theme_bw, base_size 15, dpi 300).
#   steelblue = truncation time; black = time-to-event
#   solid = subject in data; dashed = invisible (t_i < t_i^L)
#   filled steelblue circle/diamond at end of truncation bar marks in/out.
# Outputs: book/Figures/survival/fig-p1c4-left-truncation.{svg,png}

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

## --- LEFT-TRUNCATION (3 infants on age 0-365 days) -------------------------
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
save_pair(p_lt, "book/Figures/survival/fig-p1c4-left-truncation", h = 3.5, w = 8.5)
