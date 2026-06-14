# figure-p1c4-right-truncation.R
# Right-truncation schematic (3 subjects, calendar time 0-30 days). Faithful
# reproduction of the right-truncation block in code.R (~L2311-2382).
# Reproduced EXACTLY (theme_bw, base_size 15, dpi 300).
#   steelblue = truncation time (recruit -> db query)
#   black     = time-to-event (recruit -> event)
#   solid = in registry (event before query); dashed = absent from registry
#   vertical dashed line marks the database-query date.
# Outputs: book/Figures/survival/fig-p1c4-right-truncation.{svg,png}

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

## --- RIGHT-TRUNCATION (3 subjects, calendar time 0-30 days) ----------------
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
save_pair(p_rt, "book/Figures/survival/fig-p1c4-right-truncation", h = 3.5, w = 8.5)
