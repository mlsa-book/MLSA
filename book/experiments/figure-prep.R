## figure-prep.R — shared setup for MLSA per-figure scripts.
##
## Source this at the top of every figure-<chapter>-<desc>.R, then call save_fig().
## It owns the book-wide conventions so the individual figure scripts only carry
## figure-specific code:
##   - colour palettes + fixed variable->colour mappings
##   - the global ggplot theme (theme_bw, legend on the right, white facet strips)
##   - the export helper (always 600 dpi)
## See book/figure-specs.md for the rationale behind each choice.

suppressPackageStartupMessages(library(ggplot2))

## --- Colour --------------------------------------------------------------

## Okabe-Ito colourblind-safe palette (categorical scales; NOT the ggplot default)
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")

## Fixed, book-wide mappings (same variable -> same encoding everywhere)
col_complications <- c(no = "#0072B2", yes = "#D55E00")           # blue / vermillion
col_treatment     <- c(Placebo = "#0072B2", Prednisone = "#D55E00")

## Cover palette (for the rare ggplot that needs the infographic colours)
cover <- list(cream = "#f7f3e8", navy = "#18222e", gold = "#c8a800",
              cyan = "#4abcd8", blue_mid = "#5a8ab0", blue = "#6a9ac8",
              blue_light = "#9abbd8", blue_deep = "#4a7a9a")

## --- Theme ---------------------------------------------------------------

## theme_bw always; legend on the right with word titles; white facet strips.
theme_set(
  theme_bw() +
    theme(
      legend.position  = "right",
      strip.background = element_rect(fill = "white", colour = "grey30"),
      strip.text       = element_text(colour = "black")
    )
)

## Square coordinate panel for NON-faceted survival-curve plots (KM / CIF /
## transition / pseudo-value / Brier). Add with `+ theme_square_panel`.
theme_square_panel <- theme(aspect.ratio = 1)

## For FACETED plots: make each whole tile (facet strip + coordinate panel)
## square, i.e. coordinate_width = coordinate_height + strip_height. Since the
## strip sits on top, the coordinate box must be a touch wider than tall:
##   aspect.ratio = coord_h / coord_w = coord_h / (coord_h + strip_h) = 1/(1+f),
## where f = strip_height / coord_height (~0.10 for default strip text + a few
## facet rows; tune per figure). Add with `+ theme_square_tile()`.
theme_square_tile <- function(strip_frac = 0.10) theme(aspect.ratio = 1 / (1 + strip_frac))

## --- Export --------------------------------------------------------------

## Post-crop surrounding white margins of a PNG with ImageMagick (square-panel
## plots leave whitespace above/below the panel via aspect.ratio), re-adding a
## small uniform white border so the figure is not flush to its edges. fuzz = 0
## so only pure-white border pixels go — plot content (grey panel border, any
## coloured fill) is never clipped. Call directly after a bare ggsave().
trim_png <- function(path) {
  system2("convert", c(shQuote(path), "-trim", "+repage",
                       "-bordercolor", "white", "-border", "12x12",
                       shQuote(path)))
}

## Always export at 600 dpi. width/height in inches.
## trim (default TRUE) crops white margins via trim_png(); pass trim = FALSE to
## keep the raw ggsave margins.
save_fig <- function(plot, path, width, height, dpi = 600, trim = TRUE, ...) {
  ggsave(filename = path, plot = plot,
         width = width, height = height, units = "in", dpi = dpi, ...)
  if (trim) trim_png(path)
}
