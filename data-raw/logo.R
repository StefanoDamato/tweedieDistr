## tweedieDistr hex-sticker logo
## Run from the package root: source("data-raw/logo.R")
## Output: man/figures/logo.png

library(ggplot2)
library(tweedieDistr)
library(hexSticker)

# Optional exotic font via sysfonts + showtext
pkg_font <- "sans"
if (requireNamespace("sysfonts", quietly = TRUE) &&
    requireNamespace("showtext", quietly = TRUE)) {
  library(sysfonts)
  library(showtext)
  tryCatch({
    font_add_google("Exo 2", "Exo2")
    showtext_auto()
    pkg_font <- "Exo2"
  }, error = function(e) NULL)
}

# ── Parameters ----------------------------------------------------------------
# p = 1.2 → alpha = 4 > 1 → all three curves are bell-shaped (mode > 0).
# mu varies: mode shifts right, P(X=0) = exp(-lambda) decreases.
# Color mapping (left-to-right for modes, back-to-front for bar z-order):
#   gold  = mu 0.8  — leftmost mode, LARGEST  P(X=0), bar drawn behind
#   red   = mu 1.5  — middle  mode, MEDIUM  P(X=0), bar drawn in front
#   blue  = mu 2.5  — rightmost mode, SMALLEST P(X=0), bar drawn in front of all
p       <- 1.2
phi     <- 1.2
mu_vals <- c(0.8, 1.5, 2.2)
cols    <- c('#E8B86D', '#E05C5C', '#5B8DD9')   # gold (= title/border) / red / blue
gold    <- '#E8B86D'                            # hex border + title colour
bg      <- '#1B2A3B'

# ── Continuous densities ------------------------------------------------------
x_cont <- seq(0.01, 6.0, length.out = 900)
df_fam <- do.call(rbind, lapply(seq_along(mu_vals), function(i)
  data.frame(
    x   = x_cont,
    d   = dtweedie(x_cont, mean = mu_vals[i], dispersion = phi, power = p),
    grp = factor(i)
  )))

p0_vals <- sapply(mu_vals, function(mu)
  dtweedie(0, mean = mu, dispersion = phi, power = p))
ymax <- max(c(max(df_fam$d), p0_vals)) * 1.1

# Reorder factor so gold (grp "1") is drawn last → sits on top in every layer.
# Named colour values keep the mapping correct despite the level reordering.
df_fam$grp <- factor(df_fam$grp, levels = c("1", "2", "3"))

# ── Side-by-side bars with gaps, all right-edges flush at x = 0 --------------
# Same spacing / width as the original version, but shifted so the blue bar's
# right edge lands exactly at x = 0 (the discrete-continuous boundary).
# Order left-to-right: gold (largest P(X=0)) · red · blue (smallest P(X=0)).
bar_hw  <- 0.07  # half-width of each bar
centers <- -0.07 # gold, red, blue  →  right edge of blue = 0

# True stacked bar: each segment occupies its own y range (no y-overlap).
# Draw order blue → red → gold matches the curve layer order (gold on top).
# Visible from bottom to top: blue (tiny) · red (medium) · gold (tall tip).
bars_df <- data.frame(
  xmin = centers - bar_hw,
  xmax = centers + bar_hw,
  ymin = c(0,          p0_vals[3], p0_vals[2]),
  ymax = c(p0_vals[3], p0_vals[2], p0_vals[1]),
  grp  = factor(c(3, 2, 1))          # blue drawn first, gold drawn last = on top
)

# Dot at the top of each segment (same colour as its bar)
dots_df <- data.frame(
  x   = centers,
  y   = p0_vals,
  grp = factor(seq_along(mu_vals))   # 1=gold, 2=red, 3=blue
)

# ── Plot ----------------------------------------------------------------------
gg <- ggplot() +
  geom_rect(
    data = bars_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = grp),
    alpha = 1.00, colour = NA
  ) +
  # fills and lines drawn in factor level order: "2", "3", "1" → gold always on top
  geom_area(
    data = df_fam,
    aes(x, d, fill = grp, group = grp),
    alpha = 0.13, position = 'identity'
  ) +
  geom_line(
    data = df_fam,
    aes(x, d, colour = grp, group = grp),
    linewidth = 1.
  ) +
  geom_point(
    data = dots_df,
    aes(x = x, y = y, colour = grp),
    size = 1.
  ) +
  scale_fill_manual(values   = setNames(cols, c("1","2","3"))) +
  scale_colour_manual(values = setNames(cols, c("1","2","3"))) +
  scale_x_continuous(limits = c(-0.75, 6.0), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, ymax),    expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position  = 'none',
    plot.background  = element_rect(fill = 'transparent', colour = NA),
    panel.background = element_rect(fill = 'transparent', colour = NA)
  )

# ── Hex sticker --------------------------------------------------------------
dir.create('man/figures', showWarnings = FALSE, recursive = TRUE)

sticker(
  gg,
  package  = 'tweedieDistr',
  p_size   = 14,
  p_color  = gold,
  p_family = pkg_font,
  p_y      = 1.5,
  s_x      = 1.05,   # shifted right
  s_y      = 1.00,   # raised to prevent bottom clipping
  s_width  = 1.60,
  s_height = 1.00,
  h_fill   = bg,
  h_color  = gold,
  h_size   = 1.4,
  filename = 'man/figures/logo.png',
  dpi      = 320
)

