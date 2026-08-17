# Daily CRAN downloads for tweedieDistr and related fable packages
# Requires: cranlogs, dplyr, ggplot2

library(cranlogs)
library(dplyr)
library(ggplot2)
library(pkgsearch)
library(tidyr)

packages <- c("tweedieDistr", "fable.intermittent", "fable.bayesRecon")

downloads <- cran_downloads(
  packages = packages,
  from = "2026-05-27",
  to = Sys.Date()
)

# publish date of every CRAN release (including archived versions)
releases <- bind_rows(lapply(packages, function(pkg) {
  h <- cran_package_history(pkg)
  data.frame(package = pkg, version = h$Version, date = as.Date(h$date))
}))

# drop the zeros before each package's first CRAN release (version 0.1.0)
downloads <- releases %>%
  filter(version == "0.1.0") %>%
  select(package, first_release_date = date) %>%
  inner_join(downloads, by = "package") %>%
  filter(date >= first_release_date) %>%
  select(-first_release_date)

# fixed hue order, one color per package (never cycled/reassigned)
palette <- c(
  tweedieDistr        = "#2a78d6",
  fable.intermittent     = "#1baf7a",
  fable.bayesRecon           = "#eb6834"
)

p <- ggplot(downloads, aes(x = date, y = count, color = package)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = palette, name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Daily CRAN downloads",
    subtitle = paste("From", min(downloads$date), "to", Sys.Date()),
    x = NULL,
    y = "Downloads"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e1e0d9"),
    axis.text = element_text(color = "#52514e"),
    axis.title = element_text(color = "#52514e"),
    plot.title = element_text(color = "#0b0b0b", face = "bold"),
    plot.subtitle = element_text(color = "#52514e"),
    legend.position = "top"
  )

print(p)


cumulative_downloads <- downloads %>%
  group_by(package) %>%
  arrange(date) %>%
  mutate(count = cumsum(count))


release_points <- releases %>%
  inner_join(cumulative_downloads, by = c("package", "date"))

p <- ggplot(cumulative_downloads, aes(x = date, y = count, color = package)) +
  geom_line(linewidth = 0.5) +
  geom_point(
    data = release_points, aes(x = date, y = count, color = package),
    size = 6, show.legend = FALSE
  ) +
  geom_text(
    data = release_points, aes(x = date, y = count, label = version),
    color = "white", size = 2.2, fontface = "bold", show.legend = FALSE
  ) +
  scale_color_manual(values = palette, name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Cumulative CRAN downloads",
    subtitle = paste("From 2026-05-01 to", Sys.Date()),
    x = NULL,
    y = "Downloads"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e1e0d9"),
    axis.text = element_text(color = "#52514e"),
    axis.title = element_text(color = "#52514e"),
    plot.title = element_text(color = "#0b0b0b", face = "bold"),
    plot.subtitle = element_text(color = "#52514e"),
    legend.position = "top"
  )

print(p)
