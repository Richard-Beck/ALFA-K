## store theme element for consistent plotting. 
make_base_theme <- function(base_theme = "classic",
                            base_size  = 5,
                            base_family = "sans") {
  # validate and pick one of the allowed themes
  base_theme <- match.arg(base_theme,
                          choices = c("classic", "minimal", "bw", "gray","void"))
  theme_fn <- switch(base_theme,
                     classic = theme_classic,
                     minimal = theme_minimal,
                     bw      = theme_bw,
                     gray    = theme_gray,
                     void=theme_void)
  theme_fn(base_size = base_size, base_family = base_family) +
    theme(
      text              = element_text(size = base_size, family = base_family),
      axis.title        = element_text(size = base_size, family = base_family),
      axis.text         = element_text(size = base_size, family = base_family),
      legend.title      = element_text(size = base_size, family = base_family),
      legend.text       = element_text(size = base_size, family = base_family),
      strip.text        = element_text(size = base_size, family = base_family),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(3, "mm")
    )
}
