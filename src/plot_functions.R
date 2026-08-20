#library(lemon)

# Convenient function to save figs simultaneously in vectorised and non vectorised format
save_ggplot <- function(p, fp, w, h) {
  mm_to_inches <- 25.3
  w <- w / mm_to_inches
  h <- h / mm_to_inches
  ggsave(paste0(fp, ".png"), p, width = w, height = h, dpi = 400, device = "png")
  ggsave(paste0(fp, ".pdf"), p, width = w, height = h)
}

# Convenient function to change axes using lemon style
axes2lemon <- function(p, bt = "none", lt = "none") {
  p <- p + coord_capped_cart(bottom = bt, left = lt)
  return(p)
}