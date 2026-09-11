# MRCA commitment analysis
# This script analyses the trajectory along which the MRCA was evolving with respect to final trajectory

rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(ggalluvial)

# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

# PATHS
outdir <- "outputs/03_mechanism/mrca_commitment"
figdir <- "figures/03_mechanism/mrca_commitment"
dir.create(outdir, recursive = TRUE); dir.create(figdir, recursive = TRUE)
mrca_assignments_fp <- file.path(outdir, "mrca_trajectory_assignments.tsv")

# LOAD DATA
df = read_delim(mrca_assignments_fp)

trajectory_names = c("Ordering 1" = "Canonical", "Ordering 2" = "Alternative", "Ordering 3" = "Gain-enriched")
df$mrca_trajectory = trajectory_names[df$mrca_trajectory]
df$final_trajectory = trajectory_names[df$final_trajectory]

traj_levels <- c("Canonical", "Alternative", "Gain-enriched")
df$mrca_trajectory <- factor(df$mrca_trajectory, levels = traj_levels)
df$final_trajectory <- factor(df$final_trajectory, levels = traj_levels)

df_sankey <- df %>%
  group_by(mrca_trajectory) %>%
  mutate(total = sum(n),
         label = if_else(mrca_trajectory == final_trajectory, paste0(n, "/", total, " (", round(n/total*100), "%)"), ""))

traj_colors <- c(
  "Canonical" = "#0072B2",   # Deep Cerulean Blue
  "Alternative" = "#1B9E77",   # Sea Green
  "Gain-enriched" = "#D73027"   # Rich Crimson
)

# 5. Generate the Plot
p <- ggplot(df_sankey, aes(y = n, axis1 = mrca_trajectory, axis2 = final_trajectory)) +
  # Draw flows (alluvia)
  geom_alluvium(aes(fill = mrca_trajectory), width = 0.15, alpha = 0.2, knot.pos = 0.4) +
  # Draw vertical bars (strata)
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.05, alpha = 0.8, color = "white", linewidth = 0.2) +
  # Apply biological colors to flows
  scale_fill_manual(values = traj_colors) +
  # Styling
  theme_void() +
  theme(legend.position = "none")

write_tsv(df_sankey %>% dplyr::select(mrca_trajectory, final_trajectory, n), file.path(outdir, "Fig3a_source_data.tsv"))
ggsave(file.path(figdir, "Fig3a_sankey_trajectories.png"), p, height = 40/25.4, width = 40/25.4, dpi = 300)
ggsave(file.path(figdir, "Fig3a_sankey_trajectories.pdf"), p, height = 50/25.4, width = 60/25.4, dpi = 300)

