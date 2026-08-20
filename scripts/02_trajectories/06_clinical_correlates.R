# Compare clinical presentation across trajectories

rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(lemon)

# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

# Generalized function to plot donut charts for any column
plot_donut <- function(data, column_to_plot, color_palette, cohort_label = "Cohort") {
  
  # 1. Convert the input data to a standard dataframe
  df <- as.data.frame(data)
  
  # 2. Summarize the data by the specified column
  df_summary <- df %>%
    dplyr::filter(!is.na(.data[[column_to_plot]])) %>%
    group_by(.data[[column_to_plot]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      fraction = count / sum(count),
      percent_label = ifelse(fraction > 0.02, paste0(round(fraction * 100, 1), "%"), "")
    )

  # Get total N for the center text
  total_n <- sum(df_summary$count)

  # 4. Build the ggplot
  p <- ggplot(df_summary, aes(x = 2, y = fraction, fill = .data[[column_to_plot]])) +
    # The donut ring
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.25) +
    
    # The percentage text centered inside the slices
    geom_text(aes(label = percent_label), 
              position = position_stack(vjust = 0.5), 
              color = "white", size = 2.5, fontface = "bold") +
    
    # Wrap into polar coordinates
    coord_polar(theta = "y", start = 0) +
    
    # Leave the donut hole (x limits 0.5 to 2.5)
    xlim(0.5, 2.5) +
    
    # Center annotation text
    annotate("text", x = 0.5, y = 0, label = paste0(cohort_label, "\nn = ", total_n), 
             size = 3, fontface = "bold", color = "black", lineheight = 0.8) +
    
    # Apply colors
    scale_fill_manual(values = color_palette) +

    # Minimal Nature-ready theme
    theme_void() + 
    theme(legend.position = "none")

  return(p)
}



## PATHS
outdir <- "outputs/02_trajectories"
figdir <- "figures/02_trajectories/clinical"
clinfile <- "data/meta/PPCG_donors_clin_20241217.csv"
dir.create(outdir, recursive = TRUE); dir.create(figdir, recursive = TRUE) 
trajectories_fps <- list.files(outdir, pattern = "PPCG_Feb2026.*mergedseg_with_clonality.txt", full = T)

## LOAD DATA
tj_samples <- load_trajectories(trajectories_fps)
clin_data <- read_delim(clinfile, delim = ",")
tj_samples$record_id <- str_remove(tj_samples$sample, "[a-z]_DNA")
tj_samples = merge(tj_samples, clin_data)
tj_samples$gs_group = get_gs_group(tj_samples$record_id)

# Plot clinical presentation across trajectories
# Grade Group
gg_tj1 = plot_donut(tj_samples %>% dplyr::filter(trajectory == "Ordering 1"), "gs_group", color_palette = gleason_colours, cohort_label = "Canonical")
gg_tj2 = plot_donut(tj_samples %>% dplyr::filter(trajectory == "Ordering 2"), "gs_group", color_palette = gleason_colours, cohort_label = "Alternative")
gg_tj3 = plot_donut(tj_samples %>% dplyr::filter(trajectory == "Ordering 3"), "gs_group", color_palette = gleason_colours, cohort_label = "Gain-enriched")

write_tsv(tj_samples %>% dplyr::filter(trajectory == "Ordering 1") %>% dplyr::select(gs_group), file.path(outdir, "Fig2e_source_data.tsv"))
ggsave(file.path(figdir, "Fig2e_pie_grade_canonical.pdf"), gg_tj1, height = 30/25.4, width = 30/25.4, dpi = 300)
write_tsv(tj_samples %>% dplyr::filter(trajectory == "Ordering 2") %>% dplyr::select(gs_group), file.path(outdir, "Fig2f_source_data.tsv"))
ggsave(file.path(figdir, "Fig2f_pie_grade_alternative.pdf"), gg_tj2, height = 30/25.4, width = 30/25.4, dpi = 300)
write_tsv(tj_samples %>% dplyr::filter(trajectory == "Ordering 3") %>% dplyr::select(gs_group), file.path(outdir, "Fig2g_source_data.tsv"))
ggsave(file.path(figdir, "Fig2g_pie_grade_gain_enriched.pdf"), gg_tj3, height = 30/25.4, width = 30/25.4, dpi = 300)

# Significant enrichment of Grade Group 4+ in Gain-enriched
tb = table(tj_samples$gs_group, tj_samples$trajectory)
chisq.test(tb) # p-value = 1.36-e6

fisher.test(matrix(
    c(tb["4+", "Ordering 3"], sum(tb["4+", c("Ordering 1", "Ordering 2")]), 
      sum(tb[c("1", "2", "3"), "Ordering 3"]), sum(tb[c("1", "2", "3"), c("Ordering 1", "Ordering 2")])), nrow = 2, byrow = T)
)

# Early onset prostate cancer (55 years or lower)
tj_samples$eopc <- ifelse(tj_samples$age_at_tumour_collection < 55, "Early-onset", "Late-onset")
eopc_colors <- c("Early-onset" = "#5E6C85", "Late-onset" = "grey75")
write_tsv(tj_samples %>% dplyr::filter(trajectory == "Ordering 1") %>% dplyr::select(eopc), file.path(outdir, "Fig2e_eopc_source_data.tsv"))
eopc_tj1 = plot_donut(tj_samples %>% dplyr::filter(trajectory == "Ordering 1"), "eopc", color_palette = eopc_colors, cohort_label = "Canonical")
write_tsv(tj_samples %>% dplyr::filter(trajectory == "Ordering 2") %>% dplyr::select(eopc), file.path(outdir, "Fig2f_eopc_source_data.tsv"))
eopc_tj2 = plot_donut(tj_samples %>% dplyr::filter(trajectory == "Ordering 2"), "eopc", color_palette = eopc_colors, cohort_label = "Alternative")
write_tsv(tj_samples %>% dplyr::filter(trajectory == "Ordering 3") %>% dplyr::select(eopc), file.path(outdir, "Fig2g_eopc_source_data.tsv"))
eopc_tj3 = plot_donut(tj_samples %>% dplyr::filter(trajectory == "Ordering 3"), "eopc", color_palette = eopc_colors, cohort_label = "Gain-enriched")

ggsave(file.path(figdir, "Fig2e_pie_eopc_canonical.pdf"), eopc_tj1, height = 30/25.4, width = 30/25.4, dpi = 300)
ggsave(file.path(figdir, "Fig2f_pie_eopc_alternative.pdf"), eopc_tj2, height = 30/25.4, width = 30/25.4, dpi = 300)
ggsave(file.path(figdir, "Fig2g_pie_eopc_gain_enriched.pdf"), eopc_tj3, height = 30/25.4, width = 30/25.4, dpi = 300)

    
# Significant enrichment of late-onset prostate cancer in Alternative
tb_age = table(tj_samples$eopc, tj_samples$trajectory) 
chisq.test(tb_age) # 6.28e-5

fisher.test(matrix(
    c(tb_age["Late-onset", "Ordering 2"], sum(tb_age["Late-onset", c("Ordering 1", "Ordering 3")]), 
      sum(tb_age["Early-onset", "Ordering 2"]), sum(tb_age["Early-onset", c("Ordering 1", "Ordering 3")])), nrow = 2, byrow = T)
)
