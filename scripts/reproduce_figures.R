# =============================================================================
# scripts/reproduce_figures.R
#
# Standalone script to reproduce every published panel in the PPCG Evolution
# manuscript from pre-saved source-data files. Loads each *_source_data.{tsv,rds}
# from figures_data/ and writes the corresponding PDF (and PNG where the
# original used save_ggplot()) to figures/reproduced/, preserving the original
# filenames produced by the analysis scripts under scripts/.
#
# Panel order: Main Figures (Fig*) -> Extended Data (EXDF*) -> Supplementary
# (SF*), alphabetical within each block.
#
# Schematics (Fig1b, Fig2a, Fig4c, Fig5b, EXDF10a) and externally-generated
# panels (EXDF5, SF7a) are skipped; see README.md for the companion repos.
#
# To populate figures_data/ from the analysis outputs:
#   mkdir -p figures_data
#   find outputs -name "*_source_data.*" -exec cp {} figures_data/ \;
#
# Usage:
#   conda activate ppevo
#   Rscript scripts/reproduce_figures.R
# =============================================================================

rm(list = ls(all = TRUE))

# ---- 1. PACKAGES ------------------------------------------------------------
suppressPackageStartupMessages({
  library(MASS)
  library(tidyverse)
  library(patchwork)
  library(cowplot)
  library(lemon)
  library(survival)
  library(survminer)
  library(ggsankey)
  library(ggalluvial)
  library(UpSetR)
  library(ggpubr)
  library(ggh4x)
  library(RColorBrewer)
  library(scales)
  library(here)
  library(hrbrthemes)
  library(forcats)
  library(GenomicRanges)
})

# ---- 2. SHARED THEMES / PALETTES (from src/) --------------------------------
setwd(here::here())
source("src/plot_theme.R")        # theme_set, gleason_colours, trajectory_colours, ...
source("src/plot_functions.R")    # save_ggplot(), axes2lemon()

# Trajectory-named alias used by several scripts
trajectory_colors_named <- c("Canonical" = "#8DA0CB",
                             "Alternative" = "#AAF0C9",
                             "Gain-enriched" = "#C04667")

# ---- 3. PATHS ---------------------------------------------------------------
source_data_dir <- here::here("figures_data")
fig_out_dir     <- here::here("figures", "reproduced")
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(source_data_dir)) {
  stop(
    "figures_data/ not found.\n",
    "Populate it by running the analysis scripts and copying the source files:\n",
    "  mkdir -p figures_data\n",
    "  find outputs -name '*_source_data.*' -exec cp {} figures_data/ \\;"
  )
}

# Convenience accessors
sd_tsv <- function(panel) {
  read_delim(file.path(source_data_dir, paste0(panel, "_source_data.tsv")),
             delim = "\t", show_col_types = FALSE)
}
sd_rds <- function(panel) {
  readRDS(file.path(source_data_dir, paste0(panel, "_source_data.rds")))
}
fp <- function(name) file.path(fig_out_dir, name)

message("Reproduced figures will be written to: ", fig_out_dir)


# ---- 4. LOCAL HELPERS (inlined verbatim from the source scripts) ------------

## plot_donut() — from scripts/02_trajectories/06_clinical_correlates.R
plot_donut <- function(data, column_to_plot, color_palette, cohort_label = "Cohort") {
  df <- as.data.frame(data)
  df_summary <- df %>%
    dplyr::filter(!is.na(.data[[column_to_plot]])) %>%
    group_by(.data[[column_to_plot]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      fraction = count / sum(count),
      percent_label = ifelse(fraction > 0.02,
                             paste0(round(fraction * 100, 1), "%"), "")
    )
  total_n <- sum(df_summary$count)
  ggplot(df_summary, aes(x = 2, y = fraction, fill = .data[[column_to_plot]])) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.25) +
    geom_text(aes(label = percent_label),
              position = position_stack(vjust = 0.5),
              color = "white", size = 2.5, fontface = "bold") +
    coord_polar(theta = "y", start = 0) +
    xlim(0.5, 2.5) +
    annotate("text", x = 0.5, y = 0,
             label = paste0(cohort_label, "\nn = ", total_n),
             size = 3, fontface = "bold", color = "black", lineheight = 0.8) +
    scale_fill_manual(values = color_palette) +
    theme_void() +
    theme(legend.position = "none")
}

## plot_trajectory_donut() — from scripts/02_trajectories/07_ancestry.R
plot_trajectory_donut <- function(traj_table, cohort_label = "Cohort") {
  df <- as.data.frame(traj_table)
  colnames(df) <- c("trajectory", "count")
  df <- df %>%
    mutate(
      traj_name = case_when(
        trajectory == "Ordering 1" ~ "Canonical",
        trajectory == "Ordering 2" ~ "Alternative",
        trajectory == "Ordering 3" ~ "Gain-enriched",
        TRUE ~ as.character(trajectory)
      ),
      fraction = count / sum(count),
      percent_label = ifelse(fraction > 0.02,
                             paste0(round(fraction * 100, 1), "%"), "")
    ) %>%
    mutate(traj_name = factor(traj_name,
                              levels = c("Gain-enriched", "Alternative", "Canonical")))
  total_n <- sum(df$count)
  names(trajectory_colors) <- c("Canonical", "Alternative", "Gain-enriched")
  ggplot(df, aes(x = 2, y = fraction, fill = traj_name)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.25) +
    geom_text(aes(label = percent_label),
              position = position_stack(vjust = 0.5),
              color = "white", size = 2.5, fontface = "bold") +
    coord_polar(theta = "y", start = 0) +
    xlim(0.5, 2.5) +
    annotate("text", x = 0.5, y = 0,
             label = paste0(cohort_label, "\nn = ", total_n),
             size = 3, fontface = "bold", color = "black", lineheight = 0.8) +
    scale_fill_manual(values = trajectory_colors) +
    theme_void() +
    theme(legend.position = "none")
}

## plot_timeline_raincloud() — from scripts/01_landscape/04d_timing_analysis.R
##                              (Fig1c, gs_group-based timelines)
bootstrap_timeline <- function(df, b = 1000) {
  boot_timeline <- data.frame()
  df <- df %>% dplyr::filter(!is.na(event))
  for (i in 1:b) {
    idx <- sample(1:nrow(df), size = nrow(df), replace = TRUE)
    boot_sample <- df[idx, ] %>%
      group_by(event) %>%
      summarise(estimate = median(latency_rel, na.rm = TRUE),
                .groups = "drop")
    boot_timeline <- bind_rows(boot_timeline, boot_sample)
  }
  boot_timeline %>%
    group_by(event) %>%
    summarise(median = median(estimate),
              lowci = quantile(estimate, 0.025, na.rm = TRUE),
              highci = quantile(estimate, 0.975, na.rm = TRUE),
              .groups = "drop")
}

plot_timeline_raincloud_gs <- function(timeline_data, item_plot = c("MRCA", "CNA", "WGD")) {
  timeline_data <- timeline_data %>% dplyr::filter(event %in% item_plot)
  t_gs1  <- timeline_data[timeline_data$gs_group == "1", ]
  t_gs2  <- timeline_data[timeline_data$gs_group == "2", ]
  t_gs3  <- timeline_data[timeline_data$gs_group == "3", ]
  t_gs45 <- timeline_data[timeline_data$gs_group %in% c("4", "4+", "5"), ]
  timeline_plot <- rbind(
    dplyr::mutate(bootstrap_timeline(t_gs1),  ggroup = "Grade Group 1"),
    dplyr::mutate(bootstrap_timeline(t_gs2),  ggroup = "Grade Group 2"),
    dplyr::mutate(bootstrap_timeline(t_gs3),  ggroup = "Grade Group 3"),
    dplyr::mutate(bootstrap_timeline(t_gs45), ggroup = "Grade Group 4-5")
  ) %>% mutate(y_label = ggroup)
  raw_data <- timeline_data %>%
    mutate(ggroup = case_when(
      gs_group == "1" ~ "Grade Group 1",
      gs_group == "2" ~ "Grade Group 2",
      gs_group == "3" ~ "Grade Group 3",
      gs_group %in% c("4", "4+", "5") ~ "Grade Group 4-5",
      TRUE ~ NA_character_
    )) %>% filter(!is.na(ggroup))
  raw_all <- timeline_data %>% mutate(ggroup = "All")
  raw_plot_data <- bind_rows(raw_data, raw_all) %>%
    left_join(timeline_plot %>% dplyr::select(ggroup, y_label), by = "ggroup")
  order_levels <- timeline_plot$y_label
  timeline_plot$y_label <- factor(timeline_plot$y_label, levels = unique(order_levels))
  raw_plot_data$y_label <- factor(raw_plot_data$y_label, levels = unique(order_levels))
  ggplot() +
    geom_segment(data = timeline_plot,
                 aes(x = 0, xend = 1, y = y_label, yend = y_label),
                 linewidth = 2, col = "grey70", alpha = 0.2) +
    geom_jitter(data = raw_plot_data,
                aes(x = latency_rel, y = y_label, color = event),
                height = 0.15, width = 0, size = 1, alpha = 0.1) +
    geom_point(data = timeline_plot,
               aes(y = y_label, x = median, color = event),
               size = 5, stroke = 0.5, alpha = 0.8) +
    scale_colour_manual(values = c("MRCA" = "#D95F02",
                                   "CNA"  = "#9467bd",
                                   "WGD"  = "#3979bb")) +
    labs(y = "", x = "Fraction of time to diagnosis") +
    scale_y_discrete(
      limits = c("Grade Group 1", "Grade Group 2", "Grade Group 3", "Grade Group 4-5"),
      labels = c(paste0("Grade Group 1\n (n = ", n_distinct(t_gs1$sample), ")"),
                 paste0("Grade Group 2\n (n = ", n_distinct(t_gs2$sample), ")"),
                 paste0("Grade Group 3\n (n = ", n_distinct(t_gs3$sample), ")"),
                 paste0("Grade Group 4-5\n (n = ", n_distinct(t_gs45$sample), ")"))
    ) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")
}

## plot_timeline_raincloud_tj() — trajectory variant (SF12)
plot_timeline_raincloud_tj <- function(timeline_data, item_plot = c("MRCA", "CNA", "WGD")) {
  timeline_data <- timeline_data %>% dplyr::filter(event %in% item_plot)
  t_tj1 <- timeline_data[timeline_data$trajectory == "Canonical", ]
  t_tj2 <- timeline_data[timeline_data$trajectory == "Alternative", ]
  t_tj3 <- timeline_data[timeline_data$trajectory == "Gain-enriched", ]
  timeline_plot <- rbind(
    dplyr::mutate(bootstrap_timeline(t_tj1), trajectory = "Canonical"),
    dplyr::mutate(bootstrap_timeline(t_tj2), trajectory = "Alternative"),
    dplyr::mutate(bootstrap_timeline(t_tj3), trajectory = "Gain-enriched")
  ) %>% mutate(y_label = trajectory)
  raw_data <- timeline_data %>% filter(!is.na(trajectory))
  raw_all  <- timeline_data %>% mutate(trajectory = "All")
  raw_plot_data <- bind_rows(raw_data, raw_all) %>%
    left_join(timeline_plot %>% dplyr::select(trajectory, y_label), by = "trajectory")
  order_levels <- timeline_plot$y_label
  timeline_plot$y_label <- factor(timeline_plot$y_label, levels = unique(order_levels))
  raw_plot_data$y_label <- factor(raw_plot_data$y_label, levels = unique(order_levels))
  ggplot() +
    geom_segment(data = timeline_plot,
                 aes(x = 0, xend = 1, y = y_label, yend = y_label),
                 linewidth = 2, col = "grey70", alpha = 0.2) +
    geom_jitter(data = raw_plot_data,
                aes(x = latency_rel, y = y_label, color = event),
                height = 0.15, width = 0, size = 1, alpha = 0.1) +
    geom_point(data = timeline_plot,
               aes(y = y_label, x = median, color = event),
               size = 5, stroke = 0.5, alpha = 0.8) +
    scale_colour_manual(values = c("MRCA" = "#D95F02",
                                   "CNA"  = "#9467bd",
                                   "WGD"  = "#3979bb")) +
    labs(y = "", x = "Fraction of time to diagnosis", color = "Event") +
    scale_y_discrete(limits = c("Canonical", "Alternative", "Gain-enriched")) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")
}

## plot_clinical_strip() — from scripts/01_landscape/01_landscape_subclonality.R
plot_clinical_strip <- function(data, var_col, var_lab, palette) {
  ggplot(data, aes(x = smp, y = 1, fill = as.character(!!sym(var_col)))) +
    geom_tile() +
    scale_fill_manual(values = palette, na.value = "#FFFFFF", name = var_lab) +
    facet_grid(. ~ group, scales = "free_x", space = "free") +
    labs(y = var_lab) +
    theme_void() +
    theme(
      strip.text = element_blank(),
      axis.title.y = element_text(angle = 0, size = 8, hjust = 1, margin = margin(r = 5)),
      legend.position = "right",
      legend.key.size = unit(0.3, "cm"),
      legend.margin = margin(l = 0, r = 0, t = 0, b = 0),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7)
    )
}

## plot_facet() — from scripts/01_landscape/04a_sbs2age.R (SF2)
get_pval_cor <- function(x, y, test = "pearson") {
  crt <- cor.test(x, y, method = test)
  if (crt$p.value < 0.001) "p < 0.001" else as.character(round(crt$p.value, 3))
}
get_estimate_cor <- function(x, y, test = "pearson") {
  round(cor.test(x, y, method = test)$estimate, 3)
}
plot_facet <- function(df, assignment = "soft", age_var) {
  df_long <- df %>%
    pivot_longer(cols = ends_with(assignment),
                 names_to = "signature", values_to = "mutations")
  df_long$mutations <- as.numeric(df_long$mutations)
  df_long$signature <- str_replace(df_long$signature, "_", " ") %>%
    str_replace("sbs", "SBS") %>% str_replace("clock", "Clock")
  cor_data <- df_long %>%
    dplyr::group_by(signature) %>%
    dplyr::summarise(cor = get_estimate_cor(age_at_tumour_collection, mutations, "spearman"),
                     p_value = get_pval_cor(age_at_tumour_collection, mutations, "spearman")) %>%
    mutate(cor_label = paste0("rho = ", round(cor, 2), "; ", p_value))
  p <- ggplot(df_long, aes(x = age_at_tumour_collection, y = mutations)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "rlm") +
    geom_text(data = cor_data, aes(x = Inf, y = Inf, label = cor_label),
              hjust = 1.1, vjust = 1.5, size = 3, inherit.aes = FALSE) +
    facet_wrap(~signature, nrow = 1, scales = "free_y") +
    theme(strip.text = element_text(size = 10)) +
    labs(x = "Age at tumour collection", y = "Number Mutations")
  axes2lemon(p, bt = "both", lt = "both")
}

## plot_median_timing() — from scripts/01_landscape/04c_time_landmarks.R (SF4)
plot_median_timing <- function(timing_df, color) {
  t <- timing_df %>% dplyr::select(sample_id, latency, acceleration) %>%
    pivot_wider(names_from = acceleration, values_from = latency)
  cols <- colnames(t)
  smr <- t %>% as.data.frame() %>% dplyr::summarise_at(cols, median, na.rm = TRUE)
  ggplot(smr) +
    geom_rect(aes(xmin = 1 - 0.05, xmax = 1 + 0.05, ymax = `1x`,  ymin = `2.5x`), fill = color, size = .15, col = "black") +
    geom_rect(aes(xmin = 1 - 0.05, xmax = 1 + 0.05, ymax = `2.5x`, ymin = `5x`),  fill = color, size = .15, col = "black", alpha = .75) +
    geom_rect(aes(xmin = 1 - 0.05, xmax = 1 + 0.05, ymax = `5x`,   ymin = `10x`), fill = color, size = .15, col = "black", alpha = .5) +
    geom_rect(aes(xmin = 1 - 0.05, xmax = 1 + 0.05, ymax = `10x`,  ymin = `20x`), fill = color, size = .15, col = "black", alpha = .25) +
    ylim(c(0, 20)) +
    labs(x = "", y = "Latency (Years)") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank())
}

plot_cohort_diagram <- function(timing_df) {
  df <- timing_df %>%
    dplyr::filter(acceleration == "1x") %>%
    dplyr::select(sample_id, mrca_latency = latency,
                  record_id, age_at_tumour_collection,
                  relapse_ind, donor_relapse_interval, gs_grade)
  df$age_at_tumour_collection <- as.numeric(df$age_at_tumour_collection)
  df$age_mrca <- df$age_at_tumour_collection - df$mrca_latency
  df$age_mrca <- ifelse(df$age_mrca > df$age_at_tumour_collection,
                        df$age_at_tumour_collection, df$age_mrca)
  df$age_mrca <- ifelse(df$age_mrca < 0, 0, df$age_mrca)
  df$age_relapse <- ifelse(df$relapse_ind == "relapsed",
                           as.numeric(df$age_at_tumour_collection) +
                             as.numeric(df$donor_relapse_interval) / 365.25,
                           NA)
  df$facet_group <- factor(cut(1:nrow(df), breaks = 4, labels = paste("Group", 1:4)))
  df <- df[!is.na(df$gs_grade), ]
  ggplot(df, aes(y = record_id)) +
    geom_point(aes(x = age_mrca), pch = 21, alpha = .8, fill = icgc["blue"], size = 1) +
    geom_segment(aes(x = 0, xend = age_mrca, yend = record_id), col = icgc["blue"], alpha = 0.4) +
    geom_segment(aes(x = age_mrca, xend = age_at_tumour_collection, yend = record_id), col = icgc["red"], alpha = 0.4) +
    geom_point(aes(x = age_at_tumour_collection), pch = 21, alpha = .8, fill = icgc["red"], size = 1) +
    geom_point(aes(x = 0, fill = gs_grade), pch = 22, size = 1) +
    geom_point(aes(x = age_relapse), pch = 4, col = "red", alpha = .6, size = .6) +
    scale_fill_manual(values = gleason_colours) +
    facet_wrap(~ facet_group, ncol = 4, scales = "free_y") +
    theme(strip.text = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank()) +
    labs(x = "Age (Years)", y = "PPCG Patient", fill = "Grade Group")
}

## boxplot_epicmit() — from scripts/01_landscape/04e_proliferation_markers.R (SF6)
boxplot_epicmit <- function(data, xvar, yvar, xlab, ylab) {
  data %>% dplyr::filter(!is.na(!!sym(xvar))) %>%
    ggplot(aes_string(x = xvar, y = yvar)) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.1) +
    geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) +
    ggpubr::stat_compare_means(
      comparisons = list(c("1", "2"), c("2", "3"), c("3", "4+"), c("2", "4+")),
      method = "wilcox.test", size = 2) +
    labs(x = xlab, y = ylab) +
    theme(legend.position = "none")
}

## boxplot_rates() — from scripts/04_progression/02a_rate_evolution.R (SF14)
boxplot_rates <- function(df, y_str, ylab) {
  trend_stats <- df %>%
    dplyr::filter(!is.na(gs_grade)) %>%
    mutate(gs_ordered = factor(gs_grade, ordered = TRUE)) %>%
    group_by(trajectory) %>%
    summarise(p_val = DescTools::JonckheereTerpstraTest(
      x = .data[[y_str]], g = gs_ordered, alternative = "increasing")$p.value,
      .groups = 'drop') %>%
    mutate(p_label = ifelse(p_val < 0.001, "P < 0.001", sprintf("P = %.3f", p_val)),
           y_pos = max(df[[y_str]], na.rm = TRUE) * 1.5)
  box_dodge   <- position_dodge(width = 0.75)
  point_dodge <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)
  p <- ggplot(df %>% dplyr::filter(!is.na(gs_grade)),
              aes_string(y = y_str, fill = "gs_grade", x = "trajectory")) +
    geom_point(pch = 21, alpha = 0.2, position = point_dodge) +
    geom_boxplot(outlier.shape = NA, position = box_dodge, alpha = 0.5) +
    scale_y_log10() +
    scale_fill_manual(values = gleason_colours) +
    geom_text(data = trend_stats,
              aes(x = trajectory, y = y_pos, label = p_label),
              inherit.aes = FALSE, size = 2, fontface = "italic") +
    labs(x = "Trajectory", y = ylab, fill = "GS Grade") +
    scale_x_discrete(labels = c("Canonical", "Alternative", "Gain-enriched")) +
    theme(legend.position = "none")
  axes2lemon(p, bt = "both", lt = "both")
}

## plot_model() — from scripts/04_progression/02c_analysis_simulations.R (Fig4e / EXDF8)
plot_model_simulation <- function(simulations_df, model_name, gleason_colours, met_summary_data, met_progression_data) {
  names(gleason_colours) <- paste0("Grade Group ", names(gleason_colours))
  simulations_df$gg <- paste0("Grade Group ", simulations_df$gg)
  simulations_df %>% dplyr::filter(model == model_name) %>%
    dplyr::group_by(gg, year) %>%
    dplyr::summarise(first_quartile = quantile(poete, 0.25),
                     median = median(poete),
                     mean = mean(poete),
                     third_quartile = quantile(poete, 0.75),
                     .groups = "drop") %>%
    ggplot() +
    geom_step(aes(x = year, y = median, color = gg), size = .5) +
    geom_ribbon(aes(x = year, ymin = first_quartile, ymax = third_quartile, fill = gg), alpha = 0.1) +
    geom_line(data = simulations_df %>% dplyr::filter(model == model_name) %>%
                group_by(year, sample, gg) %>% summarise(poete = mean(poete), .groups = "drop"),
              aes(x = year, y = poete, group = sample, color = gg), size = 0.1, alpha = 0.1) +
    geom_hline(yintercept = median(met_progression_data$total_poete),
               color = "firebrick", linetype = "dashed", size = 0.5, alpha = 1) +
    labs(x = "Years", y = "Degree of trajectory progression") +
    scale_color_manual(values = gleason_colours, name = "Grade Group") +
    scale_fill_manual(values  = gleason_colours, name = "Grade Group") +
    xlim(-5, 15) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
    facet_wrap(~gg, scales = "free_y", nrow = 1, ncol = 4) +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    theme(legend.position = "none")
}

## plot_km() — from scripts/05_clinical_utility/02_hmf_treatment_response.R (EXDF10b-g)
plot_km <- function(data, km_fit, title, xlab, ylab, filename) {
  ggsurv <- ggsurvplot(km_fit, data = data,
    risk.table = TRUE, pval = TRUE, conf.int = FALSE,
    palette = c("darkolivegreen4", "grey50"),
    legend.title = "Evolutionary Subtype",
    title = title, xlab = xlab, ylab = ylab,
    legend = c(0.7, 0.8),
    legend.labs = c("ARPI", "Taxane"),
    ggtheme = theme_classic() + theme(
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key        = element_rect(fill = "transparent", color = NA),
      legend.margin     = margin(0, 0, 0, 0)
    ),
    size = 0.4, censor.size = 2, pval.size = 2.5, fontsize = 2.5,
    font.main = 7, font.x = 7, font.y = 7, font.tickslab = 6, font.legend = 6,
    risk.table.height = 0.2, risk.table.y.text = FALSE,
    risk.table.title = "At risk",
    tables.theme = theme_cleantable() + theme(
      plot.title = element_text(size = 7),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(0, 5.5, 0, 5.5))
  )
  pdf(fp(filename), width = 55/25.3, height = 60/25.3)
  print(ggsurv)
  dev.off()
}

## ggsurvplot wrapper for the Fig5/EXDF9 KM curves (shared styling)
plot_km_evorisk <- function(data, legend_xy = c(0.2, 0.5), filename) {
  km_fit <- survfit(Surv(time2relapse, relapse_ind) ~ evorisk, data = data)
  pdf(fp(filename), height = 50/25.3, width = 60/25.3)
  print(ggsurvplot(km_fit, data = data,
    risk.table = TRUE, pval = TRUE, conf.int = FALSE,
    xlab = "Time to metastasis (days)",
    ylab = "Metastasis-free survival",
    palette = c("#DC3220", "#005AB5"),
    legend = legend_xy,
    legend.title = "Evolutionary risk",
    legend.labs = c("High", "Low"),
    ggtheme = theme_classic() + theme(
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key        = element_rect(fill = "transparent", color = NA),
      legend.margin     = margin(0, 0, 0, 0)
    ),
    size = 0.4, censor.size = 2, pval.size = 2.5, fontsize = 2.5,
    font.main = 7, font.x = 7, font.y = 7, font.tickslab = 6, font.legend = 6,
    risk.table.height = 0.2, risk.table.y.text = FALSE,
    risk.table.title = "At risk",
    tables.theme = theme_cleantable() + theme(
      plot.title = element_text(size = 7),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(0, 5.5, 0, 5.5))
  ))
  dev.off()
}

## Palettes used inside specific panel sections (inlined verbatim)
col_genome <- "#E76A85"
cols_psa <- c("<10" = "#F0F4F8", ">10" = "#90CAF9", ">20" = "#0D47A1")
cols_stage <- c("T2" = "#F1F8E9", "T3" = "#AED581", "T4" = "#33691E")
cols_age <- c("<55" = "#F3E5F5", "55-65" = "#BA68C8", ">65" = "#4A148C")
cols_mfs <- c("FALSE" = "#F5F5F5", "TRUE" = "#B71C1C")
cols_anc <- c("EUR" = "#F0F4F8", "AFR" = "#E64B35", "AMR" = "#F39B7F",
              "EAS" = "#00A087", "SAS" = "#8491B4")

region_colours <- c("Island" = "#8dd3c7", "Shore" = "#bebada",
                    "Shelf"  = "#fb8072", "Open Sea" = "#80b1d3")
meth_colors <- c("hypo" = "#77AADD", "hyper" = "#EE8866")

gene_colours <- c(
  "TP53" = "#E41A1C", "RB1" = "#377EB8", "PTEN" = "#4DAF4A",
  "CHD1" = "#984EA3", "BRCA1" = "#FF7F00",
  "ERG" = "#A65628", "FOXP1" = "#E7298A", "LRP1B" = "#F781BF",
  "ZNF292" = "#1F78B4", "NKX3-1" = "#B2DF8A",
  "MYC" = "#800026", "PIK3CA" = "#FEB24C", "MDM4" = "#F03B20",
  "CCND1" = "#377EB8", "FOXA1" = "#1B9E77"
)

imf_colours <- c(
  "IMF1" = "#B0B0B0", "IMF2" = "#CC79A7", "IMF3" = "#0072B2",
  "IMF4" = "#DFDFDF", "IMF5" = "#D55E00", "IMF6" = "#009E73",
  "IMF7" = "#E5D8BD", "IMF8" = "#E69F00"
)


# =============================================================================
# 5. MAIN FIGURES (Fig*)  -----------------------------------------------------
# =============================================================================

# ---- Fig1a  Subclonality landscape ------------------------------------------
plot_data <- sd_tsv("Fig1a")
plot_data <- plot_data %>%
  arrange(group, desc(mut_ith)) %>%
  mutate(smp = factor(smp, levels = unique(smp)))

p_snvith <- ggplot(plot_data) +
  geom_col(aes(x = smp, y = 1), fill = "#FAFAFA", width = 1) +
  geom_col(aes(x = smp, y = mut_ith), fill = col_genome, width = 1) +
  facet_grid(. ~ group, scales = "free_x", space = "free") +
  labs(y = "Prop.\nSubclonal") +
  theme_void() +
  theme(strip.text = element_text(face = "bold", size = 12),
        axis.title.y = element_text(angle = 90, size = 9, margin = margin(r = 5)),
        legend.position = "none")

p_nsub <- ggplot(plot_data) +
  geom_col(aes(x = smp, y = nsubclones), fill = col_genome, width = 1) +
  facet_grid(. ~ group, scales = "free_x", space = "free") +
  labs(y = "#\nClones") +
  theme_void() +
  theme(strip.text = element_blank(),
        axis.title.y = element_text(angle = 90, size = 9, margin = margin(r = 5)),
        legend.position = "none")

p_clin1 <- plot_clinical_strip(plot_data, "psa_simplified", "PSA",      cols_psa)
p_clin2 <- plot_clinical_strip(plot_data, "t_stage",        "Stage",    cols_stage)
p_clin3 <- plot_clinical_strip(plot_data, "age",            "Age",      cols_age)
p_clin4 <- plot_clinical_strip(plot_data, "MFS",            "MFS",      cols_mfs)
p_clin5 <- plot_clinical_strip(plot_data, "ancestry",       "Ancestry", cols_anc)

final_plot <- p_snvith / p_nsub / p_clin1 / p_clin2 / p_clin3 / p_clin4 / p_clin5 +
  plot_layout(heights = c(4, 2, 0.5, 0.5, 0.5, 0.5, 0.5), guides = "collect") &
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

save_ggplot(final_plot, fp("Fig1a_subclonality_landscape"), w = 180, h = 65)
message("Saved: Fig1a")

# ---- Fig1c  Timelines: MRCA / CNA / WGD raincloud ---------------------------
timeline_data <- sd_tsv("Fig1c")
p <- plot_timeline_raincloud_gs(timeline_data)
p <- axes2lemon(p)
ggsave(fp("Fig1c_timelines_mrca_scna_wgd.pdf"), plot = p,
       width = 120/25.3, height = 60/25.3)
ggsave(fp("Fig1c_timelines_mrca_scna_wgd.png"), plot = p,
       width = 120/25.3, height = 80/25.3)
message("Saved: Fig1c")

# ---- Fig2b/c/d  Plackett-Luce trajectory ordering plots ---------------------
## The originals use a custom geom_stripes + violin composition rendered through
## cairo_pdf with the Roboto Condensed font. The reproduction renders an
## equivalent plot from the saved `tdata` tibble; if Roboto Condensed is not
## installed locally the default family is used (cosmetic difference only).
geom_stripes <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", ..., show.legend = NA,
                         inherit.aes = TRUE) {
  ggplot2::layer(data = data, mapping = mapping, stat = stat,
                 geom = GeomStripes, position = position,
                 show.legend = show.legend, inherit.aes = inherit.aes,
                 params = list(...))
}
GeomStripes <- ggplot2::ggproto("GeomStripes", ggplot2::Geom,
  required_aes = c("y"),
  default_aes = ggplot2::aes(ymin = -Inf, ymax = Inf,
                             odd = "#99999911", even = "#FFFFFF00",
                             alpha = NA, colour = "black",
                             linetype = "solid", size = NA),
  draw_key = ggplot2::draw_key_rect,
  draw_panel = function(data, panel_params, coord) {
    ggplot2::GeomRect$draw_panel(
      data %>%
        dplyr::mutate(x = round(.data$x),
                      xmin = .data$x - 0.5,
                      xmax = .data$x + 0.5) %>%
        dplyr::select(.data$ymin, .data$ymax, .data$xmin, .data$xmax,
                      .data$odd, .data$even, .data$alpha, .data$colour,
                      .data$linetype, .data$size) %>%
        unique() %>%
        dplyr::arrange(.data$xmin) %>%
        dplyr::mutate(.n = dplyr::row_number(),
                      fill = dplyr::if_else(.data$.n %% 2L == 1L,
                                            true = .data$odd, false = .data$even)) %>%
        dplyr::select(-.data$.n, -.data$odd, -.data$even),
      panel_params, coord)
  }
)

render_pl_trajectory <- function(panel) {
  tdata <- sd_rds(panel)
  tdata$plotpos <- -log10(tdata$value)
  tdata$plotpos <- (tdata$plotpos - min(tdata$plotpos)) + 0.11
  WGD_linepos <- -log10(median(tdata$value))
  midpoint <- median(tdata$plotpos)

  tdata2 <- tdata %>% group_by(ID) %>%
    summarise(mvalue = median(plotpos), max = max(plotpos), min = min(plotpos)) %>%
    mutate(plotpos = if_else(mvalue > midpoint, min, max), cna = 'LOH') %>%
    mutate(ID2 = str_remove(ID, ' \\[.*')) %>%
    mutate(ID2 = str_remove(ID2, '^.*: '))
  tdata2.1 <- tdata2 %>% filter(mvalue > midpoint)
  tdata2.2 <- tdata2 %>% filter(mvalue <= midpoint)
  tdata2.3 <- tdata %>% group_by(ID, freq) %>% dplyr::slice(1) %>% ungroup()
  tdata2.3$freq <- tdata2.3$freq / 1000

  event_type_colours <- c(`Mut` = "#DD8AB9", `Gain` = "#68BFA3", `LOH` = "#8E9FCA",
                          `HD` = "#F18A62", `WGD` = "#664499", `MRCA` = "#555555")
  event_type_colours <- event_type_colours[names(event_type_colours) %in% tdata$cna]
  event_type_labels <- c(`Mut` = "Mutation", `Gain` = "Gain", `LOH` = "Loss",
                         `HD` = "Homozygous Deletion", `WGD` = "WGD", `MRCA` = "MRCA")
  event_type_labels <- event_type_labels[names(event_type_labels) %in% tdata$cna]

  text_size <- 1.5
  p2 <- tdata %>%
    ggplot(aes(x = ID, y = plotpos, fill = cna)) +
    geom_stripes(odd = "#11111106", even = "#00000000") +
    geom_hline(yintercept = WGD_linepos, linetype = 'dashed',
               col = rgb(1, 1, 1, 0), linewidth = 0.25) +
    geom_violin(color = NA) +
    geom_hline(yintercept = 0, linetype = 'solid', col = 'gray50', linewidth = 0.1) +
    geom_hline(yintercept = 0.05, linetype = 'dotted', col = 'gray50', linewidth = 0.1) +
    geom_hline(yintercept = 0.1, linetype = 'dotted', col = 'gray25', linewidth = 0.1) +
    geom_bar(data = tdata2.3, aes(y = freq), width = 2/3, color = 'black',
             linewidth = 0.001, stat = "identity") +
    geom_text(data = tdata2.3,
              aes(y = -0.02, label = paste0(format(round(freq * 1000, 1), nsmall = 1), "%")),
              hjust = 1, size = text_size, family = 'Roboto Condensed') +
    stat_summary(fun = "median", geom = "crossbar",
                 width = 0.4, linewidth = 0.18,
                 colour = "black", show.legend = FALSE) +
    geom_text(data = tdata2.1, aes(label = ID2),
              vjust = 0.5, hjust = 1, nudge_y = -0.05,
              size = text_size, family = 'Roboto Condensed') +
    geom_text(data = tdata2.2, aes(label = ID2),
              vjust = 0.5, hjust = 0, nudge_y = 0.05,
              size = text_size, family = 'Roboto Condensed') +
    labs(x = "", y = "", fill = "Event", size = 'Frequency (%)') +
    scale_fill_manual(breaks = names(event_type_colours),
                      values = event_type_colours,
                      labels = event_type_labels) +
    scale_size_binned(n.breaks = 8) +
    guides(fill = guide_legend(byrow = TRUE,
                               override.aes = list(size = 2.5, pch = 21)),
           size = guide_legend(override.aes = list(fill = 'gray50', col = 'black'))) +
    hrbrthemes::theme_ipsum_rc(base_family = "Roboto Condensed",
                               ticks = FALSE, grid = FALSE,
                               axis_title_just = 'm', axis_title_size = 9) +
    theme(axis.text.y = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "null")) +
    scale_y_discrete(expand = c(0.1, 0.1)) +
    coord_flip(clip = 'off')

  # Render via cairo_pdf to mirror the original output
  tryCatch(
    {
      cairo_pdf(fp(paste0("Fig2_", substring(panel, 5, 5), "_PLPlotSlim.pdf")),
                width = 6, height = 13, family = "Roboto Condensed")
      print(p2); dev.off()
    },
    error = function(e) {
      ggsave(fp(paste0("Fig2_", substring(panel, 5, 5), "_PLPlotSlim.pdf")),
             p2, width = 6, height = 13)
    }
  )
  ggsave(fp(paste0("Fig2_", substring(panel, 5, 5), "_PLPlotSlim.png")),
         p2, width = 6, height = 13, device = "png")
}

# Loop maps panel letter b/c/d → original trajectory 1/2/3
for (panel in c("Fig2b", "Fig2c", "Fig2d")) {
  render_pl_trajectory(panel)
  message("Saved: ", panel)
}

# ---- Fig2e/f/g  Grade-Group donuts per trajectory ---------------------------
eopc_colors <- c("Early-onset" = "#5E6C85", "Late-onset" = "grey75")

tj1 <- sd_tsv("Fig2e")
gg_tj1 <- plot_donut(tj1, "gs_group", color_palette = gleason_colours, cohort_label = "Canonical")
ggsave(fp("Fig2e_pie_grade_canonical.pdf"), gg_tj1,
       height = 30/25.4, width = 30/25.4, dpi = 300)
tj1_eopc <- sd_tsv("Fig2e_eopc")
gg_eopc_tj1 <- plot_donut(tj1_eopc, "eopc", color_palette = eopc_colors, cohort_label = "Canonical")
ggsave(fp("Fig2e_pie_eopc_canonical.pdf"), gg_eopc_tj1,
       height = 30/25.4, width = 30/25.4, dpi = 300)
message("Saved: Fig2e")

tj2 <- sd_tsv("Fig2f")
gg_tj2 <- plot_donut(tj2, "gs_group", color_palette = gleason_colours, cohort_label = "Alternative")
ggsave(fp("Fig2f_pie_grade_alternative.pdf"), gg_tj2,
       height = 30/25.4, width = 30/25.4, dpi = 300)
tj2_eopc <- sd_tsv("Fig2f_eopc")
gg_eopc_tj2 <- plot_donut(tj2_eopc, "eopc", color_palette = eopc_colors, cohort_label = "Alternative")
ggsave(fp("Fig2f_pie_eopc_alternative.pdf"), gg_eopc_tj2,
       height = 30/25.4, width = 30/25.4, dpi = 300)
message("Saved: Fig2f")

tj3 <- sd_tsv("Fig2g")
gg_tj3 <- plot_donut(tj3, "gs_group", color_palette = gleason_colours, cohort_label = "Gain-enriched")
ggsave(fp("Fig2g_pie_grade_gain_enriched.pdf"), gg_tj3,
       height = 30/25.4, width = 30/25.4, dpi = 300)
tj3_eopc <- sd_tsv("Fig2g_eopc")
gg_eopc_tj3 <- plot_donut(tj3_eopc, "eopc", color_palette = eopc_colors, cohort_label = "Gain-enriched")
ggsave(fp("Fig2g_pie_eopc_gain_enriched.pdf"), gg_eopc_tj3,
       height = 30/25.4, width = 30/25.4, dpi = 300)
message("Saved: Fig2g")

# ---- Fig3a  Sankey: MRCA vs final trajectory --------------------------------
df_sankey <- sd_tsv("Fig3a") %>%
  mutate(mrca_trajectory  = factor(mrca_trajectory,
                                   levels = c("Canonical", "Alternative", "Gain-enriched")),
         final_trajectory = factor(final_trajectory,
                                   levels = c("Canonical", "Alternative", "Gain-enriched")))
traj_colors <- c("Canonical" = "#0072B2",
                 "Alternative" = "#1B9E77",
                 "Gain-enriched" = "#D73027")
p <- ggplot(df_sankey, aes(y = n, axis1 = mrca_trajectory, axis2 = final_trajectory)) +
  geom_alluvium(aes(fill = mrca_trajectory), width = 0.15, alpha = 0.2, knot.pos = 0.4) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.05, alpha = 0.8,
               color = "white", linewidth = 0.2) +
  scale_fill_manual(values = traj_colors) +
  theme_void() + theme(legend.position = "none")
ggsave(fp("Fig3a_sankey_trajectories.png"), p,
       height = 40/25.4, width = 40/25.4, dpi = 300)
ggsave(fp("Fig3a_sankey_trajectories.pdf"), p,
       height = 50/25.4, width = 60/25.4, dpi = 300)
message("Saved: Fig3a")

# ---- Fig3b  Stacked IMF activities per trajectory ---------------------------
imf_activities <- sd_tsv("Fig3b")

get_sample_order <- function(df, traj_name, target_imfs) {
  df %>% dplyr::filter(trajectory == traj_name,
                       imf_cluster %in% target_imfs) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(target_sum = sum(activity, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(target_sum) %>% dplyr::pull(sample)
}

canon_order <- get_sample_order(imf_activities, "Canonical",     c("IMF2", "IMF3"))
alt_order   <- get_sample_order(imf_activities, "Alternative",   "IMF6")
gain_order  <- get_sample_order(imf_activities, "Gain-enriched", c("IMF5", "IMF8"))

imf_names <- sort(unique(imf_activities$imf_cluster))
global_imf_fill <- imf_activities %>%
  dplyr::group_by(imf_cluster) %>%
  dplyr::summarise(total_act = sum(activity, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(total_act)) %>% dplyr::pull(imf_cluster)
imf_activities$imf_cluster <- factor(imf_activities$imf_cluster, levels = global_imf_fill)

p_tj1 <- imf_activities %>% dplyr::filter(trajectory == "Canonical") %>%
  dplyr::mutate(sample = factor(sample, levels = canon_order)) %>%
  ggplot(aes(y = sample, x = activity, fill = imf_cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = imf_colours, breaks = imf_names) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL, x = "IMF Activity") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 7, color = "black"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.line.y  = element_line(color = "black")) +
  coord_flip() + guides(fill = guide_legend(nrow = 1))

p_tj2 <- imf_activities %>% dplyr::filter(trajectory == "Alternative") %>%
  dplyr::mutate(sample = factor(sample, levels = alt_order)) %>%
  ggplot(aes(y = sample, x = activity, fill = imf_cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = imf_colours, breaks = imf_names) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL, x = NULL) +
  theme(legend.position = "none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x  = element_blank()) +
  coord_flip()

p_tj3 <- imf_activities %>% dplyr::filter(trajectory == "Gain-enriched") %>%
  dplyr::mutate(sample = factor(sample, levels = gain_order)) %>%
  ggplot(aes(y = sample, x = activity, fill = imf_cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = imf_colours, breaks = imf_names) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none", legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x  = element_blank()) +
  coord_flip()

samples_prop <- imf_activities %>%
  dplyr::select(trajectory, sample) %>% dplyr::distinct() %>%
  dplyr::group_by(trajectory) %>% dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::pull(prop, name = trajectory)

p <- cowplot::plot_grid(
  p_tj1, p_tj2, p_tj3,
  ncol = 3, align = "h",
  rel_widths = c(samples_prop["Canonical"],
                 samples_prop["Alternative"],
                 samples_prop["Gain-enriched"])
)
ggsave(fp("Fig3b_Stacked_IMF_Activities.pdf"), p,
       width = 125/25.4, height = 70/25.4, units = "in", device = "pdf")
message("Saved: Fig3b")

# ---- Fig3c  LOH segment size vs replication time ----------------------------
loh_bps_df <- sd_tsv("Fig3c")
# quantiles needed for the dashed verticals; computed from the data
rt_q  <- quantile(loh_bps_df$reptime_score, c(0.33, 0.66), na.rm = TRUE)
p <- ggplot(loh_bps_df, aes(x = reptime_score, y = size)) +
  geom_density_2d(alpha = 0.8, color = "black") +
  geom_point(size = 0.5, aes(color = gene, alpha = !is.na(gene))) +
  geom_vline(xintercept = rt_q, linetype = "dashed", color = "black") +
  facet_wrap(~trajectory) +
  scale_y_log10() +
  scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 0.05), guide = "none") +
  labs(x = "Replication timing score", y = "Size of LOH segment", color = "Gene") +
  scale_color_manual(values = gene_colours, na.value = "grey") +
  theme(legend.position = "none")
ggsave(fp("Fig3c_LOH_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)
message("Saved: Fig3c")

# ---- Fig3d  HD segment size vs replication time -----------------------------
hd_bps_df <- sd_tsv("Fig3d")
rt_q <- quantile(hd_bps_df$reptime_score, c(0.33, 0.66), na.rm = TRUE)
p <- ggplot(hd_bps_df, aes(x = reptime_score, y = size)) +
  geom_density_2d(alpha = 0.8, color = "black") +
  geom_point(size = 0.5, aes(color = gene, alpha = !is.na(gene))) +
  geom_vline(xintercept = rt_q, linetype = "dashed", color = "black") +
  facet_wrap(~trajectory) + scale_y_log10() +
  scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 0.05), guide = "none") +
  labs(x = "Replication timing score", y = "Size of HD segment", color = "Gene") +
  scale_color_manual(values = gene_colours, na.value = "grey") +
  theme(legend.position = "none")
ggsave(fp("Fig3d_HD_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)
message("Saved: Fig3d")

# ---- Fig3e  Gain segment size vs replication time ---------------------------
gain_bps_df <- sd_tsv("Fig3e")
rt_q <- quantile(gain_bps_df$reptime_score, c(0.33, 0.66), na.rm = TRUE)
p <- ggplot(gain_bps_df, aes(x = reptime_score, y = size)) +
  geom_density_2d(alpha = 0.8, color = "black") +
  geom_point(size = 0.5, aes(color = gene, alpha = !is.na(gene))) +
  geom_vline(xintercept = rt_q, linetype = "dashed", color = "black") +
  facet_wrap(~trajectory) + scale_y_log10() +
  scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 0.05), guide = "none") +
  labs(x = "Replication timing score", y = "Size of GAIN segment", color = "Gene") +
  scale_color_manual(values = gene_colours, na.value = "grey") +
  theme(legend.position = "none")
ggsave(fp("Fig3e_GAIN_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)
message("Saved: Fig3e")

# ---- Fig3f  Ancestry donuts (PPCG EUR + AFR) --------------------------------
ancestry_df <- sd_tsv("Fig3f")
ppcg_eur <- plot_trajectory_donut(
  table(ancestry_df$trajectory[ancestry_df$predicted_ancestry == "EUR"]),
  cohort_label = "PPCG European"
)
ppcg_afr <- plot_trajectory_donut(
  table(ancestry_df$trajectory[ancestry_df$predicted_ancestry == "AFR"]),
  cohort_label = "PPCG African"
)
ggsave(fp("Fig3f_pie_trajectory_EUR.pdf"), ppcg_eur,
       height = 30/25.4, width = 30/25.4, dpi = 300)
ggsave(fp("Fig3f_pie_trajectory_AFR.pdf"), ppcg_afr,
       height = 30/25.4, width = 30/25.4, dpi = 300)
message("Saved: Fig3f (EUR + AFR)")

# ---- Fig3g  HEROIC AFR validation donut -------------------------------------
heroic_tj <- sd_tsv("Fig3g")
heroic_donut <- plot_trajectory_donut(
  as.table(setNames(heroic_tj$n, heroic_tj$trajectory)),
  cohort_label = "HEROIC African"
)
ggsave(fp("Fig3g_pie_trajectory_HEROIC_AFR.pdf"), heroic_donut,
       height = 30/25.4, width = 30/25.4, dpi = 300)
message("Saved: Fig3g")

# ---- Fig4a  Trajectory distribution donut across cohorts --------------------
all_tjs <- sd_tsv("Fig4a")
all_tjs$cohort <- factor(all_tjs$cohort, levels = c("PPCG", "CombiMets", "HMF"))
plot_data <- all_tjs %>%
  group_by(cohort, trajectory) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cohort) %>%
  mutate(total = sum(count), prop = count / total)
center_data <- plot_data %>%
  dplyr::select(cohort, total) %>% distinct() %>%
  mutate(label_text = paste0("N=\n", total))
p <- ggplot(plot_data, aes(x = 2, y = prop, fill = trajectory)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
  coord_polar(theta = "y", start = 0) + xlim(0.5, 2.5) +
  facet_wrap(~ cohort) +
  geom_text(data = center_data, aes(x = 0.5, y = 0, label = label_text),
            inherit.aes = FALSE, size = 2.5, fontface = "bold", color = "#333333") +
  scale_fill_manual(values = trajectory_colors) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 7, face = "bold", margin = margin(b = 2)),
        plot.margin = margin(2, 2, 2, 2, "mm"))
ggsave(fp("Fig4a_trajectory_distribution_cohort.pdf"), plot = p,
       width = 60/25.3, height = 50/25.4, units = "in", dpi = 300)
message("Saved: Fig4a")

# ---- Fig4b  Trajectory progression boxplot by Gleason group -----------------
all_tjs <- sd_tsv("Fig4b")
trend_stats <- all_tjs %>%
  dplyr::filter(!is.na(gs_grade)) %>%
  mutate(gs_ordered = factor(gs_grade, ordered = TRUE)) %>%
  group_by(trajectory) %>%
  summarise(p_val = DescTools::JonckheereTerpstraTest(
    x = total_poete, g = gs_ordered, alternative = "increasing")$p.value,
    .groups = "drop") %>%
  mutate(p_label = ifelse(p_val < 0.001, "P < 0.001", sprintf("P = %.3f", p_val)),
         x_pos = max(all_tjs$total_poete, na.rm = TRUE) * 0.9)

box_dodge   <- position_dodge(width = 0.75)
point_dodge <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)
gc_palette <- c(gleason_colours, "Metastasis" = "#67000D")

p <- all_tjs %>% dplyr::filter(!is.na(gs_grade)) %>%
  ggplot(aes(y = trajectory, x = total_poete, fill = gs_grade)) +
  geom_point(pch = 21, alpha = 0.2, position = point_dodge) +
  geom_boxplot(outlier.shape = NA, position = box_dodge, alpha = 0.5) +
  geom_text(data = trend_stats,
            aes(x = x_pos, y = trajectory, label = p_label),
            inherit.aes = FALSE, size = 2, angle = 270, fontface = "italic") +
  scale_fill_manual(values = gc_palette) +
  theme(legend.position = "top") +
  labs(x = "Relative degree of progression along trajectory", y = "Trajectory")
p <- axes2lemon(p, lt = "both", bt = "both")
ggsave(fp("Fig4b_trajectory_progression_by_gs.pdf"), p,
       width = 90/25.3, height = 80/25.4, dpi = 300)
message("Saved: Fig4b")

# ---- Fig4d  Example simulation case (one sim per model) ---------------------
example_df <- sd_tsv("Fig4d")
sid <- unique(example_df$sample)[1]
p <- ggplot(example_df, aes(x = year, y = total_events, color = model)) +
  geom_step(size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = paste(sid, " - Simulation 1"),
       x = "Years from Diagnosis", y = "Total Accumulated Events") +
  theme(legend.position = "top",
        panel.grid.major.y = element_line(color = "grey90", linetype = "dotted")) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1"), name = "Model")
p <- axes2lemon(p, lt = "both", bt = "both")
ggsave(fp("Fig4d_trajectory_drivers_simulation_example_case_sim1.pdf"),
       plot = p, width = 100/25.3, height = 50/25.4)
message("Saved: Fig4d")

# ---- Fig4e  Simulated progression by Grade Group (uniform model) ------------
simulations_df_uniform <- sd_tsv("Fig4e")
# need per-trajectory metastatic progression data for the dashed line;
# pull from Fig4b source data (the merge result, including HMF rows)
met_progression_data <- sd_tsv("Fig4b") %>% dplyr::filter(gs_grade == "Metastasis")
met_summary_data <- met_progression_data %>% group_by(trajectory) %>%
  dplyr::summarise(met_first_quartile = quantile(total_poete, .25),
                   met_median = median(total_poete),
                   met_third_quartile = quantile(total_poete, .75),
                   .groups = "drop")
names(gleason_colours) <- c("1", "2", "3", "4-5")
p <- plot_model_simulation(simulations_df_uniform, "uniform", gleason_colours,
                           met_summary_data, met_progression_data)
ggsave(fp("Fig4e_simulated_progression_gg.pdf"), p,
       width = 160/25.3, height = 50/25.3, dpi = 300)
message("Saved: Fig4e")

# ---- Fig5a  Cox forest plot, MFS by subclonal dynamics ----------------------
met_relapse <- sd_tsv("Fig5a")
met_relapse <- met_relapse %>%
  mutate(high_pga = fct_relevel(high_pga, "Low"),
         advanced = fct_relevel(advanced, "Low"))
model <- coxph(Surv(time2relapse, relapse_ind) ~
                 ggroup + psa_cat + t_stage_simplified + age +
                 trajectory + mut_ith + nsubclones + advanced +
                 high_pga, data = met_relapse)
pdf(fp("Fig5a_mfs_by_subclonal_dynamics.pdf"), height = 100/25.3, width = 120/25.3)
print(ggforest(model, data = as.data.frame(met_relapse))); dev.off()
message("Saved: Fig5a")

# ---- Fig5c/d/e  Kaplan-Meier: evolutionary risk per Grade Group -------------
met_relapse_c <- sd_tsv("Fig5c")
plot_km_evorisk(met_relapse_c, legend_xy = c(0.2, 0.5),
                filename = "Fig5c_km_evorisk_gg2.pdf")
message("Saved: Fig5c")

met_relapse_d <- sd_tsv("Fig5d") 
plot_km_evorisk(met_relapse_d, legend_xy = c(0.65, 0.85),
                filename = "Fig5d_km_evorisk_gg3.pdf")
message("Saved: Fig5d")

met_relapse_e <- sd_tsv("Fig5e") 
plot_km_evorisk(met_relapse_e, legend_xy = c(0.65, 0.85),
                filename = "Fig5e_km_evorisk_gg4.pdf")
message("Saved: Fig5e")


# =============================================================================
# 6. EXTENDED DATA FIGURES (EXDF*)  -------------------------------------------
# =============================================================================

# ---- EXDF1ab  NRPCC vs subclones / min CCF (combined) -----------------------
df <- sd_tsv("EXDF1ab")
p1 <- ggplot(df, aes(x = nrpcc, y = as.character(nsubclones))) +
  geom_boxplot(alpha = .4) +
  geom_vline(xintercept = 10, linetype = "dashed", col = "red") +
  labs(x = "Number of reads per chromosome copy\n(NRPCC)",
       y = "Number of subclones detected")
p2 <- ggplot(df, aes(x = nrpcc, y = min_ccf)) +
  geom_point(alpha = .4) +
  geom_vline(xintercept = 10, linetype = "dashed", col = "red") +
  geom_hline(yintercept = 0.4, linetype = "dashed", col = "blue") +
  ylim(c(0, 1)) +
  labs(x = "NRPCC", y = "Minimum CCF of detected clusters")
p <- (p1 + p2) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10))
save_ggplot(p, fp("EXDF1ab_nrpcc_power_subclonality"), w = 110, h = 65)
message("Saved: EXDF1ab")

# ---- EXDF1c  NRPCC vs subclonal mutation ratio (WCC corrected) --------------
muts_wcc_df <- sd_tsv("EXDF1c")
p3 <- ggplot(muts_wcc_df,
             aes(y = nsubclonal_wcc / (nclonal_wcc + nsubclonal_wcc), x = NRPCC)) +
  geom_point(alpha = .3) +
  labs(y = "Ratio subclonal mutations", x = "NRPCC") +
  theme(legend.position = "bottom") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  ggpubr::stat_cor(method = "spearman", size = 2.5, label.x = 0, label.y = 1.0) +
  xlim(c(0, ceiling(max(muts_wcc_df$NRPCC, na.rm = TRUE)))) +
  ylim(c(0, 1))
ggsave(fp("EXDF1c_ratio_subclonal_mutations_wcc_vs_nrpcc.pdf"), p3,
       width = 50/25.3, height = 50/25.3)
message("Saved: EXDF1c")

# ---- EXDF2a  dN/dS by Grade Group, clonal vs subclonal ----------------------
dnds_df <- sd_tsv("EXDF2a")
dodge <- position_dodge2(width = .75)
p <- dnds_df %>%
  dplyr::mutate(driver_list = factor(driver_list,
                                     levels = c("PPCG", "Armenia 2018", "Wedge 2018"))) %>%
  ggplot(aes(y = group, col = timing)) +
  geom_point(aes(x = mle), size = 3, alpha = .8, position = dodge) +
  geom_linerange(aes(xmin = cilow, xmax = cihigh), position = dodge) +
  geom_vline(xintercept = 1, lty = "dashed") +
  labs(x = "dN / dS", y = "") +
  scale_y_discrete(limits = rev(c("Grade Group 4+", "Grade Group 3",
                                  "Grade Group 2", "Grade Group 1"))) +
  scale_x_log10() + facet_wrap(~driver_list, nrow = 1) +
  scale_colour_manual(values = c("Subclonal" = "#df2a55",
                                 "Clonal"    = "#3979bb",
                                 "All"       = "grey60")) +
  theme(legend.position = "top")
p <- axes2lemon(p, lt = "both")
save_ggplot(p, fp("EXDF2a_dnds_by_gleason_ppcg"), w = 120, h = 60)
message("Saved: EXDF2a")

# ---- EXDF2b  SBS signature fold-change clonal vs subclonal by GS group ------
sigtimerk <- sd_tsv("EXDF2b") %>% dplyr::mutate(grade_group = ifelse(grade_group == "4+", "4-5", grade_group))
sigcount <- dplyr::count(sigtimerk, signature)
sigprev  <- dplyr::arrange(sigcount, n) %>% dplyr::pull(signature) %>% rev() %>% unique()
sigtimerk$signature <- factor(sigtimerk$signature, levels = sigprev)
p <- sigtimerk %>%
  ggplot(aes(x = signature, y = fc_cln, col = grade_group)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = .4) +
  geom_jitter(size = 0.5, alpha = .5, position = position_dodge(width = .8)) +
  geom_boxplot(outlier.shape = NA, alpha = .7) +
  labs(y = "Fold Change", x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = .9, hjust = .9)) +
  scale_y_log10(n.breaks = 15, labels = scales::label_number(accuracy = 0.01)) +
  scale_color_manual(values = gleason_colours) +
  labs(color = "") +
  geom_text(aes(label = n, x = signature, y = 200), size = 2,
            position = position_dodge(width = .8),
            data = dplyr::count(sigtimerk, signature),
            inherit.aes = FALSE, check_overlap = FALSE) +
  theme(legend.position = 'top') +
  scale_x_discrete(guide = "axis_nested")
ggsave(fp("EXDF2b_sig_change_clonal_subclonal_by_grade_group.png"), p,
       width = 125/25.3, height = 75/25.3, dpi = 300)
ggsave(fp("EXDF2b_sig_change_clonal_subclonal_by_grade_group.pdf"), p,
       width = 125/25.3, height = 75/25.3, dpi = 300)
message("Saved: EXDF2b")

# ---- EXDF3a  UpSet plot of differentially methylated CpGs -------------------
dmcs_wide <- sd_tsv("EXDF3a") %>% as.data.frame()
set_names <- setdiff(colnames(dmcs_wide), c("cgid", "type"))
if (!"II_vs_III" %in% colnames(dmcs_wide)) dmcs_wide$II_vs_III <- 0
upset_plot <- upset(
  dmcs_wide, sets = set_names,
  keep.order = TRUE, empty.intersections = "on",
  intersections = list(list("I_vs_II"), list("I_vs_III"),
                       list("I_vs_II", "I_vs_III"), list("II_vs_III")),
  nintersects = NA,
  mainbar.y.label = "Number of Differentially Methylated CpGs",
  sets.x.label = "Total DMCs per comparison",
  text.scale = c(1.3, 1.3, 1.3, 1, 1.5, 1),
  query.legend = "bottom",
  main.bar.color = meth_colors["hypo"],
  queries = list(
    list(query = elements, params = list("type", "Hypo"),
         color = meth_colors["hypo"], active = TRUE, query.name = "Hypomethylated"),
    list(query = elements, params = list("type", "Hyper"),
         color = meth_colors["hyper"], active = TRUE, query.name = "Hypermethylated")
  )
)
pdf(fp("EXDF3a_upset_plot.pdf"), width = 100/25.3, height = 125/25.3)
print(upset_plot); dev.off()
message("Saved: EXDF3a")

# ---- EXDF3b  Methylation differences scatter --------------------------------
methylation_differences <- sd_tsv("EXDF3b")
p <- ggplot(methylation_differences, aes(x = cmp1_mean_diff, y = cmp2_mean_diff)) +
  geom_point(alpha = 0.05) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  ggpubr::stat_cor(method = "pearson") +
  xlab("Mean methylation difference\n (Canonical vs Alternative)") +
  ylab("Mean methylation difference\n (Canonical vs Gain-enriched)") +
  xlim(c(-0.4, 0.4)) + ylim(c(-0.4, 0.4))
p <- axes2lemon(p, lt = "both", bt = "both")
ggsave(fp("EXDF3b_methylation_differences_cmp1_cmp2.pdf"), p,
       width = 60/25.3, height = 60/25.3)
message("Saved: EXDF3b")

# ---- EXDF3c  DMC CGI annotation distribution --------------------------------
donut_data <- sd_tsv("EXDF3c") %>%
  mutate(cpg_island_label = factor(cpg_island_label,
                                   levels = c("Island", "Shore", "Shelf", "Open Sea")))
donut_plot <- ggplot(donut_data, aes(x = 2, y = n, fill = cpg_island_label)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.4) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +
  facet_grid(comparison_name ~ mean.diff.direction) +
  scale_fill_manual(values = region_colours, name = "CGI region") +
  theme_void() +
  theme(text = element_text(family = "Helvetica", size = 8),
        strip.text.x = element_text(face = "bold", size = 9, margin = margin(b = 5, t = 10)),
        strip.text.y = element_text(face = "bold", size = 9, angle = 270, margin = margin(l = 5, r = 10)),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 8),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))
ggsave(fp("EXDF3c_cpg_annotations.pdf"), donut_plot,
       width = 70/25.4, height = 70/25.4)
message("Saved: EXDF3c")

# ---- EXDF3d  DMC distribution across regulatory elements --------------------
reg_elements <- sd_tsv("EXDF3d")
names(meth_colors) <- c("hypo", "hyper")
make_bars <- function(data, meth_colors) {
  ggplot(data %>% mutate(subtype = fct_reorder(subtype, n)),
         aes(x = n, y = subtype, fill = mean.diff.direction)) +
    geom_col() +
    scale_fill_manual(values = meth_colors, name = "Methylation Status") +
    labs(x = "Number of DMCs", y = NULL, title = "") +
    theme(legend.position = "none")
}
p1 <- make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 3",
                                        element == "ATAC_Corces"), meth_colors)
p2 <- make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 3",
                                        element == "Prostate_PMD_Guo"), meth_colors)
p3 <- make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 3",
                                        element == "Pomerantz_Regulatory_Element"), meth_colors)
p4 <- make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 2",
                                        element == "ATAC_Corces"), meth_colors)
p5 <- make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 2",
                                        element == "Prostate_PMD_Guo"), meth_colors)
p6 <- make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 2",
                                        element == "Pomerantz_Regulatory_Element"), meth_colors)
p <- cowplot::plot_grid(p3, p2, p1, p6, p5, p4, nrow = 3, byrow = FALSE,
                        align = "hv",
                        rel_heights = c(10, 3, 3, 10, 3, 3),
                        rel_widths  = c(1, 1, 1, 1, 1, 1))
ggsave(fp("EXDF3d_dmc_regulatory_elements.pdf"), p,
       width = 150/25.3, height = 100/25.3)
message("Saved: EXDF3d")

# ---- EXDF4  GSEA bubble plot across trajectories ----------------------------
gsea_plot_df <- sd_tsv("EXDF4")
p <- ggplot(gsea_plot_df, aes(x = pathway, y = evotype)) +
  geom_point(aes(size = log10padj, fill = NES, color = is_significant), shape = 21) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(1, 5)) +
  scale_color_manual(values = c("significant" = "black", "ns" = "lightgrey")) +
  labs(x = "Pathway", y = "Trajectory", fill = "NES", size = "-log10(padj)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")
ggsave(fp("EXDF4_dea_evotypes.pdf"), plot = p,
       width = 180/25.4, height = 80/25.4)
message("Saved: EXDF4")

# ---- EXDF5  Externally generated (trajectory simulation companion repo) -----
# See README.md for the companion repository link.

# ---- EXDF6  Stacked signature activities per trajectory ---------------------
sigactivities <- sd_tsv("EXDF6")
imf_act_for_order <- sd_tsv("Fig3b")  # for canonical/alternative/gain-enriched orderings
canon_order <- get_sample_order(imf_act_for_order, "Canonical",     c("IMF2", "IMF3"))
alt_order   <- get_sample_order(imf_act_for_order, "Alternative",   "IMF6")
gain_order  <- get_sample_order(imf_act_for_order, "Gain-enriched", c("IMF5", "IMF8"))
traj_orders <- list(Canonical = canon_order,
                    Alternative = alt_order,
                    `Gain-enriched` = gain_order)

base_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#999999")
make_soft_palette <- function(n) colorRampPalette(alpha(base_cols, 0.7))(n)
make_sig_colours <- function(sigactivities) {
  sig_df <- sigactivities %>% group_by(sigtype, signature) %>%
    summarise(total = sum(activity), .groups = "drop") %>%
    arrange(sigtype, desc(total))
  split_sigs <- split(sig_df$signature, sig_df$sigtype)
  lapply(split_sigs, function(sigs) setNames(make_soft_palette(length(sigs)), sigs))
}
sig_colours <- make_sig_colours(sigactivities)

sig_order <- sigactivities %>% group_by(signature) %>%
  summarise(total = sum(activity, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(signature)
sigactivities$signature <- factor(sigactivities$signature, levels = sig_order)

samples_prop <- bind_rows(
  tibble(trajectory = "Canonical",     sample = canon_order),
  tibble(trajectory = "Alternative",   sample = alt_order),
  tibble(trajectory = "Gain-enriched", sample = gain_order)
) %>% count(trajectory, name = "n") %>%
  mutate(prop = n / sum(n)) %>% pull(prop, name = trajectory)

sigtype_levels <- c("SBS", "ID", "SV", "CN")
sig_levels_by_type <- sigactivities %>%
  group_by(sigtype, signature) %>%
  summarise(total = sum(activity), .groups = "drop") %>%
  arrange(sigtype, desc(total)) %>%
  split(.$sigtype) %>%
  lapply(function(df) as.character(df$signature))

make_sig_panel <- function(data, sigtype_name, trajectory_name, sample_order,
                           signature_order, fill_values,
                           show_y = FALSE, show_title = FALSE) {
  panel_df <- data %>%
    dplyr::filter(sigtype == sigtype_name, trajectory == trajectory_name) %>%
    dplyr::mutate(sample = factor(sample, levels = sample_order),
           signature = factor(signature, levels = signature_order)) %>%
    dplyr::select(sample, trajectory, signature, activity) %>%
    tidyr::complete(sample = factor(sample_order, levels = sample_order),
                    signature = factor(signature_order, levels = signature_order),
                    fill = list(activity = 0)) %>%
    mutate(trajectory = trajectory_name)
  ggplot(panel_df, aes(y = sample, x = activity, fill = signature)) +
    geom_bar(stat = "identity", width = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = fill_values, breaks = signature_order, drop = FALSE) +
    labs(x = "Signature activity", y = NULL,
         title = if (show_title) trajectory_name else NULL) +
    coord_flip() +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          plot.title = element_text(size = 8, hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 5),
          legend.key.width  = unit(0.3, "cm"),
          legend.key.height = unit(0.2, "cm"),
          legend.spacing.x = unit(0.1, "cm"),
          legend.spacing.y = unit(0.05, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          legend.key = element_rect(fill = NA, colour = NA),
          axis.text.y = if (show_y) element_text(size = 7, colour = "black") else element_blank(),
          axis.ticks.y = if (show_y) element_line(colour = "black") else element_blank(),
          axis.line.y  = if (show_y) element_line(colour = "black") else element_blank())
}

make_sigtype_row <- function(data, sigtype_name, traj_orders, sig_levels_by_type,
                             sig_colours, samples_prop, show_titles = FALSE) {
  p1 <- make_sig_panel(data, sigtype_name, "Canonical",
                       traj_orders[["Canonical"]],
                       sig_levels_by_type[[sigtype_name]],
                       sig_colours[[sigtype_name]], TRUE, show_titles)
  p2 <- make_sig_panel(data, sigtype_name, "Alternative",
                       traj_orders[["Alternative"]],
                       sig_levels_by_type[[sigtype_name]],
                       sig_colours[[sigtype_name]])
  p3 <- make_sig_panel(data, sigtype_name, "Gain-enriched",
                       traj_orders[["Gain-enriched"]],
                       sig_levels_by_type[[sigtype_name]],
                       sig_colours[[sigtype_name]])
  cowplot::plot_grid(p1, p2, p3, ncol = 3, align = "h",
    rel_widths = c(samples_prop["Canonical"],
                   samples_prop["Alternative"],
                   samples_prop["Gain-enriched"]))
}

sig_rows <- lapply(seq_along(sigtype_levels), function(i) {
  st <- sigtype_levels[i]
  if (!st %in% names(sig_levels_by_type)) return(NULL)
  make_sigtype_row(sigactivities, st, traj_orders, sig_levels_by_type,
                   sig_colours, samples_prop, show_titles = (i == 1))
})
sig_rows <- sig_rows[!sapply(sig_rows, is.null)]
p_sig <- cowplot::plot_grid(plotlist = sig_rows, ncol = 1, align = "v")
ggsave(fp("EXDF6_Stacked_Signature_Activities_by_Trajectory.pdf"),
       p_sig, width = 180/25.4, height = 250/25.4, units = "in", device = cairo_pdf)
message("Saved: EXDF6")

# ---- EXDF7a  SBS signature activities, clonal vs subclonal ------------------
sbs_signatures <- sd_tsv("EXDF7a")
names(trajectory_colours) <- c("Canonical", "Alternative", "Gain-enriched")
p <- sbs_signatures %>%
  mutate(timing = factor(timing, levels = c("overall", "clonal", "subclonal"))) %>%
  group_by(timing, signature) %>%
  dplyr::filter(n_distinct(trajectory) >= 2,
                sd(activity, na.rm = TRUE) > 0,
                max(activity, na.rm = TRUE) > 0) %>%
  ggplot(aes(x = signature, y = activity, fill = trajectory)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  ggpubr::stat_compare_means(aes(group = trajectory),
                             label = "p.signif", method = "kruskal.test",
                             size = 2, hide.ns = TRUE) +
  facet_wrap(~timing) +
  labs(x = "Activity", y = "Signature") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  scale_fill_manual(values = trajectory_colours)
p <- axes2lemon(p, lt = "both", bt = "both")
ggsave(fp("EXDF7a_sbs_signature_activities_clonal_subclonal.pdf"), p,
       width = 180/25.3, height = 50/25.3)
message("Saved: EXDF7a")

# ---- EXDF7b  ID signature activities, clonal vs subclonal -------------------
id_signatures <- sd_tsv("EXDF7b")
p <- id_signatures %>%
  mutate(timing = factor(timing, levels = c("overall", "clonal", "subclonal"))) %>%
  group_by(timing, signature) %>%
  dplyr::filter(n_distinct(trajectory) >= 2,
                sd(activity, na.rm = TRUE) > 0,
                max(activity, na.rm = TRUE) > 0) %>%
  ggplot(aes(x = signature, y = activity, fill = trajectory)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  ggpubr::stat_compare_means(aes(group = trajectory),
                             label = "p.signif", method = "kruskal.test",
                             size = 2, hide.ns = TRUE) +
  facet_wrap(~timing) +
  labs(x = "Activity", y = "Signature") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  scale_fill_manual(values = trajectory_colours)
p <- axes2lemon(p, lt = "both", bt = "both")
ggsave(fp("EXDF7b_id_signature_activities_clonal_subclonal.pdf"), p,
       width = 180/25.3, height = 50/25.3)
message("Saved: EXDF7b")

# ---- EXDF7c-h  Reptime/size, clonal & subclonal LOH/HD/Gain -----------------
plot_reptime_clonality <- function(panel, y_lab) {
  bps_df <- sd_tsv(panel)
  rt_q <- quantile(bps_df$reptime_score, c(0.33, 0.66), na.rm = TRUE)
  ggplot(bps_df, aes(x = reptime_score, y = size)) +
    geom_density_2d(alpha = 0.8, color = "black") +
    geom_point(size = 0.5, color = "grey", alpha = 0.05) +
    geom_vline(xintercept = rt_q, linetype = "dashed", color = "black") +
    facet_wrap(~trajectory) + scale_y_log10() +
    labs(x = "Replication timing score", y = y_lab, color = "Gene") +
    theme(legend.position = "none")
}

ggsave(fp("EXDF7c_LOH_clonal_size_vs_reptime.pdf"),
       plot_reptime_clonality("EXDF7c", "Size of clonal LOH segment"),
       width = 80/25.4, height = 60/25.3); message("Saved: EXDF7c")
ggsave(fp("EXDF7d_HD_clonal_size_vs_reptime.pdf"),
       plot_reptime_clonality("EXDF7d", "Size of clonal HD segment"),
       width = 80/25.4, height = 60/25.3); message("Saved: EXDF7d")
ggsave(fp("EXDF7e_GAIN_clonal_size_vs_reptime.pdf"),
       plot_reptime_clonality("EXDF7e", "Size of clonal GAIN segment"),
       width = 80/25.4, height = 60/25.3); message("Saved: EXDF7e")
ggsave(fp("EXDF7f_LOH_subclonal_size_vs_reptime.pdf"),
       plot_reptime_clonality("EXDF7f", "Size of subclonal LOH segment"),
       width = 80/25.4, height = 60/25.3); message("Saved: EXDF7f")
ggsave(fp("EXDF7g_HD_subclonal_size_vs_reptime.pdf"),
       plot_reptime_clonality("EXDF7g", "Size of subclonal HD segment"),
       width = 80/25.4, height = 60/25.3); message("Saved: EXDF7g")
ggsave(fp("EXDF7h_GAIN_subclonal_size_vs_reptime.pdf"),
       plot_reptime_clonality("EXDF7h", "Size of subclonal GAIN segment"),
       width = 80/25.4, height = 60/25.3); message("Saved: EXDF7h")

# ---- EXDF8  Simulated progression — punctuated + gradual models -------------
sims_alt <- sd_tsv("EXDF8")
met_progression_data <- sd_tsv("Fig4b") %>% dplyr::filter(gs_grade == "Metastasis")
met_summary_data <- met_progression_data %>% group_by(trajectory) %>%
  dplyr::summarise(met_first_quartile = quantile(total_poete, .25),
                   met_median = median(total_poete),
                   met_third_quartile = quantile(total_poete, .75),
                   .groups = "drop")
names(gleason_colours) <- c("1", "2", "3", "4-5")
p1 <- plot_model_simulation(sims_alt, "punctuated", gleason_colours,
                            met_summary_data, met_progression_data)
p2 <- plot_model_simulation(sims_alt, "gradual",   gleason_colours,
                            met_summary_data, met_progression_data)
p_combined <- cowplot::plot_grid(p1, p2, nrow = 2, labels = c("a)", "b)"))
ggsave(fp("EXDF8_simulated_progression_gg_alternative_models.pdf"), p_combined,
       width = 160/25.3, height = 100/25.3, dpi = 300)
message("Saved: EXDF8")

# ---- EXDF9a  Cox forest, MFS, SPOP wild-type subset -------------------------
met_relapse <- sd_tsv("EXDF9a") %>%
  mutate(high_pga = fct_relevel(high_pga, "Low"),
         advanced = fct_relevel(advanced, "Low"))
model <- coxph(Surv(time2relapse, relapse_ind) ~
                 ggroup + psa_cat + t_stage_simplified + age +
                 trajectory + mut_ith + nsubclones + advanced +
                 high_pga, data = met_relapse)
pdf(fp("EXDF9a_mfs_by_subclonal_dynamics_spop_wt.pdf"),
    height = 100/25.3, width = 120/25.3)
print(ggforest(model, data = as.data.frame(met_relapse))); dev.off()
message("Saved: EXDF9a")

# ---- EXDF9b Cox forest, MFS with TSG alteration burden -------------
met_relapse <- sd_tsv("EXDF9b") %>%
  mutate(high_pga = fct_relevel(high_pga, "Low"),
         advanced = fct_relevel(advanced, "Low"))
model <- coxph(Surv(time2relapse, relapse_ind) ~
              ggroup + psa_cat + t_stage_simplified + age +
              trajectory + advanced +
              high_pga + is_wgd + tsg_alt,
            data = met_relapse)
pdf(fp("EXDF9b_mfs_by_subclonal_dynamics_clinical_and_genetic.pdf"),
    height = 100/25.3, width = 120/25.3)
print(ggforest(model, data = as.data.frame(met_relapse))); dev.off()
message("Saved: EXDF9b")

# ---- EXDF9c Cox forest, MFS with clinical + genetic covariates ---------------------
met_relapse <- sd_tsv("EXDF9c") %>%
  mutate(high_pga = fct_relevel(high_pga, "Low"),
         advanced = fct_relevel(advanced, "Low"))

model <- coxph(Surv(time2relapse, relapse_ind) ~
              ggroup + psa_cat + t_stage_simplified + age +
              trajectory + advanced +
              high_pga + tp53_altered + rb1_altered + pten_altered + is_wgd,
            data = met_relapse)
pdf(fp("EXDF9c_mfs_by_subclonal_dynamics_tsg_alt.pdf"),
    height = 100/25.3, width = 120/25.3)
print(ggforest(model, data = as.data.frame(met_relapse))); dev.off()
message("Saved: EXDF9c")

# ---- EXDF9d  KM: D'Amico risk + evolutionary risk ---------------------------
met_relapse <- sd_tsv("EXDF9d")
km_fit <- survfit(Surv(time2relapse, relapse_ind) ~ damico_evorisk, data = met_relapse)
pdf(fp("EXDF9d_km_evorisk_damico.pdf"), height = 80/25.3, width = 80/25.3)
print(ggsurvplot(km_fit, data = met_relapse,
  palette = c("firebrick4", "dodgerblue4", "firebrick1",
              "deepskyblue3", "grey60", "grey60"),
  legend.title = "", risk.table = TRUE, pval = TRUE, conf.int = FALSE,
  xlab = "Days after diagnosis", ylab = "Metastasis-free survival",
  size = 0.4, censor.size = 2, pval.size = 2.5, fontsize = 2.5,
  font.main = 7, font.x = 7, font.y = 7, font.tickslab = 6, font.legend = 6,
  risk.table.height = 0.3, risk.table.y.text = FALSE, risk.table.title = "At risk",
  tables.theme = theme_cleantable() + theme(
    plot.title = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0, 5.5, 0, 5.5))
))
dev.off()
message("Saved: EXDF9d")

# ---- EXDF9e  KM: Cambridge Prognostic Group + evolutionary risk -------------
met_relapse <- sd_tsv("EXDF9e")
km_fit <- survfit(Surv(time2relapse, relapse_ind) ~ cpg_evorisk, data = met_relapse)
pdf(fp("EXDF9e_km_evorisk_cpg.pdf"), height = 80/25.3, width = 80/25.3)
print(ggsurvplot(km_fit, data = met_relapse,
  palette = c("grey80", "grey60", "coral", "deepskyblue3",
              "firebrick1", "dodgerblue3", "firebrick4", "dodgerblue4"),
  legend.title = "", risk.table = TRUE, pval = TRUE, conf.int = FALSE,
  xlab = "Days after diagnosis", ylab = "Metastasis-free survival",
  size = 0.4, censor.size = 2, pval.size = 2.5, fontsize = 2.5,
  font.main = 7, font.x = 7, font.y = 7, font.tickslab = 6, font.legend = 6,
  risk.table.height = 0.3, risk.table.y.text = FALSE, risk.table.title = "At risk",
  tables.theme = theme_cleantable() + theme(
    plot.title = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0, 5.5, 0, 5.5))
))
dev.off()
message("Saved: EXDF9e")

# ---- EXDF10b-g  TTF KM: Taxane vs ARPI, 2nd line ----------------------------
render_exdf10 <- function(panel, filename) {
  data <- sd_tsv(panel)
  km_fit <- survfit(Surv(TTF, ttf_event) ~ arm, data = data)
  plot_km(data, km_fit, "",
          "Days from treatment initiation", "Survival rate", filename)
  message("Saved: ", panel)
}
render_exdf10("EXDF10b", "EXDF10b_TTF_in_Taxane_vs_ARPI_2nd_line_O-I.pdf")
render_exdf10("EXDF10c", "EXDF10c_TTF_in_Taxane_vs_ARPI_2nd_line_O-II.pdf")
render_exdf10("EXDF10d", "EXDF10d_TTF_in_Taxane_vs_ARPI_2nd_line_O-III.pdf")
render_exdf10("EXDF10e", "EXDF10e_TTF_in_Taxane_vs_ARPI_2nd_line_O-I_non_CIN.Cluster-6.pdf")
render_exdf10("EXDF10f", "EXDF10f_TTF_in_Taxane_vs_ARPI_2nd_line_O-II_non_CIN.Cluster-6.pdf")
render_exdf10("EXDF10g", "EXDF10g_TTF_in_Taxane_vs_ARPI_2nd_line_O-III_non_CIN.Cluster-6.pdf")


# =============================================================================
# 7. SUPPLEMENTARY FIGURES (SF*)  ---------------------------------------------
# =============================================================================

# ---- SF1a  DPClust vs Phylogic CCF scatter, example PPCG0004a_DNA -----------
sample_data <- sd_tsv("SF1a")
sid <- unique(sample_data$sample)[1]
ccf_corr <- ggplot(sample_data, aes(x = truncated_ccf, y = phylogic_ccf)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgreen") +
  ggpubr::stat_cor(method = "pearson", color = "blue", size = 2) +
  labs(title = paste("Sample:", sid), x = "DPClust CCF", y = "PhylogicNDT CCF") +
  xlim(c(0, 1)) + ylim(c(0, 1))
ggsave(fp("SF1a_dpclust_vs_phylogic_ccf_scatter.pdf"), ccf_corr,
       height = 65/25.3, width = 65/25.3)
ggsave(fp("SF1a_dpclust_vs_phylogic_ccf_scatter.png"), ccf_corr,
       height = 65/25.3, width = 65/25.3)
message("Saved: SF1a")

# ---- SF1b  Histogram: Pearson CCF DPClust vs PhylogicNDT --------------------
dvp <- sd_tsv("SF1b")
ccf_hist <- ggplot(dvp, aes(x = pearson)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", alpha = 0.7) +
  geom_vline(xintercept = median(dvp$pearson, na.rm = TRUE),
             linetype = "dashed", color = "black") +
  labs(x = "CCF Pearson's correlation", y = "Number of samples") +
  xlim(c(0, 1.05))
ggsave(fp("SF1b_dpclust_vs_phylogic_ccf_pearson_histogram.pdf"), ccf_hist,
       width = 65/25.3, height = 65/25.3)
ggsave(fp("SF1b_dpclust_vs_phylogic_ccf_pearson_histogram.png"), ccf_hist,
       width = 65/25.3, height = 65/25.3)
message("Saved: SF1b")

# ---- SF1c  Heatmap: subclone count DPClust vs PhylogicNDT -------------------
heat_df <- sd_tsv("SF1c")
sf1c <- heat_df %>%
  ggplot(aes(x = n_subclones_dpclust, y = n_subclones_phylogic, fill = n)) +
  geom_tile(color = "white") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgreen") +
  geom_text(aes(label = n), color = "black", size = 2) +
  scale_fill_gradient(low = "white", high = "darkred") +
  scale_x_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
  theme(legend.position = "none") +
  labs(x = "DPClust subclones", y = "PhylogicNDT subclones", fill = "Samples")
ggsave(fp("SF1c_scatterplot_dpclust_vs_phylogic_nsubclones.pdf"), sf1c,
       width = 65/25.3, height = 65/25.3)
ggsave(fp("SF1c_scatterplot_dpclust_vs_phylogic_nsubclones.png"), sf1c,
       width = 65/25.3, height = 65/25.3)
message("Saved: SF1c")

# ---- SF1d  DPClust vs Phylogic histogram, example PPCG1060a_DNA -------------
sample_data <- sd_tsv("SF1d")
sid <- unique(sample_data$sample)[1]
ari <- NA_real_
sd_d_dpc <- ggplot(sample_data[!is.na(sample_data$cluster), ],
                   aes(x = ccf, fill = as.character(cluster))) +
  geom_histogram(alpha = 0.5, position = "identity") +
  geom_vline(aes(xintercept = cluster_ccf), linetype = "dashed", color = "black") +
  labs(title = paste("DPClust:", sid), x = "CCF", y = "Count", fill = "Cluster") +
  scale_x_continuous(breaks = seq(0, 1.5, by = 0.25), limits = c(0, 1.5)) +
  theme(legend.position = "bottom")
sd_d_phy <- ggplot(sample_data[!is.na(sample_data$cluster), ],
                   aes(x = phylogic_ccf, fill = as.character(phylogic_cluster))) +
  geom_histogram(alpha = 0.5, position = "identity") +
  geom_vline(aes(xintercept = phylogic_mean_cluster),
             linetype = "dashed", color = "black") +
  labs(title = paste("PhylogicNDT:", sid),
       x = "CCF", y = "Count", fill = "Cluster") +
  scale_x_continuous(breaks = seq(0, 1.5, by = 0.25), limits = c(0, 1.5)) +
  theme(legend.position = "bottom")
combined <- plot_grid(sd_d_dpc, sd_d_phy, ncol = 1, align = "v")
ggsave(fp("SF1d_dpclust_vs_phylogic_clusters.pdf"), combined,
       height = 90/25.3, width = 90/25.3)
ggsave(fp("SF1d_dpclust_vs_phylogic_clusters.png"), combined,
       height = 90/25.3, width = 90/25.3)
message("Saved: SF1d")

# ---- SF1e  ARI histogram across samples -------------------------------------
dvp <- sd_tsv("SF1e")
ari_hist <- ggplot(dvp, aes(x = ari)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", alpha = 0.7) +
  geom_vline(xintercept = median(dvp$ari, na.rm = TRUE),
             linetype = "dashed", color = "black") +
  labs(x = "Clustering ARI", y = "Number of samples")
ggsave(fp("SF1e_dpclust_vs_phylogic_ari_histogram.pdf"), ari_hist,
       width = 65/25.3, height = 65/25.3)
ggsave(fp("SF1e_dpclust_vs_phylogic_ari_histogram.png"), ari_hist,
       width = 65/25.3, height = 65/25.3)
message("Saved: SF1e")

# ---- SF1f  ARI histogram stratified by subclonal PGA ------------------------
dvp <- sd_tsv("SF1f")
ari_hist <- dvp %>%
  ggplot(aes(x = ari, fill = subclonal_pga_group)) +
  geom_histogram(alpha = 0.4, position = "identity") +
  geom_vline(xintercept = dvp %>% filter(subclonal_pga_group == "High subclonal PGA") %>%
               pull(ari) %>% median(na.rm = TRUE),
             linetype = "dashed", color = "salmon") +
  geom_vline(xintercept = dvp %>% filter(subclonal_pga_group == "Low subclonal PGA") %>%
               pull(ari) %>% median(na.rm = TRUE),
             linetype = "dashed", color = "lightblue") +
  scale_fill_manual(values = c("High subclonal PGA" = "salmon",
                               "Low subclonal PGA" = "lightblue")) +
  labs(x = "Clustering ARI", y = "Number of samples", fill = "Subclonal PGA group") +
  theme(legend.position = "top")
ggsave(fp("SF1f_dpclust_vs_phylogic_ari_histogram_by_subclonal_pga_group.pdf"), ari_hist,
       height = 65/25.3, width = 65/25.3)
ggsave(fp("SF1f_dpclust_vs_phylogic_ari_histogram_by_subclonal_pga_group.png"), ari_hist,
       height = 65/25.3, width = 65/25.3)
message("Saved: SF1f")

# ---- SF1g  CCF scatter, clonal CN regions -----------------------------------
sf1g <- sd_tsv("SF1g")
ccf_clonal <- ggplot(sf1g, aes(x = truncated_ccf, y = phylogic_ccf)) +
  geom_point(alpha = 0.01) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgreen") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  ggpubr::stat_cor(method = "spearman", color = "blue", size = 2) +
  ggtitle("Clonal CN regions") +
  labs(x = "DPClust CCF", y = "PhylogicNDT CCF") +
  xlim(c(0, 1)) + ylim(c(0, 1))
ggsave(fp("SF1g_dpclust_vs_phylogic_ccf_scatter_clonal_cn.png"), ccf_clonal,
       height = 65/25.3, width = 65/25.3)
ggsave(fp("SF1g_dpclust_vs_phylogic_ccf_scatter_clonal_cn.pdf"), ccf_clonal,
       height = 65/25.3, width = 65/25.3)
message("Saved: SF1g")

# ---- SF1h  CCF scatter, subclonal CN regions --------------------------------
sf1h <- sd_tsv("SF1h")
ccf_sub <- ggplot(sf1h, aes(x = truncated_ccf, y = phylogic_ccf)) +
  geom_point(alpha = 0.01) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgreen") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  ggpubr::stat_cor(method = "spearman", color = "blue", size = 2) +
  labs(x = "DPClust CCF", y = "PhylogicNDT CCF") +
  ggtitle("Subclonal CN regions") +
  xlim(c(0, 1)) + ylim(c(0, 1))
ggsave(fp("SF1h_dpclust_vs_phylogic_ccf_scatter_subclonal_cn.png"), ccf_sub,
       height = 65/25.3, width = 65/25.3)
ggsave(fp("SF1h_dpclust_vs_phylogic_ccf_scatter_subclonal_cn.pdf"), ccf_sub,
       height = 65/25.3, width = 65/25.3)
message("Saved: SF1h")

# ---- SF2  Clock signature (SBS1/5/40) burden vs age -------------------------
df <- sd_tsv("SF2")
p1 <- plot_facet(df, "soft", "age_at_tumour_collection")
p2 <- plot_facet(df, "hard", "age_at_tumour_collection")
p <- (p1 / p2) + plot_annotation(tag_levels = "a")
save_ggplot(p, fp("SF2_clock_sbs2age"), w = 200, h = 100)
message("Saved: SF2")

# ---- SF3  MRCA latency cohort diagram ---------------------------------------
pt_rate_mrca_estimates <- sd_tsv("SF3")
p <- plot_cohort_diagram(pt_rate_mrca_estimates)
ggsave(fp("SF3_mrca_latency_cohort_diagram.pdf"), p,
       width = 180/25.3, height = 180/25.3)
message("Saved: SF3")

# ---- SF4  Global vs patient-specific clock rate -----------------------------
global_rate_mrca <- sd_tsv("SF4")
pt_rate_mrca     <- sd_tsv("SF3")
p1 <- merge(
  pt_rate_mrca %>% dplyr::filter(acceleration == "1x") %>%
    dplyr::rename(pt_latency = latency) %>% dplyr::select(sample_id, pt_latency),
  global_rate_mrca %>% dplyr::filter(acceleration == "1x") %>%
    dplyr::rename(global_latency = latency) %>% dplyr::select(sample_id, global_latency),
  by = "sample_id") %>%
  ggplot(aes(x = pt_latency, y = global_latency)) +
  geom_point(alpha = 0.5) +
  labs(y = "MRCA latency global clock rate",
       x = "MRCA latency patient specific rate") +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
  ggpubr::stat_cor(method = "spearman") +
  xlim(c(0, 70)) + ylim(c(0, 70))
p2 <- plot_median_timing(global_rate_mrca, "#E61C2A") +
  labs(y = "Median MRCA latency (Years)")
p3 <- plot_cohort_diagram(global_rate_mrca)
design <- "A##BBBB#
           CCCCCCCC
           CCCCCCCC
           CCCCCCCC"
p <- p2 + p1 + p3 + plot_annotation(tag_levels = 'a') +
  plot_layout(design = design)
save_ggplot(p, fp("SF4_global_vs_patient_rate"), w = 180, h = 250)
message("Saved: SF4")

# ---- SF5  MRCA / CNA latency distributions by Grade Group -------------------
sf5 <- sd_rds("SF5")
t_mrca <- sf5$mrca
earliest_cnas <- sf5$earliest_cnas

p1 <- t_mrca %>% dplyr::filter(!is.na(grade_group)) %>%
  ggplot(aes(x = grade_group, y = latency, col = country)) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.1) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) +
  ggpubr::stat_compare_means(
    comparisons = list(c("1", "2"), c("2", "3"), c("3", "4+"), c("2", "4+")),
    method = "wilcox.test", label = "p.signif", size = 2) +
  facet_wrap(~acceleration) +
  labs(x = "Grade Group", y = "MRCA Latency (Years)", fill = "Country")
p1 <- axes2lemon(p1)

p2 <- ggplot(earliest_cnas %>% dplyr::filter(!is.na(grade_group)),
             aes(y = latency, x = grade_group)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.15) +
  ggpubr::stat_compare_means(
    comparisons = list(c("1", "2"), c("2", "3"), c("3", "4+"), c("2", "4+")),
    method = "wilcox.test", label = "p.signif", size = 2) +
  labs(x = "Grade Group", y = "Latency from first\n timed gain to MRCA")
p2 <- axes2lemon(p2) + theme(legend.position = "none")

top_row    <- cowplot::plot_grid(p1, NULL, nrow = 1, rel_widths = c(9, 0))
bottom_row <- cowplot::plot_grid(p2, NULL, nrow = 1, rel_widths = c(3, 6))
p <- cowplot::plot_grid(top_row, bottom_row, nrow = 2,
                        labels = c("a", "b"), label_size = 10,
                        rel_heights = c(2, 1), rel_widths = c(3, 1))
ggsave(fp("SF5_latency_by_grades.pdf"), p,
       width = 120/25.3, height = 120/25.3)
message("Saved: SF5")

# ---- SF6  Proliferation markers by Grade Group ------------------------------
sf6 <- sd_rds("SF6")
methylation_data <- sf6$methylation
gsea_rank        <- sf6$gsea_rank
p1 <- boxplot_epicmit(methylation_data, "grade_group", "epiCMIT",
                      "Grade Group", "epiCMIT")
p2 <- boxplot_epicmit(methylation_data, "grade_group", "epiCMIT.hyper",
                      "Grade Group", "epiCMIT hypermethylation")
p3 <- boxplot_epicmit(methylation_data, "grade_group", "epiCMIT.hypo",
                      "Grade Group", "epiCMIT hypomethylation")
rank <- gsea_rank$stat; names(rank) <- gsea_rank$gene
# GSEA enrichment plots require fgsea + the original pathway gene sets. If
# fgsea is unavailable, the bottom row is blank.
p4 <- tryCatch({
  pathways_fp <- "data/meta/pathways/pathways.tsv"
  if (!file.exists(pathways_fp))
    stop("pathways file not available")
  pathways <- read_delim(pathways_fp, show_col_types = FALSE) %>%
    dplyr::filter(str_detect(pathway, "HALLMARK") & str_detect(pathway, "G2M|E2F"))
  ps <- split(pathways$gene, pathways$pathway)
  fgsea::plotEnrichment(ps[["HALLMARK_G2M_CHECKPOINT"]], rank) +
    labs(title = "G2M Checkpoint")
}, error = function(e) ggplot() + theme_void() +
     labs(title = "G2M Checkpoint\n(pathway file missing)"))
p5 <- tryCatch({
  pathways_fp <- "data/meta/pathways/pathways.tsv"
  pathways <- read_delim(pathways_fp, show_col_types = FALSE) %>%
    dplyr::filter(str_detect(pathway, "HALLMARK") & str_detect(pathway, "G2M|E2F"))
  ps <- split(pathways$gene, pathways$pathway)
  fgsea::plotEnrichment(ps[["HALLMARK_E2F_TARGETS"]], rank) +
    labs(title = "E2F Targets")
}, error = function(e) ggplot() + theme_void() +
     labs(title = "E2F Targets\n(pathway file missing)"))
p <- plot_grid(p1, p2, p3, p4, p5, NULL, nrow = 2, ncol = 3,
               rel_widths = c(2, 2, 2, 3, 3, 0),
               rel_heights = c(1, 1, 1, 1, 1, 1),
               labels = c("a)", "b)", "c)", "d)", "e)", ""))
ggsave(fp("SF6_proliferation_markers_by_grade_group.png"), p,
       width = 180/25.3, height = 120/25.3, dpi = 300)
ggsave(fp("SF6_proliferation_markers_by_grade_group.pdf"), p,
       width = 180/25.3, height = 120/25.3)
message("Saved: SF6")

# ---- SF7a  Externally generated (trajectory-derivation methodology) ---------
# See README.md for the companion repository link.

# ---- SF7b  Sankey: Woodcock et al. vs PPCG trajectory -----------------------
df_sankey <- sd_tsv("SF7b")
names(trajectory_colors) <- c("Canonical", "Alternative", "Gain-enriched")
df_sankey$ppcg_trajectory <- factor(df_sankey$ppcg_trajectory,
                                    levels = c("Canonical", "Alternative", "Gain-enriched"))
p <- ggplot(df_sankey, aes(y = n, axis1 = woodcock_trajectory, axis2 = ppcg_trajectory)) +
  geom_alluvium(aes(fill = ppcg_trajectory), width = 0.15, alpha = 0.2, knot.pos = 0.4) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.05, alpha = 0.8,
               color = "white", linewidth = 0.2) +
  scale_fill_manual(values = trajectory_colors) +
  theme_void() + theme(legend.position = "none")
ggsave(fp("SF7b_sankey_woodcock_ppcg.pdf"), p,
       height = 50/25.4, width = 50/25.4, dpi = 300)
message("Saved: SF7b")

# ---- SF7c  Trajectory counts by country -------------------------------------
country_counts <- sd_tsv("SF7c")
p <- ggplot(country_counts, aes(x = country, y = n, fill = trajectory)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Country of Origin", y = "Number of Samples", fill = "Trajectory") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = trajectory_colors)
ggsave(fp("SF7c_trajectory_counts_by_country.pdf"), p,
       width = 70/25.3, height = 50/25.3, dpi = 300)
message("Saved: SF7c")

# ---- SF8 (amps)  Genome-wide GISTIC amplification scores --------------------
amps <- sd_tsv("SF8_amps")
p <- amps %>%
  ggplot(aes(x = Start, y = G.score, color = trajectory)) +
  geom_line(alpha = 0.5) +
  scale_color_manual(values = trajectory_colors) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  facet_wrap(~Chromosome, scales = "free_x", nrow = 1)
p <- axes2lemon(p)
ggsave(fp("SF8_gistic_amps.pdf"), p,
       width = 100/25.3, height = 60/25.3, dpi = 300)
message("Saved: SF8 (amps)")

# ---- SF8 (dels)  Genome-wide GISTIC deletion scores -------------------------
dels <- sd_tsv("SF8_dels")
p <- dels %>%
  ggplot(aes(x = Start, y = G.score, color = trajectory)) +
  geom_line(alpha = 0.5) +
  scale_color_manual(values = trajectory_colors) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  facet_wrap(~Chromosome, scales = "free_x", nrow = 1)
p <- axes2lemon(p)
ggsave(fp("SF8_gistic_dels.pdf"), p,
       width = 100/25.3, height = 60/25.3, dpi = 300)
message("Saved: SF8 (dels)")

# ---- SF9  Per-chromosome-arm GISTIC deletion scores -------------------------
## Reproduced without the gene-name annotations used in the original (gene
## coordinates are not in the source data). Reviewers can recompute these via
## biomaRt — see scripts/03_mechanism/02_compare_gistic.R lines 140-200.
sf9 <- sd_tsv("SF9")
chr_sizes_fp <- "data/ext/cytoBand_hg19.txt"
if (file.exists(chr_sizes_fp)) {
  chrs <- data.table::fread(chr_sizes_fp) %>%
    dplyr::select(V1, V2, V3, V4) %>%
    dplyr::rename(chr = V1, start = V2, end = V3, cytoband = V4) %>%
    dplyr::mutate(chr_arm = paste0(chr, ifelse(grepl("p", cytoband), "p", "q"))) %>%
    dplyr::mutate(chr_arm = gsub("chr", "", chr_arm))
  chr_arms <- chrs %>% group_by(chr, chr_arm) %>%
    summarise(start = min(start), end = max(end), .groups = "drop")
  relevant_loss_arms <- c("1p","1q","2q","3p","4q","5q","6q","8p","10q",
                          "11q","12p","12q","13q","16q","17p","17q","22q")
  lp <- lapply(relevant_loss_arms, function(arm) {
    row <- chr_arms[chr_arms$chr_arm == arm, ]
    if (nrow(row) == 0) return(NULL)
    chr_n <- gsub("chr", "", row$chr)
    sub <- sf9 %>% dplyr::filter(Chromosome == chr_n,
                                 Start <= row$end, End >= row$start) %>%
      dplyr::mutate(coordinate = (Start + End) / (2 * 1e6))
    p = ggplot(sub, aes(x = coordinate, y = G.score)) +
      geom_line(alpha = 0.5, aes(color = trajectory)) +
      scale_color_manual(values = trajectory_colors) +
      labs(title = paste0("Deletions in ", arm),
           x = "Coordinate (Mb)", y = "GISTIC score") +
      theme(legend.position = "none") 
    p = axes2lemon(p)
    return(p)
  })
  lp <- lp[!sapply(lp, is.null)]
  p <- cowplot::plot_grid(plotlist = lp, ncol = 3)
  ggsave(fp("SF9_Del_chr_arms_gistic.pdf"), p,
         width = 180/25.3, height = 250/25.3, dpi = 300)
  message("Saved: SF9")
} else {
  message("Skipped SF9: cytoband file ", chr_sizes_fp, " not available")
}

# ---- SF10  Per-chromosome-arm GISTIC amplification scores -------------------
sf10 <- sd_tsv("SF10")
if (file.exists(chr_sizes_fp)) {
  relevant_gain_arms <- c("1q","2p","3q","5p","7p","7q","8p","8q","9p","9q",
                          "10q","11q","12q","14q","16p","16q","17q")
  lp <- lapply(relevant_gain_arms, function(arm) {
    row <- chr_arms[chr_arms$chr_arm == arm, ]
    if (nrow(row) == 0) return(NULL)
    chr_n <- gsub("chr", "", row$chr)
    sub <- sf10 %>% dplyr::filter(Chromosome == chr_n,
                                  Start <= row$end, End >= row$start) %>%
      dplyr::mutate(coordinate = (Start + End) / (2 * 1e6))
    p = ggplot(sub, aes(x = coordinate, y = G.score)) +
      geom_line(alpha = 0.5, aes(color = trajectory)) +
      scale_color_manual(values = trajectory_colors) +
      labs(title = paste0("Amplifications in ", arm),
           x = "Coordinate (Mb)", y = "GISTIC score") +
      theme(legend.position = "none")
    p = axes2lemon(p)
    return(p)
  })
  lp <- lp[!sapply(lp, is.null)]
  p <- cowplot::plot_grid(plotlist = lp, ncol = 3)
  ggsave(fp("SF10_Amp_chr_arms_gistic.pdf"), p,
         width = 180/25.3, height = 250/25.3, dpi = 300)
  message("Saved: SF10")
} else {
  message("Skipped SF10: cytoband file ", chr_sizes_fp, " not available")
}

# ---- SF11a  TCGA molecular subtype donuts per trajectory --------------------
donut_data <- sd_tsv("SF11a")
p <- ggplot(donut_data, aes(x = 2, y = proportion, fill = tcga_subtype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  facet_wrap(~ trajectory) +
  scale_fill_brewer(palette = "Set2") +
  xlim(0.5, 2.5) +
  theme_void(base_size = 9, base_family = "sans") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold", margin = margin(b = 10)),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)) +
  labs(fill = "TCGA subtype") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             title.position = "top", title.hjust = 0.5)) +
  geom_text(
    data = donut_data %>% group_by(trajectory) %>%
      summarise(sample_count = sum(count), .groups = "drop"),
    aes(x = 0.5, y = 0, label = paste0("N=", sample_count)),
    inherit.aes = FALSE, size = 3, fontface = "italic", color = "grey30")
ggsave(fp("SF11a_evotype_overlap_tcga_subtypes_donut.pdf"), p,
       width = 225/45.3, height = 100/25.3)
message("Saved: SF11a")

# ---- SF11b  You et al. expression subtypes per trajectory -------------------
exp_donut <- sd_tsv("SF11b")
p <- ggplot(exp_donut, aes(x = 2, y = proportion, fill = expression_subtype_you)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  facet_wrap(~ trajectory) +
  scale_fill_brewer(palette = "Set2") +
  xlim(0.5, 2.5) +
  theme_void(base_size = 9, base_family = "sans") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold", margin = margin(b = 10)),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)) +
  labs(fill = "You expression subtype") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             title.position = "top", title.hjust = 0.5)) +
  geom_text(
    data = exp_donut %>% group_by(trajectory) %>%
      summarise(sample_count = sum(count), .groups = "drop"),
    aes(x = 0.5, y = 0, label = paste0("N=", sample_count)),
    inherit.aes = FALSE, size = 3, fontface = "italic", color = "grey30")
ggsave(fp("SF11b_evotype_overlap_you_subtypes.pdf"), p,
       width = 100/25.3, height = 100/25.3)
message("Saved: SF11b")

# ---- SF12  Trajectory timelines (multi-panel) -------------------------------
timeline_data <- sd_tsv("SF12")
p3 <- plot_timeline_raincloud_tj(timeline_data)
p3 <- axes2lemon(p3)
# The original SF12 also includes p1 (MRCA latency per trajectory) and p2 (CNA
# latency to MRCA per trajectory) panels. With only the saved timeline_data
# we have everything needed to reproduce these panels; we derive them inline
# from `timeline_data`.
mrca_df <- timeline_data %>% dplyr::filter(event == "MRCA") %>%
  dplyr::mutate(latency = age_sampling - latency_rel * age_sampling)
p1 <- ggplot(mrca_df, aes(x = trajectory, y = latency)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.15) +
  ggpubr::stat_compare_means(
    comparisons = list(c("Canonical", "Alternative"),
                       c("Canonical", "Gain-enriched"),
                       c("Alternative", "Gain-enriched")),
    method = "wilcox.test", label = "p.signif", size = 2) +
  scale_x_discrete(limits = c("Canonical", "Alternative", "Gain-enriched")) +
  labs(x = "Trajectory", y = "MRCA Latency (Years)")
p1 <- axes2lemon(p1) + theme(legend.position = "none")

scna_df <- timeline_data %>% dplyr::filter(event == "CNA") %>%
  dplyr::mutate(latency = age_sampling - latency_rel * age_sampling)
p2 <- ggplot(scna_df, aes(y = latency, x = trajectory)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.15) +
  ggpubr::stat_compare_means(
    comparisons = list(c("Canonical", "Alternative"),
                       c("Alternative", "Gain-enriched"),
                       c("Canonical", "Gain-enriched")),
    method = "wilcox.test", label = "p.signif", size = 2) +
  labs(x = "Trajectory", y = "Latency from first\n timed gain to MRCA")
p2 <- axes2lemon(p2) + theme(legend.position = "none")

top_row    <- cowplot::plot_grid(p3, NULL, nrow = 1,
                                  rel_widths = c(13, 0), labels = c("a", ""))
bottom_row <- cowplot::plot_grid(NULL, p1, NULL, p2, NULL, nrow = 1,
                                  rel_widths = c(1, 5, 1, 5, 1),
                                  labels = c("", "b", "", "c", ""))
p <- cowplot::plot_grid(top_row, bottom_row, nrow = 2, label_size = 10,
                        rel_heights = c(1, 1), rel_widths = c(1, 1))
ggsave(fp("SF12_latency_by_trajectories.pdf"), p,
       width = 155/25.3, height = 120/25.3)
message("Saved: SF12")

# ---- SF13  Sankey: CombiMets primary -> metastasis trajectory transitions ---
pair_tj <- sd_tsv("SF13")
sankey_data <- pair_tj %>% ggsankey::make_long(s1_tj, s2_tj)
p <- ggplot(sankey_data,
            aes(x = x, next_x = next_x, node = node, next_node = next_node,
                fill = factor(node))) +
  geom_sankey(flow.alpha = 0.5, node.color = "black") +
  geom_sankey_label(aes(label = node), size = 3, color = "white",
                    data = sankey_data %>% filter(x %in% c("s1_tj", "s2_tj"))) +
  theme_minimal() +
  labs(title = "Trajectory changes from primary tumour to matched metastasis",
       x = "Sample Timepoint", y = "", fill = "Trajectory") +
  scale_x_discrete(labels = c("s1_tj" = "Primary", "s2_tj" = "Metastasis")) +
  scale_fill_manual(values = c("1" = "#8DA0CB",
                               "2" = "#AAF0C9",
                               "3" = "#C04667")) +
  theme(legend.position = "none")
ggsave(fp("SF13_tj_transitions_sankey.pdf"), p,
       width = 150/25.3, height = 100/25.3)
message("Saved: SF13")

# ---- SF14  Evolutionary rates by Gleason score ------------------------------
evo_metrics_rates <- sd_tsv("SF14") %>% dplyr::mutate(gs_grade = ifelse(gs_grade == "4+", "4-5", gs_grade))
p1 <- boxplot_rates(evo_metrics_rates, "clonal_pga_per_year", "Clonal PGA per Year")
p2 <- boxplot_rates(evo_metrics_rates, "subclonal_pga_per_year", "Subclonal PGA per Year")
p3 <- boxplot_rates(evo_metrics_rates, "clonal_poteo_per_year",
                    "Yearly clonal rate \nof trajectory progression")
p4 <- boxplot_rates(evo_metrics_rates, "subclonal_ndrivers_per_year",
                    "Yearly subclonal rate \nof trajectory progression")
p <- cowplot::plot_grid(p1, p2, p3, p4, nrow = 2,
                        labels = c("a)", "b)", "c)", "d)"))
ggsave(fp("SF14_evolutionary_rates_by_gs.pdf"), p,
       width = 130/25.4, height = 130/25.4, units = "in")
message("Saved: SF14")

# ---- SF15  Rate forest plots (3 Cox models) ---------------------------------
met_relapse <- sd_tsv("SF15") %>% as.data.frame()
p1 <- ggforest(coxph(Surv(time2relapse, relapse_ind) ~ imminent,
                     data = met_relapse[met_relapse$gs_group == "2", ]),
               main = "")
p2 <- ggforest(coxph(Surv(time2relapse, relapse_ind) ~
                       imminent + gs_group + psa + age + t_stage,
                     data = met_relapse[met_relapse$gs_group != "1", ]),
               main = "")
p3 <- ggforest(coxph(Surv(time2relapse, relapse_ind) ~
                       median_time_to_met + gs_group + psa + t_stage,
                     data = met_relapse[met_relapse$gs_group != "1", ]),
               main = "")
combined_forest <- cowplot::plot_grid(
  p3, p2, p1,
  rel_heights = c(5, 5, 2),
  ncol = 1,
  labels = c("a)", "b)", "c)"),
  label_size = 9,
  align = "v", axis = "l"
)
pdf(fp("SF15_rate_forest_plots.pdf"), width = 120/25.3, height = 180/25.3)
print(combined_forest); dev.off()
message("Saved: SF15")

# =============================================================================
message("\nAll reproducible panels written to: ", fig_out_dir)
message("Schematics (Fig1b, Fig2a, Fig4c, Fig5b, EXDF10a) are not generated by code.")
message("External panels (EXDF5, SF7a) are produced by companion repos — see README.md.")
