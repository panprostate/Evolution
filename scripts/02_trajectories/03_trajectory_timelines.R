# Analysis of PPCG MRCA, WGD and SCNA timing estimates in relation to diagnosis
# Across each of the three different evolutionary trajectories

rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(patchwork)
library(lemon)

# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

# Combine MRCA, CNA, and WGD timing estimates into a single long-format data frame
prepare_timeline_data <- function(time_cnas, t_mrca, t_wgd) {
  # Prepare WGD data
  t_wgd <- t_wgd %>%
    dplyr::filter(acceleration == "1x") %>%
    dplyr::mutate(event = "WGD") %>%
    dplyr::mutate(wgd_latency = latency + mrca_latency) %>%  # Total latency from diagnosis
    dplyr::mutate(age_sampling = get_age(sample_id), gs_group = get_gs_group(sample_id)) %>% 
    dplyr::mutate(rt_event_pt_rate = age_sampling - wgd_latency)  

  # Prepare MRCA data
  t_mrca <- t_mrca %>%
    dplyr::filter(acceleration == "1x") %>%
    dplyr::mutate(event = "MRCA") %>%
    dplyr::mutate(age_sampling = get_age(sample_id), gs_group = get_gs_group(sample_id)) %>% 
    dplyr::mutate(rt_event_pt_rate = age_sampling - latency)  # Convert to real time

  # Prepare CNA data
  time_cnas <- time_cnas %>%
    dplyr::mutate(event = "CNA", rt_event_pt_rate = rt_cn_pt_rate) %>% 
    dplyr::mutate(age_sampling = get_age(sample_id), gs_group = get_gs_group(sample_id))

  # Combine all events into a single data frame
  # Note: We will not consider SCNA events in WGD samples, as timable events in WGD samples
  # occur as part of the WGD event and are not independent events that can be timed separately
  timeline_data <- bind_rows(
    time_cnas %>% dplyr::select(sample = sample_id, rt_event_pt_rate, event, trajectory, age_sampling) %>% dplyr::filter(!sample %in% t_wgd$sample_id),
    t_mrca %>% dplyr::select(sample = sample_id, rt_event_pt_rate, event, trajectory, age_sampling),
    t_wgd %>% dplyr::select(sample = sample_id, rt_event_pt_rate, event, trajectory, age_sampling)
  )

  # Normalize by age
  timeline_data <- timeline_data %>%
    dplyr::mutate(
      latency_rel = rt_event_pt_rate / age_sampling,
      age_normalized = age_sampling / age_sampling
    )

  return(timeline_data)
}

# Refactor bootstrap function to handle long-format data
bootstrap_timeline <- function(df, b = 1000) {
  boot_timeline <- data.frame()
  df <- df %>% dplyr::filter(!is.na(event))
  for (i in 1:b) {
    idx <- sample(1:nrow(df), size = nrow(df), replace = TRUE)
    boot_sample <- df[idx, ] %>%
      group_by(event) %>%
      summarise(
        estimate = median(latency_rel, na.rm = TRUE),
        .groups = "drop"
      )
    boot_timeline <- bind_rows(boot_timeline, boot_sample)
  }

  # Summarize bootstrap results
  boot_summary <- boot_timeline %>%
    group_by(event) %>%
    summarise(
      median = median(estimate),
      lowci = quantile(estimate, 0.025, na.rm = TRUE),
      highci = quantile(estimate, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  return(boot_summary)
}

plot_timeline_raincloud <- function(timeline_data, item_plot = c("MRCA", "CNA", "WGD")) {
  # Filter data for selected events
  timeline_data <- timeline_data %>% dplyr::filter(event %in% item_plot)

  # Bootstrap estimates for each group
  t_all <- timeline_data
  t_tj1 <- timeline_data[timeline_data$trajectory == "Canonical",]
  t_tj2 <- timeline_data[timeline_data$trajectory == "Alternative",]
  t_tj3 <- timeline_data[timeline_data$trajectory == "Gain-enriched",]

  timeline_plot <- rbind(
    dplyr::mutate(bootstrap_timeline(t_tj1), trajectory = "Canonical"),
    dplyr::mutate(bootstrap_timeline(t_tj2), trajectory = "Alternative"),
    dplyr::mutate(bootstrap_timeline(t_tj3), trajectory = "Gain-enriched")
  )

  # Create final label with (n=...)
  timeline_plot <- timeline_plot %>% 
    mutate(y_label = trajectory)

  # Prepare raw data for plotting
  raw_data <- timeline_data %>% filter(!is.na(trajectory))

  raw_all <- timeline_data %>% mutate(trajectory = "All")
  raw_plot_data <- bind_rows(raw_data, raw_all)

  label_map <- timeline_plot %>% dplyr::select(trajectory, y_label)
  raw_plot_data <- left_join(raw_plot_data, label_map, by = "trajectory")

  order_levels <- timeline_plot$y_label
  timeline_plot$y_label <- factor(timeline_plot$y_label, levels = unique(order_levels))
  raw_plot_data$y_label <- factor(raw_plot_data$y_label, levels = unique(order_levels))

  # Plot
  p <- ggplot() +
    geom_segment(data = timeline_plot, aes(x = 0, xend = 1, y = y_label, yend = y_label), 
                 linewidth = 2, col = "grey70", alpha = 0.2)

  p <- p +
    geom_jitter(data = raw_plot_data, aes(x = latency_rel, y = y_label, color = event), 
                height = 0.15, width = 0, size = 1, alpha = 0.1) +
    # geom_segment(data = timeline_plot, aes(y = y_label, yend = y_label, x = lowci, xend = highci, color = event), 
    #               linewidth = 2, alpha = 0.25) +
    geom_point(data = timeline_plot, aes(y = y_label, x = median, color = event), 
                size = 5, stroke = 0.5, alpha = 0.8)

  p <- p + scale_colour_manual(values = c("MRCA" = "#D95F02", "CNA" = "#9467bd", "WGD" = "#3979bb"))
  p <- p + labs(
    y = "", 
    x = "Fraction of time to diagnosis", 
    color = "Event"
  ) +
    scale_y_discrete(
      limits = c("Canonical", "Alternative", "Gain-enriched"), 
      labels = c("Canonical", "Alternative", "Gain-enriched")
    ) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")

  return(p)
}

# PATHS
outdir <- "outputs/02_trajectories/trajectory_timeline"; figdir <- "figures/02_trajectories/trajectory_timeline"
timedir <- "outputs/01_landscape/real_timing"
dir.create(figdir)
dir.create(outdir)
trajectories_fps <- list.files("outputs/02_trajectories/", pattern = "PPCG_Feb2026_.*mergedseg_with_clonality.txt", full = TRUE)
clinfile <- "data/meta/PPCG_donors_clin_20241217.csv"
wgd_time_fp <- file.path(timedir, "latency_estimates_wgd_patient_rate_summary.tsv")
mrca_time_fp <- file.path(timedir, "latency_estimates_mrca_patient_rate_summary.tsv")
cna_time_fp <- file.path(timedir, "time_cnas.rds")

### LOAD DATA AND WRANGLE DATA --------------------------------------
t_mrca <- read_delim(mrca_time_fp, delim = "\t") %>% dplyr::filter(acceleration == "1x")
t_wgd <- read_delim(wgd_time_fp, delim = "\t") %>% dplyr::filter(acceleration == "1x")
time_cnas <- readRDS(cna_time_fp) %>% dplyr::rename(sample_id = sample)
tj_samples <- load_trajectories(trajectories_fps) %>% dplyr::rename(sample_id = sample)

trajectory_dict <- c("Ordering 1" = "Canonical", "Ordering 2" = "Alternative", "Ordering 3" = "Gain-enriched")
tj_samples$trajectory <- trajectory_dict[tj_samples$trajectory]

# Add trajectory information
t_mrca <- merge(t_mrca, tj_samples)
t_wgd <- merge(t_wgd, tj_samples)
time_cnas <- merge(time_cnas, tj_samples)

# ANALYSIS -----------------------------------------------------------
# MRCA timing
trajectories = unique(t_mrca$trajectory)
for (tj in trajectories){
    print(paste0("Trajectory: ", tj))
    boot_mrca = boot_median(t_mrca %>% dplyr::filter(acceleration == "1x", trajectory == tj) %>% dplyr::pull(latency), n = 1000)
    median(boot_mrca); quantile(boot_mrca, 0.025); quantile(boot_mrca, 0.975)
    print(paste0("MRCA: ", round(median(boot_mrca), 2), " (", round(quantile(boot_mrca, 0.025), 2), "-", round(quantile(boot_mrca, 0.975), 2), ")"))
}

p1 = t_mrca %>%
  # dplyr::mutate(acceleration = paste0("Rate of acceleration: ", acceleration)) %>%
  ggplot(aes(x = trajectory, y = latency)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.2), alpha = 0.15) +  
  ggpubr::stat_compare_means(
    comparisons = list(c("Canonical", "Alternative"), c("Canonical", "Gain-enriched"), c("Alternative", "Gain-enriched")),
    method = "wilcox.test",
    label = "p.signif", 
    size = 2
  ) + 
  facet_wrap(~acceleration) +
  scale_x_discrete(limits = c("Canonical", "Alternative", "Gain-enriched")) +
  labs(x = "Grade Group", y = "MRCA Latency (Years)") 

p1 = axes2lemon(p1) + theme(legend.position = "none")


# 2. WGD timing
# Frequency of WGD across trajectories
tj_samples$is_wgd = tj_samples$sample_id %in% get_wgd()
table(tj_samples$trajectory, tj_samples$is_wgd) 
# Canonical: 34/453; # Alternative: 12/130; Gain-enriched: 2/47
fisher.test(table(tj_samples$trajectory, tj_samples$is_wgd)) # p = 0.54

table(t_wgd$trajectory) # Canonical: 32; Alternative: 12; Gain-enriched: 2

for (tj in trajectories){
    print(paste0("Trajectory: ", tj))
    boot_wgd = boot_median(t_wgd %>% dplyr::filter(acceleration == "1x", trajectory == tj) %>% dplyr::mutate(latency = latency + mrca_latency) %>% dplyr::pull(latency), n = 1000)
    median(boot_wgd, na.rm = T); quantile(boot_wgd, 0.025, na.rm = TRUE); quantile(boot_wgd, 0.975, na.rm = TRUE)
    print(paste0("WGD: ", round(median(boot_wgd, na.rm = T), 2), " (", round(quantile(boot_wgd, 0.025, na.rm = T), 2), "-", round(quantile(boot_wgd, 0.975, na.rm = T), 2), ")"))
}

# 3. Time of earliest SCNA event relative to diagnosis

# Filtering of timed SCNA segments
# 1. Filter out segments with less than 10 mutations (only reliable timing calls)
# 2. Filter fits with mtime_cn == 1 or mtime_cn == 0 as those are more likely to represent poor fits
time_cnas$total_muts <- time_cnas$clock_early + time_cnas$clock_late + time_cnas$clock_na
time_cnas <- time_cnas[time_cnas$total_muts >= 10,]
time_cnas <- time_cnas[time_cnas$mtime_cn > 0 & time_cnas$mtime_cn < 1,]

# Bootstrap estimates of earliest SCNA time
earliest_cnas = time_cnas %>% dplyr::filter(!sample_id %in% get_wgd()) %>% dplyr::group_by(sample_id, trajectory) %>%
    dplyr::mutate(latency = rt_mrca_pt_rate - rt_cn_pt_rate) %>%
    dplyr::arrange(desc(latency)) %>% 
    # We just need to take the earliest-occurring event of each type, as WGD and MRCA are single per sample
    # We don't get the ties as we only care about earliest occurring
    slice_min(order_by = latency, n = 1, with_ties = FALSE) %>%
    ungroup() 

for (tj in trajectories){
    print(paste0("Trajectory: ", tj))
    boot_scna = boot_median(earliest_cnas %>% dplyr::filter(trajectory == tj) %>% dplyr::mutate(latency = get_age(sample_id) - rt_cn_pt_rate) %>% dplyr::pull(latency), n = 1000)
    median(boot_scna); quantile(boot_scna, 0.025); quantile(boot_scna, 0.975)
    print(paste0("SCNA: ", round(median(boot_scna), 2), " (", round(quantile(boot_scna, 0.025), 2), "-", round(quantile(boot_scna, 0.975), 2), ")"))
}

p2 = ggplot(earliest_cnas, aes(y = latency, x = trajectory)) + 
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.2), alpha = 0.15) +  
  ggpubr::stat_compare_means(
    comparisons = list(c("Canonical", "Alternative"), c("Alternative", "Gain-enriched"), c("Canonical", "Gain-enriched")),
    method = "wilcox.test",
    label = "p.signif", 
    size = 2
  ) + 
  labs(x = "Trajectory", y = "Latency from first\n timed gain to MRCA")

p2 = axes2lemon(p2) + theme(legend.position = "none")

# top_row = cowplot::plot_grid(p1, NULL, nrow = 1, rel_widths = c(9, 0))
# bottom_row = cowplot::plot_grid(p2, NULL, nrow = 1, rel_widths = c(3, 6))

# PLOT TIMELINES ------------------------------------------------------
timeline_data <- prepare_timeline_data(time_cnas, t_mrca, t_wgd)

# For SCNAs, we consider the earliest occurring SCNA event per sample
timeline_data <- timeline_data %>%
    dplyr::filter(!is.na(latency_rel)) %>% 
    group_by(sample, event) %>%
    # We just need to take the earliest-occurring event of each type, as WGD and MRCA are single per sample
    # We don't get the ties as we only care about earliest occurring
    slice_min(order_by = latency_rel, n = 1, with_ties = FALSE) %>%
    ungroup()


# Generate plot
p3 <- plot_timeline_raincloud(timeline_data)
p3 <- axes2lemon(p3)

top_row = cowplot::plot_grid(p3, NULL, nrow = 1, rel_widths = c(13, 0), labels = c("a", ""))
bottom_row = cowplot::plot_grid(NULL, p1, NULL, p2, NULL, nrow = 1, rel_widths = c(1, 5, 1, 5, 1), labels = c("", "b", "", "c", ""))
p = cowplot::plot_grid(top_row, bottom_row, nrow = 2, label_size = 10, rel_heights = c(1, 1), rel_widths = c(1, 1))

write_tsv(
  timeline_data %>% dplyr::select(trajectory, event, latency_rel, age_sampling),
  file.path(outdir, "SF12_source_data.tsv")
)
ggsave(file.path(figdir, "SF12_latency_by_trajectories.pdf"), width = 155/25.3, height = 120/25.3)


