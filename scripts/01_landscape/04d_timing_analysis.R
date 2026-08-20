# Analysis of PPCG MRCA, WGD and SCNA timing estimates in relation to diagnosis
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
    dplyr::mutate(age_sampling = get_age(sample), gs_group = get_gs_group(sample))

  # Combine all events into a single data frame
  # Note: We will not consider SCNA events in WGD samples, as timable events in WGD samples
  # occur as part of the WGD event and are not independent events that can be timed separately
  timeline_data <- bind_rows(
    time_cnas %>% dplyr::select(sample, rt_event_pt_rate, event, gs_group, age_sampling) %>% dplyr::filter(!sample %in% t_wgd$sample_id),
    t_mrca %>% dplyr::select(sample = sample_id, rt_event_pt_rate, event, gs_group, age_sampling),
    t_wgd %>% dplyr::select(sample = sample_id, rt_event_pt_rate, event, gs_group, age_sampling)
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
  t_gs12 <- timeline_data[timeline_data$gs_group %in% c("1", "2"),]
  t_gs1 <- timeline_data[timeline_data$gs_group == "1",]
  t_gs2 <- timeline_data[timeline_data$gs_group == "2",]
  t_gs3  <- timeline_data[timeline_data$gs_group == "3",]
  t_gs45 <- timeline_data[timeline_data$gs_group %in% c("4", "4+", "5"),]

  timeline_plot <- rbind(
    dplyr::mutate(bootstrap_timeline(t_gs1), ggroup = "Grade Group 1"),
    dplyr::mutate(bootstrap_timeline(t_gs2), ggroup = "Grade Group 2"),
    dplyr::mutate(bootstrap_timeline(t_gs3), ggroup = "Grade Group 3"),
    dplyr::mutate(bootstrap_timeline(t_gs45), ggroup = "Grade Group 4-5")
  )

  # Create final label with (n=...)
  timeline_plot <- timeline_plot %>% 
    mutate(y_label = ggroup)

  # Prepare raw data for plotting
  raw_data <- timeline_data %>%
    mutate(ggroup = case_when(
      gs_group == "1" ~ "Grade Group 1",
      gs_group == "2" ~ "Grade Group 2",
      gs_group == "3" ~ "Grade Group 3",
      gs_group %in% c("4", "4+", "5") ~ "Grade Group 4-5",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(ggroup))

  raw_all <- timeline_data %>% mutate(ggroup = "All")
  raw_plot_data <- bind_rows(raw_data, raw_all)

  label_map <- timeline_plot %>% dplyr::select(ggroup, y_label)
  raw_plot_data <- left_join(raw_plot_data, label_map, by = "ggroup")

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
    x = "Fraction of time to diagnosis"
  ) +
    scale_y_discrete(
      limits = c("Grade Group 1", "Grade Group 2", "Grade Group 3", "Grade Group 4-5"), 
      labels = c(paste0("Grade Group 1\n (n = ", n_distinct(t_gs1$sample), ")"), paste0("Grade Group 2\n (n = ", n_distinct(t_gs2$sample), ")"), paste0("Grade Group 3\n (n = ", n_distinct(t_gs3$sample), ")"), paste0("Grade Group 4-5\n (n = ", n_distinct(t_gs45$sample), ")"))  
    ) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")

  return(p)
}

# PATHS
outdir <- "outputs/01_landscape/real_timing/"; figdir <- "figures/01_landscape/real_timing/"
dir.create(figdir)
dir.create(outdir)
clinfile <- "data/meta/PPCG_donors_clin_20241217.csv"
wgd_time_fp <- file.path(outdir, "latency_estimates_wgd_patient_rate_summary.tsv")
mrca_time_fp <- file.path(outdir, "latency_estimates_mrca_patient_rate_summary.tsv")
cna_time_fp <- file.path(outdir, "time_cnas.rds")

### LOAD DATA AND WRANGLE DATA --------------------------------------
t_mrca <- read_delim(mrca_time_fp, delim = "\t")
t_wgd <- read_delim(wgd_time_fp, delim = "\t")
time_cnas <- readRDS(cna_time_fp)

# ANALYSIS ----------------------------------------------------------
# 1. MRCA timing and Grade Group
t_mrca$grade_group <- get_gs_group(t_mrca$sample_id, clinfile)
t_mrca$acceleration <- factor(t_mrca$acceleration, levels = paste0(c(1, 2.5, 5, 7.5, 10, 20), "x"))
t_mrca$country <- get_country(t_mrca$sample_id, clinfile)

boot_mrca = boot_median(t_mrca %>% dplyr::filter(acceleration == "1x") %>% dplyr::pull(latency), n = 1000)
median(boot_mrca); quantile(boot_mrca, 0.025); quantile(boot_mrca, 0.975)

t_mrca %>% dplyr::filter(acceleration == "1x") %>% group_by(grade_group) %>% summarise(median_latency = median(latency, na.rm = TRUE), n = n())

p1 = t_mrca %>%
  dplyr::filter(!is.na(grade_group)) %>%
  # dplyr::mutate(acceleration = paste0("Rate of acceleration: ", acceleration)) %>%
  ggplot(aes(x = grade_group, y = latency, col = country)) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.1) +  
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) + 
  ggpubr::stat_compare_means(
    comparisons = list(c("1", "2"), c("2", "3"), c("3", "4+"), c("2", "4+")),
    method = "wilcox.test",
    label = "p.signif", 
    size = 2
  ) + 
  facet_wrap(~acceleration) +
  labs(x = "Grade Group", y = "MRCA Latency (Years)", fill = "Country")

p1 = axes2lemon(p1)

t_mrca$is_low_grade = t_mrca$grade_group %in% c("1", "2")
lm(latency ~ is_low_grade + country, data = t_mrca %>% dplyr::filter(acceleration == "1x")) %>% summary()

# 2. WGD timing
# Bootstrap estimates of WGD time
boot_wgd = boot_median(t_wgd %>% dplyr::filter(acceleration == "1x") %>% dplyr::mutate(latency = latency + mrca_latency) %>% dplyr::pull(latency), n = 1000)
median(boot_wgd); quantile(boot_wgd, 0.025); quantile(boot_wgd, 0.975)

# 3. Time of earliest SCNA event relative to diagnosis

# Filtering of timed SCNA segments
# 1. Filter out segments with less than 10 mutations (only reliable timing calls)
# 2. Filter fits with mtime_cn == 1 or mtime_cn == 0 as those are more likely to represent poor fits
time_cnas$total_muts <- time_cnas$clock_early + time_cnas$clock_late + time_cnas$clock_na
time_cnas <- time_cnas[time_cnas$total_muts >= 10,]
time_cnas <- time_cnas[time_cnas$mtime_cn > 0 & time_cnas$mtime_cn < 1,]

# Bootstrap estimates of earliest SCNA time
earliest_cnas = time_cnas %>% dplyr::filter(!sample %in% get_wgd()) %>% dplyr::group_by(sample) %>%
    dplyr::mutate(latency = rt_mrca_pt_rate - rt_cn_pt_rate) %>%
    dplyr::arrange(desc(latency)) %>% 
    # We just need to take the earliest-occurring event of each type, as WGD and MRCA are single per sample
    # We don't get the ties as we only care about earliest occurring
    slice_min(order_by = latency, n = 1, with_ties = FALSE) %>%
    ungroup() 

boot_scna = boot_median(earliest_cnas %>% dplyr::mutate(latency = get_age(sample) - rt_cn_pt_rate) %>% dplyr::pull(latency), n = 1000)
median(boot_scna); quantile(boot_scna, 0.025); quantile(boot_scna, 0.975)

earliest_cnas$grade_group = get_gs_group(earliest_cnas$sample)
p2 = ggplot(earliest_cnas %>% dplyr::filter(!is.na(grade_group)), aes(y = latency, x = grade_group)) + 
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.2), alpha = 0.15) +  
  ggpubr::stat_compare_means(
    comparisons = list(c("1", "2"), c("2", "3"), c("3", "4+"), c("2", "4+")),
    method = "wilcox.test",
    label = "p.signif", 
    size = 2
  ) + 
  labs(x = "Grade Group", y = "Latency from first\n timed gain to MRCA")

p2 = axes2lemon(p2) + theme(legend.position = "none")

top_row = cowplot::plot_grid(p1, NULL, nrow = 1, rel_widths = c(9, 0))
bottom_row = cowplot::plot_grid(p2, NULL, nrow = 1, rel_widths = c(3, 6))
p = cowplot::plot_grid(top_row, bottom_row, nrow = 2, labels = c("a", "b"), label_size = 10, rel_heights = c(2, 1), rel_widths = c(3, 1))

saveRDS(list(
  mrca          = t_mrca %>% dplyr::select(grade_group, latency, country, acceleration),
  earliest_cnas = earliest_cnas %>% dplyr::select(grade_group, latency)
), file.path(outdir, "SF5_source_data.rds"))
ggsave(file.path(figdir, "SF5_latency_by_grades.pdf"), width = 120/25.3, height = 120/25.3)

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
p <- plot_timeline_raincloud(timeline_data)
p <- axes2lemon(p) 

# Save plot
write_tsv(
  timeline_data %>% dplyr::select(event, gs_group, latency_rel),
  file.path(outdir, "Fig1c_source_data.tsv")
)
ggsave(file.path(figdir, "Fig1c_timelines_mrca_scna_wgd.pdf"), plot = p, width = 120/25.3, height = 60/25.3)
ggsave(file.path(figdir, "Fig1c_timelines_mrca_scna_wgd.png"), plot = p, width = 120/25.3, height = 80/25.3)
