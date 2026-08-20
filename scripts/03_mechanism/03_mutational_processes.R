# Analysis of mutational processes across evolutionary trajectories

rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(RColorBrewer)
library(lemon)
library(colorspace)
library(scales)
library(cowplot)
library(grid)
library(qs)
library(foreach)
library(doMC)
library(ggpubr)

NCORES=10
registerDoMC(NCORES)

# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

# Coming from plot_theme.R
names(trajectory_colours) = c("Canonical", "Alternative", "Gain-enriched")


load_sigtimer <- function(fp, tj_samples, sigtype, sig_names, method = "nnls", add_zero = TRUE){
  sigtimer_obj = qs::qread(fp)
  df = sigtimer_obj[[method]]@estimate %>% dplyr::filter(info == "mean") %>% dplyr::filter(signature != "error")
  sid = extract_ppcg_id(fp, full = F)
  trajectory = tj_samples %>% dplyr::filter(sample == sid) %>% dplyr::pull(trajectory) %>% unique()

  if (length(trajectory) == 0|length(sid) == 0|nrow(df) == 0) {
    return(data.frame())
  }

  out = data.frame(
    sample = sid,
    trajectory = trajectory,
    signature = df$signature,
    activity = df$prop,
    sigtype = sigtype, 
    timing = ifelse(df$epoch == "epoch_1", "clonal", "subclonal")
  )

  # add 0 for missing signatures including clonal and subclonal
  if (add_zero) {
    out = out %>% 
     tidyr::complete(
      sample, trajectory, sigtype, signature = sig_names, 
      timing = c("clonal", "subclonal"), fill = list(activity = 0)) %>% 
      dplyr::select(sample, trajectory, sigtype, signature, activity, timing)
  }
}

get_sample_order <- function(df, traj_name, target_imfs) {
  df %>%
    dplyr::filter(trajectory == traj_name, imf_cluster %in% target_imfs) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(target_sum = sum(activity, na.rm = TRUE)) %>%
    dplyr::arrange(target_sum) %>% # Ascending so highest is at the top of the y-axis
    dplyr::pull(sample)
}


base_cols <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

make_soft_palette <- function(n) {
  colorRampPalette(alpha(base_cols, 0.7))(n)
}

make_sig_colours <- function(sigactivities) {
  
  sig_df <- sigactivities %>%
    group_by(sigtype, signature) %>%
    summarise(total = sum(activity), .groups = "drop") %>%
    arrange(sigtype, desc(total))
  
  split_sigs <- split(sig_df$signature, sig_df$sigtype)
  
  sig_colours <- lapply(split_sigs, function(sigs) {
    cols <- make_soft_palette(length(sigs))
    setNames(cols, sigs)
  })
  
  sig_colours
}

make_sig_panel <- function(data,
                           sigtype_name,
                           trajectory_name,
                           sample_order,
                           signature_order,
                           fill_values,
                           show_y = FALSE,
                           show_title = FALSE) {

  panel_df <- data %>%
    filter(sigtype == sigtype_name, trajectory == trajectory_name) %>%
    mutate(
      sample = factor(sample, levels = sample_order),
      signature = factor(signature, levels = signature_order)
    ) %>%
    select(sample, trajectory, signature, activity) %>%
    tidyr::complete(
      sample = factor(sample_order, levels = sample_order),
      signature = factor(signature_order, levels = signature_order),
      fill = list(activity = 0)
    ) %>%
    mutate(trajectory = trajectory_name)

  ggplot(panel_df, aes(y = sample, x = activity, fill = signature)) +
    geom_bar(stat = "identity", width = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = fill_values, breaks = signature_order, drop = FALSE) +
    labs(
      x = "Signature activity",
      y = NULL,
      title = if (show_title) trajectory_name else NULL
    ) +
    coord_flip() +
    theme_bw(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(size = 8, face = "plain", hjust = 0.5),

      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
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
      axis.line.y  = if (show_y) element_line(colour = "black") else element_blank()
    )
}

make_sigtype_row <- function(data,
                             sigtype_name,
                             traj_orders,
                             sig_levels_by_type,
                             sig_colours,
                             samples_prop,
                             show_titles = FALSE) {

  p1 <- make_sig_panel(data, sigtype_name, "Canonical",
                       traj_orders[["Canonical"]],
                       sig_levels_by_type[[sigtype_name]],
                       sig_colours[[sigtype_name]],
                       TRUE, show_titles)

  p2 <- make_sig_panel(data, sigtype_name, "Alternative",
                       traj_orders[["Alternative"]],
                       sig_levels_by_type[[sigtype_name]],
                       sig_colours[[sigtype_name]])

  p3 <- make_sig_panel(data, sigtype_name, "Gain-enriched",
                       traj_orders[["Gain-enriched"]],
                       sig_levels_by_type[[sigtype_name]],
                       sig_colours[[sigtype_name]])

  row_plot <- cowplot::plot_grid(
    p1, p2, p3,
    ncol = 3,
    align = "h",
    rel_widths = c(
      samples_prop["Canonical"],
      samples_prop["Alternative"],
      samples_prop["Gain-enriched"]
    )
  )
}



# PATHS
outdir <- "outputs/03_mechanism/mutational_processes"
figdir <- "figures/03_mechanism/mutational_processes"
dir.create(outdir, recursive = TRUE); dir.create(figdir, recursive = TRUE)
trajectories_fps <- list.files("outputs/02_trajectories/", pattern = "PPCG_Feb2026.*mergedseg_with_clonality.txt", full = TRUE)
imf_activities_fps <- "data/ext/Activities_Clusters_Classification_PPCG_CIN.rds"
sigtimer_dir <- "outputs/00_preprocessing/sigtimer/"
sbs_timing_fps <- list.files(file.path(sigtimer_dir, "sbs", "clonal_subclonal"), pattern = ".*qs", full = TRUE, recursive = TRUE)
id_timing_fps <- list.files(file.path(sigtimer_dir, "id", "clonal_subclonal"), pattern = ".*qs", full = TRUE, recursive = TRUE)

# LOAD DATA
tj_samples = load_trajectories(trajectories_fps)
tj_samples = tj_samples %>% dplyr::mutate(trajectory = case_when(
    trajectory == "Ordering 1" ~ "Canonical",
    trajectory == "Ordering 2" ~ "Alternative",
    trajectory == "Ordering 3" ~ "Gain-enriched",
    TRUE ~ trajectory
))

imf_activities = readRDS(imf_activities_fps) %>% dplyr::rename(sample = WGS_AssayID)
# Add trajectory information on IMF activities
imf_activities = merge(tj_samples, imf_activities)
colnames(imf_activities) = str_replace(colnames(imf_activities), "CIN.Cluster-", "IMF")
imfs = colnames(imf_activities)[grepl("IMF", colnames(imf_activities))]
imf_activities = imf_activities %>% tidyr::pivot_longer(cols=all_of(imfs), names_to='imf_cluster', values_to='activity') 

# Load individual signature activities per sample
sigactivities = readRDS("data/ext/Sample_per_Signature_Activities_PPCG_dCIN.rds") 
signatures = colnames(sigactivities)
# Merge with trajectory information
sigactivities = merge(tj_samples, sigactivities %>% rownames_to_column("sample")) %>% tidyr::pivot_longer(cols=all_of(signatures), names_to='signature', values_to='activity')
sigactivities = sigactivities %>% dplyr::mutate(sigtype = case_when(
        str_detect(signature, "SBS") ~ "SBS",
        str_detect(signature, "ID") ~ "ID",
        str_detect(signature, "SV") ~ "SV",
        str_detect(signature, "CX") ~ "CN",
    )) 

# PLOT IMF activities across trajectories
canon_order <- get_sample_order(imf_activities, "Canonical", c("IMF2", "IMF3"))
imf_names <- sort(unique(imf_activities$imf_cluster))
imf_colours <- c(
  "IMF1" = "#B0B0B0", # Muted Gray (Background)
  "IMF2" = "#CC79A7", # Vivid Purple (Gain-enriched)
  "IMF3" = "#0072B2", # Deep Blue (Canonical)
  "IMF4" = "#DFDFDF", # Light Gray (Background)
  "IMF5" = "#D55E00", # Green (Alternative)
  "IMF6" = "#009E73", # Deep Vermilion/Red (Gain-enriched)
  "IMF7" = "#E5D8BD", # Soft Beige (Background)
  "IMF8" = "#E69F00"  # Vivid Orange (Canonical)
)

# Order IMF clusters by overall activity across all samples
global_imf_fill <- imf_activities %>%
  dplyr::group_by(imf_cluster) %>%
  dplyr::summarise(total_act = sum(activity, na.rm = TRUE)) %>%
  dplyr::arrange(desc(total_act)) %>% # Highest overall activity first
  dplyr::pull(imf_cluster)

imf_activities <- imf_activities %>%
  dplyr::mutate(imf_cluster = factor(imf_cluster, levels = global_imf_fill))


p_tj1 <- imf_activities %>% 
  dplyr::filter(trajectory == "Canonical") %>% 
  dplyr::mutate(sample = factor(sample, levels = canon_order)) %>%
  ggplot(aes(y = sample, x = activity, fill = imf_cluster)) +
  geom_bar(stat = "identity", width = 1) + # width = 1 removes white lines between bars
  scale_fill_manual(values = imf_colours, breaks = imf_names) + 
  scale_x_continuous(expand = c(0, 0)) + # Flushes bars directly against the axis
  labs(y = NULL, x = "IMF Activity") +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  coord_flip() +
  guides(fill = guide_legend(nrow = 1)) # Forces legend into a single, clean line

# Middle Row: Alternative (Order by IMF6)
alt_order <- get_sample_order(imf_activities, "Alternative", "IMF6")

p_tj2 <- imf_activities %>% 
  dplyr::filter(trajectory == "Alternative") %>% 
  dplyr::mutate(sample = factor(sample, levels = alt_order)) %>%
  ggplot(aes(y = sample, x = activity, fill = imf_cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = imf_colours, breaks = imf_names) + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL, x = NULL) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) + 
  coord_flip()

# Bottom Row: Gain-enriched (Order by IMF5 + IMF8)
gain_order <- get_sample_order(imf_activities, "Gain-enriched", c("IMF5", "IMF8"))

p_tj3 <- imf_activities %>% 
  dplyr::filter(trajectory == "Gain-enriched") %>% 
  dplyr::mutate(sample = factor(sample, levels = gain_order)) %>%
  ggplot(aes(y = sample, x = activity, fill = imf_cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = imf_colours, breaks = imf_names) + 
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = "none", 
    legend.title = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) + 
  coord_flip()

# --- 5. ASSEMBLE GRID ---
samples_prop <- tj_samples %>% 
  group_by(trajectory) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(prop = n / sum(n)) %>% 
  dplyr::pull(prop, name = trajectory)

p <- cowplot::plot_grid(
  p_tj1, p_tj2, p_tj3, 
  ncol = 3, 
  align = "h", 
  rel_widths = c(
    samples_prop["Canonical"], # Adding a slight buffer to top row for the legend
    samples_prop["Alternative"], 
    samples_prop["Gain-enriched"] # Adding buffer to bottom row for the x-axis text
  )
)

# Export at standard Nature double-column width (183 mm)
write_tsv(imf_activities, file.path(outdir, "Fig3b_source_data.tsv"))
ggsave(file.path(figdir, "Fig3b_Stacked_IMF_Activities.pdf"), p, width = 125/25.4, height = 70/25.4, units = "in", device = "pdf")

# BOXPLOT OF IMF activities across trajectories
p = imf_activities %>% ggplot(aes(x = trajectory, y = activity, fill = trajectory)) +
    geom_boxplot() +
    labs(x = "Trajectory", y = "IMF Activity") +
    scale_fill_manual(values = imf_colours) +
    theme(legend.position = "top") + 
    facet_wrap(~imf_cluster) +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("Canonical", "Alternative"), c("Canonical", "Gain-enriched"), c("Alternative", "Gain-enriched")))

ggsave(file.path(figdir, "imf_activities_boxplot.png"), p, width = 10, height = 12)

# GET TABLE WITH MEDIAN IMF ACTIVITIES AND WILCOXON P-VALUES FOR HIGHER VS THE OTHER TRAJECTORIES
imf_summary <- imf_activities %>%
    group_by(trajectory, imf_cluster) %>%
    summarise(median_activity = mean(activity, na.rm = TRUE)) %>%
    pivot_wider(names_from = trajectory, values_from = median_activity)

wilcox_results <- imf_activities %>%
    group_by(imf_cluster) %>%
    summarise(
        p_canon_vs_alt = wilcox.test(activity[trajectory == "Canonical"], activity[trajectory == "Alternative"])$p.value,
        p_canon_vs_gain = wilcox.test(activity[trajectory == "Canonical"], activity[trajectory == "Gain-enriched"])$p.value,
        p_alt_vs_gain = wilcox.test(activity[trajectory == "Alternative"], activity[trajectory == "Gain-enriched"])$p.value, 
        p_canon_vs_all = wilcox.test(activity[trajectory == "Canonical"], activity[trajectory != "Canonical"])$p.value,
        p_alt_vs_all = wilcox.test(activity[trajectory == "Alternative"], activity[trajectory != "Alternative"])$p.value,
        p_gain_vs_all = wilcox.test(activity[trajectory == "Gain-enriched"], activity[trajectory != "Gain-enriched"])$p.value
    )

imf_summary <- imf_summary %>%
    left_join(wilcox_results, by = "imf_cluster") %>% 
    dplyr::select(imf_cluster, Canonical, Alternative, Gain_enriched = `Gain-enriched`, p_canon_vs_alt, p_canon_vs_gain, p_alt_vs_gain, p_canon_vs_all, p_alt_vs_all, p_gain_vs_all) %>% 
    dplyr::arrange(as.character(imf_cluster))

write_delim(imf_summary, file.path(outdir, "ST_IMF_activities_summary.csv"), delim = ",")

# PLOT INDIVIDUAL SIGNATURE ACTIVITIES
traj_orders <- list(
  Canonical     = canon_order,
  Alternative   = alt_order,
  `Gain-enriched` = gain_order
)

sig_colours <- make_sig_colours(sigactivities)

sig_order <- sigactivities %>%
  group_by(signature) %>%
  summarise(total = sum(activity, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  pull(signature)

sigactivities <- sigactivities %>%
  mutate(signature = factor(signature, levels = sig_order))

samples_prop <- bind_rows(
  tibble(trajectory = "Canonical", sample = canon_order),
  tibble(trajectory = "Alternative", sample = alt_order),
  tibble(trajectory = "Gain-enriched", sample = gain_order)
) %>%
  count(trajectory, name = "n") %>%
  mutate(prop = n / sum(n)) %>%
  pull(prop, name = trajectory)

sigtype_levels <- c("SBS", "ID", "SV", "CN")

sig_levels_by_type <- sigactivities %>%
  group_by(sigtype, signature) %>%
  summarise(total = sum(activity), .groups = "drop") %>%
  arrange(sigtype, desc(total)) %>%
  split(.$sigtype) %>%
  lapply(function(df) df$signature)

sig_rows <- lapply(seq_along(sigtype_levels), function(i) {
  st <- sigtype_levels[i]
  if (!st %in% names(sig_levels_by_type)) return(NULL)
  
  make_sigtype_row(
    sigactivities,
    st,
    traj_orders,
    sig_levels_by_type,
    sig_colours,
    samples_prop,
    show_titles = (i == 1)
  )
})

sig_rows <- sig_rows[!sapply(sig_rows, is.null)]

p_sig <- cowplot::plot_grid(
  plotlist = sig_rows,
  ncol = 1,
  align = "v"
)

write_tsv(sigactivities, file.path(outdir, "EXDF6_source_data.tsv"))
ggsave(
  file.path(figdir, "EXDF6_Stacked_Signature_Activities_by_Trajectory.pdf"),
  p_sig,
  width = 180/25.4,
  height = 250/25.4,
  units = "in",
  device = cairo_pdf
)

# GET TABLE WITH MEDIAN SIGNATURE ACTIVITIES AND WILCOXON P-VALUES FOR HIGHER VS THE OTHER TRAJECTORIES
sig_summary <- sigactivities %>% 
    group_by(trajectory, signature, sigtype) %>%
    summarise(median_activity = mean(activity, na.rm = TRUE)) %>%
    pivot_wider(names_from = trajectory, values_from = median_activity)

wilcox_results_sig <- sigactivities %>%
    group_by(signature) %>%
    summarise(
        p_canon_vs_alt = wilcox.test(activity[trajectory == "Canonical"], activity[trajectory == "Alternative"])$p.value,
        p_canon_vs_gain = wilcox.test(activity[trajectory == "Canonical"], activity[trajectory == "Gain-enriched"])$p.value,
        p_alt_vs_gain = wilcox.test(activity[trajectory == "Alternative"], activity[trajectory == "Gain-enriched"])$p.value
    )

sig_summary <- sig_summary %>%
    left_join(wilcox_results_sig, by = "signature") %>% 
    dplyr::select(signature, sigtype, Canonical, Alternative, Gain_enriched = `Gain-enriched`, p_canon_vs_alt, p_canon_vs_gain, p_alt_vs_gain) %>% 
    dplyr::arrange(sigtype)

write_delim(sig_summary, file.path(outdir, "ST_signature_activities_summary.csv"), delim = ",")


# ANALYSIS OF SIGNATURE ACTIVITIES IN CLONAL VS SUBCLONAL PERIOD

# A. SBS SIGNATURES
timed_sbs = foreach(fp = sbs_timing_fps, .combine = rbind) %dopar% {
  load_sigtimer(fp, tj_samples, "SBS", sig_names = unique(sigactivities %>% dplyr::filter(sigtype == "SBS") %>% pull(signature) %>% as.character()))
}

sbs_signatures = rbind(
  sigactivities %>% dplyr::filter(sigtype == "SBS") %>% dplyr::select(sample, trajectory, signature, activity) %>% dplyr::mutate(timing = "overall"),
  timed_sbs %>% dplyr::select(sample, trajectory, signature, activity, timing) 
) %>% dplyr::filter(signature %in% unique(sigactivities %>% dplyr::filter(sigtype == "SBS") %>% pull(signature) %>% as.character())) 

p = sbs_signatures %>% 
  mutate(timing = factor(timing, levels = c("overall", "clonal", "subclonal"))) %>% 
  group_by(timing, signature) %>%
  # Filter out signatures that don't have at least 2 trajectories in that facet
  dplyr::filter(
    n_distinct(trajectory) >= 2, 
    sd(activity, na.rm = TRUE) > 0,
    max(activity, na.rm = TRUE) > 0
  ) %>%
ggplot(aes(x = signature, y = activity, fill = trajectory)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  ggpubr::stat_compare_means( 
    aes(group = trajectory), 
    label = "p.signif", 
    method = "kruskal.test",
    size = 2, hide.ns = TRUE
  ) +
  facet_wrap(~timing) +
  labs(x = "Activity", y = "Signature") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) + 
  scale_fill_manual(values = trajectory_colours) 

p = axes2lemon(p, lt = "both", bt = "both")

write_tsv(
  sbs_signatures %>% dplyr::select(timing, trajectory, signature, activity),
  file.path(outdir, "EXDF7a_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF7a_sbs_signature_activities_clonal_subclonal.pdf"), width = 180/25.3, height = 50/25.3)


# B. ID SIGNATURES
timed_id = foreach(fp = id_timing_fps, .combine = rbind) %dopar% {
  load_sigtimer(fp, tj_samples, "ID", sig_names = unique(sigactivities %>% dplyr::filter(sigtype == "ID") %>% pull(signature) %>% as.character()))
}

id_signatures = rbind(
  sigactivities %>% dplyr::filter(sigtype == "ID") %>% dplyr::select(sample, trajectory, signature, activity) %>% dplyr::mutate(timing = "overall"),
  timed_id %>% dplyr::select(sample, trajectory, signature, activity, timing) 
) %>% dplyr::filter(signature %in% unique(sigactivities %>% dplyr::filter(sigtype == "ID") %>% pull(signature) %>% as.character()))


p = id_signatures %>% 
  mutate(timing = factor(timing, levels = c("overall", "clonal", "subclonal"))) %>% 
  group_by(timing, signature) %>%
  # Filter out signatures that don't have at least 2 trajectories in that facet
  dplyr::filter(
    n_distinct(trajectory) >= 2, 
    sd(activity, na.rm = TRUE) > 0,
    max(activity, na.rm = TRUE) > 0
  ) %>%
ggplot(aes(x = signature, y = activity, fill = trajectory)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  ggpubr::stat_compare_means( 
    aes(group = trajectory), 
    label = "p.signif", 
    method = "kruskal.test",
    size = 2, hide.ns = TRUE
  ) +
  facet_wrap(~timing) +
  labs(x = "Activity", y = "Signature") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) + 
  scale_fill_manual(values = trajectory_colours) 

p = axes2lemon(p, lt = "both", bt = "both")

write_tsv(
  id_signatures %>% dplyr::select(timing, trajectory, signature, activity),
  file.path(outdir, "EXDF7b_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF7b_id_signature_activities_clonal_subclonal.pdf"), width = 180/25.3, height = 50/25.3)
