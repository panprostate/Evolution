# Compare the prevalence of evolutionary subtypes across cohorts and degree of trajectory progression

rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(lemon)
library(ggsankey)

# FUNCTIONS
source("src/trajectory_drivers.R")
source("src/plot_theme.R")
source("src/plot_functions.R")
source("src/utils.R")

# PATHS
ppcg_trajectory_fps <- list.files("outputs/02_trajectories", pattern = "PPCG_Feb2026.*mergedseg_with_clonality.txt", full = T)
hmf_trajectory_fps <-  list.files("outputs/04_progression/", pattern = "HMF_PPCG.*mergedseg_with_clonality.txt", full = T)
hmf_trajectory_fps <- c(hmf_trajectory_fps[grepl("canonical", hmf_trajectory_fps)], 
                           hmf_trajectory_fps[grepl("alternative", hmf_trajectory_fps)],
                           hmf_trajectory_fps[grepl("gain_enriched", hmf_trajectory_fps)])
combimets_trajectory_fps <- list.files("outputs/04_progression/combimets_trajectories", pattern = "CombiMets.*mergedseg_with_clonality.txt", full = T)
combimets_trajectory_fps <- c(combimets_trajectory_fps[grepl("canonical", combimets_trajectory_fps)], 
                        combimets_trajectory_fps[grepl("alternative", combimets_trajectory_fps)],
                        combimets_trajectory_fps[grepl("gain_enriched", combimets_trajectory_fps)])
combimets_pairs_fp <- "outputs/04_progression/combimets_trajectories/combimets_pairs_trajectories.tsv"

ppcg_cna_paths <- list.files("data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/Subclonal_SCNA", pattern = "Subclonal_SCNA.txt", full.names = T, recursive = T)
hmf_cna_paths <- list.files("data/controlled/hmf/bb", pattern = ".tsv", full.names = T, recursive = T)


outdir <- "outputs/04_progression/trajectory_progression"
figdir <- "figures/04_progression/trajectory_progression"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# Load trajectories
ppcg_tjs = load_trajectories(ppcg_trajectory_fps)
hmf_tjs = load_trajectories(hmf_trajectory_fps, parse_ppcg_id = F)
combi_tjs = load_trajectories(combimets_trajectory_fps, parse_ppcg_id = T)
# for this analysis, we only want to include combimets metastatic samples
combi_tjs = combi_tjs[!combi_tjs$sample %in% get_met_smps(qc=F),]

all_tjs = rbind(
    data.frame(sample = ppcg_tjs$sample, trajectory = ppcg_tjs$trajectory, cohort = "PPCG"),
    data.frame(sample = hmf_tjs$sample, trajectory = hmf_tjs$trajectory, cohort = "HMF"),
    data.frame(sample = combi_tjs$sample, trajectory = combi_tjs$trajectory, cohort = "CombiMets")
)
all_tjs$cohort <- factor(all_tjs$cohort, levels = c("PPCG", "CombiMets", "HMF"))

fisher.test(table(all_tjs$cohort != "PPCG", all_tjs$trajectory == "Ordering 3"))

# Dougnut plot of trajectory distribution by cohort
plot_data <- all_tjs %>%
  group_by(cohort, trajectory) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(cohort) %>%
  mutate(
    total = sum(count),
    prop = count / total
  )

center_data <- plot_data %>%
  dplyr::select(cohort, total) %>%
  distinct() %>%
  mutate(label_text = paste0("N=\n", total))

p <- ggplot(plot_data, aes(x = 2, y = prop, fill = trajectory)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) + 
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) + 
  facet_wrap(~ cohort) +
  geom_text(data = center_data, aes(x = 0.5, y = 0, label = label_text),
            inherit.aes = FALSE, size = 2.5, fontface = "bold", color = "#333333") +
  scale_fill_manual(values = trajectory_colors) +
  theme_void() + 
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 7, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

write_tsv(all_tjs %>% dplyr::select(trajectory, cohort), file.path(outdir, "Fig4a_source_data.tsv"))
ggsave(file.path(figdir, "Fig4a_trajectory_distribution_cohort.pdf"), plot = p, width = 60/25.3, height = 50/25.4, units = "in", dpi = 300)

# CombiMets trajectory transitions
pair_tj = read_delim(combimets_pairs_fp, delim = "\t")
sankey_data <- pair_tj %>%
  make_long(s1_tj, s2_tj)

p = ggplot(sankey_data,
       aes(x = x, 
           next_x = next_x,
           node = node, 
           next_node = next_node,
           fill = factor(node))) +

  geom_sankey(flow.alpha = 0.5, node.color = "black") +

  geom_sankey_label(
    aes(label = node),
    size = 3,
    color = "white",
    data = sankey_data %>% filter(x %in% c("s1_tj", "s2_tj"))
  ) +

  theme_minimal() +
  labs(
    title = "Trajectory changes from primary tumour to matched metastasis",
    x = "Sample Timepoint",
    y = "",
    fill = "Trajectory"
  ) + 
  scale_x_discrete(
    labels = c("s1_tj" = "Primary", "s2_tj" = "Metastasis")
  ) + 
  scale_fill_manual(values = 
  c("1" = "#8DA0CB", 
  "2" = "#AAF0C9", 
  "3" = "#C04667")) +
  theme(legend.position = "none")

write_tsv(pair_tj %>% dplyr::select(s1_tj, s2_tj), file.path(outdir, "SF13_source_data.tsv"))
ggsave(file.path(figdir, "SF13_tj_transitions_sankey.pdf"), plot = p, width = 150/25.3, height = 100/25.3)

# Calculate progression
ppcg_progression = get_trajectory_drivers(ppcg_trajectory_fps, ppcg_trajectory_fps) %>% dplyr::mutate(sample = extract_ppcg_id(Tumour_Name, full = F))
ppcg_progression$gs_grade = get_gs_group(ppcg_progression$sample)
hmf_progression = get_trajectory_drivers(hmf_trajectory_fps, ppcg_trajectory_fps) %>% dplyr::mutate(sample = Tumour_Name)
hmf_progression$gs_grade = "Metastasis"
combi_progression = get_trajectory_drivers(combimets_trajectory_fps, ppcg_trajectory_fps) %>% dplyr::mutate(sample = Tumour_Name)
combi_progression$gs_grade = "Metastasis"

all_tjs = merge(all_tjs, rbind(ppcg_progression, hmf_progression, combi_progression), by = "sample")

dodge = position_dodge2(width = 0.5)
gleason_colours = c(gleason_colours, "Metastasis" = "#67000D")

# Calculate Jonckheere-Terpstra test for trend across GS grades within each trajectory
metric = "total_poete"
trend_stats <- all_tjs %>%
    dplyr::filter(!is.na(gs_grade)) %>%
    mutate(gs_ordered = factor(gs_grade, ordered = TRUE)) %>%
    group_by(trajectory) %>%
    summarise(
    # Compute the Jonckheere-Terpstra test using DescTools
    p_val = DescTools::JonckheereTerpstraTest(
        x = .data[[metric]], 
        g = gs_ordered, 
        alternative = "increasing" 
    )$p.value,
    .groups = 'drop'
    ) %>%
    mutate(
    # Format P-values for Nature style
    p_label = ifelse(p_val < 0.001, "P < 0.001", sprintf("P = %.3f", p_val)),
    # Dynamic y-position for the plot
    x_pos = max(all_tjs[[metric]], na.rm = TRUE) * 0.9
    )

box_dodge <- position_dodge(width = 0.75)
point_dodge <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)

p = all_tjs %>% dplyr::filter(!is.na(gs_grade)) %>% 
    ggplot(aes_string(y = "trajectory", x = metric, fill = "gs_grade")) + 
    geom_point(pch = 21, alpha = 0.2, position = point_dodge) + 
    geom_boxplot(outlier.shape = NA, position = box_dodge, alpha = 0.5) + 
    # scale_x_log10(limits = c(0.01, 1)) +
    geom_text(data = trend_stats, 
                aes(x = x_pos, y = trajectory, label = p_label), 
                inherit.aes = FALSE, 
                size = 2, angle = 270, 
                fontface = "italic") +
    scale_fill_manual(values = gleason_colours) + 
    theme(legend.position = "top") + 
    labs(x = "Relative degree of progression along trajectory", y = "Trajectory")

p = axes2lemon(p, lt = "both", bt = "both")

write_tsv(
  all_tjs %>% dplyr::filter(!is.na(gs_grade)) %>% dplyr::select(trajectory, total_poete, gs_grade),
  file.path(outdir, "Fig4b_source_data.tsv")
)
ggsave(file.path(figdir, "Fig4b_trajectory_progression_by_gs.pdf"), width = 90/25.3, height = 80/25.4, dpi = 300)

all_tjs = all_tjs %>% dplyr::mutate(ggnum = ifelse(gs_grade == "Metastasis", 5, as.numeric(gsub("4\\+", "4", gs_grade))))

table(all_tjs$ggnum)


all_tjs$patient = str_remove(all_tjs$sample, "[a-z]_DNA")

# Check that the p-value is under 0.001 when using LME and accounting for patients with matched primary and met
summary(nlme::lme(total_poete ~ ggnum, random = ~1|patient, data = all_tjs %>% dplyr::filter(!is.na(ggnum) & trajectory == "Ordering 1"))) # < 0.0001
summary(nlme::lme(total_poete ~ ggnum, random = ~1|patient, data = all_tjs %>% dplyr::filter(!is.na(ggnum) & trajectory == "Ordering 2"))) # < 0.0001
summary(nlme::lme(total_poete ~ ggnum, random = ~1|patient, data = all_tjs %>% dplyr::filter(!is.na(ggnum) & trajectory == "Ordering 3"))) # < 0.0001


# Check that p-value is significant when excluding combimets with Jonckheere-Terpstra test
all_tjs_no_combi = all_tjs %>% dplyr::filter(cohort != "CombiMets")
trend_stats_no_combi <- all_tjs_no_combi %>%
    dplyr::filter(!is.na(gs_grade)) %>%
    mutate(gs_ordered = factor(gs_grade, ordered = TRUE)) %>%
    group_by(trajectory) %>%
    summarise(
    # Compute the Jonckheere-Terpstra test using DescTools
    p_val = DescTools::JonckheereTerpstraTest(
        x = .data[[metric]], 
        g = gs_ordered, 
        alternative = "increasing" 
    )$p.value,
    .groups = 'drop'
    ) 


write_delim(all_tjs, file.path(outdir, "trajectory_progression_data.tsv"), delim = "\t")
