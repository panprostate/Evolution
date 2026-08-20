# Estimate the "speed" of tumour evolution by dividing number of alterations
# by period spanned within each time period (i.e. clonal or subclonal)
# Schematic:
# Z-----------------------------M------------S---D
# where Z is zygote, M is MRCA, S is subclone and D is diagnosis.
# Time between Z and M is reflected by number of clock clonal mutations
# Time between M and S is reflected by number of clock subclonal mutations 
# Time between S and D can't be estimated as mutations in this period have low VAF
# By normalising number of drivers in trajectory, number of SCNAs and 
# number of non-clock mutations in those time periods, we can estimate time-averaged "speed" of evolution

rm(list = ls(all = TRUE))

# PACKAGES
library(data.table)
library(tidyverse)
library(doMC)
library(foreach)
library(DescTools)
library(lemon)
library(cowplot)

# FUNCTIONS
source("src/utils.R")
source("src/plot_theme.R")
source("src/plot_functions.R")
source("src/trajectory_drivers.R")

boxplot_rates <- function(df, y_str, ylab){
    dodge = position_dodge2(width = 0.2)
    trend_stats <- df %>%
        dplyr::filter(!is.na(gs_grade)) %>%
        mutate(gs_ordered = factor(gs_grade, ordered = TRUE)) %>%
        group_by(trajectory) %>%
        summarise(
        # Compute the Jonckheere-Terpstra test using DescTools
        p_val = DescTools::JonckheereTerpstraTest(
            x = .data[[y_str]], 
            g = gs_ordered, 
            alternative = "increasing" 
        )$p.value,
        .groups = 'drop'
        ) %>%
        mutate(
        # Format P-values for Nature style
        p_label = ifelse(p_val < 0.001, "P < 0.001", sprintf("P = %.3f", p_val)),
        # Dynamic y-position for the plot
        y_pos = max(df[[y_str]], na.rm = TRUE) * 1.5 
        )

    box_dodge <- position_dodge(width = 0.75)
    point_dodge <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)

    p = ggplot(df %>% dplyr::filter(!is.na(gs_grade)), aes_string(y = y_str, fill = "gs_grade", x = "trajectory")) + 
        geom_point(pch = 21, alpha = 0.2, position = point_dodge) + 
        geom_boxplot(outlier.shape = NA, position = box_dodge, alpha = 0.5) + 
        scale_y_log10() + 
        # from plot_theme.R
        scale_fill_manual(values = gleason_colours) + 
        # Add exact trend P-values calculated previously via JT
        geom_text(data = trend_stats, 
                    aes(x = trajectory, y = y_pos, label = p_label), 
                    inherit.aes = FALSE, 
                    size = 2, 
                    fontface = "italic") + 
        labs(x = "Trajectory", y = ylab, fill = "GS Grade") +
        scale_x_discrete(labels = c("Canonical", "Alternative", "Gain-enriched")) + 
        theme(legend.position = "none")
    
    p = axes2lemon(p, bt = "both", lt = "both")

    return(p)
}

# PATHS
outdir = "outputs/04_progression/rate"
figdir = "figures/04_progression/rate"
dir.create(outdir); dir.create(figdir)
evo_metrics_fp <- "outputs/00_preprocessing/evo_metrics.tsv"
trajectories_fps <- list.files("outputs/02_trajectories", pattern = "PPCG_Feb2026.*mergedseg_with_clonality.txt", full = T)
mrca_time_fp <- "outputs/01_landscape/real_timing/latency_estimates_mrca_patient_rate_summary.tsv"

# LOAD DATA
tj_samples = load_trajectories(trajectories_fps)
lookup_tj = c("Ordering 1" = "Canonical", "Ordering 2" = "Alternative", "Ordering 3" = "Gain-enriched")
tj_samples$trajectory = lookup_tj[tj_samples$trajectory]
mrca_time = read_delim(mrca_time_fp, delim = "\t") %>% dplyr::filter(acceleration == "1x") %>% dplyr::select(smp = sample_id, latency)
evo_metrics = read_delim(evo_metrics_fp, delim = "\t") %>% dplyr::mutate(smp = paste0(smp, "_DNA"))
evo_metrics = merge(evo_metrics, mrca_time) %>% dplyr::rename(time_subclonal = latency)
evo_metrics$record_id = extract_ppcg_pt(evo_metrics$smp)

## ANALYSIS
# 1. Get real-time in clonal and subclonal periods
evo_metrics$age = get_age(evo_metrics$smp)
evo_metrics$time_clonal = evo_metrics$age - evo_metrics$time_subclonal

# 2. Gather number of Mbs altered per year
evo_metrics$clonal_pga_per_year = evo_metrics$clonal_pga / evo_metrics$time_clonal
evo_metrics$subclonal_pga_per_year = evo_metrics$subclonal_pga / evo_metrics$time_subclonal

# 3. Gather number of trajectory-specific drivers per year
evo_metrics = merge(evo_metrics, tj_samples %>% dplyr::rename(smp = sample))
trajectory_drivers = get_trajectory_drivers(trajectories_fps, trajectories_fps) %>% rownames_to_column("smp") %>% dplyr::mutate(smp = extract_ppcg_id(smp, full = F))
evo_metrics = merge(evo_metrics, trajectory_drivers, by = "smp")

# calculate rates per year
evo_metrics$subclonal_ndrivers_per_year = evo_metrics$subclonal_events / evo_metrics$time_subclonal
evo_metrics$clonal_ndrivers_per_year = evo_metrics$clonal_events / evo_metrics$time_clonal
evo_metrics$clonal_poteo_per_year = evo_metrics$clonal_poteo / evo_metrics$time_clonal
evo_metrics$subclonal_poteo_per_year = evo_metrics$subclonal_poteo / evo_metrics$time_subclonal

# 4. Plot yearly rates per GS
evo_metrics$gs_grade = get_gs_group(evo_metrics$smp)
write_delim(evo_metrics, file.path(outdir, "evo_metrics_rates.tsv"), delim = "\t")
p1 = boxplot_rates(evo_metrics %>% dplyr::filter(nsubclones > 0), "clonal_pga_per_year", "Clonal PGA per Year")
p2 = boxplot_rates(evo_metrics %>% dplyr::filter(nsubclones > 0), "subclonal_pga_per_year", "Subclonal PGA per Year")
p3 = boxplot_rates(evo_metrics %>% dplyr::filter(nsubclones > 0), "clonal_poteo_per_year", "Yearly clonal rate \nof trajectory progression")
p4 = boxplot_rates(evo_metrics %>% dplyr::filter(nsubclones > 0), "subclonal_ndrivers_per_year", "Yearly subclonal rate \nof trajectory progression")

cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = c("a)", "b)", "c)", "d)"))
write_tsv(
  evo_metrics %>% dplyr::filter(nsubclones > 0) %>%
    dplyr::select(trajectory, gs_grade,
                  clonal_pga_per_year, subclonal_pga_per_year,
                  clonal_poteo_per_year, subclonal_ndrivers_per_year),
  file.path(outdir, "SF14_source_data.tsv")
)
ggsave(file.path(figdir, "SF14_evolutionary_rates_by_gs.pdf"), width = 130/25.4, height = 130/25.4, units = "in")
