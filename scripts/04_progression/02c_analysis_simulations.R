# Analysis results of simulated trajectories of tumour evolution after diagnosis
rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(lemon)
library(survival)
library(survminer)
library(cowplot)

# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

plot_model <- function(simulations_df, model_name, gleason_colours){
        
    # just create additional repeated rows for each grade group so that we can plot
    met_summ = met_summary_data %>% group_by(trajectory) %>% 
        tidyr::expand_grid(gg = c("1", "2", "3", "4-5")) %>% 
        ungroup() %>% 
        group_by(gg) %>% 
        mutate(
            met_first_quartile = met_first_quartile,
            met_median = met_median,
            met_third_quartile = met_third_quartile
        )

    names(gleason_colours) = paste0("Grade Group ", names(gleason_colours))
    simulations_df$gg = paste0("Grade Group ", simulations_df$gg)
    p = simulations_df %>% dplyr::filter(model == model_name) %>% dplyr::group_by(gg, year) %>% 
        dplyr::summarise(
            first_quartile = quantile(poete, 0.25),
            median = median(poete),
            mean = mean(poete),
            third_quartile = quantile(poete, 0.75)
        ) %>% 
        ggplot() +
        geom_step(aes(x = year, y = median, color = gg), size = .5) +
        geom_ribbon(aes(x = year, ymin = first_quartile, ymax = third_quartile, fill = gg), alpha = 0.1) +
        geom_line(data = simulations_df %>% dplyr::filter(model == model_name) %>% group_by(year, sample, gg) %>% summarise(poete = mean(poete)), aes(x = year, y = poete, group = sample, color = gg), size = 0.1, alpha = 0.1) +
        geom_hline(yintercept = median(met_progression_data$total_poete), color = "firebrick", linetype = "dashed", size = 0.5, alpha = 1) +
        labs(x = "Years", y = "Degree of trajectory progression") + 
        scale_color_manual(values = gleason_colours, name = "Grade Group") +
        scale_fill_manual(values = gleason_colours, name = "Grade Group") +
        xlim(-5, 15) + 
        geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
        # facet_wrap(~trajectory, scales = "free_y") +
        facet_wrap(~gg, scales = "free_y", nrow = 1, ncol = 4) +
        scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) + 
        theme(legend.position = "none") 
    
    return(p)
}

# PATHS
outdir = "outputs/04_progression/simulation/"
figdir = "figures/04_progression/simulation/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

clinfile <- "data/meta/PPCG_donors_clin_20241217.csv"
simulations_fp = file.path(outdir, "trajectory_drivers_simulation_results.tsv")
trajectory_progression_fp = file.path("outputs/04_progression/trajectory_progression", "trajectory_progression_data.tsv")
trajectories_fps = list.files("outputs/02_trajectories", pattern = "PPCG_Feb2026.*mergedseg_with_clonality.txt", full.names = T)
evo_metrics_fp = file.path("outputs/04_progression/rate", "evo_metrics_rates.tsv")

# LOAD DATA
evo_metrics = read_delim(evo_metrics_fp, delim = "\t")
simulations_df = read_delim(simulations_fp, delim = "\t")
trajectory_progression_data = read_delim(trajectory_progression_fp, delim = "\t")

# Calculate POETE for each sample in simulations_df
ms_list = lapply(trajectories_fps, function(x) read.table(x, header = T))
ms_list = lapply(ms_list, function(ms) {
    cbind(ms, `cna_id`=apply(ms, 1, function(ms_row) { paste(ms_row["CNA"], ms_row["ID"], sep='_') }))
})
expected_num_of_events = sapply(ms_list, function(ms) {
    # table of cna_ids gets number of samples to have each events_per_sample
    # dividing this by total number of samples gets frequency of each event
    # then we get sum of these frequencies to calculate
    # the expected number of trajectory events per sample
    sum(table(ms$cna_id) / length(unique(ms$Tumour_Name)))
})
names(expected_num_of_events) = c("Canonical", "Alternative", "Gain-enriched")

simulations_df$poete = simulations_df$total_events / expected_num_of_events[simulations_df$trajectory]

# Add landmark metastasis progression metrics for each trajectory
met_progression_data = trajectory_progression_data %>% dplyr::filter(cohort == "HMF") 
met_summary_data = met_progression_data %>% group_by(trajectory) %>% 
    dplyr::summarise(
        met_first_quartile = quantile(total_poete, .25), 
        met_median = median(total_poete),
        met_third_quartile = quantile(total_poete, .75)
    )

sid = "PPCG0006a_DNA"
example_df = simulations_df %>% dplyr::filter(sample == sid, model != "accelerating")
example_summary = example_df %>% group_by(year, model) %>% summarise(
    median_events = median(total_events),
    q25 = quantile(total_events, 0.25),
    q75 = quantile(total_events, 0.75)
)

# Plot just one simulation (last one)
p = ggplot(example_df %>% dplyr::filter(sim == 100), aes(x = year, y = total_events, color = model)) + 
    geom_step(size = 1) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    labs(
        title = paste(sid, " - Simulation 1"),
        x = "Years from Diagnosis",
        y = "Total Accumulated Events"
    ) +
    theme(
        legend.position = "top",
        panel.grid.major.y = element_line(color = "grey90", linetype = "dotted")
    ) + 
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1"), name = "Model")

p = axes2lemon(p, lt = "both", bt = "both")

write_tsv(
  example_df %>% dplyr::filter(sim == 100) %>% dplyr::select(sample, year, total_events, model),
  file.path(outdir, "Fig4d_source_data.tsv")
)
ggsave(file.path(figdir, "Fig4d_trajectory_drivers_simulation_example_case_sim1.pdf"), plot = p, width = 100/25.3, height = 50/25.4)

# Plot median and 25%-75% range for simulations across all GGs
# gleason_colours from plot_theme.R
names(gleason_colours) = c("1", "2", "3", "4-5")

# Main plot: uniform model
p = plot_model(simulations_df, "uniform", gleason_colours)
write_tsv(
  simulations_df %>% dplyr::filter(model == "uniform") %>% dplyr::select(trajectory, sample, gg, year, model, poete),
  file.path(outdir, "Fig4e_source_data.tsv")
)
ggsave(file.path(figdir, "Fig4e_simulated_progression_gg.pdf"), width = 160/25.3, height = 50/25.3, dpi = 300)

# Supplementary: punctuated and gradually accelerating models
p1 = plot_model(simulations_df, "punctuated", gleason_colours)
p2 = plot_model(simulations_df, "gradual", gleason_colours)

# Combine plots
p_combined = cowplot::plot_grid(p1, p2, nrow = 2, labels = c("a)", "b)"))
write_tsv(
  simulations_df %>% dplyr::filter(model %in% c("punctuated", "gradual")) %>% dplyr::select(trajectory, sample, gg, year, model, poete),
  file.path(outdir, "EXDF8_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF8_simulated_progression_gg_alternative_models.pdf"), width = 160/25.3, height = 100/25.3, dpi = 300)

# Estimate time at which trajectory progression reaches median metastatic progression level in each GG
sims_median = simulations_df %>% dplyr::group_by(model, gg, year) %>% dplyr::summarise(poete = median(poete, na.rm = T))

sims_median$met_poete = median(met_progression_data$total_poete)
sims_median$higher_met = sims_median$poete >= sims_median$met_poete
times_to_met_all = sims_median %>% dplyr::group_by(model, gg) %>% 
    dplyr::summarise(
        median_time_to_met = min(year[higher_met])
    ) 

# Calculate across all samples
sample_sims_df = simulations_df %>% dplyr::group_by(model, gg, trajectory, sample, year) %>% dplyr::summarise(poete = median(poete, na.rm = T))
dic_median = met_summary_data$met_median
names(dic_median) = c("Canonical", "Alternative", "Gain-enriched")
sample_sims_df$met_poete = dic_median[sample_sims_df$trajectory]
sample_sims_df$higher_met = sample_sims_df$poete >= sample_sims_df$met_poete
times_to_met_sample = sample_sims_df %>% dplyr::group_by(model, gg, trajectory, sample) %>% 
    dplyr::summarise(
        median_time_to_met = min(year[higher_met])
    )

# Estimate how many samples reach metastatic progression levels within final time of simulations (15 years post diagnosis)
times_to_met_sample %>% group_by(model, gg) %>% dplyr::summarise(mean(median_time_to_met < 15.1))


##### SURVIVAL ANALYSIS BY RATE
imminent_samples = times_to_met_sample %>% dplyr::filter(model == "uniform", median_time_to_met < 10) %>% dplyr::pull(sample)
no_imminent_samples = times_to_met_sample %>% dplyr::filter(model == "uniform", median_time_to_met >= 10) %>% dplyr::pull(sample)
met_relapse <- get_relapse_data(clinfile, evo_metrics$record_id, type = "metastasis") %>% dplyr::rename(age_sample = age) %>% merge(evo_metrics, by = "record_id") %>% dplyr::rename(sample = smp) 
met_relapse$record_id <- extract_ppcg_pt(met_relapse$sample)

met_relapse$imminent = ifelse(met_relapse$sample %in% imminent_samples, "Imminent", ifelse(met_relapse$sample %in% no_imminent_samples, "No Imminent", NA))
met_relapse$median_time_to_met = pmin(15, times_to_met_sample$median_time_to_met[match(met_relapse$sample, times_to_met_sample$sample)])
met_relapse$gs_group <- get_gs_group(met_relapse$record_id)
met_relapse$t_stage = str_remove(met_relapse$t_stage, "[a-z]$")



p1 <- ggforest(coxph(Surv(time2relapse, relapse_ind) ~ imminent, data = met_relapse[met_relapse$gs_group == "2" ,]), main = "")
p2 <- ggforest(coxph(Surv(time2relapse, relapse_ind) ~ imminent + gs_group + psa + age + t_stage, data = met_relapse[met_relapse$gs_group != "1",]), main = "")
p3 <- ggforest(coxph(Surv(time2relapse, relapse_ind) ~ median_time_to_met + gs_group + psa + t_stage, data = met_relapse[met_relapse$gs_group != "1",]), main = "")
combined_forest <- cowplot::plot_grid(
  p3, p2, p1, 
  rel_heights = c(5, 5, 2), # Adjust heights to give more space to the top plot
  ncol = 1,                 # Stack vertically to preserve text space
  labels = c("a)", "b)", "c)"),   # Nature-style panel labels
  label_size = 9,
  align = "v",              # Vertically align the plot panels
  axis = "l"                # Align by the left axis
)

write_tsv(
  met_relapse %>% dplyr::select(time2relapse, relapse_ind, imminent, median_time_to_met, gs_group, psa, age, t_stage),
  file.path(outdir, "SF15_source_data.tsv")
)
pdf(file.path(figdir, "SF15_rate_forest_plots.pdf"), width = 120/25.3, height = 180/25.3)
print(combined_forest)
dev.off()

