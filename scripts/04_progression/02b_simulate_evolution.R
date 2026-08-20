# Simulate trajectory evolution after diagnosis
rm(list = ls(all = TRUE))

# PACKAGES
library(foreach)
library(doMC)
library(data.table)
library(tidyverse)

# FUNCTIONS
source("src/utils.R")
source("src/helpers_simulate_evolution.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

# PATHS
outdir = "outputs/04_progression/simulation"
figdir = "figures/04_progression/simulation"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
evo_metrics_fp = file.path("outputs/04_progression/rate", "evo_metrics_rates.tsv")

# LOAD DATA
evo_metrics = read_delim(evo_metrics_fp, delim = "\t")

# evo_metrics$gs_grade = ifelse(
#     evo_metrics$gs_grade %in% c("1", "2"), "1-2",
#     ifelse(evo_metrics$gs_grade == "4+", "4-5",
#            ifelse(evo_metrics$gs_grade == "3", "3", NA))
# )
evo_metrics$gs_grade = ifelse(evo_metrics$gs_grade == "4+", "4-5", evo_metrics$gs_grade)

# PARAMETERS FOR SIMULATION
start_latency = 5
future_horizon = 15
dt = 0.1
# Simulate for samples where observed prolonged latency to go back in time
evo_metrics = evo_metrics %>% dplyr::filter(time_subclonal >= start_latency)

set.seed(42)
n_sims = 100
registerDoMC(cores = 20)

models_to_run <- c("uniform", "gradual", "punctuated", "accelerating")

simulations_df <- foreach(i = 1:nrow(evo_metrics), .combine = rbind) %dopar% {

    sample_id <- evo_metrics$smp[i]
    gg <- evo_metrics$gs_grade[i]
    trajectory <- evo_metrics$trajectory[i]

    clonal_events <- evo_metrics$clonal_events[i]
    subclonal_events <- evo_metrics$subclonal_events[i]
    time_subclonal <- evo_metrics$time_subclonal[i]

    if (any(is.na(c(gg, trajectory, clonal_events, subclonal_events, time_subclonal)))) {
        return(NULL)
    }

    sim_list <- list()

    for (sim in seq_len(n_sims)) {
        model_param_list <- list(
            uniform = list(),
            gradual = list(shape = 1.2), # < 1 means slowing down over time
            punctuated = list(
                # Randomly assign between 1 and 5 bursts before diagnosis
                n_bursts = sample(1:5, 1), 
                # How sharp the steps are (e.g., 10 is very step-like, 2 is more wavy)
                burst_sharpness = 7 
            ),
            accelerating = list(shape = 2.0) # > 1 means accelerating over time; this is very extreme acceleration
        )

        for (model_name in models_to_run) {
            sim_df <- simulate_events_monotonic(
                clonal_events = clonal_events,
                subclonal_events = subclonal_events,
                time_subclonal = time_subclonal,
                years = seq(-5, 15, by = 0.1),
                future_horizon = future_horizon,
                model = model_name,
                model_params = model_param_list[[model_name]]
            )

            sim_df$gg <- gg
            sim_df$trajectory <- trajectory
            sim_df$sample <- sample_id
            sim_df$sim <- sim
            sim_df$model <- model_name

            sim_list[[length(sim_list) + 1]] <- sim_df
        }
    }

    data.table::rbindlist(sim_list)
} %>% as.data.frame()

write_delim(simulations_df, file.path(outdir, "trajectory_drivers_simulation_results.tsv"), delim = "\t")
