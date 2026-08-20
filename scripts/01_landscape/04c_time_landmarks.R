# Time MRCA latency (latency to last subclone) and WGD latency in PPCG

rm(list = ls(all = TRUE))

## PACKAGES
library(readr)
library(dplyr)
library(ggplot2)
library(lemon)
library(patchwork)
library(foreach)
library(doMC)
library(tidyverse)
library(MASS)
library(data.table)
library(VariantAnnotation)

## PATHS
outdir <- "outputs/01_landscape/real_timing"
figdir <- "figures/01_landscape/real_timing"
mutation_timer_dir <- "data/processed/mutation_timer/snvs"
dir.create(outdir, recursive = T)
dir.create(figdir, recursive = T)
# NB: first entry is subclonal, second is clonal
sbs_assign_fp <- "data/processed/SNV_signature_annotations__2024-06-21/SNVs_all_signature_20240621.tsv"
clonal_branches_fp <- file.path(outdir, "branch_deamination_dense.rds")
# NB: first entry is pre-WGD, second is post-WGD
wgd_branches_fp <- file.path(outdir, "wgd_deamination.rds")
vcf_fs <- list.files(mutation_timer_dir, pattern = ".*vcf", full = T)
bb_fs <- list.files(mutation_timer_dir, pattern = ".*bb.rds", full = T)
purdir <- "data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/purity_ploidy"

CORES = 20
registerDoMC(cores = CORES)

## FUNCTIONS
source("src/utils.R")
source("src/real_timing.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

plot_median_timing <- function(timing_df, color){
    t = timing_df %>% dplyr::select(sample_id, latency, acceleration) %>% pivot_wider(names_from = acceleration, values_from = latency)
    cols = colnames(t)
    smr = t %>% as.data.frame() %>% dplyr::summarise_at(cols, median, na.rm = T)
    p <- ggplot(smr) + 
        geom_rect(aes(xmin = 1 - 0.05, xmax = 1 + 0.05, ymax = `1x`, ymin = `2.5x`), fill = color, size = .15, col = "black") + 
        geom_rect(aes(xmin = 1 - 0.05, xmax = 1 + 0.05, ymax = `2.5x`, ymin = `5x`), fill = color, size = .15, col = "black", alpha = .75) + 
        geom_rect(aes(xmin = 1 - 0.05, xmax = 1 + 0.05, ymax = `5x`, ymin = `10x`), fill = color, size = .15, col = "black", alpha = .5) + 
        geom_rect(aes(xmin = 1 - 0.05, xmax = 1 + 0.05, ymax = `10x`, ymin = `20x`), fill = color, size = .15, col = "black", alpha = .25) + 
        ylim(c(0,20)) + 
        labs(x = "", y = "Latency (Years)") + 
        theme(
            axis.text.x = element_blank(),  
            axis.ticks.x = element_blank(), 
            axis.line.x = element_blank()   
            )

    # p <- axes2lemon(p, lt = "both", bt = "both")
    return(p)
}

plot_cohort_diagram <- function(timing_df){
    df <- timing_df %>%
        dplyr::filter(acceleration == "1x") %>%
        dplyr::select(sample_id, mrca_latency = latency,
                      record_id, age_at_tumour_collection,
                      relapse_ind, donor_relapse_interval, gs_grade)
    df$age_at_tumour_collection <- as.numeric(df$age_at_tumour_collection)
    df$age_mrca <- df$age_at_tumour_collection - df$mrca_latency
    df$age_mrca <- ifelse(df$age_mrca > df$age_at_tumour_collection, df$age_at_tumour_collection, df$age_mrca)
    df$age_mrca <- ifelse(df$age_mrca < 0, 0, df$age_mrca)

    df$age_relapse <- ifelse(df$relapse_ind == "relapsed",
                             as.numeric(df$age_at_tumour_collection) + as.numeric(df$donor_relapse_interval) / 365.25,
                             NA)

    # Create an artificial variable to split the plot into four groups
    df$facet_group <- factor(cut(1:nrow(df), breaks = 4, labels = paste("Group", 1:4)))

    df <- df[!is.na(df$gs_grade),]

    p <- ggplot(df, aes(y = record_id)) +
        geom_point(aes(x = age_mrca), pch = 21, alpha = .8, fill = icgc["blue"], size = 1) +
        geom_segment(aes(x = 0, xend = age_mrca, yend = record_id), col = icgc["blue"], alpha = 0.4) +
        geom_segment(aes(x = age_mrca, xend = age_at_tumour_collection, yend = record_id), col = icgc["red"], alpha = 0.4) +
        geom_point(aes(x = age_at_tumour_collection), pch = 21, alpha = .8, fill = icgc["red"], size = 1) +
        geom_point(aes(x = 0, fill = gs_grade), pch = 22, size = 1) +
        geom_point(aes(x = age_relapse), pch = 4, col = "red", alpha = .6, size = .6) +
        scale_fill_manual(values = gleason_colours) +
        facet_wrap(~ facet_group, ncol = 4, scales = "free_y") +
        theme(
            strip.text = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.title = element_blank()) +
        labs(x = "Age (Years)", y = "PPCG Patient", fill = "Grade Group")

    return(p)
}

summarise_latency <- function(list_latencies, accelerations){
    df = data.table::rbindlist(lapply(list_latencies, function(x) as.data.frame(x) %>% rownames_to_column("sample_id")), idcol = "iter")

    return(data.table::rbindlist(lapply(accelerations, function(a){
        df %>% dplyr::mutate(latency = .data[[paste0(a, "x")]]) %>% group_by(sample_id) %>% 
            dplyr::summarise(
                latency = median(latency, na.rm = TRUE), 
                low_ci = quantile(latency, 0.025, na.rm = TRUE),
                up_ci = quantile(latency, 0.975, na.rm = TRUE)
            ) %>% dplyr::mutate(acceleration = paste0(a, "x")) 
    })))
}

## LOAD DATA
prim_smps <- get_prim_smps(qc = T)
wgd_smps <- dplyr::intersect(get_wgd(qc = T), prim_smps)
clonal_branches <- t(simplify2array(readRDS(clonal_branches_fp)))
wgd_branches <- t(simplify2array(readRDS(wgd_branches_fp)))
clin_data <- read_delim("data/meta/PPCG_donors_clin_20241217.csv")
accelerations <- c(1, 2.5, 5, 7.5, 10, 20)
boot = 1000
sbs_a = read_delim(sbs_assign_fp)

# Estimate MRCA latency
# Get average clock-rate per year
clock_mutations = rowSums(clonal_branches)
ages = get_age(names(clock_mutations))
age_lm = rlm(clock_mutations ~ 0 + ages) # 28.3 mutations per year
clock_rate = coefficients(age_lm)[1]
t_mrca = mol2time(rownames(clonal_branches), clonal_branches[,2], clonal_branches[,1], clock_rate, accelerations = accelerations, event = "MRCA", boot = boot)
saveRDS(t_mrca, file.path(outdir, "bootstrap_absolute_latency_mrca_global_rate.rds"))

global_rate_mrca_estimates = summarise_latency(t_mrca, accelerations)
write_delim(global_rate_mrca_estimates, file.path(outdir, "latency_estimates_mrca_global_rate_summary2.tsv"), delim = "\t")

# Per-patient rate
t_mrca_pt_rate = mol2time_patient_rate(rownames(clonal_branches), clonal_branches[,2], clonal_branches[,1], accelerations = accelerations, event = "MRCA", boot = boot)
saveRDS(t_mrca_pt_rate, file.path(outdir, "bootstrap_absolute_latency_mrca_patient_rate.rds"))

pt_rate_mrca_estimates = summarise_latency(t_mrca_pt_rate, accelerations)
write_delim(pt_rate_mrca_estimates, file.path(outdir, "latency_estimates_mrca_patient_rate_summary.tsv"), delim = "\t")

# Add clinical variables to timing estimates so source data files are self-contained
clin_for_merge <- clin_data %>%
  dplyr::select(record_id, age_at_tumour_collection, relapse_ind, donor_relapse_interval) %>%
  dplyr::mutate(age_at_tumour_collection = as.numeric(age_at_tumour_collection))

pt_rate_mrca_estimates <- pt_rate_mrca_estimates %>%
  dplyr::mutate(record_id = smp2pt(sample_id)) %>%
  dplyr::left_join(clin_for_merge, by = "record_id") %>%
  dplyr::mutate(gs_grade = get_gs_group(record_id))

global_rate_mrca_estimates <- global_rate_mrca_estimates %>%
  dplyr::mutate(record_id = smp2pt(sample_id)) %>%
  dplyr::left_join(clin_for_merge, by = "record_id") %>%
  dplyr::mutate(gs_grade = get_gs_group(record_id))

# median MRCA latency and 95%CI interval for median
medians = pt_rate_mrca_estimates %>% 
    dplyr::group_by(acceleration) %>% 
    dplyr::summarise(median_latency = median(latency, na.rm = T))

# 10.4 years with 1x rate acceleration; 3.34 with 5x rate

# categorical classification of latency to MRCA
# early if latency > 10 years
# late if latency < 10 years
thr = medians$median_latency[medians$acceleration == "1x"]

pt_rate_mrca_estimates %>% dplyr::filter(acceleration == "1x") %>% 
    dplyr::mutate(
        record_id = smp2pt(sample_id),
        mrca_subgroups = ifelse(latency >= thr, "late", "early")
    ) %>% dplyr::select(record_id, mrca_subgroups) %>% 
    write_delim(file.path(outdir, "is_late_mrca.tsv"), delim = "\t")

# Plot rate estimates for patient-specific rate
p <- plot_cohort_diagram(pt_rate_mrca_estimates)

write_tsv(pt_rate_mrca_estimates, file.path(outdir, "SF3_source_data.tsv"))
ggsave(file.path(figdir, "SF3_mrca_latency_cohort_diagram.pdf"), p, width = 180/25.3, height = 180/25.3)

# Compare timing estimates between patient-specific and global clock rates
p1 = merge(
    pt_rate_mrca_estimates %>% dplyr::filter(acceleration == "1x") %>% dplyr::rename(pt_latency = latency) %>% dplyr::select(sample_id, pt_latency), 
    global_rate_mrca_estimates %>% dplyr::filter(acceleration == "1x") %>% dplyr::rename(global_latency = latency) %>% dplyr::select(sample_id, global_latency), 
by = c("sample_id")) %>% 
    ggplot(aes(x = pt_latency, y = global_latency)) + 
    geom_point(alpha = 0.5) + 
    labs(y = "MRCA latency global clock rate", x = "MRCA latency patient specific rate") + 
    geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") + 
    ggpubr::stat_cor(method = "spearman") + 
    xlim(c(0,70)) + ylim(c(0,70))

p2 <- plot_median_timing(global_rate_mrca_estimates, "#E61C2A") + labs(y = "Median MRCA latency (Years)")

p3 <- plot_cohort_diagram(global_rate_mrca_estimates)

design <- "A##BBBB#
           CCCCCCCC
           CCCCCCCC
           CCCCCCCC"

p <- p2 + p1 + p3 + 
    plot_annotation(tag_levels = 'a') + 
    plot_layout(design = design)

write_tsv(global_rate_mrca_estimates, file.path(outdir, "SF4_source_data.tsv"))
save_ggplot(p, file.path(figdir, "SF4_global_vs_patient_rate"), w = 180, h = 250)

# Estimate WGD latency to MRCA
t_wgd_pt_rate = mol2time_patient_rate(rownames(wgd_branches), wgd_branches[,1], wgd_branches[,2], clonal_branches[rownames(wgd_branches),1], accelerations, event = "WGD", boot = boot)
saveRDS(t_wgd_pt_rate, file.path(outdir, "bootstrap_absolute_latency_wgd_patient_rate.rds"))

pt_rate_wgd_estimates = summarise_latency(t_wgd_pt_rate, accelerations)
pt_rate_wgd_estimates$age = get_age(pt_rate_wgd_estimates$sample_id)
pt_rate_wgd_estimates = merge(pt_rate_wgd_estimates, pt_rate_mrca_estimates %>% dplyr::rename(mrca_latency = latency) %>% dplyr::select(sample_id, mrca_latency, acceleration), by = c("sample_id", "acceleration")) %>% data.table::setkey(NULL)

pt_rate_wgd_estimates %>% dplyr::filter(acceleration == "1x") %>% {median(1 - .$latency / .$age, na.rm = T)} # 0.67
pt_rate_wgd_estimates %>% dplyr::filter(acceleration == "1x") %>% {median(.$latency, na.rm = TRUE)} # 21.6 years

pt_rate_wgd_estimates %>% dplyr::filter(acceleration == "1x") %>% {sort((.$latency + .$mrca_latency) / .$age)} 

write_delim(pt_rate_wgd_estimates, file.path(outdir, "latency_estimates_wgd_patient_rate_summary.tsv"), delim = "\t")

# Time CNA events

# 1) Get the number of clock mutations before and after each timable segment
# and calculate molecular time with clock mutations using MutationTimeR formula
time_segments <- rbindlist(foreach(smp=unique(pt_rate_mrca_estimates$sample)) %dopar% {
    print(smp)
    pull_cn_early_late(smp, bb_fs, vcf_fs, sbs_a)
}) %>% as.data.frame()

saveRDS(time_segments, file.path(outdir, "time_segments.rds"))

# 2) Add the MRCA, global and per-patient rate to table
time_segments <- readRDS(file.path(outdir, "time_segments.rds"))

mrca_time_pt_rate <- pt_rate_mrca_estimates %>% dplyr::filter(acceleration == "1x") %>% dplyr::select(sample_id, latency) %>% 
    dplyr::rename(lt_mrca_pt_rate = latency, sample = sample_id) %>% 
    dplyr::mutate(age = get_age(sample))
mrca_time_pt_rate$mrca_rt_pt_rate = mrca_time_pt_rate$age - mrca_time_pt_rate$lt_mrca_pt_rate

time_segments <- left_join(time_segments, mrca_time_pt_rate)

mrca_time_global_rate <- global_rate_mrca_estimates %>% dplyr::filter(acceleration == "1x") %>% dplyr::select(sample_id, latency) %>% 
    dplyr::rename(lt_mrca_global_rate = latency, sample = sample_id) %>% 
    dplyr::mutate(age = get_age(sample))
mrca_time_global_rate$rt_mrca_global_rate = mrca_time_global_rate$age - mrca_time_global_rate$lt_mrca_global_rate

time_segments <- left_join(time_segments, mrca_time_global_rate)

saveRDS(time_segments, file.path(outdir, "time_segments_with_mrca.rds"))

# 3) Merge segments that are adjacent to each other with same CN type
# and get weighted average
time_segments <- readRDS(file.path(outdir, "time_segments_with_mrca.rds"))

time_segments <- time_segments[!is.na(time_segments$type),]
time_segments$type <- as.character(time_segments$type)

# takes just a bit to run; ~2 mins
time_cnas <- merge_segments(time_segments, 1e6)

# transform to real-time by dividing by time at which MRCA happened
# values higher than 1 should be 1
time_cnas$mtime_cn <- ifelse(time_cnas$mtime_cn > 1, 1, time_cnas$mtime_cn)
time_cnas$rt_cn_pt_rate <- time_cnas$mtime_cn * time_cnas$rt_mrca_pt_rate

saveRDS(time_cnas, file.path(outdir, "time_cnas.rds"))
