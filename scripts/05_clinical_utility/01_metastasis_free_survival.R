#' Metastasis-free survival analysis for the PPCG Evolution study
#' 
#' This script performs a metastasis-free survival analysis across the PPCG cohort using evolutionary metrics
#' 
#' 

rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(survival)
library(survminer)
library(GenomicRanges)
library(foreach)
library(doMC)

CORES = 20
# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")
source("src/trajectory_utils.R")

get_mutations = function(gene_name, driver_genes.gr, bb.grs, driver_calls_df, gene_type = "tsg", biallelic = FALSE){
    # Filter to gene of interest
    stopifnot(gene_name %in% driver_genes.gr$gene)
    stopifnot(gene_type %in% c("tsg", "og"))
    gene.gr = driver_genes.gr[driver_genes.gr$gene == gene_name]
    # Get driver calls in this gene
    gene_mutants = driver_calls_df[driver_calls_df$symbol == gene_name,] %>% dplyr::select(-symbol) %>% {names(.)[. != "."]} # NB: . means no alteration; otherwise is the functional annotation of driver mutation
    # Get CN status of the gene in each sample
    is_loh = unlist(foreach(bb.gr = bb.grs) %dopar% {
        sum(width(pintersect(bb.gr[bb.gr$is_loh], gene.gr))) > 0
    })
    is_hd = unlist(foreach(bb.gr = bb.grs) %dopar% {
        sum(width(pintersect(bb.gr[bb.gr$is_hd], gene.gr))) > 0
    })
    is_gain = unlist(foreach(bb.gr = bb.grs) %dopar% {
        sum(width(pintersect(bb.gr[bb.gr$is_gain], gene.gr))) > 0
    })

    if (gene_type == "og"){
        altered_samples = union(names(bb.grs)[is_gain], gene_mutants)
    } 
    if (gene_type == "tsg"){
        if (biallelic){
            altered_samples = union(names(bb.grs)[is_hd], dplyr::intersect(gene_mutants, names(bb.grs)[is_loh]))
        } else {
            altered_samples = union(names(bb.grs)[is_loh], union(gene_mutants, names(bb.grs)[is_hd]))
        }
    }
    return(altered_samples)
}

preload_cna = function(cna_path, only_clonal = FALSE){
    # Read SCNA file
    bb_data <- read_delim(cna_path, delim = "\t") %>% dplyr::filter(chr %in% c(1:22), !is.na(startpos), !is.na(nMin1_A))
    # Calculate ploidy
    widths <- bb_data$endpos - bb_data$startpos
    bb_data$avg_cn <- ifelse(
        is.na(bb_data$frac2_A),
        bb_data$nMin1_A + bb_data$nMaj1_A,
        bb_data$frac1_A * (bb_data$nMin1_A + bb_data$nMaj1_A) + bb_data$frac2_A * (bb_data$nMin2_A + bb_data$nMaj2_A)
    )
    ploidy <- sum(widths * bb_data$avg_cn) / sum(widths)
    bb_data$rel_cn = bb_data$avg_cn / ploidy
    # Calculate proportion of genome with minor copy number 0
    if (only_clonal){
        bb_data$is_loh = bb_data$frac1_A > 0.9 & bb_data$nMin1_A == 0
    } else {
        bb_data$is_loh = ifelse(
            # only one subclone, so check minor copy number equal to 0 in this case
            is.na(bb_data$frac2_A), bb_data$nMin1_A == 0, 
            # multiple subclones, check if for any of them minor copy number equal to 0
            bb_data$nMin1_A == 0 | bb_data$nMin2_A == 0) # only consider clonal LOH
    }
    prop_loh = sum(widths[bb_data$is_loh]) / sum(widths)
    # Calculate proportion of genome with HD
    if (only_clonal){
        bb_data$is_hd = bb_data$frac1_A > 0.9 & bb_data$nMin1_A == 0 & bb_data$nMaj1_A == 0
    } else {
        bb_data$is_hd = ifelse(
            # only one subclone, so check both alleles equal to 0 in this case
            is.na(bb_data$frac2_A), bb_data$nMin1_A == 0 & bb_data$nMaj1_A == 0,
            # multiple subclones, check if for any of them both alleles equal to 0
            (bb_data$nMin1_A == 0 & bb_data$nMaj1_A == 0) | (bb_data$nMin2_A == 0 & bb_data$nMaj2_A == 0)
        )
    }
    
    prop_hd = sum(widths[bb_data$is_hd]) / sum(widths)
    # Calculate proportion of genome with Gain (total copy number >= 3 if diploid; >= 5 if WGD)
    # WGD calls - should be carried out as per PCAWG (WGD if ploidy < 2.9 - [proportion of the genome with minor copy number 0])
    is_wgd = 2.9 < (ploidy + prop_loh)
    gain_threshold = ifelse(is_wgd, 5, 3)
    if (only_clonal) {
        bb_data$is_gain = 
            (bb_data$nMin1_A + bb_data$nMaj1_A) >= gain_threshold & bb_data$frac1_A > 0.9 | 
            (bb_data$nMin1_A + bb_data$nMaj1_A) >= gain_threshold & (bb_data$nMin2_A + bb_data$nMaj2_A) >= gain_threshold
    } else {
        bb_data$is_gain = ifelse(
            # only one subclone
            is.na(bb_data$frac2_A), (bb_data$nMin1_A + bb_data$nMaj1_A) >= gain_threshold,
            # multiple subclones
            ((bb_data$nMin1_A + bb_data$nMaj1_A) >= gain_threshold) | ((bb_data$nMin2_A + bb_data$nMaj2_A) >= gain_threshold)
        )
    }
    prop_gain = sum(widths[bb_data$is_gain]) / sum(widths)

    # Estimate fraction of the genome altered
    normal_allele_cn = ifelse(is_wgd, 2, 1)
    bb_data$is_altered = !is.na(bb_data$frac2_A) | bb_data$nMin1_A != normal_allele_cn | bb_data$nMaj1_A != normal_allele_cn
    pga = sum(widths[bb_data$is_altered]) / sum(widths)
    bb_data$pga = pga
    bb.gr = GRanges(seqnames = bb_data$chr, IRanges(start = bb_data$startpos, end = bb_data$endpos), is_loh = bb_data$is_loh, is_hd = bb_data$is_hd, is_gain = bb_data$is_gain, is_altered = bb_data$is_altered, sample = extract_ppcg_id(cna_path, full = FALSE))
    return(bb.gr)
}

registerDoMC(cores = CORES)
# PATHS
outdir <- "outputs/05_clinical_utility"
figdir <- "figures/05_clinical_utility"
dir.create(outdir, recursive = TRUE); dir.create(figdir, recursive = TRUE)
clinfile <- "data/meta/PPCG_donors_clin_20241217.csv"
evo_metrics_fp <- "outputs/00_preprocessing/evo_metrics.tsv"
trajectories_fps <- list.files("outputs/02_trajectories/", pattern = "PPCG_Feb2026.*mergedseg_with_clonality.txt", full = TRUE)
mrca_timing_fp <- "outputs/01_landscape/real_timing/latency_estimates_mrca_patient_rate_summary.tsv"
cnadir <- "data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/Subclonal_SCNA"
cna_fs <- list.files(cnadir, full.names = TRUE)
drivers_fp <- "data/raw/Somatic_drivers/DRIVER_ANNOTATIONS__2023-09-25_release4__SNV+indel__filteredVariants__driverGenes.txt"
drivers_annotation <- "data/meta/drivers/prostate_ppcg.tsv"

# LOAD DATA
evo_metrics = read_delim(evo_metrics_fp, delim = "\t") 
evo_metrics$record_id = extract_ppcg_pt(evo_metrics$smp) 
evo_metrics$sample = paste0(evo_metrics$smp, "_DNA")
mrca_timing = read_delim(mrca_timing_fp, delim = "\t") %>% dplyr::filter(acceleration == "1x") %>% dplyr::select(sample = sample_id, mrca_latency = latency)
tj_samples = load_trajectories(trajectories_fps)
tj_metrics = get_trajectory_metrics(trajectories_fps, ppcg_trajectory_fps = trajectories_fps) %>% dplyr::mutate(sample = extract_ppcg_id(Tumour_Name, full = F))
met_relapse = get_relapse_data(clinfile, extract_ppcg_pt(tj_samples$sample), type = "metastasis") 
drivers = read_delim(drivers_fp, delim = "\t") %>% dplyr::select(-n_mut_samples, -gene_source)
driver_genes.gr = read_delim(drivers_annotation) %>% {GRanges(seqnames = .$chr, IRanges(start = .$start, end = .$end), gene = .$gene)}

# precompute GenomicRanges with LOH, HD and Gains to query CN status of drivers
cache_bb = file.path(outdir, ".bb_grs.rds")
if ( file.exists(cache_bb) ){
    bb.grs = readRDS(cache_bb)
} else {
    bb.grs = foreach(cna_path = cna_fs) %dopar% {
        preload_cna(cna_path, only_clonal = TRUE)
    }
    names(bb.grs) = extract_ppcg_id(cna_fs, full = FALSE)
    saveRDS(bb.grs, cache_bb)
}

# Merge all the information with clinical data
tj_data = merge(tj_samples, tj_metrics, by = "sample")
evo_metrics = merge(evo_metrics, tj_data, by = "sample")
evo_metrics = left_join(evo_metrics, mrca_timing, by = "sample")
met_relapse = merge(met_relapse, evo_metrics, by = "record_id") # 632 primary samples

# ANALYSIS
# 1. Clean clinical covariates
met_relapse$psa_cat <- cut(met_relapse$psa, 
    breaks = c(-Inf, 10, 20, Inf), 
    labels = c("0-10", "10-20", "+20")
)
met_relapse$t_stage_simplified <- str_remove(met_relapse$t_stage, "[a-z]$")
met_relapse$ggroup = ifelse(met_relapse$grade_group %in% c(1, 2), "1-2", ifelse(met_relapse$grade_group %in% c(4, 5), "4+", met_relapse$grade_group))
met_relapse$age_cat <- ifelse(met_relapse$age < 55, "<55", ifelse(met_relapse$age < 65, "55-65", "65+"))

# 2. Define thresholds for the different metrics and binarise
thr_ith = median(met_relapse$mut_ith, na.rm = TRUE); met_relapse$high_ith = ifelse(met_relapse$mut_ith > thr_ith, "High", "Low")
thr_pga = median(met_relapse$total_pga, na.rm = TRUE); met_relapse$high_pga = ifelse(met_relapse$total_pga > thr_pga, "High", "Low")
thr_poete = median(met_relapse$total_poete, na.rm = TRUE); met_relapse$advanced = ifelse(met_relapse$total_poete > thr_poete, "Advanced", "Non-advanced")
met_relapse$advanced = ifelse(met_relapse$total_poete < quantile(met_relapse$total_poete, 0.33, na.rm = TRUE), "Low", ifelse(met_relapse$total_poete > quantile(met_relapse$total_poete, 0.66, na.rm = TRUE), "High", "Intermediate"))
met_relapse$advanced = ifelse(met_relapse$total_poete > quantile(met_relapse$total_poete, 0.66, na.rm = TRUE), "High", "Low")
thr_mrca = median(met_relapse$mrca_latency, na.rm = TRUE); met_relapse$short_mrca_latency = ifelse(met_relapse$mrca_latency < thr_mrca, "Short", "Long")

# 3. Get alternative genetic metrics (TBC)
met_relapse$is_wgd = met_relapse$sample %in% get_wgd()
tp53_mutants = get_mutations("TP53", driver_genes.gr, bb.grs, drivers, gene_type = "tsg", biallelic = FALSE)
spop_mutants = get_mutations("SPOP", driver_genes.gr, bb.grs, drivers, gene_type = "tsg", biallelic = FALSE)
rb1_mutants = get_mutations("RB1", driver_genes.gr, bb.grs, drivers, gene_type = "tsg", biallelic = FALSE)
pten_mutants = get_mutations("PTEN", driver_genes.gr, bb.grs, drivers, gene_type = "tsg", biallelic = FALSE)


met_relapse$tp53_altered = ifelse(met_relapse$sample %in% tp53_mutants, "Altered", "Wild-type")
met_relapse$spop_altered = ifelse(met_relapse$sample %in% spop_mutants, "Altered", "Wild-type")
met_relapse$rb1_altered = ifelse(met_relapse$sample %in% rb1_mutants, "Altered", "Wild-type")
met_relapse$pten_altered = ifelse(met_relapse$sample %in% pten_mutants, "Altered", "Wild-type")

met_relapse$tsg_alt = (met_relapse$tp53_altered == "Altered") + (met_relapse$rb1_altered == "Altered") + (met_relapse$pten_altered == "Altered")
# met_relapse$tsg_alt = ifelse(met_relapse$tsg_alt >= 2, ">=2", ifelse(met_relapse$tsg_alt == 1, "1", "0"))

# Set reference levels for the different variables
met_relapse = met_relapse %>% dplyr::mutate(
    high_pga = fct_relevel(high_pga, "Low"),
    advanced = fct_relevel(advanced, "Low")
)

# 4. Run survival Cox model for each metric, adjusting for clinical covariates
model = coxph(Surv(time2relapse, relapse_ind) ~ 
        ggroup + psa_cat + t_stage_simplified + age + # clinical features
        trajectory + mut_ith + nsubclones + advanced + # target evolutionary features
        high_pga # confounding genetic features
        ,data = met_relapse)

summary(model)

coefficients = summary(model)$coefficients

write_tsv(
  met_relapse %>% dplyr::select(time2relapse, relapse_ind, ggroup, psa_cat, t_stage_simplified,
                                 age, trajectory, advanced, mut_ith, nsubclones, high_pga),
  file.path(outdir, "Fig5a_source_data.tsv")
)
pdf(file.path(figdir, "Fig5a_mfs_by_subclonal_dynamics.pdf"), height = 100/25.3, width = 120/25.3)
ggforest(model, data = met_relapse)
dev.off()

# Association with better outcomes for Alternative trajectory extends to SPOP WT cases
model = coxph(Surv(time2relapse, relapse_ind) ~ 
        ggroup + psa_cat + t_stage_simplified + age + # clinical features
        trajectory + mut_ith + nsubclones + advanced + # target evolutionary features
        high_pga # confounding genetic features
        ,data = met_relapse[met_relapse$spop_altered != "Altered",])

summary(model)

write_tsv(
  met_relapse %>% dplyr::filter(spop_altered != "Altered") %>%
    dplyr::select(time2relapse, relapse_ind, ggroup, psa_cat, t_stage_simplified,
                  age, trajectory, mut_ith, nsubclones, advanced, high_pga),
  file.path(outdir, "EXDF9a_source_data.tsv")
)
pdf(file.path(figdir, "EXDF9a_mfs_by_subclonal_dynamics_spop_wt.pdf"), height = 100/25.3, width = 120/25.3)
ggforest(model, data = met_relapse[met_relapse$spop_altered != "Altered",])
dev.off()

# Certify that evolutionary risk as a marker complements alterations in determining risk
model = coxph(Surv(time2relapse, relapse_ind) ~ 
        ggroup + psa_cat + t_stage_simplified  + age + # clinical features
        trajectory + advanced + # target evolutionary features
        high_pga + tp53_altered + rb1_altered + pten_altered + is_wgd # confounding genetic features
        ,data = met_relapse)

summary(model)

write_tsv(
  met_relapse %>% dplyr::select(time2relapse, relapse_ind, ggroup, psa_cat, t_stage_simplified, age, 
                                trajectory, advanced, 
                                high_pga, tp53_altered, rb1_altered, pten_altered, is_wgd),
  file.path(outdir, "EXDF9c_source_data.tsv")
)
pdf(file.path(figdir, "EXDF9c_mfs_by_subclonal_dynamics_clinical_and_genetic.pdf"), height = 100/25.3, width = 120/25.3)
ggforest(model, data = met_relapse)
dev.off()

model = coxph(Surv(time2relapse, relapse_ind) ~
        ggroup + psa_cat + t_stage_simplified  + age + # clinical features
        trajectory + advanced + # target evolutionary features
        high_pga + is_wgd + tsg_alt # confounding genetic features
        ,data = met_relapse)
        
summary(model)

write_tsv(
  met_relapse %>% dplyr::select(time2relapse, relapse_ind, ggroup, psa_cat, t_stage_simplified, age, 
                                trajectory, advanced, 
                                high_pga, is_wgd, tsg_alt),
  file.path(outdir, "EXDF9b_source_data.tsv")
)
pdf(file.path(figdir, "EXDF9b_mfs_by_subclonal_dynamics_tsg_alt.pdf"), height = 100/25.3, width = 120/25.3)
ggforest(model, data = met_relapse)
dev.off()

# 5. Kaplan-Meier curves by Grade Groups with evolutionary risk classifier
# met_relapse$advanced = met_relapse$total_poteo > median(met_relapse$total_poteo, na.rm = TRUE)
met_relapse = met_relapse %>% dplyr::mutate(evorisk = case_when(
    trajectory == "Ordering 1" & advanced == "High" ~ "High", 
    trajectory == "Ordering 1" & advanced != "High" ~ "Low",
    trajectory == "Ordering 2" ~ "Low", 
    trajectory == "Ordering 3" ~ "High"
))

model = coxph(Surv(time2relapse, relapse_ind) ~ 
        ggroup + psa_cat + t_stage_simplified + age + # clinical features
        evorisk + # target evolutionary features
        high_pga + tsg_alt # confounding genetic features
        ,data = met_relapse)

summary(model)

# A. GRADE GROUP 2
write_tsv(
  met_relapse %>% dplyr::filter(grade_group == "2") %>% dplyr::select(
    time2relapse, relapse_ind, evorisk),
  file.path(outdir, "Fig5c_source_data.tsv")
)
pdf(file.path(figdir, "Fig5c_km_evorisk_gg2.pdf"),  height = 50/25.3, width = 60/25.3)
km_fit = survfit(Surv(time2relapse, relapse_ind) ~ evorisk, 
    data = met_relapse[met_relapse$grade_group == "2",]
)
ggsurvplot(km_fit, data = met_relapse[met_relapse$grade_group == "2",],
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    xlab = "Time to metastasis (days)", 
    ylab = "Metastasis-free survival",
    palette = c("#DC3220", "#005AB5"),

    # Legend management: Move it INSIDE the plot to save precious canvas space
    legend = c(0.2, 0.5),
    legend.title = "Evolutionary risk",
    legend.labs = c("High", "Low"),

    # Make the legend background perfectly transparent so it doesn't block the lines
    ggtheme = theme_classic() + theme(
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.margin = margin(0, 0, 0, 0)
    ),

    # --- SHRINK ALL ELEMENTS FOR 60x80mm ---
    size = 0.4,                  # Thinner survival lines (default is too thick for small plots)
    censor.size = 2,             # Smaller censor ticks
    pval.size = 2.5,             # P-value text size (in mm)
    fontsize = 2.5,              # Risk table inner text size (in mm)
    
    # Scale down all thematic fonts (size in points: Nature uses ~6-8pt)
    font.main = 7,
    font.x = 7,
    font.y = 7,
    font.tickslab = 6,
    font.legend = 6,
    
    # --- FIX 3: RISK TABLE MANAGEMENT ---
    risk.table.height = 0.2,    # Use 30% of the vertical space for the table
    risk.table.y.text = FALSE,   # CRITICAL: Replaces y-axis text with color bars to save horizontal space
    risk.table.title = "At risk",
    
    # Clean up the table theme
    tables.theme = theme_cleantable() + 
        theme(
            plot.title = element_text(size = 7), # Table title size
            axis.text.x = element_blank(),       # Remove duplicate x-axis labels from table
            axis.title.y = element_blank(), 
            # This explicitly anchors the left and right margins to align with the plot above it
            plot.margin = margin(0, 5.5, 0, 5.5)
        )
)
dev.off()

model = coxph(Surv(time2relapse, relapse_ind) ~ evorisk +  high_pga + psa_cat + t_stage_simplified, data = met_relapse[met_relapse$grade_group == "2",])
summary(model)

model = coxph(Surv(time2relapse, relapse_ind) ~ evorisk, data = met_relapse[met_relapse$grade_group == "2",])
summary(model)

pdf(file.path(figdir, "mfs_by_subclonal_dynamics_gg2_clinical.pdf"), height = 50/25.3, width = 100/25.3)
ggforest(model, data = met_relapse[met_relapse$grade_group == "2",])
dev.off()

model = coxph(Surv(time2relapse, relapse_ind) ~ trajectory + advanced +  high_pga + psa_cat + t_stage_simplified, data = met_relapse[met_relapse$grade_group == "2",])
summary(model)

# B. GRADE GROUP 3
write_tsv(
  met_relapse %>% dplyr::filter(grade_group == "3") %>% dplyr::select(time2relapse, relapse_ind, evorisk),
  file.path(outdir, "Fig5d_source_data.tsv")
)
pdf(file.path(figdir, "Fig5d_km_evorisk_gg3.pdf"),  height = 50/25.3, width = 60/25.3)
km_fit = survfit(Surv(time2relapse, relapse_ind) ~ evorisk, 
    data = met_relapse[met_relapse$grade_group == "3",]
)
ggsurvplot(km_fit, data = met_relapse[met_relapse$grade_group == "3",],
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    xlab = "Time to metastasis (days)", 
    ylab = "Metastasis-free survival",
    palette = c("#DC3220", "#005AB5"),

    # Legend management: Move it INSIDE the plot to save precious canvas space
    legend = c(0.65, 0.85),
    legend.title = "Evolutionary risk",
    legend.labs = c("High", "Low"),

    # Make the legend background perfectly transparent so it doesn't block the lines
    ggtheme = theme_classic() + theme(
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.margin = margin(0, 0, 0, 0)
    ),

    # --- SHRINK ALL ELEMENTS FOR 60x80mm ---
    size = 0.4,                  # Thinner survival lines (default is too thick for small plots)
    censor.size = 2,             # Smaller censor ticks
    pval.size = 2.5,             # P-value text size (in mm)
    fontsize = 2.5,              # Risk table inner text size (in mm)
    
    # Scale down all thematic fonts (size in points: Nature uses ~6-8pt)
    font.main = 7,
    font.x = 7,
    font.y = 7,
    font.tickslab = 6,
    font.legend = 6,
    
    # --- FIX 3: RISK TABLE MANAGEMENT ---
    risk.table.height = 0.2,    # Use 30% of the vertical space for the table
    risk.table.y.text = FALSE,   # CRITICAL: Replaces y-axis text with color bars to save horizontal space
    risk.table.title = "At risk",
    
    # Clean up the table theme
    tables.theme = theme_cleantable() + 
        theme(
            plot.title = element_text(size = 7), # Table title size
            axis.text.x = element_blank(),       # Remove duplicate x-axis labels from table
            axis.title.y = element_blank(), 
            # This explicitly anchors the left and right margins to align with the plot above it
            plot.margin = margin(0, 5.5, 0, 5.5)
        )
)
dev.off()

model = coxph(Surv(time2relapse, relapse_ind) ~ evorisk + high_pga + psa_cat + t_stage_simplified, data = met_relapse[met_relapse$grade_group == "3",])
summary(model)

model = coxph(Surv(time2relapse, relapse_ind) ~ evorisk, data = met_relapse[met_relapse$grade_group == "3",])
summary(model)

pdf(file.path(figdir, "mfs_by_subclonal_dynamics_gg3_clinical.pdf"), height = 50/25.3, width = 100/25.3)
ggforest(model, data = met_relapse[met_relapse$grade_group == "3",])
dev.off()

model = coxph(Surv(time2relapse, relapse_ind) ~  trajectory + advanced + high_pga + psa_cat + t_stage_simplified, data = met_relapse[met_relapse$grade_group == "3",])
summary(model)

# GRADE GROUP 4+
write_tsv(
  met_relapse %>% dplyr::filter(grade_group %in% c("4", "5")) %>% dplyr::select(time2relapse, relapse_ind, evorisk),
  file.path(outdir, "Fig5e_source_data.tsv")
)
pdf(file.path(figdir, "Fig5e_km_evorisk_gg4.pdf"),  height = 50/25.3, width = 60/25.3)
km_fit = survfit(Surv(time2relapse, relapse_ind) ~ evorisk, 
    data = met_relapse[met_relapse$grade_group %in% c("4", "5"),]
)
ggsurvplot(km_fit, data = met_relapse[met_relapse$grade_group %in% c("4", "5"),],
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    xlab = "Time to metastasis (days)", 
    ylab = "Metastasis-free survival",
    palette = c("#DC3220", "#005AB5"),

    # Legend management: Move it INSIDE the plot to save precious canvas space
    legend = c(0.65, 0.85),
    legend.title = "Evolutionary risk",
    legend.labs = c("High", "Low"),

    # Make the legend background perfectly transparent so it doesn't block the lines
    ggtheme = theme_classic() + theme(
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.margin = margin(0, 0, 0, 0)
    ),

    # --- SHRINK ALL ELEMENTS FOR 60x80mm ---
    size = 0.4,                  # Thinner survival lines (default is too thick for small plots)
    censor.size = 2,             # Smaller censor ticks
    pval.size = 2.5,             # P-value text size (in mm)
    fontsize = 2.5,              # Risk table inner text size (in mm)
    
    # Scale down all thematic fonts (size in points: Nature uses ~6-8pt)
    font.main = 7,
    font.x = 7,
    font.y = 7,
    font.tickslab = 6,
    font.legend = 6,
    
    # --- FIX 3: RISK TABLE MANAGEMENT ---
    risk.table.height = 0.2,    # Use 30% of the vertical space for the table
    risk.table.y.text = FALSE,   # CRITICAL: Replaces y-axis text with color bars to save horizontal space
    risk.table.title = "At risk",
    
    # Clean up the table theme
    tables.theme = theme_cleantable() + 
        theme(
            plot.title = element_text(size = 7), # Table title size
            axis.text.x = element_blank(),       # Remove duplicate x-axis labels from table
            axis.title.y = element_blank(), 
            # This explicitly anchors the left and right margins to align with the plot above it
            plot.margin = margin(0, 5.5, 0, 5.5)
        )
)
dev.off()

model = coxph(Surv(time2relapse, relapse_ind) ~ evorisk + high_pga + psa_cat + t_stage_simplified, data = met_relapse[met_relapse$grade_group %in% c("4", "5"),])
summary(model)

model = coxph(Surv(time2relapse, relapse_ind) ~ evorisk, data = met_relapse[met_relapse$grade_group %in% c("4", "5"),])
summary(model)

pdf(file.path(figdir, "mfs_by_subclonal_dynamics_gg4_5_clinical.pdf"), height = 50/25.3, width = 100/25.3)
ggforest(model, data = met_relapse[met_relapse$grade_group %in% c("4", "5"),])
dev.off()

model = coxph(Surv(time2relapse, relapse_ind) ~ trajectory + advanced + high_pga + psa_cat + t_stage_simplified, data = met_relapse[met_relapse$grade_group %in% c("4", "5"),])
summary(model)

# 6. Analysis in other prognostic schemes
# A. D'Amico risk groups
met_relapse$damico_evorisk = ifelse(met_relapse$damico_risk != "low", paste0(met_relapse$damico_risk, "_Evo", met_relapse$evorisk), met_relapse$damico_risk)

km_fit = survfit(Surv(time2relapse, relapse_ind) ~ damico_evorisk, 
    data = met_relapse
)

write_tsv(
  met_relapse %>% dplyr::select(time2relapse, relapse_ind, damico_evorisk),
  file.path(outdir, "EXDF9d_source_data.tsv")
)
pdf(file.path(figdir, "EXDF9d_km_evorisk_damico.pdf"),  height = 80/25.3, width = 80/25.3)
ggsurvplot(km_fit, data = met_relapse,
    palette = c("firebrick4", "dodgerblue4", "firebrick1", "deepskyblue3", "grey60", "grey60"),
    legend.title = "",
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    xlab = "Days after diagnosis", 
    ylab = "Metastasis-free survival", 
     # --- SHRINK ALL ELEMENTS FOR 60x80mm ---
    size = 0.4,                  # Thinner survival lines (default is too thick for small plots)
    censor.size = 2,             # Smaller censor ticks
    pval.size = 2.5,             # P-value text size (in mm)
    fontsize = 2.5,              # Risk table inner text size (in mm)
    
    # Scale down all thematic fonts (size in points: Nature uses ~6-8pt)
    font.main = 7,
    font.x = 7,
    font.y = 7,
    font.tickslab = 6,
    font.legend = 6,
    
    # --- FIX 3: RISK TABLE MANAGEMENT ---
    risk.table.height = 0.3,    # Use 30% of the vertical space for the table
    risk.table.y.text = FALSE,   # CRITICAL: Replaces y-axis text with color bars to save horizontal space
    risk.table.title = "At risk",
    
    # Clean up the table theme
    tables.theme = theme_cleantable() + 
        theme(
            plot.title = element_text(size = 7), # Table title size
            axis.text.x = element_blank(),       # Remove duplicate x-axis labels from table
            axis.title.y = element_blank(), 
            # This explicitly anchors the left and right margins to align with the plot above it
            plot.margin = margin(0, 5.5, 0, 5.5)
        )
)
dev.off()

model = coxph(Surv(time2relapse, relapse_ind) ~ damico_evorisk, data = met_relapse[met_relapse$damico_risk == "intermediate",])
summary(model)

model = coxph(Surv(time2relapse, relapse_ind) ~ damico_evorisk, data = met_relapse[met_relapse$damico_risk == "high",])
summary(model)

# B. Cambridge Prognostic Groups
met_relapse$cpg = sapply(1:nrow(met_relapse), function(i){
    get_cambridge_prognostic_groups(met_relapse$grade_group[i], met_relapse$psa[i], met_relapse$t_stage[i])
})

met_relapse$cpg_evorisk = ifelse(met_relapse$cpg %in% c("CPG3", "CPG4", "CPG5"), paste0(met_relapse$cpg, "_Evo", met_relapse$evorisk), met_relapse$cpg)

km_fit = survfit(Surv(time2relapse, relapse_ind) ~ cpg_evorisk, 
    data = met_relapse
)

write_tsv(
  met_relapse %>% dplyr::select(time2relapse, relapse_ind, cpg_evorisk),
  file.path(outdir, "EXDF9e_source_data.tsv")
)
pdf(file.path(figdir, "EXDF9e_km_evorisk_cpg.pdf"),  height = 80/25.3, width = 80/25.3)
ggsurvplot(km_fit, data = met_relapse,
    palette = c("grey80", "grey60", "coral", "deepskyblue3", "firebrick1", "dodgerblue3", "firebrick4", "dodgerblue4"),
    legend.title = "",
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    xlab = "Days after diagnosis", 
    ylab = "Metastasis-free survival", 
     # --- SHRINK ALL ELEMENTS FOR 60x80mm ---
    size = 0.4,                  # Thinner survival lines (default is too thick for small plots)
    censor.size = 2,             # Smaller censor ticks
    pval.size = 2.5,             # P-value text size (in mm)
    fontsize = 2.5,              # Risk table inner text size (in mm)
    
    # Scale down all thematic fonts (size in points: Nature uses ~6-8pt)
    font.main = 7,
    font.x = 7,
    font.y = 7,
    font.tickslab = 6,
    font.legend = 6,
    
    # --- FIX 3: RISK TABLE MANAGEMENT ---
    risk.table.height = 0.3,    # Use 30% of the vertical space for the table
    risk.table.y.text = FALSE,   # CRITICAL: Replaces y-axis text with color bars to save horizontal space
    risk.table.title = "At risk",
    
    # Clean up the table theme
    tables.theme = theme_cleantable() + 
        theme(
            plot.title = element_text(size = 7), # Table title size
            axis.text.x = element_blank(),       # Remove duplicate x-axis labels from table
            axis.title.y = element_blank(), 
            # This explicitly anchors the left and right margins to align with the plot above it
            plot.margin = margin(0, 5.5, 0, 5.5)
        )
)
dev.off()

model = coxph(Surv(time2relapse, relapse_ind) ~ evorisk, data = met_relapse[met_relapse$cpg == "CPG3",])
summary(model)

model = coxph(Surv(time2relapse, relapse_ind) ~ evorisk, data = met_relapse[met_relapse$cpg == "CPG4",])
summary(model)

model = coxph(Surv(time2relapse, relapse_ind) ~ evorisk, data = met_relapse[met_relapse$cpg == "CPG5",])
summary(model)
