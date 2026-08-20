# Pull summary evolutionary metrics from DPClust and Battenberg files for all samples

library(readr)
library(dplyr)
library(stringr)
library(survival)
source("src/utils.R")
source("src/dpclust_parsing_functions.R")

## PATHS
outdir <- "outputs/00_preprocessing"
fs::dir_create(outdir, recurse = T)
dpclust_dir <- "data/processed/PPCG_SNV_clustering/"
batt_dir <- "data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/Subclonal_SCNA"
clinfile <- "data/meta/PPCG_donors_clin_20241217.csv"
anatomical_loc_file <- "data/meta/Sample_Donor_Tissue_Origin_2024.csv"

## LOAD DATA
clin_meta <- read_delim(clinfile, ",")
anatomical_loc = read_delim(anatomical_loc_file, ",")

## GET EVO METRICS
# total number of subclones (assumption -> subclones is i) not highest CCF DPClust, ii) CCF < 0.85)
ccf_subclone <- 0.9
dpc_fps <- list.files(dpclust_dir, pattern = "bestClusterInfo.txt", recursive = T, full = T)
nsubclones <- sapply(dpc_fps, function(dpc_fp) {
    get_nsubclones(dpc_fp, ccf_subclone)
})

# mutation ITH
mut_ith <- sapply(dpc_fps, function(dpc_fp) {
    get_mut_ith(dpc_fp, ccf_subclone)
})

# total subclonal mutation burden
subclonal_burden <- sapply(dpc_fps, function(dpc_fp) {
    get_subclone_size(dpc_fp, ccf_subclone)
})

# largest subclonal expansion
largest_subclone <- sapply(dpc_fps, function(dpc_fp) {
    get_largest_subclone(dpc_fp, ccf_subclone)
})

sum_ccf_subclone <- sapply(dpc_fps, function(dpc_fp) {
    get_sum_ccf_subclones(dpc_fp, ccf_subclone)
})

# PGA, clonal and subclonal
batt_fps <- list.files(batt_dir, pattern = "SCNA.txt", full = T)
wgd_samps <- get_wgd(qc = FALSE)

# subclonal PGA
subclonal_pga <- sapply(batt_fps, function(batt_fp) {
    print(batt_fp)
    get_subclonal_pga(batt_fp, wgd_samps = wgd_samps)
})

# clonal PGA
clonal_pga <- sapply(batt_fps, function(batt_fp) {
    get_clonal_pga(batt_fp, wgd_samps = wgd_samps)
})

# Number of SCNAs
nscnas <- sapply(batt_fps, function(batt_fp) {
    get_nscna(batt_fp)
})

# total pga
total_pga <- subclonal_pga + clonal_pga

## store everything into dataframe
pga_df <- data.frame(
    clonal_pga = clonal_pga, total_pga = total_pga, subclonal_pga = subclonal_pga,
    nscnas = nscnas,
    smp = str_extract(basename(names(clonal_pga)), "[^_]+")
)

write_delim(pga_df, file.path(outdir, "pga_metrics.tsv"), delim = "\t")

evo_df <- data.frame(
    nsubclones = nsubclones, mut_ith = mut_ith,
    subclonal_burden = subclonal_burden,
    largest_subclone = largest_subclone,
    ccf_sum_subclone = sum_ccf_subclone,
    smp = str_extract(basename(names(nsubclones)), "[^_]+")
)

evo_df = full_join(pga_df, evo_df, by = "smp")

evo_df = dplyr::select(evo_df, smp, nsubclones, mut_ith, subclonal_burden, largest_subclone, ccf_sum_subclone, clonal_pga, subclonal_pga, nscnas, total_pga)

write_delim(evo_df, file.path(outdir, "evo_metrics.tsv"), delim = "\t")

# Filter to primary tumour samples and add clinical information
evo_df <- evo_df[evo_df$smp %in% str_remove(get_prim_smps(), "_DNA"), ]
evo_df$group <- get_gs_group(evo_df$smp)
evo_df$record_id = gsub("[a-z]$", "", evo_df$smp)

evo_df = left_join(evo_df, clin_meta %>% dplyr::select(record_id, psa_at_tumour_collection, path_t_stage, age_at_tumour_collection, path_t_stage, mets_ind, predicted_ancestry, donor_interval_of_last_followup), by = "record_id")
write_delim(evo_df, file.path(outdir, "ST2.tsv"), delim = "\t")