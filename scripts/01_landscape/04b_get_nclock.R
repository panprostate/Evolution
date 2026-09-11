# Get number of clock mutations (SBS1+SBS5+SBS40) on each period
# i.e. pre and post MRCA and pre and post WGD
# Normalised by ploidy
# Clonal mutations: normalised by time-averaged ploidy back to diploid (i.e. 2/effGenome)
# Subclonal mutations: normalised by observed final ploidy (i.e. 2/get_ploidy(purdir, id))
# Pre-WGD mutations: assumed to predominantly occur in diploid genome
# Post-WGD mutations: assumed to predominantly occur in 

# LIBRARIES
library(foreach)
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(doParallel)
library(ggplot2)
library(survival)
library(VariantAnnotation)
library(MutationTimeR)
source("src/utils.R")
source("src/real_timing.R")
source("src/mutation_timer.R")

CORES=10
doParallel::registerDoParallel(cores=CORES)

## PATHS
outdir <- "outputs/01_landscape/real_timing"
fs::dir_create(outdir, recurse = T)

mutation_timer_dir <- "data/processed/mutation_timer/snvs"
ccfdir <- "data/processed/PPCG_SNV_clustering"
purdir <- "data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/purity_ploidy"
sqmet_fp <- "data/qc/Sanger_QC_Metrics_01_June_2020.txt"
sbs_assign_fp <- "data/processed/SNV_signature_annotations__2024-06-21/SNVs_all_signature_20240621.tsv"

vcf_fs <- list.files(mutation_timer_dir, pattern = ".*vcf", full = T)
bb_fs <- list.files(mutation_timer_dir, pattern = ".*bb.rds", full = T)
dpc_fps <- list.files(ccfdir, pattern = "bestClusterInfo.txt", recursive = T, full = T)

# Find total number of clock mutations on each evo period
# (i.e. pre and post WGD and pre and post MRCA)

## LOAD DATA
cvg = read_delim(sqmet_fp)
sbs_a = read_delim(sbs_assign_fp)

## CALCULATE TOTAL DEAMINATION BURDEN IN ALL PRIMARY SAMPLES
## ONLY ONE SAMPLE PER PATIENT
prim_smps = get_prim_smps()
is_mmrd = get_mmrd(sbs_a)
prim_smps = prim_smps[!extract_ppcg_pt(prim_smps) %in% is_mmrd]

branch_deam_dense = foreach(smp = prim_smps) %dopar% {
    print(smp)
    vcf_mtr_fp = vcf_fs[str_detect(vcf_fs, smp)]
    if (length(vcf_mtr_fp) != 1){return(c(NA, NA))} else {return(get_branches_dense(vcf_mtr_fp, sbs_a, purdir, cvg, TRUE))}
}

names(branch_deam_dense) <- prim_smps
saveRDS(branch_deam_dense, file.path(outdir, "branch_deamination_dense.rds"))

branch_deam_linear = foreach(smp = prim_smps) %dopar% {
    print(smp)
    vcf_mtr_fp = vcf_fs[str_detect(vcf_fs, smp)]
    if (length(vcf_mtr_fp) != 1){return(c(NA, NA))} else {return(get_branches_linear(vcf_mtr_fp, sbs_a, purdir, cvg, TRUE))}
}

names(branch_deam_linear) <- prim_smps
saveRDS(branch_deam_linear, file.path(outdir, "branch_deamination_linear.rds"))

wgd_smps = get_wgd()
wgd_smps = wgd_smps[wgd_smps %in% prim_smps]
# NB: Estimates will provide number of mutations in clonal and subclonal branches
clock_wgd_muts = foreach(id = wgd_smps) %dopar% {
    print(id)
    return(time_clock_wgd(vcf_fs, bb_fs, sbs_a, purdir, ccfdir, id))
}

names(clock_wgd_muts) <- wgd_smps

saveRDS(clock_wgd_muts, file.path(outdir, "wgd_deamination.rds"))
