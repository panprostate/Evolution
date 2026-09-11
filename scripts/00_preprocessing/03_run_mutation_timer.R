# Run Mutation TimeR

source("src/utils.R")
source("src/mutation_timer.R")
options(bitmapType='cairo')
outdir <- fs::dir_create("data/processed/mutation_timer/")

# Datasets
snvdir <- "data/raw/Somatic_variants/SNVs/Filtered_SNV_VCFs_20_April_2021/"
indeldir <- "data/raw/Somatic_variants/InDels/Tier_1/"
# ccfdir <- "data/raw/pptech_exchange/Working_Groups/PP_Evo/SNV_clonality_DPClust_results_16_02_2023/"
ccfdir <- "data/processed/PPCG_SNV_clustering/"
purdir <- "data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/purity_ploidy"
cnadir <- "data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/Subclonal_SCNA/"
#wgdfile <- "data/raw/data_releases/WGS_Data_Release/Somatic_variants/WGD/WGD_consensus_estimates_25_May_2021.txt"

# 1) Run MutationTimeR
pipeline_muttime(
    snvdir = snvdir,
    indeldir = indeldir,
    ccfdir = ccfdir,
    purdir = purdir,
    cnadir = cnadir,
    wgdfile = NULL,
    outdir = outdir,
    ncores = 10
)
