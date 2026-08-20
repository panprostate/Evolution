# Winner's curse correction is applied to scale number of mutations in a subclone and its CCF across the different analyses performed in the study
# (e.g. real-timing, trajectories, etc.)
# Here, we compare the number of subclonal mutations with Winner's Curse Correction (WCC) versus without performing any correction
# to explicitly determine its impact

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)   
library(VariantAnnotation)    
library(data.table)
library(foreach)     
library(doMC)        
library(ggpubr)  
library(cowplot)
library(lemon)
library(patchwork)


# FUNCTIONS ---------------------------------------------------------------
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")
source("src/spoilsport.R")
source("src/dpclust_parsing_functions.R")

# PATHS ---------------------------------------------------------------

ccfdir <- "data/processed/PPCG_SNV_clustering/"
purdir <- "data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/purity_ploidy/"
sqmet_fp <- "data/qc/Sanger_QC_Metrics_01_June_2020.txt"
ccf_fps <- list.files(ccfdir, pattern = "bestClusterInfo.txt", recursive = T, full = T)

## LOAD DATA
cvg = read_delim(sqmet_fp)
outdir <- "outputs/00_preprocessing/wcc/"; dir.create(outdir, recursive = TRUE)
figdir <- "figures/00_preprocessing/wcc"; dir.create(figdir, recursive = TRUE)

smps <- c(get_prim_smps(), get_met_smps())

# CCF to consider a cluster as clonal
clonal_threshold = 0.9

# Register parallele backend with 10 cores
registerDoMC(cores = 10)
muts_wcc_df = foreach(sid = smps) %dopar% {
    print(sid)
    # Load all the information for the sample
    ccf_fp = ccf_fps[str_detect(ccf_fps, sid)]
    if (length(ccf_fp) != 1){
      return(data.frame(
          sample = sid,
          nclonal = NA,
          nsubclonal = NA,
          nclonal_wcc = NA,
          nsubclonal_wcc = NA
      ))
    }

    ccf = read_delim(ccf_fp)
    pur = get_purity(purdir, sid)
    plo = get_ploidy(purdir, sid)
    cov = round(cvg[["Tumour.total.depth"]][cvg$PPCG_Sample_ID == sid])
    # SpoilSport requires cellular prevelances, convert CCFs into CPs
    ccf$CP = ccf$location * pur
    ccf$is_clonal = ccf$location >= clonal_threshold
    wcc_out = wcc(pur, plo, cov, 3, ccf$CP)
    
    # convert SpoilSport output columns back to corrected CCF and number of mutations
    pcCCF = wcc_out$tp / pur
    pc_nmuts = ccf$no.of.mutations*wcc_out$sf
    nclonal_wcc = round(sum(pc_nmuts[ccf$is_clonal]))
    nsubclonal_wcc = round(sum(pc_nmuts[!ccf$is_clonal]))
    
    return(data.frame(
        sample = sid,
        nclonal = sum(ccf$no.of.mutations[ccf$is_clonal]),
        nsubclonal = sum(ccf$no.of.mutations[!ccf$is_clonal]),
        nclonal_wcc = nclonal_wcc,
        nsubclonal_wcc = nsubclonal_wcc, 
        purity = pur, 
        ploidy = plo
    ))
} %>% data.table::rbindlist()

write_delim(muts_wcc_df, file.path(outdir, "muts_winner_curse_correction.tsv"), delim = "\t")

# Plot changes in number mutations upon WCC correction
muts_wcc_df$NRPCC = get_nrpcc(muts_wcc_df$sample)
# Filter to samples passing QC
muts_wcc_df = muts_wcc_df %>% dplyr::filter(sample %in% qc_pass())

# 1 - CLONAL MUTATIONS
# plot number of clonal and subclonal mutations with and without WCC
p1 = ggplot(muts_wcc_df, aes(x = log10(nclonal), y = log10(nclonal_wcc))) +
    geom_point(alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Number clonal mutations\n(log10)",
         y = "Number of clonal mutations\n(WCC corrected, log10)") +
    theme(legend.position = "bottom")

p1 = axes2lemon(p1, lt = "both", bt = "both")
ggsave(file.path(figdir, "clonal_mutations_wcc_vs_no_wcc.pdf"), width = 50/25.3, height = 50/25.3)

# 2- SUBCLONAL MUTATIONS
# plot number of clonal and subclonal mutations with and without WCC
p2 = ggplot(muts_wcc_df, aes(x = log10(nsubclonal), y = log10(nsubclonal_wcc))) +
    geom_point(alpha = .3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Number of subclonal mutations\n(log10)",       
         y = "Number of subclonal mutations\n(WCC corrected, log10)") +
    theme(legend.position = "bottom")

p2 = axes2lemon(p2, lt = "both", bt = "both")
ggsave(file.path(figdir, "subclonal_mutations_wcc_vs_no_wcc.pdf"), width = 50/25.3, height = 50/25.3)

## COMPARE AGAINST NRPCC
muts_wcc_df = muts_wcc_df %>% filter(!is.na(nclonal_wcc) & !is.na(nsubclonal_wcc) & nsubclonal_wcc > 0 & nclonal_wcc > 0)

max_nrpcc = ceiling(max(muts_wcc_df$NRPCC))
lims = c(0, max_nrpcc)

p3 = ggplot(muts_wcc_df %>% dplyr::filter(NRPCC >= 10), aes(y = nsubclonal_wcc/(nclonal_wcc + nsubclonal_wcc), x = NRPCC)) +
    geom_point(alpha = .3) +
    labs(y = "Ratio subclonal mutations",       
         x = "NRPCC") +
    theme(legend.position = "bottom") + 
    geom_smooth(method = "lm", color = "blue", se = FALSE) + 
    stat_cor(method = "spearman", size = 2.5, label.x = 0, label.y = 1.0) + 
    xlim(lims) + ylim(c(0,1))

# EXDF1c - Relationship between NRPCC and ratio of subclonal mutations
write_tsv(
  muts_wcc_df %>% dplyr::filter(NRPCC >= 10) %>%
    dplyr::select(NRPCC, nclonal_wcc, nsubclonal_wcc),
  file.path(outdir, "EXDF1c_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF1c_ratio_subclonal_mutations_wcc_vs_nrpcc.pdf"), width = 50/25.3, height = 50/25.3)

p = p1 + p2 + p3 + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 10, face = "bold"))

ggsave(file.path(figdir, "number_mutations_wcc_vs_no_wcc.pdf"), p, width = 150/25.3, height = 50/25.3)

# Average correction
mean(muts_wcc_df$nsubclonal_wcc/muts_wcc_df$nsubclonal)
