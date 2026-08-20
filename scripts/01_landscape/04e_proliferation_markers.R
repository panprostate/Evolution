# Analysis of markers of proliferation across Gleason groups

rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(cowplot)
library(qs)
library(DESeq2)
library(edgeR)
library(fgsea)
library(msigdbr)
library(GSVA)

# PATHS
outdir <- "outputs/01_landscape/proliferation_markers"
figdir <- "figures/01_landscape/proliferation_markers"
dir.create(outdir); dir.create(figdir)
clinfile <- "data/meta/PPCG_donors_clin_20241217.csv"
methylation_fp <- "data/processed/PPCG.annotation.1870.samples.BEN.TUM.MET.April2025.csv"
mrca_time_fp <- file.path("outputs/01_landscape/real_timing", "latency_estimates_mrca_patient_rate_summary.tsv")
rna_fp <- "data/raw/RNA/PPCG_RNAseq_RUVIIIPRPS_Supervised.qs"
pathways_fp <- "data/meta/pathways/pathways.tsv"

# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")
source("src/rna_functions.R")

boxplot_epicmit <- function(data, xvar, yvar, xlab, ylab){

    p = data %>% dplyr::filter(!is.na(!!sym(xvar))) %>% ggplot(aes_string(x = xvar, y = yvar)) +
            geom_point(position = position_jitter(width = 0.2), alpha = 0.1) +  
            geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) + 
            ggpubr::stat_compare_means(
                comparisons = list(c("1", "2"), c("2", "3"), c("3", "4+"), c("2", "4+")),
                method = "wilcox.test", 
                size = 2
            ) +
            labs(x = xlab, y = ylab) +
            theme(legend.position = "none")
}

# LOAD DATA
t_mrca = read_delim(mrca_time_fp, delim = "\t") %>% dplyr::filter(acceleration == "1x")
t_mrca$grade_group <- get_gs_group(t_mrca$sample_id, clinfile)
methylation_data = read_delim(methylation_fp, delim = ",")
rna_counts = qread(rna_fp)  # Load RNA counts (SummarizedExperiment or similar)
pathways = read_delim(pathways_fp)  # Load pathway gene sets
# Filter to proliferation-related pathways from HALLMARK collection (G2M and E2F targets)
pathways = pathways %>% dplyr::filter(str_detect(pathway, "HALLMARK") & str_detect(pathway, "G2M|E2F"))  # Filter to relevant pathways

# Epigenetic mitotic clocks: epiCMIT
# Combine information
methylation_data$sample_id = methylation_data$Matching_WGS_Assay
methylation_data = merge(methylation_data, t_mrca)

p1 = boxplot_epicmit(methylation_data, "grade_group", "epiCMIT", "Grade Group", "epiCMIT")
p2 = boxplot_epicmit(methylation_data, "grade_group", "epiCMIT.hyper", "Grade Group", "epiCMIT hypermethylation")
p3 = boxplot_epicmit(methylation_data, "grade_group", "epiCMIT.hypo", "Grade Group", "epiCMIT hypomethylation")

p = plot_grid(p1, p2, p3, nrow = 1)
ggsave(file.path(figdir, "epiCMIT_by_grade_group.png"), p, width = 120/25.3, height = 65/25.3, dpi = 300)
ggsave(file.path(figdir, "epiCMIT_by_grade_group.pdf"), p, width = 120/25.3, height = 65/25.3)

# Differential expression analysis by Grade Group
# Prepare pathways as a named list of gene vectors for GSVA
pathways = lapply(unique(pathways$pathway), function(pathway){
  return(pathways$gene[pathways$pathway == pathway])
}) %>% setNames(unique(pathways$pathway))

# Get raw counts and filter to samples with trajectory data
idx = rna_counts@colData[["Matching_WGS_Sample.ID.ADD.NA."]] %in% get_prim_smps()

countData = assays(rna_counts)[["Raw.count"]][,idx]  # Raw counts matrix for selected samples
countData = expression_filter(countData)            # Filter out lowly expressed genes

# Combine colData with RUVIII.PRPS adjustment factors
colData = cbind(
  rna_counts@colData[idx,], 
  rna_counts@metadata[["RUVIII.PRPS.AdjustedForPurity.Supervised.W"]][idx,]
)

# Define technical covariates (RUV factors)
technical_covariates = colnames(rna_counts@metadata[["RUVIII.PRPS.AdjustedForPurity.Supervised.W"]])

# Add Grade Group information to ColData
colData$grade_group = get_gs_group(colData$Matching_WGS_Sample.ID.ADD.NA., clinfile, asnum = T, collapse_gs4 = F)

# filter NA values in grade group
na_ggroup = is.na(colData$grade_group)
colData = colData[!na_ggroup, ]
countData = countData[, !na_ggroup]

# Run DESeq2 analysis comparing Grade Groups as numerical variable controlling for technical covariates
deseq_grade_group = run_deseq2("grade_group", countData, colData, technical_covariates, control_gleason_group = F, return_comparison = FALSE)
results = results(deseq_grade_group, name = "grade_group")
gsea_grade = run_gsea(results, pathways)

# E2F targets: p-value = 0.03, NES = 1.21
# G2M checkpoint: p-value = 0.04, NES = 1.20

rank = results$stat; names(rank) = rownames(results)
p4 = plotEnrichment(pathways[["HALLMARK_G2M_CHECKPOINT"]], rank) + labs(title = "G2M Checkpoint")
p5 = plotEnrichment(pathways[["HALLMARK_E2F_TARGETS"]], rank) + labs(title = "E2F Targets")

p = plot_grid(p1, p2, p3, p4, p5, NULL, nrow = 2, ncol = 3, rel_widths = c(2,2,2, 3, 3, 0), rel_heights = c(1,1, 1, 1, 1, 1), labels = c("a)", "b)", "c)", "d)", "e)", ""))

saveRDS(
  list(
    methylation = methylation_data %>%
      dplyr::select(grade_group, epiCMIT, epiCMIT.hyper, epiCMIT.hypo),
    gsea_rank = data.frame(gene = names(rank), stat = rank)
  ),
  file.path(outdir, "SF6_source_data.rds")
)
ggsave(file.path(figdir, "SF6_proliferation_markers_by_grade_group.png"), p, width = 180/25.3, height = 120/25.3, dpi = 300)
ggsave(file.path(figdir, "SF6_proliferation_markers_by_grade_group.pdf"), p, width = 180/25.3, height = 120/25.3)
