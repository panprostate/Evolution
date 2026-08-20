### Calculate ssGSEA scores in PPCG 

rm(list = ls(all = TRUE))  # Clear environment

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)   # Data manipulation and piping
library(GSVA)        # Gene Set Variation Analysis (ssGSEA)
library(msigdbr)     # MSigDB gene sets
library(DESeq2)      # Differential expression analysis
library(qs)          # Fast serialization
library(edgeR)       # CPM calculation for filtering
library(fgsea)       # Fast GSEA

# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")
source("src/rna_functions.R")


# PATHS -------------------------------------------------------------------
outdir <- "outputs/02_trajectories/rnaseq"
figdir <- "figures/02_trajectories/rnaseq"
dir.create(outdir); dir.create(figdir)

# Filepaths for RNA counts, pathways, and trajectory data
rna_fp <- "data/raw/RNA/PPCG_RNAseq_RUVIIIPRPS_Supervised.qs"
pathways_fp <- "data/meta/pathways/pathways.tsv"
trajectories_fps <- list.files("outputs/02_trajectories", pattern = "PPCG_Feb2026.*mergedseg_with_clonality.txt", full = T)

## LOAD DATA
rna_counts = qread(rna_fp)  # Load RNA counts (SummarizedExperiment or similar)
pathways = read_delim(pathways_fp)  # Load pathway gene sets
pathways = pathways %>% dplyr::filter(str_detect(pathway, "HALLMARK")|
                                       str_detect(pathway, "SETLUR") |
                                       str_detect(pathway, "HIERONYMOUS"))  # Filter to relevant pathways
                                       
tj_samples = load_trajectories(trajectories_fps) %>% 
  mutate(patient = extract_ppcg_pt(sample))  %>% # Load trajectory info 
  mutate(grade_group = get_gs_group(sample))

# Prepare pathways as a named list of gene vectors for GSVA
pathways = lapply(unique(pathways$pathway), function(pathway){
  return(pathways$gene[pathways$pathway == pathway])
}) %>% setNames(unique(pathways$pathway))

# BUILD DESeq2 DESIGN
# Get raw counts and filter to samples with trajectory data
idx = rna_counts@colData[["Matching_WGS_Sample.ID.ADD.NA."]] %in% tj_samples$sample

countData = assays(rna_counts)[["Raw.count"]][,idx]  # Raw counts matrix for selected samples
countData = expression_filter(countData)            # Filter out lowly expressed genes

# Combine colData with RUVIII.PRPS adjustment factors
colData = cbind(
  rna_counts@colData[idx,], 
  rna_counts@metadata[["RUVIII.PRPS.AdjustedForPurity.Supervised.W"]][idx,]
)

# Define technical covariates (RUV factors)
technical_covariates = colnames(rna_counts@metadata[["RUVIII.PRPS.AdjustedForPurity.Supervised.W"]])

# Add trajectory variables to colData
m = match(colData[["Matching_WGS_Sample.ID.ADD.NA."]], tj_samples$sample)
stopifnot(all(colData[["Matching_WGS_Sample.ID.ADD.NA."]] == tj_samples$sample[m]))  # Confirm matching samples

colData$trajectory = tj_samples$trajectory[m]
colData$grade_group = tj_samples$grade_group[m]
colData$is_tj1 = colData$trajectory == "Ordering 1"
colData$is_tj2 = colData$trajectory == "Ordering 2"
colData$is_tj3 = colData$trajectory == "Ordering 3"

# RUN differential expression analysis for each trajectory vs others
# TJ1 vs ALL
deseq_tj1_vs_all = run_deseq2("is_tj1", countData, colData, technical_covariates)
saveRDS(deseq_tj1_vs_all, file.path(outdir, "deseq2_tj1_vs_all.rds"))
gsea_tj1_vs_all = run_gsea(deseq_tj1_vs_all, pathways)
saveRDS(gsea_tj1_vs_all, file.path(outdir, "gsea_tj1_vs_all.rds"))

# TJ2 vs ALL
deseq_tj2_vs_all = run_deseq2("is_tj2", countData, colData, technical_covariates)
saveRDS(deseq_tj2_vs_all, file.path(outdir, "deseq2_tj2_vs_all.rds"))
gsea_tj2_vs_all = run_gsea(deseq_tj2_vs_all, pathways)
saveRDS(gsea_tj2_vs_all, file.path(outdir, "gsea_tj2_vs_all.rds"))

# TJ3 vs ALL
deseq_tj3_vs_all = run_deseq2("is_tj3", countData, colData, technical_covariates)
saveRDS(deseq_tj3_vs_all, file.path(outdir, "deseq2_tj3_vs_all.rds"))
gsea_tj3_vs_all = run_gsea(deseq_tj3_vs_all, pathways)
saveRDS(gsea_tj3_vs_all, file.path(outdir, "gsea_tj3_vs_all.rds"))


# Add labels
gsea_tj1_vs_all$evotype <- "TJ1"
gsea_tj2_vs_all$evotype <- "TJ2"
gsea_tj3_vs_all$evotype <- "TJ3"

# Combine into one long data frame
gsea_all <- bind_rows(gsea_tj1_vs_all, gsea_tj2_vs_all, gsea_tj3_vs_all)

gsea_plot_df <- gsea_all %>%
  dplyr::select(pathway, NES, padj, evotype) %>%
  filter(!grepl("biscut", pathway)) %>%  # Filter out biscut pathways
  mutate(is_significant = ifelse(padj < 0.05, "significant", "ns")) %>% 
  mutate(log10padj = -log10(padj + 1e-10))

gsea_plot_df$pathway <- str_replace(gsea_plot_df$pathway, "HALLMARK_", "H_")

p <- ggplot(gsea_plot_df, aes(x = pathway, y = evotype)) +
  geom_point(aes(size = log10padj, fill = NES, color = is_significant), shape = 21) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(1, 5)) +
  scale_color_manual(values = c("significant" = "black", "ns" = "lightgrey")) +
  labs(x = "Pathway",
       y = "Trajectory",
       fill = "NES",
       size = "-log10(padj)") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  
    legend.position = "none"
  )

write_tsv(gsea_plot_df, file.path(outdir, "EXDF4_source_data.tsv"))
ggsave(
  file.path(figdir, "EXDF4_dea_evotypes.pdf"),
  plot = p,
  width = 180 / 25.4,
  height = 80 / 25.4
)


p <- ggplot(gsea_plot_df, aes(x = pathway, y = evotype)) +
  geom_point(aes(size = log10padj, fill = NES, color = is_significant), shape = 21) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(1, 5)) +
  scale_color_manual(values = c("significant" = "black", "ns" = "lightgrey")) +
  labs(x = "Pathway",
       y = "Trajectory",
       fill = "NES",
       size = "-log10(padj)") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  
    legend.position = "bottom"
  )

ggsave(
  file.path(figdir, "dea_evotypes_legend.pdf"),
  plot = p,
  width = 180 / 25.4,  
  height = 150 / 25.4  
)