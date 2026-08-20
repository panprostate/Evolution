#' Run differential methylation analysis
#' conda activate rnbeads
#' 

rm(list = ls(all = TRUE))

# PACKAGES
library(qs)
library(dplyr)
library(readr)
library(RnBeads)
library(RnBeads.hg19)
library(hexbin)
library(ggrastr)
library(lemon)
library(UpSetR)
library(forcats)
library(cowplot)
library(tidyr)

source("src/utils.R")
source("src/methylation_helpers.R")
source("src/dmc.R")
source("src/plot_theme.R")
source("src/plot_functions.R")

# Define input files
f_annotation <- "/storage/scratch01/users/nemensha/ppcg_methylation/ext/PPCG_betas_anno.ATAC.Pom.PMD.PMI.confounder_correlation_from_ct11.1488.hypoxia.GE.meth.correlation.928.samples.RUV1.2.1.3_k10.RDS"
f_methylation <- "/storage/scratch01/users/nemensha/ppcg_methylation/ext/rnbSet_unnormalized/"
trajectories_fps <- list.files("outputs/02_trajectories/", pattern = "PPCG_Feb2026_.*mergedseg_with_clonality.txt", full = TRUE)
outdir <- "outputs/02_trajectories/methylation"
figdir <- "figures/02_trajectories/methylation"
dir.create(outdir, showWarnings = FALSE); dir.create(figdir, showWarnings = FALSE);

# Ensure output directory exists
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# 1. RUN DIFFERENTIAL METHYLATION ANALYSIS WITH RnBeads
# -----------------------------------------------------------------------

# Setup RnBeads options
RnBeads::rnb.options(disk.dump.big.matrices = TRUE)
RnBeads::rnb.options(region.types = character(0))

# Load RnBeads dataset and metadata
rnb_set = RnBeads::load.rnb.set(f_methylation)

# Load CpG annotations
cpg_meta = readRDS(f_annotation)

# Load data for differential analysis
tj_df = load_trajectories(trajectories_fps)

# Filter RnBeads dataset to samples in PPCG dataset and samples that pass methylation QC
bool_drop_samples <- !(rnb_set@pheno$Matching_WGS_Sample %in% tj_df$sample)
bool_drop_methylation_qc <- find_methylation_qc_fail_samples(rnb_set@pheno)
rnb_set <- RnBeads::remove.samples(
  rnb_set,
  bool_drop_samples | bool_drop_methylation_qc
)

# Annotate the metadata
sample_ids <- rnb_set@pheno$Matching_WGS_Sample
rnb_set@pheno <- dplyr::left_join(rnb_set@pheno, dplyr::select(tj_df, Matching_WGS_Sample = sample, Ordering = trajectory))
stopifnot(
  all(sample_ids == rnb_set@pheno$Matching_WGS_Sample)
)

# Get RnBeads settings
set_rnbeads_analysis_options(
  analysis_name = "PPEVO",
  identifiers_column = "Matching_WGS_Sample"
)

# rnb_result <- rnb.run.analysis(data.source = rnb_set, data.type = "rnb.set", dir.reports = file.path(outdir, "rnb_dmc"))
# qs2::qs_save(rnb_result, fs::path(outdir, "rnb_result.qs2"))
rnb_result = qs2::qs_read(file.path(outdir, "rnb_result.qs2"))

# ------------------------------------------------------------------------
# 2. UPSET PLOT WITH DIFFERENTIALLY METHYLATED REGIONS
# ------------------------------------------------------------------------
# comparison 1: Canonical (O-I) vs Alternative (O-II); 
# comparison 2: Canonical (O-I) vs Gain-enriched (O-III); 
# comparison 3: Alternative (O-II) vs Gain-enriched (O-III);
dm1 <- data.table::fread(file.path(outdir, "rnb_dmc", "differential_methylation_data", "diffMethTable_site_cmp1.csv"))
dm2 <- data.table::fread(file.path(outdir, "rnb_dmc", "differential_methylation_data", "diffMethTable_site_cmp2.csv"))
dm3 <- data.table::fread(file.path(outdir, "rnb_dmc", "differential_methylation_data", "diffMethTable_site_cmp3.csv"))

# Find significant DMCs for each comparison: abs mean difference > 0.2 and q-value < 0.05
dmc1 <- dm_filter(dm1)
dmc2 <- dm_filter(dm2)
dmc3 <- dm_filter(dm3)

# UPSET PLOT WITH DIFFERENTIALLY METHYLATED REGIONS

# Extract cgids based on methylation direction
dmcs = rbind(dmc1 %>% dplyr::mutate(type = ifelse(mean.diff < 0, "Hypo", "Hyper")) %>% dplyr::select(cgid, type) %>% dplyr::mutate(comparison = "I_vs_II"),
      dmc2 %>% dplyr::mutate(type = ifelse(mean.diff < 0, "Hypo", "Hyper")) %>% dplyr::select(cgid, type) %>% dplyr::mutate(comparison = "I_vs_III"),
      dmc3 %>% dplyr::mutate(type = ifelse(mean.diff < 0, "Hypo", "Hyper")) %>% dplyr::select(cgid, type) %>% dplyr::mutate(comparison = "II_vs_III")
)

dmcs_wide = dmcs %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols = c(cgid, type),       # Keep ID and metadata (Hypo/Hyper) as rows
    names_from = comparison,       # The sets we are intersecting
    values_from = present, 
    values_fill = list(present = 0) # Fill missing combinations with 0
  ) %>%
  as.data.frame() # CRITICAL: UpSetR requires a base data.frame, not a tibble/data.table

if (nrow(dmc3) == 0){
    # if no dmcs for last comparison, manually add the column with 0s
    dmcs_wide$II_vs_III <- 0
}

# Identify the columns that represent our sets (all except cgid and type)
set_names <- setdiff(colnames(dmcs_wide), c("cgid", "type"))

meth_colors <- c("hypo" = "#77AADD", "hyper" = "#EE8866")

upset_plot = upset(
  dmcs_wide,
  sets = set_names,
  keep.order = TRUE,
  empty.intersections = "on", 
  intersections = list(
    list("I_vs_II"),
    list("I_vs_III"),
    list("I_vs_II", "I_vs_III"),
    list("II_vs_III")
  ),
  nintersects = NA, # Override the default top 40 cutoff so the 0-count isn't dropped
  # Typography and labels
  mainbar.y.label = "Number of Differentially Methylated CpGs",
  sets.x.label = "Total DMCs per comparison",
  text.scale = c(1.3, 1.3, 1.3, 1, 1.5, 1), 
  query.legend = "bottom",
  main.bar.color = meth_colors["hypo"], 
  # The Queries
  queries = list(
    list(
      query = elements, 
      params = list("type", "Hypo"), 
      color = meth_colors["hypo"], 
      active = TRUE, 
      query.name = "Hypomethylated"
    ),
    list(
      query = elements, 
      params = list("type", "Hyper"), 
      color = meth_colors["hyper"], 
      active = TRUE, 
      query.name = "Hypermethylated"
    )
  )
)

write_tsv(dmcs_wide, file.path(outdir, "EXDF3a_source_data.tsv"))
pdf(file.path(figdir, "EXDF3a_upset_plot.pdf"), width = 100/25.3, height = 125/25.3)
print(upset_plot)
dev.off()

# ------------------------------------------------------------------------
# 3. CORRELATION OF DIFFERENTIAL METHYLATION BY EVOTYPES
# ------------------------------------------------------------------------
all(dm1$cgid == dm2$cgid) & all(dm1$cgid == dm3$cgid)
methylation_differences = data.frame(
    cgid = dm1$cgid,
    cmp1_mean_diff = dm1$mean.diff,
    cmp2_mean_diff = dm2$mean.diff,
    cmp3_mean_diff = dm3$mean.diff
)

p = ggplot(methylation_differences, aes(x=cmp1_mean_diff, y=cmp2_mean_diff)) +
    rasterise(geom_point(alpha = 0.05)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="blue") +
    geom_vline(xintercept=0, linetype="dashed", color="grey40") +
    ggpubr::stat_cor(method = "pearson") +
    xlab("Mean methylation difference\n (Canonical vs Alternative)") +
    ylab("Mean methylation difference\n (Canonical vs Gain-enriched)") + 
    xlim(c(-0.4, 0.4)) + ylim(c(-0.4, 0.4)) 

p = axes2lemon(p, lt = "both", bt = "both")

write_tsv(
  methylation_differences %>% dplyr::select(cmp1_mean_diff, cmp2_mean_diff),
  file.path(outdir, "EXDF3b_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF3b_methylation_differences_cmp1_cmp2.pdf"), p, width = 60/25.3, height = 60/25.3)

# ------------------------------------------------------------------------
# 4. DIFFERENTIALLY METHYLATED REGIONS IN DISTINCT CPG DENSE CONTEXTS
# ------------------------------------------------------------------------
annotation <- cpg_meta %>% dplyr::mutate(cgid = IlmnID)

# Assert that all comparisons share CpG sites
records_identical <- all(dm1$cgid == dm2$cgid) & all(dm1$cgid == dm3$cgid)
stopifnot(records_identical)

# Add annotations and summarise for each differential methylation cpg set
ds1 <- dm_filter(dm1) %>% dmc_annotation_reporter(annotation)
ds2 <- dm_filter(dm2) %>% dmc_annotation_reporter(annotation)
ds3 <- dm_filter(dm3) %>% dmc_annotation_reporter(annotation)

# Write files to output directory. This writes to secondary_analysis with cgi, tss, etc. as calculated by annotation_reporter
write_annotation_summary(ds1, outdir)
write_annotation_summary(ds2, outdir)
write_annotation_summary(ds3, outdir)

cgi_data <- read_dmc_files(outdir, "cgi")

region_colours <- c(
  "Island"   = "#8dd3c7", # Soft Teal
  "Shore"    = "#bebada", # Soft Lavender
  "Shelf"    = "#fb8072", # Soft Coral
  "Open Sea" = "#80b1d3"  # Soft Steel Blue
)

donut_data <- cgi_data %>%
  mutate(
    cpg_island_label = factor(
      cpg_island_label, 
      levels = c("Island", "Shore", "Shelf", "Open Sea")
    )
  )

donut_plot <- ggplot(
  donut_data, 
  # Map 'perc' to Y to determine slice sizes
  # Map a hardcoded X value (e.g., 2) to place the ring
  aes(x = 2, y = n, fill = cpg_island_label)
) +
  # Draw the stacked bars with a crisp white border for separation
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.4) +
  # Wrap the Y-axis into a circle
  coord_polar(theta = "y", start = 0) +  
  xlim(0.5, 2.5) +
  # Create a compact grid: Regions across columns, Comparisons down rows
  facet_grid(comparison_name ~ mean.diff.direction) +
  # Apply our shared colors
  scale_fill_manual(values = region_colours, name = "CGI region") +
  # Strip away ALL background elements, axes, and grids
  theme_void() +
  # Add back carefully styled text for the facet labels
  theme(
    text = element_text(family = "Helvetica", size = 8),
    strip.text.x = element_text(
      face = "bold", size = 9, margin = margin(b = 5, t = 10)
    ),
    strip.text.y = element_text(
      face = "bold", size = 9, angle = 270, margin = margin(l = 5, r = 10)
    ),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

# Render the plot
write_tsv(donut_data, file.path(outdir, "EXDF3c_source_data.tsv"))
ggsave(donut_plot, filename = file.path(figdir, "EXDF3c_cpg_annotations.pdf"), width = 70/25.4, height = 70/25.4)

# ------------------------------------------------------------------------
# 5. DIFFERENTIALLY METHYLATED REGIONS IN REGULATORY REGIONS
# ------------------------------------------------------------------------
atac_data <- read_dmc_files(outdir, "atac") %>% 
    # keep only DMCs overlapping ATAC
    dplyr::filter(Corces_ATAC == "ATAC_Corces") %>% 
    dplyr::mutate(subtype = "ATAC peaks") %>% 
    dplyr::select(element = Corces_ATAC, subtype, mean.diff.direction, n, comparison_name)

pmd_data <-  read_dmc_files(outdir, "pmd") %>%
  # keep only DMCs overlapping PMDs
  dplyr::filter(Prostate_PMD_Guo == "Prostate_PMD_Guo") %>%
  dplyr::mutate(subtype = "PMDs") %>%
  dplyr::select(element = Prostate_PMD_Guo, subtype, mean.diff.direction, n, comparison_name)

pomerantz_data <- read_dmc_files(outdir, "chromhmm") %>%
  # keep only DMCs overlapping Pomerantz regulatory elements
  dplyr::filter(ChromHMM_Pomerantz != "") %>%
  dplyr::mutate(element = "Pomerantz_Regulatory_Element") %>%
  dplyr::select(element, subtype = ChromHMM_Pomerantz, mean.diff.direction, n, comparison_name)

reg_elements = rbind(atac_data, pmd_data, pomerantz_data)

make_bars <- function(data, meth_colors){
    ggplot(data %>% mutate(subtype = fct_reorder(subtype, n)), aes(x = n, y = subtype, fill = mean.diff.direction)) +
        geom_col() +
        scale_fill_manual(values = meth_colors, name = "Methylation Status") +
        labs(
            x = "Number of DMCs",
            y = NULL,
            title = ""
        ) + 
        theme(legend.position = "none")
}

names(meth_colors) <- c("hypo", "hyper")
p1 = make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 3", element == "ATAC_Corces"), meth_colors)
p2 = make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 3", element == "Prostate_PMD_Guo"), meth_colors)
p3 = make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 3", element == "Pomerantz_Regulatory_Element"), meth_colors)
p4 = make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 2", element == "ATAC_Corces"), meth_colors)
p5 = make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 2", element == "Prostate_PMD_Guo"), meth_colors)
p6 = make_bars(reg_elements %>% filter(comparison_name == "Ordering 1_vs_Ordering 2", element == "Pomerantz_Regulatory_Element"), meth_colors)

p = cowplot::plot_grid(
    p3, p2, p1, p6, p5, p4, nrow = 3, byrow = FALSE, align = "hv", 
    rel_heights = c(10, 3, 3, 10, 3, 3), rel_widths = c(1, 1, 1, 1, 1, 1)
)

write_tsv(reg_elements, file.path(outdir, "EXDF3d_source_data.tsv"))
ggsave(p, filename = file.path(figdir, "EXDF3d_dmc_regulatory_elements.pdf"), width = 150/25.3, height = 100/25.3)
