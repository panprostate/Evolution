# Analysis of number of subclones and proportion clonal mutations per Gleason

rm(list = ls(all = TRUE))


# PACKAGES
library(readr)
library(dplyr)
library(stringr)
library(patchwork)
library(ggsignif)
library(nlme)
library(DescTools)
library(lemon)

## PATHS
outdir <- "outputs/01_landscape/subclonality"
figdir <- "figures/01_landscape/subclonality"
dir.create(outdir); dir.create(figdir)
evo_metrics_fp <- "outputs/00_preprocessing/evo_metrics.tsv"
clinfile <- "data/meta/PPCG_donors_clin_20241217.csv"

## FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

## LOAD DATA
evo_metrics <- read_delim(evo_metrics_fp)
clin_data <- read_delim(clinfile)

## Filter to QC pass samples and classify by Grade Group / met status
smps <- get_prim_smps()
evo_metrics <- evo_metrics[evo_metrics$smp %in% str_remove(smps, "_DNA"), ]
# number of subclones in primary tumour samples
table((evo_metrics$nsubclones > 0)[evo_metrics$smp %in% str_remove(get_prim_smps(), "_DNA")])

evo_metrics$group <- get_gs_group(evo_metrics$smp)

## ANALYSES

evo_metrics$patient_id <- extract_ppcg_pt(evo_metrics$smp)
evo_metrics$group <- factor(
  evo_metrics$group,
  levels = c("1", "2", "3", "4+"),
  ordered = TRUE
)

# Convert to numeric for trend test
evo_metrics$group_n <- as.numeric(evo_metrics$group)

# Fit mixed-effects model
model <- nlme::lme(
  mut_ith ~ group_n,                 # fixed effect: ordinal progression
  random = ~1 | patient_id,          # random intercept per patient
  data = evo_metrics[!is.na(evo_metrics$group),],
  method = "REML"
)

summary(model)
beta1 <- summary(model)$tTable["group_n", "Value"]
pval <- summary(model)$tTable["group_n", "p-value"] / 2  # one-sided

# Jonckheere-Terpstra test
JonckheereTerpstraTest(evo_metrics$mut_ith, evo_metrics$group_n) # p = 5.1e-10

# Number of samples with subclonal expansions per Gleason
tb <- table(evo_metrics$nsubclones > 0, evo_metrics$group)

# Plot percentage of subclonal alterations
p <- ggplot(evo_metrics[!is.na(evo_metrics$group),], aes(x = group, y = mut_ith*100)) + 
    geom_boxplot() + 
    geom_point(size = 1, alpha = .3, position = position_dodge2(width = 0.2)) + 
    labs(x = "", y = "Subclonal mutations (%)") + 
    ylim(c(0, 100))
  
p <- axes2lemon(p)

save_ggplot(p, file.path(figdir, "mutation_ith_by_group"), w = 50, h = 50)

## MAIN PLOT
# 1. CLEAN & ORDER DATA

# Add clinical metadata on top
evo_metrics = evo_metrics[evo_metrics$group != "Metastasis",]
evo_metrics$record_id = gsub("[a-z]$", "", evo_metrics$smp)
evo_metrics = left_join(evo_metrics, clin_data %>% dplyr::select(-group), by = "record_id")

evo_metrics$psa_simplified = ifelse(
    evo_metrics$psa_at_tumour_collection < 10, "<10", 
    ifelse(evo_metrics$psa_at_tumour_collection > 20, ">20", ifelse(
        evo_metrics$psa_at_tumour_collection >= 10 & evo_metrics$psa_at_tumour_collection <= 20, ">10", NA
    )))
    
evo_metrics$t_stage = gsub("[a-z]", "", evo_metrics$path_t_stage)
evo_metrics = mutate(evo_metrics, age = case_when(
    age_at_tumour_collection < 55 ~ "<55", 
    age_at_tumour_collection > 65 ~ ">65", 
    age_at_tumour_collection >= 55 & age_at_tumour_collection <= 65 ~ "55-65", 
    TRUE ~ NA
))

evo_metrics$MFS = evo_metrics$mets_ind == "mets"
evo_metrics$ancestry = evo_metrics$predicted_ancestry

# Filter data with complete group info for display
plot_data <- evo_metrics[!is.na(evo_metrics$group),]

write_delim(plot_data, file.path(outdir, "plot_data_subclonal_landscape.tsv"), delim = "\t")

# Create explicit ordering: Sort by Group, then by Mutation ITH (descending)
plot_data <- plot_data %>%
  arrange(group, desc(mut_ith)) %>%
  mutate(smp = factor(smp, levels = unique(smp))) # Lock factor levels

# 2. DEFINE COLORS (Nature-style palettes)
# Define specific palettes for clinical features to ensure distinct visual separation
col_genome <- "#E76A85" 

# PSA: Light Grey-Blue -> Deep Blue
cols_psa <- c(
  "<10" = "#F0F4F8",  # Very light grey-blue
  ">10" = "#90CAF9",  # Medium blue
  ">20" = "#0D47A1"   # Deep navy
)

# Stage: Light Grey-Green -> Deep Green
cols_stage <- c(
  "T2" = "#F1F8E9",   # Very light grey-green
  "T3" = "#AED581",   # Medium green
  "T4" = "#33691E"    # Deep forest green
)

# Age: Light Grey-Purple -> Deep Purple
# (Assuming bins are <55, 55-65, >65)
cols_age <- c(
  "<55"   = "#F3E5F5", # Very light lavender
  "55-65" = "#BA68C8", # Medium orchid
  ">65"   = "#4A148C"  # Deep purple
)

# MFS: Binary (Grey vs. Signal Red)
cols_mfs <- c(
  "FALSE" = "#F5F5F5", # Light Grey (Noise)
  "TRUE"  = "#B71C1C"  # Deep Red (Signal)
)

# Ancestry: Distinct Qualitative (NPG style)
cols_anc <- c(
  "EUR" = "#F0F4F8", # Very light grey-blue (The dominant background color)
  "AFR" = "#E64B35",  # Muted Red
  "AMR" = "#F39B7F",  # Muted Orange/Salmon
  "EAS" = "#00A087",  # Deep Teal
  "SAS" = "#8491B4"   # Periwinkle/Muted Purple
)

# 2. PLOT GENERATION
# ---------------------------------------------------------

# Top: Proportion Subclonal (Re-using genomic color)
p_snvith <- ggplot(plot_data) + 
    geom_col(aes(x = smp, y = 1), fill = "#FAFAFA", width = 1) + # Lightest grey bg
    geom_col(aes(x = smp, y = mut_ith), fill = col_genome, width = 1) + 
    facet_grid(.~group, scales = "free_x", space = "free") + 
    labs(y = "Prop.\nSubclonal") + 
    theme_void() + 
    theme(
        strip.text = element_text(face = "bold", size = 12),
        axis.title.y = element_text(angle = 90, size = 9, margin = margin(r=5)),
        legend.position = "none"
    )

# Middle: Number of Clones (NOW MATCHING COLOR)
p_nsub <- ggplot(plot_data) + 
    # Optional: Log scale often looks better for counts, but linear is fine too
    geom_col(aes(x = smp, y = nsubclones), fill = col_genome, width = 1) + 
    facet_grid(.~group, scales = "free_x", space = "free") + 
    labs(y = "#\nClones") + 
    theme_void() + 
    theme(
        strip.text = element_blank(),
        axis.title.y = element_text(angle = 90, size = 9, margin = margin(r=5)),
        legend.position = "none"
    )

# Clinical Strips (Helper Function)
plot_clinical_strip <- function(data, var_col, var_lab, palette) {
  ggplot(data, aes(x = smp, y = 1, fill = as.character(!!sym(var_col)))) +
    geom_tile() +
    scale_fill_manual(values = palette, na.value = "#FFFFFF", name = var_lab) +
    facet_grid(.~group, scales = "free_x", space = "free") +
    labs(y = var_lab) +
    theme_void() +
    theme(
        strip.text = element_blank(),
        axis.title.y = element_text(angle = 0, size = 8, hjust = 1, margin = margin(r=5)), 
        legend.position = "right", 
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(l=0, r=0, t=0, b=0),
        legend.title = element_text(size = 8, face="bold"),
        legend.text = element_text(size = 7)
    )
}

# Generate Strips with new colors
p_clin1 <- plot_clinical_strip(plot_data, "psa_simplified", "PSA", cols_psa)
p_clin2 <- plot_clinical_strip(plot_data, "t_stage", "Stage", cols_stage)
p_clin3 <- plot_clinical_strip(plot_data, "age", "Age", cols_age)
p_clin4 <- plot_clinical_strip(plot_data, "MFS", "MFS", cols_mfs)
p_clin5 <- plot_clinical_strip(plot_data, "ancestry", "Ancestry", cols_anc)

# 3. ASSEMBLY
# ---------------------------------------------------------
final_plot <- p_snvith / p_nsub / p_clin1 / p_clin2 / p_clin3 / p_clin4 / p_clin5 + 
  plot_layout(heights = c(4, 2, 0.5, 0.5, 0.5, 0.5, 0.5), guides = "collect") & # Collect gathers all legends to the side
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

write_tsv(
  plot_data %>% dplyr::select(smp, group, mut_ith, nsubclones,
                               psa_simplified, t_stage, age, MFS, ancestry),
  file.path(outdir, "Fig1a_source_data.tsv")
)
save_ggplot(final_plot, file.path(figdir, "Fig1a_subclonality_landscape"), w = 180, h = 65)

