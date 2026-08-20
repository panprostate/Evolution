#' Analysis of TCGA Molecular Subtype and You et al Subtypes overlap with PPCG Evolution Trajectories
#' Then quickly map to Woodcock et al. trajectories and map across countries
#'
#' This script performs comprehensive analysis of the overlap between TCGA prostate 
#' cancer molecular subtypes and expression subtypes with the PPCG evolution trajectory classifications (evotypes).
#' The analysis implements the molecular taxonomy described in "The Molecular Taxonomy 
#' of Prostate Cancer" (Cell 2015) and integrates expression subtypes from You et al. 
#' (Cancer Research 2016).
#' 
#' Key Analyses:
#' 1. Classification into 8 TCGA molecular subtypes:
#'    - ETS fusions (TMPRSS2-ERG, -ETV4, -ETV1, -FLI1)
#'    - Key mutations (SPOP, FOXA1, IDH1)
#'    - Other/unclassified tumors
#' 2. Integration with RNA expression subtypes from You et al.
#' 3. Visualization of subtype distributions across evolution trajectories
#' 
#' Input files required:
#'   - PPCG trajectory classifications
#'   - ETS fusion status calls
#'   - Somatic driver mutation annotations
#'   - RNA expression subtype classifications
#' 
#' Outputs:
#'   - Annotated trajectory data with molecular subtypes
#'   - Donut plots showing subtype distributions by evotype
#'   - Statistical summaries of subtype overlaps

# Clear workspace for clean analysis environment
rm(list = ls(all = TRUE))

# =============================================================================
# PACKAGE DEPENDENCIES
# =============================================================================

library(tidyverse)      # Data manipulation and visualization
library(RColorBrewer)   # Color palettes for plots
library(ggalluvial)    # Alluvial plots for subtype overlaps

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Load custom utility functions and plotting themes
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

#' Find Mutated Samples for a Specific Gene
#' 
#' Identifies PPCG samples that carry mutations in a specified gene from
#' the driver mutations dataset.
#' 
#' @param driver_mutations_data data.frame containing mutation data
#' @param gene_symbol Character. Gene symbol to search for mutations
#' @return Character vector of PPCG sample IDs with mutations in the specified gene
find_mutated_samples <- function(driver_mutations_data, gene_symbol) {
    # Extract PPCG sample columns from the dataset
    ppcg_sample_columns <- colnames(driver_mutations_data) %>% 
        {.[grepl("PPCG", .)]}
    
    # Find samples with mutations (indicated by non-"." values)
    mutation_indices <- driver_mutations_data[driver_mutations_data$symbol == gene_symbol, ppcg_sample_columns] %>%
        {. != "."}
    
    return(ppcg_sample_columns[mutation_indices])
}

#' Find Samples with Specific ETS Fusions
#' 
#' Identifies samples carrying a specific ETS fusion from the fusion status data.
#' 
#' @param ets_fusion_data data.frame containing ETS fusion status
#' @param fusion_gene Character. Gene involved in the fusion (e.g., "ERG", "ETV4")
#' @return Character vector of PPCG donor IDs with the specified fusion
find_fusion_samples <- function(ets_fusion_data, fusion_gene) {
    return(ets_fusion_data$PPCG_Donor_ID[ets_fusion_data[[fusion_gene]]])
}

#' Classify Samples into TCGA Molecular Subtypes
#' 
#' Implements the 8-subtype classification system from "The Molecular Taxonomy 
#' of Prostate Cancer" (Cell 2015). Classifies samples into one of:
#' 1. TMPRSS2-ERG fusions
#' 2. TMPRSS2-ETV4 fusions  
#' 3. TMPRSS2-ETV1 fusions
#' 4. TMPRSS2-FLI1 fusions
#' 5. SPOP mutations
#' 6. FOXA1 mutations
#' 7. IDH1 mutations
#' 8. Other/unclassified
#' 
#' @param ets_fusion_data data.frame containing ETS fusion status calls
#' @param driver_mutations_data data.frame containing driver mutation annotations
#' @return List containing vectors of sample IDs for each molecular subtype
classify_tcga_molecular_subtypes <- function(ets_fusion_data, driver_mutations_data) {
    
    # ETS fusion subtypes - identify samples with specific fusion genes
    erg_fusion_samples <- find_fusion_samples(ets_fusion_data, "ERG")
    etv4_fusion_samples <- find_fusion_samples(ets_fusion_data, "ETV4")
    etv1_fusion_samples <- find_fusion_samples(ets_fusion_data, "ETV1")
    fli1_fusion_samples <- find_fusion_samples(ets_fusion_data, "FLI1")

    # Key driver mutation subtypes - identify samples with functional mutations
    spop_mutated_samples <- find_mutated_samples(driver_mutations_data, "SPOP") %>% 
        extract_ppcg_pt()
    foxa1_mutated_samples <- find_mutated_samples(driver_mutations_data, "FOXA1") %>% 
        extract_ppcg_pt()
    idh1_mutated_samples <- find_mutated_samples(driver_mutations_data, "IDH1") %>% 
        extract_ppcg_pt()

    # Return structured list of molecular subtypes
    return(list(
        erg_fusions = erg_fusion_samples, 
        etv4_fusions = etv4_fusion_samples, 
        etv1_fusions = etv1_fusion_samples, 
        fli1_fusions = fli1_fusion_samples, 
        spop_mutants = spop_mutated_samples, 
        foxa1_mutants = foxa1_mutated_samples, 
        idh1_mutants = idh1_mutated_samples
    ))
}

#' Assess Overlap Between Molecular Subtypes
#' 
#' Creates a matrix showing the number of overlapping samples between
#' different molecular subtypes. This is used as a quality control check
#' to ensure subtypes are mutually exclusive (minimal overlaps expected).
#' 
#' @param molecular_subtypes_list List containing sample vectors for each subtype
#' @return Matrix showing pairwise overlap counts between subtypes
assess_subtype_overlaps <- function(molecular_subtypes_list) {
    num_subtypes <- length(molecular_subtypes_list)
    overlap_matrix <- matrix(0, nrow = num_subtypes, ncol = num_subtypes)
    rownames(overlap_matrix) <- colnames(overlap_matrix) <- names(molecular_subtypes_list)

    # Calculate pairwise overlaps between all subtype combinations
    for (i in seq_along(molecular_subtypes_list)) {
        for (j in seq_along(molecular_subtypes_list)) {
            if (i != j) {
                overlap_matrix[i, j] <- length(intersect(
                    molecular_subtypes_list[[i]], 
                    molecular_subtypes_list[[j]]
                ))
            }
        }
    }
    return(overlap_matrix)
}

#' Annotate Trajectory Samples with TCGA Molecular Subtypes
#' 
#' Adds molecular subtype annotations to the trajectory samples dataset
#' based on the TCGA classification scheme. Samples are assigned to the
#' first matching subtype, with remaining samples classified as "Other".
#' 
#' @param trajectory_samples data.frame containing trajectory sample data
#' @param molecular_subtypes_list List containing sample vectors for each subtype
#' @return data.frame with added tcga_subtype column
annotate_molecular_subtypes <- function(trajectory_samples, molecular_subtypes_list) {
    trajectory_samples_annotated <- trajectory_samples %>% 
        mutate(patient_id = extract_ppcg_pt(sample)) %>% 
        mutate(tcga_subtype = case_when(
            patient_id %in% molecular_subtypes_list$erg_fusions ~ "ERG fusion", 
            patient_id %in% molecular_subtypes_list$etv4_fusions ~ "ETV4 fusion", 
            patient_id %in% molecular_subtypes_list$etv1_fusions ~ "ETV1 fusion", 
            patient_id %in% molecular_subtypes_list$fli1_fusions ~ "FLI1 fusion", 
            patient_id %in% molecular_subtypes_list$spop_mutants ~ "SPOP mutation", 
            patient_id %in% molecular_subtypes_list$foxa1_mutants ~ "FOXA1 mutation", 
            patient_id %in% molecular_subtypes_list$idh1_mutants ~ "IDH1 mutation", 
            TRUE ~ "Other"
        ))
    return(trajectory_samples_annotated)
}


# =============================================================================
# FILE PATH CONFIGURATION
# =============================================================================

# Define output directories
OUTPUT_DIR <- "outputs/02_trajectories"
FIGURE_DIR <- "figures/02_trajectories/subtypes"
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Define input file paths
trajectories_filepaths <- list.files(OUTPUT_DIR, pattern = "PPCG_Feb2026.*mergedseg_with_clonality.txt", full.names = TRUE)
ets_fusions_filepath <- "data/meta/ets_fusion_status/ets_calls.tsv"
expression_subtypes_filepath <- "data/processed/You_classes_ppcg.csv"
woodcock_filepath <- "data/meta/woodcock_subtypes.csv"

# Driver mutations file contains only mutations with functional impact
driver_mutations_filepath <- "data/raw/Somatic_drivers/DRIVER_ANNOTATIONS__2023-09-25_release4__SNV+indel__filteredVariants__driverGenes.txt"

# =============================================================================
# DATA LOADING AND PREPROCESSING
# =============================================================================
# Load trajectory classifications and molecular annotation data
trajectory_samples <- load_trajectories(trajectories_filepaths)
tj_lookup <- c("Ordering 1" = "Canonical", "Ordering 2" = "Alternative", "Ordering 3" = "Gain-enriched")
trajectory_samples$trajectory <- tj_lookup[trajectory_samples$trajectory]
driver_mutations_data <- read_delim(driver_mutations_filepath)
ets_fusion_data <- read_delim(ets_fusions_filepath, delim = "\t")
expression_subtypes_data <- read_delim(expression_subtypes_filepath, delim = ",")
woodcock_subtypes <- read_delim(woodcock_filepath, delim = ",")

# =============================================================================
# TCGA MOLECULAR SUBTYPE CLASSIFICATION
# =============================================================================

# Classify samples into TCGA molecular subtypes
tcga_molecular_subtypes <- classify_tcga_molecular_subtypes(ets_fusion_data, driver_mutations_data)

# Quality control: assess overlaps between subtypes (expecting minimal overlaps except SPOP-FOXA1)
subtype_overlap_matrix <- assess_subtype_overlaps(tcga_molecular_subtypes)
print("Overlap matrix between TCGA molecular subtypes:")
print(subtype_overlap_matrix)

# Add TCGA molecular subtype annotations to trajectory data
trajectory_samples_with_subtypes <- annotate_molecular_subtypes(trajectory_samples, tcga_molecular_subtypes)

# =============================================================================
# TCGA SUBTYPE DISTRIBUTION ANALYSIS BY EVOTYPE
# =============================================================================

# Calculate proportions of each TCGA subtype within each evolutionary trajectory
subtype_proportions_by_evotype <- trajectory_samples_with_subtypes %>% 
    group_by(trajectory) %>% 
    dplyr::count(tcga_subtype) %>%
    dplyr::mutate(proportion = n / sum(n)) %>% 
    arrange(trajectory, desc(proportion))

print("TCGA subtype proportions by evotype:")
print(subtype_proportions_by_evotype)
# Key findings from the analysis:
# O-I: 73.7% mapping to ETS fusions (68% ERG); 25.4% as "Other" (n = 115)
# O-II: 29.2% mapping to SPOP mutation; 53.1% (n = 69) as "Other"
# O-III: 26.2% mapping to ETS fusions (20.4% ERG); 65.3% as "Other" (n = 32)

# =============================================================================
# DONUT PLOT VISUALIZATION OF SUBTYPE DISTRIBUTIONS
# =============================================================================

# Prepare data for donut plot visualization
donut_plot_data <- trajectory_samples_with_subtypes %>%
    group_by(trajectory, tcga_subtype) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(trajectory) %>%
    mutate(proportion = count / sum(count)) %>%
    # Arrange is crucial for correct label positioning in stacked plots
    arrange(trajectory, desc(tcga_subtype)) %>% 
    mutate(
        # Calculate position for labels (center of the slice)
        label_y_position = cumsum(proportion) - 0.5 * proportion,
        # Create a clean label string (e.g., "15%")
        label_text = scales::percent(proportion, accuracy = 1)
    )
  
# Create donut plot showing TCGA subtype distributions by evolutionary trajectory
donut_plot <- ggplot(donut_plot_data, aes(x = 2, y = proportion, fill = tcga_subtype)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    facet_wrap(~ trajectory) +
    scale_fill_brewer(palette = "Set2") +
    # Limit x-axis to create the hole in the middle of the donut
    xlim(0.5, 2.5) + 
    theme_void(base_size = 9, base_family = "sans") + 
    theme(
        legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold", margin = margin(b = 10)),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9, face = "bold"),
        plot.margin = margin(10, 10, 10, 10) 
    ) + 
    labs(fill = "TCGA subtype") +
    guides(fill = guide_legend(
        nrow = 2,                 # Divide legend elements into 2 rows
        byrow = TRUE,             # Arrange elements by row
        title.position = "top",   # Place title above legend boxes
        title.hjust = 0.5         # Center the legend title
    )) + 
    # Optional: Add sample size (N=...) in the center hole
    geom_text(
        data = donut_plot_data %>% 
            group_by(trajectory) %>% 
            summarise(sample_count = sum(count)), 
        aes(x = 0.5, y = 0, label = paste0("N=", sample_count)), 
        inherit.aes = FALSE, 
        size = 3, 
        fontface = "italic", 
        color = "grey30"
    )

# Save donut plot
write_tsv(
  donut_plot_data %>% dplyr::select(trajectory, tcga_subtype, count, proportion),
  file.path(OUTPUT_DIR, "SF11a_source_data.tsv")
)
ggsave(
    file.path(FIGURE_DIR, "SF11a_evotype_overlap_tcga_subtypes_donut.pdf"),
    donut_plot,
    width = 225/45.3,
    height = 100/25.3
)


# =============================================================================
# RNA EXPRESSION SUBTYPE ANALYSIS (YOU ET AL. 2016)
# =============================================================================

# Integrate expression subtypes from You et al. (Cancer Research 2016)
# calculated by PPCG RNA working group
trajectory_samples_with_all_subtypes <- trajectory_samples_with_subtypes %>%
    mutate(
        patient_id = extract_ppcg_pt(sample),
        expression_subtype_you = expression_subtypes_data$You_class[
            match(patient_id, extract_ppcg_pt(expression_subtypes_data$Assay_ID))
        ]
    )

# Save comprehensive dataset with all molecular annotations
write_delim(
    trajectory_samples_with_all_subtypes, 
    file.path(OUTPUT_DIR, "PPCG_trajectories_with_tcga_subtypes_and_you_subtypes.txt"), 
    delim = "\t"
)

# Calculate proportions of You expression subtypes within each evolutionary trajectory
expression_subtype_proportions_by_evotype <- trajectory_samples_with_all_subtypes %>% 
    dplyr::filter(!is.na(expression_subtype_you)) %>%  
    group_by(trajectory) %>% 
    dplyr::count(expression_subtype_you) %>%
    dplyr::mutate(proportion = n / sum(n)) %>% 
    arrange(trajectory, desc(proportion))

print("You expression subtype proportions by evotype:")
print(expression_subtype_proportions_by_evotype)

# Key findings from the expression subtype analysis:
# O-I: 50.8% mapping to You subtype PCS1; 25.6% mapping to You subtype PCS2; 23.6% mapping to You subtype PCS3
# O-II: 46.8% mapping to You subtype PCS1; 19.0% mapping to You subtype PCS2; 34.2% mapping to You subtype PCS3
# O-III: 76.0% mapping to You subtype PCS1; 24.0% mapping to You subtype PCS3

# =============================================================================
# EXPRESSION SUBTYPE VISUALIZATION
# =============================================================================

# Create donut plot showing You expression subtype distributions by evolutionary trajectory
# Similar to TCGA subtype donut plot, but using expression_subtype_you instead of tcga_subtype
expression_donut_plot_data <- trajectory_samples_with_all_subtypes %>% 
    dplyr::filter(!is.na(expression_subtype_you)) %>%  
    dplyr::mutate(expression_subtype_you = factor(expression_subtype_you, levels = c("NMF1", "NMF2", "NMF3"), labels = c("PCS1", "PCS2", "PCS3"))) %>%
    group_by(trajectory, expression_subtype_you) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(trajectory) %>%
    mutate(proportion = count / sum(count)) %>%
    arrange(trajectory, desc(expression_subtype_you)) %>% 
    mutate(
        label_y_position = cumsum(proportion) - 0.5 * proportion,
        label_text = scales::percent(proportion, accuracy = 1)
    )

expression_subtype_plot <- ggplot(expression_donut_plot_data, aes(x = 2, y = proportion, fill = expression_subtype_you)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    facet_wrap(~ trajectory) +
    scale_fill_brewer(palette = "Set2") +
    xlim(0.5, 2.5) +
    theme_void(base_size = 9, base_family = "sans") +
    theme(
        legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold", margin = margin(b = 10)),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(fill = "You expression subtype") +
    guides(fill = guide_legend(
        nrow = 2,
        byrow = TRUE,
        title.position = "top",
        title.hjust = 0.5
    )) +
    geom_text(
        data = expression_donut_plot_data %>%   
            group_by(trajectory) %>% 
            summarise(sample_count = sum(count)), 
        aes(x = 0.5, y = 0, label = paste0("N=", sample_count)), 
        inherit.aes = FALSE, 
        size = 3, 
        fontface = "italic", 
        color = "grey30"
    )


# Save expression subtype plot
write_tsv(
  expression_donut_plot_data %>% dplyr::select(trajectory, expression_subtype_you, count, proportion),
  file.path(OUTPUT_DIR, "SF11b_source_data.tsv")
)
ggsave(
    file.path(FIGURE_DIR, "SF11b_evotype_overlap_you_subtypes.pdf"),
    expression_subtype_plot,
    width = 100 / 25.3,
    height = 100 / 25.3
)


# =============================================================================
# WOODCOCK SUBTYPES
# =============================================================================

woodcock_subtypes = read_delim(woodcock_filepath)
woodcock_subtypes = merge(woodcock_subtypes %>% dplyr::rename(woodcock_trajectory = Ordering), trajectory_samples %>% dplyr::rename(PPCG_Sample_ID = sample, ppcg_trajectory = trajectory))

# Sankey diagram
names(trajectory_colors) = c("Canonical", "Alternative", "Gain-enriched")
df_sankey = woodcock_subtypes %>% 
    group_by(woodcock_trajectory, ppcg_trajectory) %>% 
    summarise(n = n(), .groups = "drop")

df_sankey$ppcg_trajectory = factor(df_sankey$ppcg_trajectory, levels = c("Canonical", "Alternative", "Gain-enriched"))
p <- ggplot(df_sankey, aes(y = n, axis1 = woodcock_trajectory, axis2 = ppcg_trajectory)) +
  # Draw flows (alluvia)
  geom_alluvium(aes(fill = ppcg_trajectory), width = 0.15, alpha = 0.2, knot.pos = 0.4) +
  # Draw vertical bars (strata)
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.05, alpha = 0.8, color = "white", linewidth = 0.2) +
  # Apply biological colors to flows
  scale_fill_manual(values = trajectory_colors) +
  # Styling
  theme_void() +
  theme(legend.position = "none")

write_tsv(df_sankey, file.path(OUTPUT_DIR, "SF7b_source_data.tsv"))
ggsave(file.path(FIGURE_DIR, "SF7b_sankey_woodcock_ppcg.pdf"), p, height = 50/25.4, width = 50/25.4, dpi = 300)

table(woodcock_subtypes$woodcock_trajectory, woodcock_subtypes$ppcg_trajectory)

# 74/78 canonical samples in woodcock map to our canonical
# 38/49 alternative samples in woodcock map to our alternative

# =============================================================================
# COUNTRY OF ORIGIN
# =============================================================================
trajectory_samples$country = get_country(trajectory_samples$sample)

# Bar plot of trajectory count per country
country_counts = trajectory_samples %>% group_by(country, trajectory) %>% summarise(n = n(), .groups = "drop")
p = ggplot(country_counts, aes(x = country, y = n, fill = trajectory)) + 
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Country of Origin", y = "Number of Samples", fill = "Trajectory") +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = trajectory_colors)

write_tsv(country_counts, file.path(OUTPUT_DIR, "SF7c_source_data.tsv"))
ggsave(file.path(FIGURE_DIR, "SF7c_trajectory_counts_by_country.pdf"), p, width = 70/25.3, height = 50/25.3, dpi = 300)
