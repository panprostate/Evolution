#' Quality Control Analysis of DPClust Clusters for PPCG Evolution Study
#' 
#' This script performs comprehensive quality control analysis of mutation clusters
#' derived from DPClust in the PPCG Evolution study. The analysis includes:
#' 
#' 1. Comparison between DPClust and PhylogicNDT clustering results: CCF correlation and ARI calcuation
#' 2. Assessment of cancer cell fraction (CCF) correlations between methods
#' 3. Analysis of variant allele frequency (VAF) distributions by copy number status
#' 4. Quality control metrics for sequencing depth and purity adjustments
#' 
#' Input files required:
#'   - Cellularity/Ploidy files from Battenberg
#'   - DPClust cluster assignments and cluster definitions
#'   - Cancer Cell Fraction (CCF) input files
#'   - PhylogicNDT clustering results (optional)
#' 
#' Outputs:
#'   - Individual sample QC plots and summaries
#'   - Cross-method comparison statistics
#'   - Aggregate quality control visualizations


# Clear workspace for clean analysis environment
rm(list = ls(all = TRUE))

# =============================================================================
# PACKAGE DEPENDENCIES
# =============================================================================

# Load required packages for data manipulation, visualization, and analysis
library(dplyr)           # Data manipulation and transformation
library(VariantAnnotation) # Genomic variant annotation
library(readr)           # Fast data reading/writing
library(cowplot)         # Publication-quality plot arrangements
library(data.table)      # High-performance data operations
library(doMC)            # Parallel computing backend
library(foreach)         # Parallel loops
library(ggplot2)         # Grammar of graphics visualization
library(lemon)           # Extended ggplot2 functionality
library(mclust)          # Model-based clustering (for ARI calculation)
library(patchwork)       # Plot composition

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================

# Configure parallel processing
PARALLEL_CORES <- 20
registerDoMC(PARALLEL_CORES)

# Define output directories
OUTPUT_DIR <- "outputs/00_preprocessing/cluster_qc/"
FIGURE_DIR <- "figures/00_preprocessing/cluster_qc/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# Load custom utility functions and plotting themes
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Parse DPClust Output and Merge with PhylogicNDT CCF Data
#' 
#' Combines DPClust cluster assignments, cluster definitions, and PhylogicNDT CCF data into
#' a comprehensive dataset with mutation annotations, copy number information,
#' and clonality classifications.
#' 
#' @param dpclust_assignments_filepath Path to DPClust assignment file (.bed)
#' @param dpclust_clusters_filepath Path to DPClust cluster definitions (.txt)
#' @param ccf_filepath Path to CCF input file (.txt)
#' @param clonal_threshold Numeric. CCF threshold for clonal classification (default: 0.9)
#' @return data.frame with merged DPClust and CCF information
parse_dpclust_data <- function(dpclust_assignments_filepath, dpclust_clusters_filepath, 
                              ccf_filepath, clonal_threshold = 0.9) {
    
    # Load and prepare input datasets
    ccf_data <- read_delim(ccf_filepath) %>% 
        dplyr::mutate(mutation_id = paste0("chr", chr, ":", end))
    
    dpclust_assignments <- read_delim(dpclust_assignments_filepath) %>% 
        dplyr::mutate(mutation_id = paste0("chr", chr, ":", end))
    
    dpclust_clusters <- read_delim(dpclust_clusters_filepath)

    # Merge CCF data with cluster assignments
    merged_dpclust_data <- merge(ccf_data, 
                                dpclust_assignments[, c("cluster", "mutation_id")], 
                                by = "mutation_id")
    
    # Create cluster CCF lookup table
    cluster_ccf_lookup <- dpclust_clusters$location
    names(cluster_ccf_lookup) <- as.character(dpclust_clusters$cluster.no)

    # Add cluster CCF and determine clonality status
    merged_dpclust_data$cluster_ccf <- cluster_ccf_lookup[as.character(merged_dpclust_data$cluster)]
    merged_dpclust_data$clonal <- ifelse(merged_dpclust_data$cluster_ccf >= clonal_threshold, TRUE, FALSE)

    # Calculate average copy number metrics for subclonal copy number regions
    merged_dpclust_data$avg_total_cn <- ifelse(
        !is.na(merged_dpclust_data$frac2), 
        merged_dpclust_data$frac1 * (merged_dpclust_data$nMaj1 + merged_dpclust_data$nMin1) + 
        merged_dpclust_data$frac2 * (merged_dpclust_data$nMaj2 + merged_dpclust_data$nMin2), 
        merged_dpclust_data$nMaj1 + merged_dpclust_data$nMin1
    )

    merged_dpclust_data$avg_minor_cn <- ifelse(
        !is.na(merged_dpclust_data$frac2), 
        merged_dpclust_data$frac1 * (merged_dpclust_data$nMin1) + 
        merged_dpclust_data$frac2 * (merged_dpclust_data$nMin2), 
        merged_dpclust_data$nMin1
    )

    merged_dpclust_data$avg_major_cn <- ifelse(
        !is.na(merged_dpclust_data$frac2), 
        merged_dpclust_data$frac1 * (merged_dpclust_data$nMaj1) + 
        merged_dpclust_data$frac2 * (merged_dpclust_data$nMaj2), 
        merged_dpclust_data$nMaj1
    )

    # Calculate sequencing depth and VAF metrics
    merged_dpclust_data$sequencing_depth <- merged_dpclust_data$mut.count + merged_dpclust_data$WT.count
    merged_dpclust_data$VAF <- merged_dpclust_data$mut.count / merged_dpclust_data$sequencing_depth
    
    # Normalize depth by sample median for quality control
    merged_dpclust_data$normalized_depth <- merged_dpclust_data$sequencing_depth / 
        median(merged_dpclust_data$sequencing_depth, na.rm = TRUE)
    
    return(merged_dpclust_data)
}

#' Adjust VAF by Tumor Purity and Copy Number
#' 
#' Adjusts variant allele frequencies to account for tumor purity and copy number
#' alterations. The adjusted VAF represents the fraction of tumor cells carrying
#' the mutation.
#' 
#' @param vafs Numeric vector of variant allele frequencies
#' @param purity Numeric. Tumor purity (cellularity) value
#' @param copy_number_values Numeric vector of total copy number values
#' @return Numeric vector of purity-adjusted VAFs
adjust_vaf_for_purity <- function(vafs, purity, copy_number_values) {
    # Internal function to adjust individual VAF values
    adjust_single_vaf <- function(vaf, purity, avg_total_cn) {
        return(vaf / (purity * avg_total_cn / (purity * avg_total_cn + (1 - purity) * 2)))
    }
    
    # Apply adjustment to all VAF values
    return(sapply(1:length(vafs), function(i) {
        adjust_single_vaf(vafs[i], purity, copy_number_values[i])
    }))
}

#' Process Complete Sample Dataset
#' 
#' Integrates all data sources for a single sample including DPClust results,
#' PhylogicNDT clustering (if available), purity/ploidy information, and
#' calculates adjusted VAF values.
#' 
#' @param sample_id Character. Sample identifier
#' @param purity_filepaths Character vector. Paths to purity files
#' @param ccf_filepaths Character vector. Paths to CCF files
#' @param dpclust_clusters_filepaths Character vector. Paths to DPClust cluster files
#' @param dpclust_assignments_filepaths Character vector. Paths to DPClust assignment files
#' @param phylogic_filepaths Character vector. Paths to PhylogicNDT files
#' @return data.frame Complete sample dataset with all integrated information
process_sample_data <- function(sample_id, purity_filepaths, ccf_filepaths, 
                               dpclust_clusters_filepaths, dpclust_assignments_filepaths, 
                               phylogic_filepaths) {
    
    # Extract filepaths for current sample
    current_purity_filepath <- purity_filepaths[grep(sample_id, purity_filepaths)]
    current_ccf_filepath <- ccf_filepaths[grep(sample_id, ccf_filepaths)]
    current_dpclust_clusters_filepath <- dpclust_clusters_filepaths[grep(sample_id, dpclust_clusters_filepaths)]
    current_dpclust_assignments_filepath <- dpclust_assignments_filepaths[grep(sample_id, dpclust_assignments_filepaths)]

    # Validate that all required files are present and unique
    if (any(c(
            length(current_purity_filepath), 
            length(current_ccf_filepath), 
            length(current_dpclust_clusters_filepath), 
            length(current_dpclust_assignments_filepath)
        ) == 0)) {
            return(data.frame())
    }

    # Parse DPClust information: CCF, VAF, CN status, cluster assignment, clonality
    sample_dpclust_data <- parse_dpclust_data(
        dpclust_assignments_filepath = current_dpclust_assignments_filepath,
        dpclust_clusters_filepath = current_dpclust_clusters_filepath,
        ccf_filepath = current_ccf_filepath
    )

    # Load PhylogicNDT information if available
    current_phylogic_filepath <- phylogic_filepaths[grep(sample_id, phylogic_filepaths)]
    if (length(current_phylogic_filepath) == 1) {
        phylogic_data <- read_delim(current_phylogic_filepath) %>%
            dplyr::mutate(mutation_id = paste0(Chromosome, ":", Start_position)) %>%
            dplyr::select(
                mutation_id, 
                phylogic_ccf = preDP_ccf_mean, 
                phylogic_cluster = Cluster_Assignment, 
                phylogic_mean_cluster = clust_ccf_mean
            )
        sample_dpclust_data <- merge(sample_dpclust_data, phylogic_data, by = "mutation_id", all.x = TRUE)
    } else {
        # Add empty columns if PhylogicNDT data not available
        sample_dpclust_data$phylogic_ccf <- NA
        sample_dpclust_data$phylogic_cluster <- NA
        sample_dpclust_data$phylogic_mean_cluster <- NA
    }
    
    # Load purity and ploidy information
    sample_dpclust_data$purity <- read_purity(current_purity_filepath)
    sample_dpclust_data$ploidy <- read_ploidy(current_purity_filepath)
    
    # Adjust VAF by purity and copy number
    sample_dpclust_data$adjusted_vaf <- adjust_vaf_for_purity(
        sample_dpclust_data$VAF, 
        sample_dpclust_data$purity[1], 
        sample_dpclust_data$avg_total_cn
    )
    
    # Add sample ID for identification
    sample_dpclust_data$sample <- sample_id

    return(sample_dpclust_data)
}

#' Calculate Adjusted Rand Index (ARI)
#' 
#' Computes the Adjusted Rand Index between two clustering assignments to
#' quantify clustering agreement.
#' 
#' @param data_frame data.frame containing clustering assignments
#' @param clusters_column1 Character. Column name for first clustering method (default: "cluster")
#' @param clusters_column2 Character. Column name for second clustering method (default: "phylogic_cluster")
#' @return Numeric. ARI value between -1 and 1
calculate_ari <- function(data_frame, clusters_column1 = "cluster", clusters_column2 = "phylogic_cluster") {
    return(
        adjustedRandIndex(
            data_frame[[clusters_column1]], 
            data_frame[[clusters_column2]]
        )
    )
}

# =============================================================================
# FILE PATH CONFIGURATION
# =============================================================================

# Define input file paths for all required data sources
purity_filepaths <- list.files(
    "data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/purity_ploidy/",
    pattern = "*.txt", 
    full.names = TRUE
)

dpclust_clusters_filepaths <- list.files(
    "data/processed/PPCG_SNV_clustering", 
    pattern = "bestClusterInfo.txt$", 
    full.names = TRUE, 
    recursive = TRUE
)

dpclust_assignments_filepaths <- list.files(
    "data/processed/PPCG_SNV_clustering/", 
    pattern = "bestConsensusAssignments.bed$", 
    full.names = TRUE, 
    recursive = TRUE
)

ccf_filepaths <- list.files(
    "data/processed/PPCG_SNV_clustering/DPClust_inputs", 
    pattern = "txt$", 
    full.names = TRUE
)

phylogic_filepaths <- list.files(
    "data/processed/phylogicndt/phylogic_clustering", 
    pattern = "*mut_ccfs.txt$", 
    full.names = TRUE
)

# Load sample list
qc_pass_samples <- qc_pass()

# =============================================================================
# DATA PROCESSING AND INTEGRATION
# =============================================================================
# Process DPClust information for all samples into a combined dataframe
# Note: The commented code shows the original parallel processing approach
# dpclust_data = foreach(sample_id = qc_pass_samples, .combine = rbind) %dopar% {
#     process_sample_data(
#         sample_id = sample_id, 
#         purity_filepaths, 
#         ccf_filepaths, 
#         dpclust_clusters_filepaths, 
#         dpclust_assignments_filepaths, 
#         phylogic_filepaths
#     )
# }

# Rename subclonal fraction column to CCF for clarity
# dpclust_data = dpclust_data %>% dplyr::rename(ccf = cluster_ccf) 

# Save processed data for future use
# saveRDS(dpclust_data, file.path(OUTPUT_DIR, "info_dpclust_phylogic_all_samples.rds"))

# Load pre-processed data (using existing workflow)
dpclust_data = readRDS(file.path(OUTPUT_DIR, "info_dpclust_phylogic_all_samples.rds"))

dpclust_data$truncated_ccf <- ifelse(dpclust_data$ccf > 1, 1, dpclust_data$ccf)
dpclust_data$copy_number_status <- ifelse(
    round(dpclust_data$avg_total_cn) != dpclust_data$avg_total_cn, 
    "Subclonal CN", 
    "Clonal CN"
)

# =============================================================================
# DPCLUST VS PHYLOGICNDT COMPARISON ANALYSIS
# =============================================================================

# Compare DPClust and PhylogicNDT clustering results for each sample
individual_sample_figure_dir = file.path(FIGURE_DIR, "individual_samples")
dir.create(individual_sample_figure_dir, recursive = TRUE, showWarnings = FALSE)

dpclust_vs_phylogic_comparison = foreach(sample_id = qc_pass_samples) %dopar% {
    # Filter data for current sample with valid clustering assignments
    sample_data = dpclust_data %>% 
        dplyr::filter(sample == sample_id & !is.na(phylogic_ccf) & !is.na(cluster))
    
    # Group clonal clusters together for ARI calculation (arbitrary label 101)
    dpclust_clonal_clusters = sample_data$cluster[sample_data$clonal]
    if (length(dpclust_clonal_clusters) > 2) {
        sample_data$cluster = ifelse(
            sample_data$cluster %in% dpclust_clonal_clusters, 
            101, 
            sample_data$cluster
        )
    }
    
    # Calculate clustering agreement using Adjusted Rand Index
    clustering_ari = calculate_ari(sample_data)
    
    # Create histogram plots comparing CCF distributions between methods
    dpclust_ccf_plot = ggplot(
        sample_data[!is.na(sample_data$cluster),], 
        aes(x = ccf, fill = as.character(cluster))
    ) +
        geom_histogram(alpha = 0.5, position = "identity") +
        geom_vline(aes(xintercept = cluster_ccf), linetype = "dashed", color = "black") +
        labs(
            title = paste("DPClust:", sample_id, "ARI: ", round(clustering_ari, 2)), 
            x = "CCF", 
            y = "Count", 
            fill = "Cluster"
        ) + 
        scale_x_continuous(breaks = seq(0, 1.5, by = 0.25), limits = c(0, 1.5)) +
        theme(legend.position = "bottom")

    phylogic_ccf_plot = ggplot(
        sample_data[!is.na(sample_data$cluster),], 
        aes(x = phylogic_ccf, fill = as.character(phylogic_cluster))
    ) + 
        geom_histogram(alpha = 0.5, position = "identity") +
        geom_vline(aes(xintercept = phylogic_mean_cluster), linetype = "dashed", color = "black") +
        labs(
            title = paste("PhylogicNDT:", sample_id), 
            x = "CCF", 
            y = "Count", 
            fill = "Cluster"
        ) + 
        scale_x_continuous(breaks = seq(0, 1.5, by = 0.25), limits = c(0, 1.5)) +
        theme(legend.position = "bottom")
    
    # Combine and save histogram plots
    combined_histogram_plot = plot_grid(dpclust_ccf_plot, phylogic_ccf_plot, ncol = 1, align = "v")
    if (sample_id == "PPCG1060a_DNA") {
        write_tsv(
          sample_data %>% dplyr::select(sample, cluster, ccf, cluster_ccf,
                                        phylogic_ccf, phylogic_cluster, phylogic_mean_cluster),
          file.path(OUTPUT_DIR, "SF1d_source_data.tsv")
        )
    }
    ggsave(
        file.path(individual_sample_figure_dir, paste0(sample_id, "_dpclust_vs_phylogic_clusters.png")),
        combined_histogram_plot, 
        height = 90/25.3, 
        width = 90/25.3
    )
    ggsave(
        file.path(individual_sample_figure_dir, paste0(sample_id, "_dpclust_vs_phylogic_clusters.pdf")), 
        combined_histogram_plot, 
        height = 90/25.3, 
        width = 90/25.3
    )

    # Calculate CCF correlation between methods
    ccf_pearson_correlation = cor(sample_data$truncated_ccf, sample_data$phylogic_ccf, method = "pearson")

    # Create CCF correlation scatterplot
    ccf_correlation_plot = ggplot(sample_data, aes(x = truncated_ccf, y = phylogic_ccf)) + 
        geom_point(alpha = 0.2) +
        geom_smooth(method = "lm", color = "blue", se = FALSE) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgreen") +
        ggpubr::stat_cor(method = "pearson", color = "blue", size = 2) +
        labs(title = paste("Sample:", sample_id), x = "DPClust CCF", y = "PhylogicNDT CCF") + 
        xlim(c(0, 1)) + 
        ylim(c(0, 1))

    if (sample_id == "PPCG0004a_DNA") {
        write_tsv(
          sample_data %>% dplyr::select(sample, truncated_ccf, phylogic_ccf),
          file.path(OUTPUT_DIR, "SF1a_source_data.tsv")
        )
    }
    ggsave(
        file.path(individual_sample_figure_dir, paste0(sample_id, "_dpclust_vs_phylogic_ccf_scatter.png")),
        ccf_correlation_plot,
        height = 65/25.3, 
        width = 65/25.3
    )

    ggsave(
        file.path(individual_sample_figure_dir, paste0(sample_id, "_dpclust_vs_phylogic_ccf_scatter.pdf")), 
        ccf_correlation_plot, 
        height = 65/25.3, 
        width = 65/25.3
    )

    # Compare number of subclones between methods
    n_subclones_dpclust = length(unique(sample_data$cluster[!sample_data$clonal]))
    
    # Determine clonal clusters in PhylogicNDT data
    phylogic_clonal_clusters = sample_data$phylogic_cluster[sample_data$phylogic_mean_cluster >= 0.9]
    if (length(phylogic_clonal_clusters) == 0) {
        # If no clusters above threshold, take highest CCF cluster
        phylogic_clonal_clusters = sample_data %>% 
            arrange(desc(phylogic_mean_cluster)) %>% 
            head(n = 1) %>% 
            dplyr::pull(phylogic_cluster)
    }
    
    sample_data$phylogic_clonal = sample_data$phylogic_cluster %in% phylogic_clonal_clusters
    n_subclones_phylogic = length(unique(sample_data$phylogic_cluster[!sample_data$phylogic_clonal])) 

    # Compare number of subclonal mutations
    n_subclonal_mutations_dpclust = sum(!sample_data$clonal)
    n_subclonal_mutations_phylogic = sum(!sample_data$phylogic_clonal, na.rm = TRUE)

    # Compile summary statistics
    sample_summary_stats = data.frame(
        sample = sample_id,
        ari = clustering_ari,
        pearson = ccf_pearson_correlation,
        n_subclones_dpclust = n_subclones_dpclust,
        n_subclones_phylogic = n_subclones_phylogic,
        n_subclonal_mutations_dpclust = n_subclonal_mutations_dpclust,
        n_subclonal_mutations_phylogic = n_subclonal_mutations_phylogic
    )
    
    # Save individual sample summary
    write_delim(
        sample_summary_stats, 
        file.path(individual_sample_figure_dir, paste0(sample_id, "_dpclust_vs_phylogic_summary.tsv")), 
        delim = "\t"
    )
    
    return(sample_summary_stats)
} %>% data.table::rbindlist()

# =============================================================================
# AGGREGATE QUALITY CONTROL VISUALIZATIONS
# =============================================================================

median(dpclust_vs_phylogic_comparison$ari, na.rm = TRUE)
median(dpclust_vs_phylogic_comparison$pearson, na.rm = TRUE)

# Plot distribution of clustering ARI across samples
clustering_ari_histogram = ggplot(dpclust_vs_phylogic_comparison, aes(x = ari)) + 
    geom_histogram(binwidth = 0.1, fill = "lightblue", alpha = 0.7) + 
    geom_vline(xintercept = median(dpclust_vs_phylogic_comparison$ari, na.rm = TRUE), linetype = "dashed", color = "black") +
    labs(x = "Clustering ARI", y = "Number of samples")

write_tsv(dpclust_vs_phylogic_comparison %>% dplyr::select(ari), file.path(OUTPUT_DIR, "SF1e_source_data.tsv"))
ggsave(file.path(FIGURE_DIR, "SF1e_dpclust_vs_phylogic_ari_histogram.png"), clustering_ari_histogram, width = 65/25.3, height = 65/25.3)
ggsave(file.path(FIGURE_DIR, "SF1e_dpclust_vs_phylogic_ari_histogram.pdf"), clustering_ari_histogram, width = 65/25.3, height = 65/25.3)

# Plot distribution of correlation between CCFs across samples
ccf_correlation_histogram = ggplot(dpclust_vs_phylogic_comparison, aes(x = pearson)) + 
    geom_histogram(binwidth = 0.1, fill = "lightblue", alpha = 0.7) + 
    geom_vline(xintercept = median(dpclust_vs_phylogic_comparison$pearson, na.rm = TRUE), linetype = "dashed", color = "black") +
    labs(x = "CCF Pearson's correlation", y = "Number of samples") + 
    xlim(c(0, 1.05))

write_tsv(dpclust_vs_phylogic_comparison %>% dplyr::select(pearson), file.path(OUTPUT_DIR, "SF1b_source_data.tsv"))
ggsave(file.path(FIGURE_DIR, "SF1b_dpclust_vs_phylogic_ccf_pearson_histogram.png"), ccf_correlation_histogram, width = 65/25.3, height = 65/25.3)
ggsave(file.path(FIGURE_DIR, "SF1b_dpclust_vs_phylogic_ccf_pearson_histogram.pdf"), ccf_correlation_histogram, width = 65/25.3, height = 65/25.3)

# Plot number of subclones across samples
subclones_comparison_heatmap = dpclust_vs_phylogic_comparison %>%
  dplyr::count(n_subclones_dpclust, n_subclones_phylogic) %>% 
  ggplot(aes(x = n_subclones_dpclust, y = n_subclones_phylogic, fill = n)) +
  geom_tile(color = "white") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgreen") +
  geom_text(aes(label = n), color = "black", size = 2) +
  scale_fill_gradient(low = "white", high = "darkred") +
  scale_x_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
  theme(legend.position = "none") +
  labs(x = "DPClust subclones", y = "PhylogicNDT subclones", fill = "Samples")

write_tsv(
  dpclust_vs_phylogic_comparison %>% dplyr::count(n_subclones_dpclust, n_subclones_phylogic),
  file.path(OUTPUT_DIR, "SF1c_source_data.tsv")
)
ggsave(file.path(FIGURE_DIR, "SF1c_scatterplot_dpclust_vs_phylogic_nsubclones.png"), subclones_comparison_heatmap, width = 65/25.3, height = 65/25.3)
ggsave(file.path(FIGURE_DIR, "SF1c_scatterplot_dpclust_vs_phylogic_nsubclones.pdf"), subclones_comparison_heatmap, width = 65/25.3, height = 65/25.3)

# Plot number of subclonal mutations across samples
subclonal_mutations_comparison = ggplot(dpclust_vs_phylogic_comparison, aes(x = log10(n_subclonal_mutations_dpclust+1), y = log10(n_subclonal_mutations_phylogic+1))) + 
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgreen") +
    labs(x = "Subclonal mutations\n(DPClust)", y = "Subclonal mutations\n(PhylogicNDT)") + 
    ggpubr::stat_cor(method = "spearman", color = "blue", size = 2)

ggsave(file.path(FIGURE_DIR, "scatterplot_dpclust_vs_phylogic_nsubclonal_mutations.png"), subclonal_mutations_comparison, width = 65/25.3, height = 65/25.3)

# Create combined summary plot
phylogic_comparison_summary = plot_grid(
    clustering_ari_histogram, 
    ccf_correlation_histogram, 
    subclones_comparison_heatmap, 
    subclonal_mutations_comparison, 
    ncol = 2
)

ggsave(file.path(FIGURE_DIR, "dpclust_vs_phylogic_summary.png"), phylogic_comparison_summary, height = 120/25.3, width = 120/25.3)

# =============================================================================
# EVOLUTIONARY METRICS INTEGRATION AND ANALYSIS
# =============================================================================

# Load evolutionary metrics for further analysis
evolutionary_metrics = read_delim("outputs/00_preprocessing/evo_metrics.tsv") %>% 
    dplyr::mutate(sample = paste0(smp, "_DNA"))

# Merge with comparison results to analyze relationship between clustering quality and evolution
dpclust_vs_phylogic_with_evolution = merge(evolutionary_metrics, dpclust_vs_phylogic_comparison, by = "sample") 

# Display samples with lowest clustering agreement
arrange(dpclust_vs_phylogic_with_evolution, ari) %>% 
    dplyr::select(sample, ari, subclonal_pga)

# Analyze relationship between clustering quality and subclonal proportion of genome altered (PGA)
ari_vs_subclonal_pga_plot = ggplot(
    dpclust_vs_phylogic_with_evolution %>% 
        dplyr::filter(n_subclones_dpclust > 0 & n_subclones_phylogic > 0), 
    aes(x = ari, y = subclonal_pga)
) + 
    geom_point() + 
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    ggpubr::stat_cor(method = "spearman", color = "blue") +
    labs(x = "Clustering ARI", y = "Subclonal PGA")

ggsave(file.path(FIGURE_DIR, "dpclust_vs_phylogic_ari_vs_subclonal_pga.png"), ari_vs_subclonal_pga_plot, height = 65/25.3, width = 65/25.3)

# Statistical test for correlation between ARI and subclonal PGA
cor.test(dpclust_vs_phylogic_with_evolution$ari, dpclust_vs_phylogic_with_evolution$subclonal_pga, method = "spearman")

# Stratified analysis by subclonal PGA levels
ari_histogram_all_samples = ggplot(dpclust_vs_phylogic_with_evolution, aes(x = ari)) + 
    geom_histogram(alpha = 0.7, fill = "lightblue") +
    geom_vline(xintercept = median(dpclust_vs_phylogic_with_evolution$ari, na.rm = TRUE), linetype = "dashed", color = "palegreen") +
    labs(x = "Clustering ARI", y = "Number of samples") +
    theme(legend.position = "bottom")

high_subclonal_pga_threshold = quantile(dpclust_vs_phylogic_with_evolution$subclonal_pga, na.rm = TRUE, .66)
dpclust_vs_phylogic_with_evolution = dpclust_vs_phylogic_with_evolution %>% 
    dplyr::mutate(subclonal_pga_group = ifelse(subclonal_pga > high_subclonal_pga_threshold, "High subclonal PGA", "Low subclonal PGA")) 

ari_histogram = dpclust_vs_phylogic_with_evolution %>%
    ggplot(aes(x = ari, fill = subclonal_pga_group)) +
    geom_histogram(alpha = 0.4, position = "identity") +
    geom_vline(xintercept = dpclust_vs_phylogic_with_evolution %>% filter(subclonal_pga_group == "High subclonal PGA") %>% pull(ari) %>% median(na.rm = TRUE), linetype = "dashed", color = "salmon") +
    geom_vline(xintercept = dpclust_vs_phylogic_with_evolution %>% filter(subclonal_pga_group == "Low subclonal PGA") %>% pull(ari) %>% median(na.rm = TRUE), linetype = "dashed", color = "lightblue") +
    scale_fill_manual(values = c("High subclonal PGA" = "salmon", "Low subclonal PGA" = "lightblue")) +
    labs(x = "Clustering ARI", y = "Number of samples", fill = "Subclonal PGA group") +
    theme(legend.position = "top")

write_tsv(dpclust_vs_phylogic_with_evolution %>% dplyr::select(ari, subclonal_pga_group), file.path(OUTPUT_DIR, "SF1f_source_data.tsv"))
ggsave(file.path(FIGURE_DIR, "SF1f_dpclust_vs_phylogic_ari_histogram_by_subclonal_pga_group.png"), ari_histogram, height = 65/25.3, width = 65/25.3)
ggsave(file.path(FIGURE_DIR, "SF1f_dpclust_vs_phylogic_ari_histogram_by_subclonal_pga_group.pdf"), ari_histogram, height = 65/25.3, width = 65/25.3)

# High subclonal PGA samples (top third)
high_subclonal_pga_samples = dpclust_vs_phylogic_with_evolution[
    dpclust_vs_phylogic_with_evolution$subclonal_pga > high_subclonal_pga_threshold & 
    dpclust_vs_phylogic_with_evolution$n_subclones_dpclust > 0 & 
    dpclust_vs_phylogic_with_evolution$n_subclones_phylogic > 0,
]

ari_histogram_high_pga = ggplot(high_subclonal_pga_samples, aes(x = ari)) + 
    geom_histogram(alpha = 0.7, fill = "lightblue") +
    geom_vline(
        xintercept = median(high_subclonal_pga_samples$ari, na.rm = TRUE), 
        linetype = "dashed", color = "palegreen"
    ) +
    labs(x = "Clustering ARI\n(High subclonal PGA)", y = "Number of samples") +
    theme(legend.position = "bottom")

# Low subclonal PGA samples (bottom two-thirds)
low_subclonal_pga_samples = dpclust_vs_phylogic_with_evolution[
    dpclust_vs_phylogic_with_evolution$subclonal_pga <= high_subclonal_pga_threshold & 
    dpclust_vs_phylogic_with_evolution$n_subclones_dpclust > 0 & 
    dpclust_vs_phylogic_with_evolution$n_subclones_phylogic > 0,
]

ari_histogram_low_pga = ggplot(low_subclonal_pga_samples, aes(x = ari)) + 
    geom_histogram(alpha = 0.7, fill = "lightblue") +
    geom_vline(
        xintercept = median(low_subclonal_pga_samples$ari, na.rm = TRUE), 
        linetype = "dashed", color = "palegreen"
    ) +
    labs(x = "Clustering ARI\n(Low subclonal PGA)", y = "Number of samples") +
    theme(legend.position = "bottom")

# Correlation in regions with subclonal CNA
ccf_scatter_subclonal = ggplot(dpclust_data %>% dplyr::filter(copy_number_status == "Subclonal CN"), aes(x = truncated_ccf, y = phylogic_ccf)) + 
    geom_point(alpha = 0.01) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgreen") +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    ggpubr::stat_cor(method = "spearman", color = "blue", size = 2) +
    labs(x = "DPClust CCF", y = "PhylogicNDT CCF") + 
    ggtitle("Subclonal CN regions") +
    xlim(c(0, 1)) +
    ylim(c(0, 1))

ccf_scatter_clonal = ggplot(dpclust_data %>% dplyr::filter(copy_number_status == "Clonal CN"), aes(x = truncated_ccf, y = phylogic_ccf)) + 
    geom_point(alpha = 0.01) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgreen") +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    ggpubr::stat_cor(method = "spearman", color = "blue", size = 2) +
    ggtitle("Clonal CN regions") +
    labs(x = "DPClust CCF", y = "PhylogicNDT CCF") + 
    xlim(c(0, 1)) +
    ylim(c(0, 1))

write_tsv(
  dpclust_data %>% dplyr::filter(copy_number_status == "Clonal CN") %>%
    dplyr::select(truncated_ccf, phylogic_ccf, copy_number_status),
  file.path(OUTPUT_DIR, "SF1g_source_data.tsv")
)
ggsave(file.path(FIGURE_DIR, "SF1g_dpclust_vs_phylogic_ccf_scatter_clonal_cn.png"), ccf_scatter_clonal, height = 65/25.3, width = 65/25.3)
write_tsv(
  dpclust_data %>% dplyr::filter(copy_number_status == "Subclonal CN") %>%
    dplyr::select(truncated_ccf, phylogic_ccf, copy_number_status),
  file.path(OUTPUT_DIR, "SF1h_source_data.tsv")
)
ggsave(file.path(FIGURE_DIR, "SF1h_dpclust_vs_phylogic_ccf_scatter_subclonal_cn.png"), ccf_scatter_subclonal, height = 65/25.3, width = 65/25.3)

# Combine stratified histograms
ari_by_subclonal_pga_plot = cowplot::plot_grid(
    ari_histogram_all_samples, 
    ari_histogram_high_pga, 
    ari_histogram_low_pga, 
    ccf_scatter_clonal, 
    ccf_scatter_subclonal, 
    NULL, 
    ncol = 3,
    nrow = 2, 
    labels = c("a)", "b)", "c)", "d)", "e)", ""), 
    label_size = 9
)


ggsave(
    file.path(FIGURE_DIR, "dpclust_vs_phylogic_ari_histogram_by_subclonal_pga.png"), 
    ari_by_subclonal_pga_plot, 
    height = 120/25.3, 
    width = 180/25.3
)

# =============================================================================
# SEQUENCING DEPTH AND VAF QUALITY CONTROL ANALYSIS
# =============================================================================

# Analyse relationship between purity and mutation CCF
ccf_vs_purity_plot = ggplot(dpclust_data, aes(x = ccf, y = purity)) + 
    geom_point(alpha = 0.01) + 
    ggpubr::stat_cor(method = "spearman", color = "blue", size = 2) +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    labs(y = "Tumor Purity", x = "CCF") + 
    xlim(c(0,1)) + 
    ylim(c(0,1))

# Analyse relationship between purity and mutation VAF
vaf_vs_purity_plot = ggplot(dpclust_data, aes(x = VAF, y = purity)) + 
    geom_point(alpha = 0.01) + 
    ggpubr::stat_cor(method = "spearman", color = "blue", size = 2) +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    labs(y = "Tumor Purity", x = "VAF") + 
    xlim(c(0,1)) + 
    ylim(c(0,1))

# Analyze relationship between mutation CCF and sequencing depth
ccf_vs_depth_plot = ggplot(dpclust_data, aes(x = ccf, y = sequencing_depth)) + 
    geom_point(alpha = 0.01) + 
    ggpubr::stat_cor(method = "spearman", color = "blue", size = 2) +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    labs(x = "CCF", y = "Sequencing Depth") + 
    xlim(c(0,1)) + 
    ylim(c(0,300))

# Analyze relationship between mutation VAF and sequencing depth
vaf_vs_depth_plot = ggplot(dpclust_data, aes(x = VAF, y = sequencing_depth)) + 
    geom_point(alpha = 0.01) + 
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    ggpubr::stat_cor(method = "spearman", color = "blue", size = 2) + 
    labs(x = "VAF", y = "Sequencing Depth") + 
    ylim(c(0,300)) + 
    xlim(c(0,1)) 


# Combine depth quality control plots
qc_summary = plot_grid(
    vaf_vs_depth_plot, 
    vaf_vs_purity_plot, 
    ccf_vs_purity_plot, 
    ccf_vs_depth_plot, 
    ncol = 4, 
    labels = c("a)", "b)", "c)", "d)"), 
    label_size = 10
)

ggsave(file.path(FIGURE_DIR, "cluster_qc_summary.png"), qc_summary, height = 65/25.3, width = 160/25.3)

# =============================================================================
# COPY NUMBER AND VAF VALIDATION ANALYSIS
# =============================================================================

# Demonstrate that clustering has properly accounted for copy number status, purity, and ploidy
# by showing expected VAF patterns for different karyotypes

# Create expected VAF plots for different copy number states in clonal mutations
dpclust_data$karyotype = paste0(round(dpclust_data$avg_major_cn), ":", round(dpclust_data$avg_minor_cn))

# Clonal mutations in diploid regions (2:2) - expected peaks at 0.5 (heterozygous) and 0.25 (homozygous)
vaf_2_2_plot = dpclust_data %>% 
    dplyr::filter(karyotype == "2:2" & clonal) %>%
    dplyr::filter(copy_number_status == "Clonal CN") %>%
    ggplot(aes(x = adjusted_vaf)) + 
    geom_histogram(binwidth = 0.02, fill = "lightblue") + 
    labs(x = "Adjusted VAF", y = "Count") + 
    ggtitle("Copy Number: 2:2") +
    xlim(c(0,1.5)) + 
    # Expected peaks at 0.5 (heterozygous) and 0.25 (homozygous)
    geom_vline(xintercept = c(0.5, 0.25), linetype = "dashed", color = "green", size = 1)

# Clonal mutations in copy number loss (2:1) - expected peaks at 0.66 (on major allele) and 0.33 (on minor allele)  
vaf_2_1_plot = dpclust_data %>% 
    dplyr::filter(karyotype == "2:1" & clonal) %>%
    dplyr::filter(copy_number_status == "Clonal CN") %>%
    ggplot(aes(x = adjusted_vaf)) + 
    geom_histogram(binwidth = 0.02, fill = "lightblue") + 
    labs(x = "Adjusted VAF", y = "Count") + 
    ggtitle("Copy Number: 2:1") +
    xlim(c(0,1.5)) + 
    # Expected peaks at 0.66 (major allele) and 0.33 (minor allele)
    geom_vline(xintercept = c(0.66, 0.33), linetype = "dashed", color = "green", size = 1)

# Clonal mutations in copy number loss (1:1) - expected peak at 0.5 (heterozygous in diploid state)
vaf_1_1_plot = dpclust_data %>% 
    dplyr::filter(karyotype == "1:1" & clonal) %>%
    dplyr::filter(copy_number_status == "Clonal CN") %>%
    ggplot(aes(x = adjusted_vaf)) + 
    geom_histogram(binwidth = 0.02, fill = "lightblue") + 
    labs(x = "Adjusted VAF", y = "Count") + 
    ggtitle("Copy Number: 1:1") +
    xlim(c(0,1.5)) + 
    # Expected peak at 0.5
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "green", size = 1)

# Clonal mutations in copy-neutral LOH (2:0) - expected peaks at 1.0 (both copies) and 0.5 (one copy)
vaf_2_0_plot = dpclust_data %>% 
    dplyr::filter(karyotype == "2:0" & clonal) %>%
    dplyr::filter(copy_number_status == "Clonal CN") %>%
    ggplot(aes(x = adjusted_vaf)) + 
    geom_histogram(binwidth = 0.02, fill = "lightblue") + 
    labs(x = "Adjusted VAF", y = "Count") + 
    ggtitle("Copy Number: 2:0") +
    xlim(c(0,1.5)) + 
    # Expected peaks at 1.0 (both copies) and 0.5 (one copy)
    geom_vline(xintercept = c(1.0, 0.5), linetype = "dashed", color = "green", size = 1)

# Clonal mutations in hemizygous deletion (1:0) - expected peak at 1.0 (single remaining copy)
vaf_1_0_plot = dpclust_data %>% 
    dplyr::filter(karyotype == "1:0" & clonal) %>%
    dplyr::filter(copy_number_status == "Clonal CN") %>%
    ggplot(aes(x = adjusted_vaf)) + 
    geom_histogram(binwidth = 0.02, fill = "lightblue") + 
    labs(x = "Adjusted VAF", y = "Count") + 
    ggtitle("Copy Number: 1:0") +
    xlim(c(0,1.5)) + 
    # Expected peak at 1.0
    geom_vline(xintercept = 1.0, linetype = "dashed", color = "green", size = 1)

# Combine all copy number VAF validation plots
copy_number_vaf_validation_plot = plot_grid(
    vaf_2_2_plot, 
    vaf_2_1_plot, 
    vaf_1_1_plot, 
    vaf_2_0_plot, 
    vaf_1_0_plot,
    ncol = 2
)

ggsave(file.path(FIGURE_DIR, "cluster_qc_cn_vaf_fits.png"), copy_number_vaf_validation_plot, height = 100/25.3, width = 110/25.3)


