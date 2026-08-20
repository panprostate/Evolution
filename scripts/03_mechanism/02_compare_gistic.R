# Evaluate GISTIC results for the three different trajectories
# To identify potential targets of positive selection for recurrent SCNAs

rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(lemon)
library(data.table)
library(GenomicRanges)
library(biomaRt)

# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

# PATHS
figdir <- "figures/03_mechanism/compare_gistic/"
outdir <- "outputs/03_mechanism/compare_gistic/"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
events_fps <- list.files("outputs/02_trajectories/", pattern = ".*_PLPlot_data.rds", full = TRUE, recursive = TRUE)
trajectories_fps <- list.files("outputs/02_trajectories/", pattern = "PPCG_Feb2026_.*mergedseg_with_clonality.txt", full.names = TRUE, recursive = TRUE)
gistic_outdir <- "scripts/03_mechanism/gistic/gistic_output"
gistic_scores_fp <- list.files(gistic_outdir, pattern = "scores.gistic", full.names = TRUE, recursive = TRUE)
chr_sizes_fp <- "data/ext/cytoBand_hg19.txt"
driver_gr_fp <- "data/ext/driver_gene_gr_hg19.Rdata"

# LOAD DATA
# Save chromosome coordinates in a data frame for plotting purposes
chrs <- data.table::fread(chr_sizes_fp) %>% 
    dplyr::select(V1, V2, V3, V4) %>% 
    dplyr::rename(chr = V1, start = V2, end = V3, cytoband = V4) %>%
    dplyr::mutate(chr_arm = paste0(chr, ifelse(grepl("p", cytoband), "p", "q"))) %>% 
    dplyr::mutate(chr_arm = gsub("chr", "", chr_arm)) 

chr_arms <- chrs %>% group_by(chr, chr_arm) %>% summarise(start = min(start), end = max(end)) %>% ungroup()

# Save GISTIC scores in a data frame, with an additional column for the trajectory
gistic_scores <- lapply(1:length(gistic_scores_fp), function(i) {
    scores = read.table(gistic_scores_fp[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    scores$trajectory <- basename(dirname(gistic_scores_fp[i]))
    return(scores)
}) %>% data.table::rbindlist()

gistic_scores <- gistic_scores %>% dplyr::filter(trajectory != "all")

# Load driver gene coordinates as a GRanges object for annotation purposes
load(driver_gr_fp)

# Gather all recurrent SCNAs across trajectories, and annotate type of SCNA and trajectory
trajectory_events <- lapply(trajectories_fps, function(fp){
    events = read.table(fp) %>% dplyr::filter(CNA != "dMut")
    events$trajectory = paste0("ordering_", str_extract(basename(fp), "DN\\d+")) %>% str_remove("DN")
    events = events %>% dplyr::mutate(ID = paste0(CNA, "_", ID))
    events = events %>% dplyr::filter(!duplicated(ID))    
    events$formatted_id = as.factor(events$ID)
    levels(events$formatted_id) = sapply(levels(events$formatted_id), function(id) {
        text_to_add = ""
        #if (substr(id, 1, 3) == "chr") {
        if (length(grep("_chr", id)) > 0) {
            locus_split = strsplit(id, "_")[[1]]
            #locus_gr = GRanges(seqnames=gsub(".*chr", "", locus_split[1]), ranges=IRanges(start=as.numeric(locus_split[2]), end=as.numeric(locus_split[3])))
            locus_gr = GRanges(seqnames=gsub(".*chr", "", locus_split[2]), ranges=IRanges(start=as.numeric(locus_split[3]), end=as.numeric(locus_split[4])))
            ol = suppressWarnings(findOverlaps(locus_gr, driver_gr))
            if (length(ol) > 0) {
            text_to_add = paste0(" (", paste(driver_gr$symbol[subjectHits(ol)], collapse=", "), ")")
            }
        }
        paste0(id, text_to_add)
    })
    return(events)
}) %>% data.table::rbindlist()

# Trajectory colors loaded from plot_theme.R
# trajectory_colors <- c("grey30", "#1B9E77", "#D95F02", "#7570B3")
names(trajectory_colors) <- c("Canonical", "Alternative", "Gain-enriched")

lookup <- c("all" = "All", "ordering_1" = "Canonical", "ordering_2" = "Alternative", "ordering_3" = "Gain-enriched")
gistic_scores$trajectory <- lookup[gistic_scores$trajectory]
trajectory_events$trajectory <- lookup[trajectory_events$trajectory]

write_delim(gistic_scores, file.path(outdir, "ST_GISTIC_by_trajectory.tsv"), delim = "\t")

# ANALYSIS -----------------------------------------------------------
# Plot genome-wide distribution of GISTIC scores for each trajectory 
p = gistic_scores %>% 
    dplyr::filter(Type == "Amp") %>% 
    ggplot(aes(x = Start, y = G.score, color = trajectory)) +
    geom_line(alpha = 0.5) + 
    scale_color_manual(values = trajectory_colors) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +
    facet_wrap(~Chromosome, scales = "free_x", nrow = 1)

p = axes2lemon(p)

write_tsv(
  gistic_scores %>% dplyr::filter(Type == "Amp") %>% dplyr::select(Chromosome, Start, G.score, trajectory),
  file.path(outdir, "SF8_amps_source_data.tsv")
)
ggsave(file.path(figdir, "SF8_gistic_amps.pdf"), width = 100/25.3, height = 60/25.3, dpi = 300)


p = gistic_scores %>% 
    dplyr::filter(Type == "Del") %>% 
    ggplot(aes(x = Start, y = G.score, color = trajectory)) +
    geom_line(alpha = 0.5) + 
    scale_color_manual(values = trajectory_colors) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +
    facet_wrap(~Chromosome, scales = "free_x", nrow = 1)

p = axes2lemon(p)

write_tsv(
  gistic_scores %>% dplyr::filter(Type == "Del") %>% dplyr::select(Chromosome, Start, G.score, trajectory),
  file.path(outdir, "SF8_dels_source_data.tsv")
)
ggsave(file.path(figdir, "SF8_gistic_dels.pdf"), width = 100/25.3, height = 60/25.3, dpi = 300)

# PEAKS IN RECURRENT EVENTS ACROSS TRAJECTORIES
# event_set = trajectory_events %>% dplyr::filter(!duplicated(formatted_id))

# lp = lapply(1:nrow(event_set), function(i){
#     label = event_set$formatted_id[i]
#     cna = event_set$CNA[i]
#     chr = event_set$chr[i]
#     start = event_set$start[i]
#     end = event_set$end[i]

#     if (cna == "dGain") {
#         gistic_subset = gistic_scores %>% dplyr::filter(Type == "Amp", Chromosome == chr, Start <= end, End >= start)
#     } else {
#         gistic_subset = gistic_scores %>% dplyr::filter(Type == "Del", Chromosome == chr, Start <= end, End >= start)
#     }
#     p = gistic_subset %>% 
#         ggplot(aes(x = Start, y = G.score, color = trajectory)) +
#         geom_line(alpha = 0.5) + 
#         scale_color_manual(values = trajectory_colors) +
#         labs(title = label)
#     p = axes2lemon(p)
#     ggsave(file.path(figdir, paste0(label, "_gistic.pdf")), width = 100/25.3, height = 60/25.3, dpi = 300)
# })


potential_targets <- data.frame(
    genes = c("MCM7","MDM4", "LRP1B", "FOXP1", "CHD1", "ZNF292", "ETV1", "NKX3-1", "MYC", "NCOA2", "CDKN2A", "PTEN", "ZBTB16", "CCND1", "CDKN1B", "BRCA2", "RB1", "FOXA1", "ZFHX3", "TP53", "TMPRSS2", "ERG", "BRCA1", "HIXIM1", "ETV6", "ZNRF3", "NCOR2", "NCOR1"),
    arms = c("7q", "1q", "2q", "3p", "5q", "6q", "7p", "8p", "8q", "8q", "9p", "10q", "11q", "11q", "12p", "13q", "13q", "14q", "16q", "17p", "21q", "21q", "17q", "17q", "12p", "22q", "12q", "17q"),
    gene_type = c("OG", "OG", "TSG", "TSG", "TSG", "TSG", "TSG", "TSG", "OG", "OG", "TSG", "TSG", "TSG", "OG", "TSG", "TSG", "TSG", "OG", "TSG", "TSG", "TSG", "TSG", "TSG", "TSG", "TSG", "TSG", "TSG", "TSG")
)

# Uncomment to rerun
# mart <- biomaRt::useEnsembl(biomart = "ensembl", 
#                     dataset = "hsapiens_gene_ensembl", 
#                     GRCh = 37)

# gene_coords <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters = "hgnc_symbol", values = potential_targets$genes, mart = mart) %>% 
#     dplyr::select(hgnc_symbol, chromosome_name, start_position, end_position) %>% 
#     dplyr::mutate(coordinate = (start_position + end_position)/2)

# gene_coords$gene_type = potential_targets$gene_type[match(gene_coords$hgnc_symbol, potential_targets$genes)]

# # Save into cache
# saveRDS(gene_coords, file = file.path(outdir, ".potential_target_gene_coords.rds"))
gene_coords = readRDS(file.path(outdir, ".potential_target_gene_coords.rds"))


# PEAKS ACROSS CHROMOSOME ARMS
# Deletions
lp_deletions = lapply(1:nrow(chr_arms), function(i){
    arm = chr_arms$chr_arm[i]
    chr = gsub("chr", "", chr_arms$chr[i])
    start = chr_arms$start[i]
    end = chr_arms$end[i]

    genes = gene_coords %>% dplyr::filter(chromosome_name == chr, coordinate >= start, coordinate <= end) %>% dplyr::mutate(coordinate = coordinate/1e6)

    # if unclear target, check for alternative drivers in the region
    # if (nrow(genes) == 0) {
    #     genes = as.data.frame(driver_gr) %>% dplyr::mutate(coordinate = (start + end)/2) %>% dplyr::rename(chromosome_name = seqnames, hgnc_symbol = symbol) %>% dplyr::select(coordinate, chromosome_name, hgnc_symbol)
    #     genes = genes %>% dplyr::filter(chromosome_name == chr, coordinate <= end, coordinate >= start)
    # }

    # Find putative drivers in the chromosome arm, and annotate them in the plot
    putative_drivers = as.data.frame(driver_gr) %>% dplyr::mutate(coordinate = (start + end)/2) %>% dplyr::rename(chromosome_name = seqnames, hgnc_symbol = symbol) %>% dplyr::select(coordinate, chromosome_name, hgnc_symbol) 
    putative_drivers = putative_drivers %>% dplyr::filter(chromosome_name == chr, coordinate <= end, coordinate >= start) %>% dplyr::filter(!hgnc_symbol %in% genes$hgnc_symbol) %>% dplyr::mutate(coordinate = coordinate/1e6)

    # Plot deletions
    gistic_subset = gistic_scores %>% dplyr::filter(Type == "Del", Chromosome == chr, Start <= end, End >= start) %>% dplyr::mutate(coordinate = (Start + End)/(2*1e6))
    p = gistic_subset %>% 
        ggplot(aes(x = coordinate, y = G.score)) +
        geom_line(alpha = 0.5, aes(color = trajectory)) + 
        scale_color_manual(values = trajectory_colors) +
        labs(title = paste0("Deletions in ", arm), x = "Coordinate (Mb)", y = "GISTIC score") + 
        theme(legend.position = "none")
    
    p = p + 
        geom_vline(data = genes %>% dplyr::filter(gene_type == "TSG"), aes(xintercept = coordinate), linetype = "dashed", color = "#005A8E") +
        geom_text(data = genes %>% dplyr::filter(gene_type == "TSG"), aes(x = coordinate, y = max(gistic_subset$G.score)*0.8, label = hgnc_symbol), angle = 90, vjust = -0.5, size = 2, color = "#005A8E") + 
        geom_vline(data = putative_drivers, aes(xintercept = coordinate), linetype = "dotted", color = "grey30") + 
        geom_text(data = putative_drivers, aes(x = coordinate, y = max(gistic_subset$G.score, na.rm = TRUE)*0.8, label = hgnc_symbol), angle = 90, vjust = -0.5, size = 1, color = "grey30")
    
    p = axes2lemon(p)
    ggsave(file.path(figdir, paste0("Del_", arm, "_gistic.pdf")), width = 100/25.3, height = 60/25.3, dpi = 300)
    return(p)
})

names(lp_deletions) <- chr_arms$chr_arm

# Amplifications
lp_gains = lapply(1:nrow(chr_arms), function(i){
    arm = chr_arms$chr_arm[i]
    chr = gsub("chr", "", chr_arms$chr[i])
    start = chr_arms$start[i]
    end = chr_arms$end[i]

    genes = gene_coords %>% dplyr::filter(chromosome_name == chr, coordinate >= start, coordinate <= end) %>% dplyr::mutate(coordinate = coordinate/1e6)

    # if unclear target, check for alternative drivers in the region
    # if (nrow(genes) == 0) {
    #     genes = as.data.frame(driver_gr) %>% dplyr::mutate(coordinate = (start + end)/2) %>% dplyr::rename(chromosome_name = seqnames, hgnc_symbol = symbol) %>% dplyr::select(coordinate, chromosome_name, hgnc_symbol)
    #     genes = genes %>% dplyr::filter(chromosome_name == chr, coordinate <= end, coordinate >= start)
    # }

    # Find putative drivers in the chromosome arm, and annotate them in the plot
    putative_drivers = as.data.frame(driver_gr) %>% dplyr::mutate(coordinate = (start + end)/2) %>% dplyr::rename(chromosome_name = seqnames, hgnc_symbol = symbol) %>% dplyr::select(coordinate, chromosome_name, hgnc_symbol) 
    putative_drivers = putative_drivers %>% dplyr::filter(chromosome_name == chr, coordinate <= end, coordinate >= start) %>% dplyr::filter(!hgnc_symbol %in% genes$hgnc_symbol) %>% dplyr::mutate(coordinate = coordinate/1e6)

    # Plot amplifications
    gistic_subset = gistic_scores %>% dplyr::filter(Type == "Amp", Chromosome == chr, Start <= chr_arms$end[i], End >= chr_arms$start[i]) %>% dplyr::mutate(coordinate = (Start + End)/(2*1e6))
    p = gistic_subset %>% 
        ggplot(aes(x = coordinate, y = G.score)) +
        geom_line(alpha = 0.5, aes(color = trajectory)) + 
        scale_color_manual(values = trajectory_colors) +
        labs(title = paste0("Amplifications in ", arm), x = "Coordinate (Mb)", y = "GISTIC score") + 
        theme(legend.position = "none")

    p = p + 
        geom_vline(data = genes %>% dplyr::filter(gene_type == "OG"), aes(xintercept = coordinate), linetype = "dashed", color = "#991F36") +
        geom_text(data = genes %>% dplyr::filter(gene_type == "OG"), aes(x = coordinate, y = max(gistic_subset$G.score)*0.9, label = hgnc_symbol), angle = 90, vjust = -0.5, size = 2, color = "#991F36") + 
        geom_vline(data = putative_drivers, aes(xintercept = coordinate), linetype = "dotted", color = "grey30") + 
        geom_text(data = putative_drivers, aes(x = coordinate, y = max(gistic_subset$G.score, na.rm = TRUE)*0.8, label = hgnc_symbol), angle = 90, vjust = -0.5, size = 1, color = "grey30")
    
    p = axes2lemon(p)

    ggsave(file.path(figdir, paste0("Amp_", arm, "_gistic.pdf")), width = 100/25.3, height = 60/25.3, dpi = 300)
    return(p)
})

names(lp_gains) <- chr_arms$chr_arm

relevant_loss_arms = c("1p", "1q", "2q", "3p", "4q", "5q", "6q", "8p", "10q", "11q", "12p", "12q", "13q", "16q", "17p", "17q", "22q")
p = cowplot::plot_grid(plotlist = lp_deletions[relevant_loss_arms], ncol = 3)
write_tsv(
  gistic_scores %>% dplyr::filter(Type == "Del") %>% dplyr::select(Chromosome, Start, End, G.score, trajectory),
  file.path(outdir, "SF9_source_data.tsv")
)
ggsave(file.path(figdir, "SF9_Del_chr_arms_gistic.pdf"), width = 180/25.3, height = 250/25.3, dpi = 300)

relevant_gain_arms = c("1q", "2p", "3q", "5p", "7p", "7q", "8p", "8q", "9p", "9q", "10q", "11q", "12q", "14q", "16p", "16q", "17q")
p = cowplot::plot_grid(plotlist = lp_gains[relevant_gain_arms], ncol = 3)
write_tsv(
  gistic_scores %>% dplyr::filter(Type == "Amp") %>% dplyr::select(Chromosome, Start, End, G.score, trajectory),
  file.path(outdir, "SF10_source_data.tsv")
)
ggsave(file.path(figdir, "SF10_Amp_chr_arms_gistic.pdf"), width = 180/25.3, height = 250/25.3, dpi = 300)
