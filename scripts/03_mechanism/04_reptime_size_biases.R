# Analysis of size and replication time for LOH, HD and Gains across each evolutionary subtype
# Highlighting specific genes relevant for prostate cancer evolution

rm(list = ls(all = TRUE))

# LIBRARIES
library(tidyverse)
library(ggpubr)
library(rtracklayer)
library(GenomicRanges)
library(foreach)
library(doMC)
# library(ggrastr)
CORES = 20
registerDoMC(cores = CORES)

# FUNCTIONS
source("src/utils.R")
source("src/plot_theme.R")
source("src/plot_functions.R")

filter_alterations <- function(bb_scna, type_alteration, clonality){
    # Calculate ploidy 
    widths <- bb_scna$endpos - bb_scna$startpos
    bb_scna$avg_cn <- ifelse(
        is.na(bb_scna$frac2_A),
        bb_scna$nMin1_A + bb_scna$nMaj1_A,
        bb_scna$frac1_A * (bb_scna$nMin1_A + bb_scna$nMaj1_A) + bb_scna$frac2_A * (bb_scna$nMin2_A + bb_scna$nMaj2_A)
    )
    ploidy <- sum(widths * bb_scna$avg_cn) / sum(widths)

    # Estimate WGD status
    # Calculate proportion of genome with minor copy number 0
    bb_scna$is_loh = ifelse(is.na(bb_scna$frac2_A), bb_scna$nMin1_A == 0, bb_scna$nMin1_A == 0 | bb_scna$nMin2_A == 0) # only consider clonal LOH
    prop_loh = sum(widths[bb_scna$is_loh]) / sum(widths)
    is_wgd = 2.9 < (ploidy + prop_loh)

    # LOH alterations
    if (type_alteration == "LOH"){
        if (clonality == "clonal"){
            bb_scna = bb_scna %>% dplyr::filter(
                # Clonal LOH in either subpopulation
                (frac1_A > 0.9 & nMin1_A == 0 & nMaj1_A > 0) | 
                (frac1_A < 0.1 & nMin2_A == 0 & nMaj2_A > 0) | 
                # Clonal LOH triggered by subclonal HD
                (frac1_A < 0.9 & nMin1_A == 0 & nMin2_A == 0)
            )
        } else if (clonality == "subclonal"){
            bb_scna = bb_scna %>% dplyr::filter(
                # Subclonal LOH
                (frac1_A <= 0.9 & nMin1_A == 0 & nMin2_A > 0) | 
                (frac1_A <= 0.9 & nMin1_A > 0 & nMin2_A == 0)
            )
        } else {
            bb_scna = bb_scna %>% dplyr::filter(
                # LOH in any of the two alleles
                nMin1_A == 0 | nMin2_A == 0
            )
        }
    }


    # Homozygous deletions
    if (type_alteration == "HD"){
        if (clonality == "clonal"){
            bb_scna = bb_scna %>% dplyr::filter(
                # Clonal HD
                (frac1_A > 0.9 & nMin1_A == 0 & nMaj1_A == 0) | (frac1_A < 0.1 & nMin2_A == 0 & nMaj2_A == 0)
            )
        } else if (clonality == "subclonal"){
            bb_scna = bb_scna %>% dplyr::filter(
                # Subclonal HD
                (frac1_A < 0.9 & nMin1_A == 0 & nMaj1_A == 0 & nMaj2_A > 0) | (frac1_A < 0.9 & nMin2_A == 0 & nMaj2_A == 0 & nMaj1_A > 0)
            )
        } else {
            bb_scna = bb_scna %>% dplyr::filter(
                # HD in any of the two alleles
                (nMin1_A == 0 & nMaj1_A == 0) | (nMin2_A == 0 & nMaj2_A == 0)
            )
        }
    }

    # Copy-number gains
    gain_threshold = ifelse(is_wgd, 5, 3)

    if (type_alteration == "Gain"){
        if (clonality == "clonal"){
            bb_scna = bb_scna %>% dplyr::filter(
                # Clonal gain
                frac1_A > 0.9 & (nMin1_A + nMaj1_A) >= gain_threshold | 
                frac2_A > 0.9 & (nMin2_A + nMaj2_A) >= gain_threshold |
                # Subclonal gain where both subclones have gain
                (frac1_A < 0.9 & (nMin1_A + nMaj1_A) >= gain_threshold & (nMin2_A + nMaj2_A) >= gain_threshold)
            )
        } else if (clonality == "subclonal"){
            bb_scna = bb_scna %>% dplyr::filter(
                # Subclonal gain where either clone has gain
                frac1_A < 0.9 & 
                    (((nMin1_A + nMaj1_A) >= gain_threshold) & ((nMin2_A + nMaj2_A) < gain_threshold)) | 
                    (((nMin1_A + nMaj1_A) < gain_threshold) & ((nMin2_A + nMaj2_A) >= gain_threshold))
            )
        } else {
            bb_scna = bb_scna %>% dplyr::filter(
                # Gain in any of the two subclones
                (nMin1_A + nMaj1_A) >= gain_threshold | (nMin2_A + nMaj2_A) >= gain_threshold
            )
        }
    }

    # Any type of alteration
    normal_allele_cn = ifelse(is_wgd, 2, 1)
    if (type_alteration == "all"){
        if (clonality == "clonal"){
            bb_scna = bb_scna %>% dplyr::filter(
                (frac1_A > 0.9 | frac1_A < 0.1) & ((nMin1_A != normal_allele_cn) | (nMaj1_A != normal_allele_cn))
            )
        } else if (clonality == "subclonal"){
            bb_scna = bb_scna %>% dplyr::filter(
                frac1_A < 0.9 & frac1_A > 0.1 & ((nMin1_A != normal_allele_cn | nMaj1_A != normal_allele_cn) | (nMin2_A != normal_allele_cn | nMaj2_A != normal_allele_cn))
            )
        } else {
            bb_scna = bb_scna %>% dplyr::filter(
                (nMin1_A != normal_allele_cn | nMaj1_A != normal_allele_cn) | (nMin2_A != normal_allele_cn | nMaj2_A != normal_allele_cn)
            )
        }
    }
    return(bb_scna)
}

extract_breakpoints <- function(bb_scna, type_alteration, clonality){     
    # Filter start and end of chromosomes, which are not real breakpoints
    # And set up to width = 1 to only focus on real breakpoint    
    bb_scna = bb_scna %>%
        group_by(chr) %>%
        arrange(startpos) %>%
        # For the first entry in bb_scna for each chromosome, modify startpos to be the same value as endpos
        mutate(startpos = ifelse(row_number() == 1, endpos - 1, startpos)) %>%
        # For the last entry in bb_scna for each chromosome, modify endpos to be the same value as startpos
        mutate(endpos = ifelse(row_number() == n(), startpos + 1, endpos)) %>%
        ungroup() %>% 
        arrange(as.numeric(chr), startpos) %>% 
        mutate(chr = paste0("chr", chr))
    
    # Filter out artificial breakpoints when only one alteration present in chromosome
    bb_scna = bb_scna %>%
        group_by(chr) %>%
        dplyr::filter(n() > 1) %>%
        ungroup()

    # Filter alterations by type and clonality
    bb_scna = filter_alterations(bb_scna, type_alteration, clonality)
    
    # Return as GRanges
    bkps.gr = GRanges(
        seqnames = bb_scna$chr, 
        ranges = IRanges(start = bb_scna$startpos, end = bb_scna$endpos)
    )

    metadata(bkps.gr) = bb_scna %>% dplyr::select(-chr, -startpos, -endpos)
    
    return(bkps.gr)
}

add_gene_names <- function(bkps.gr, driver_genes.gr){
    # Add gene names to breakpoints GRanges
    overlaps = findOverlaps(bkps.gr, driver_genes.gr)
    bkps.gr$gene = NA
    bkps.gr$gene[queryHits(overlaps)] = driver_genes.gr$gene[subjectHits(overlaps)]
    return(bkps.gr)
}

add_replication_time <- function(bkps.gr, reptime.gr){
    start_bkps = GRanges(seqnames = seqnames(bkps.gr), ranges = IRanges(start = start(bkps.gr), end = start(bkps.gr) + 1))
    end_bkps = GRanges(seqnames = seqnames(bkps.gr), ranges = IRanges(start = end(bkps.gr) - 1, end = end(bkps.gr)))

    # Get replication time scores for start and end of breakpoints
    # taking into account potential for some breakpoints not overlapping with rep time data
    start_overlaps = findOverlaps(start_bkps, reptime.gr)
    end_overlaps = findOverlaps(end_bkps, reptime.gr)
    bkps.gr$reptime_score_start = NA
    bkps.gr$reptime_score_end = NA
    bkps.gr$reptime_score_start[queryHits(start_overlaps)] = reptime.gr$score[subjectHits(start_overlaps)]
    bkps.gr$reptime_score_end[queryHits(end_overlaps)] = reptime.gr$score[subjectHits(end_overlaps)]
    bkps.gr$reptime_score = rowMeans(cbind(bkps.gr$reptime_score_start, bkps.gr$reptime_score_end), na.rm = TRUE)
    return(bkps.gr)
}

# PATHS
outdir = "outputs/03_mechanism/mutational_processes"
figdir = "figures/03_mechanism/mutational_processes"
dir.create(outdir, recursive = TRUE); dir.create(figdir, recursive = TRUE)
trajectories_fps <- list.files("outputs/02_trajectories/", pattern = "PPCG_Feb2026.*mergedseg_with_clonality.txt", full = TRUE)
bbdir = "data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/Subclonal_SCNA"
lncap_reptime_fp = "data/ext/GSE98730_LNCaP_WA.bw"
cosmic_genes_fp = "data/meta/drivers/cancergenes_cosmic_cgc81.tsv"
ppcg_genes_fp = "data/meta/drivers/prostate_ppcg.tsv"
discriminative_genes_fp = "outputs/02_trajectories/most_discriminative_events.Rdata"
# LOAD DATA
tj_samples = load_trajectories(trajectories_fps)
reptime.gr = import(lncap_reptime_fp, format = "BigWig")
reptime.gr$class = ifelse(
    reptime.gr$score < quantile(reptime.gr$score, 0.33), "Late",
    ifelse(reptime.gr$score < quantile(reptime.gr$score, 0.66), "Middle", "Early")
)

driver_genes = rbind(
    read_delim(cosmic_genes_fp), read_delim(ppcg_genes_fp)
) %>% dplyr::filter(!duplicated(gene))

driver_genes = rbind(
    driver_genes, 
    data.frame(
        gene = c("MYC", "CHD1", "LRP1B", "ZNF292"), 
        chr = c(8, 5, 2, 6), 
        start = c(128747680,98189689, 140231423, 90513573), 
        end = c(128755197, 98264711, 142131701, 90587093)
    )
)

driver_genes.gr = GRanges(
    seqnames = paste0("chr", driver_genes$chr),
    ranges = IRanges(start = driver_genes$start, end = driver_genes$end),
    gene = driver_genes$gene
)

# 1 - ANALYSIS OF ALL LOH SEGMENTS
gene_colours = c(

    # TSGs - plot with losses
    "TP53" = "#E41A1C", "RB1" = "#377EB8", "PTEN" = "#4DAF4A", "CHD1" = "#984EA3", "BRCA1" = "#FF7F00", 
    "ERG" = "#A65628", "FOXP1" = "#E7298A", "LRP1B" = "#F781BF", "ZNF292" = "#e6ab02", "NKX3-1" = "#B2DF8A", 

    # Oncogenes - plot with gains
    "MYC" = "#800026", "PIK3CA" = "#FEB24C", "MDM4" = "#F03B20", "CCND1" = "#377EB8", "FOXA1" = "#1B9E77"
)

# Define TSG genes of interest
load(discriminative_genes_fp)
# Top enriched in Canonical

event_or_fdr %>% as.data.frame() %>% rownames_to_column("event") %>% dplyr::filter(str_detect(event, "HD|LOH")) %>% dplyr::filter(Can_OR > 1.5) %>% dplyr::arrange(Can_FDR) %>% head()
# TMPRSS2-ERG, TP53, BRCA1, PTEN, FOXP1, ZNRF3

event_or_fdr %>% as.data.frame() %>% rownames_to_column("event") %>% dplyr::filter(str_detect(event, "HD|LOH")) %>% dplyr::filter(Alt_OR > 1.5) %>% dplyr::arrange(Alt_FDR) %>% head()
# CHD1, CHD1, LRP1B, ZNF292, chr4 and chr1 losses with no known targets

# Most discriminative events: ERG, TP53, BRCA1, PTEN, FOXP1, CHD1, LRP1B, ZNF292
# Plus recurrent events that occur in all trajectories: NKX3-1 and RB1

loh_targets = c("ERG", "TP53", "BRCA1", "PTEN", "FOXP1", "CHD1", "LRP1B", "ZNF292", "RB1", "NKX3-1")
loh_targets.gr = driver_genes.gr[driver_genes.gr$gene %in% loh_targets]

# Gather breakpoints for all samples and annotate with gene names and replication time
loh_bps = foreach(i = 1:nrow(tj_samples)) %dopar% {
    sid = tj_samples$sample[i]
    evotype = tj_samples$trajectory[i]
    bb_fp = list.files(bbdir, pattern = sid, full.names = TRUE)
    bb_scna = read_delim(bb_fp) %>% dplyr::filter(chr %in% 1:22)
    bps_df = data.frame()

    for (timing in c("clonal", "subclonal")){
        loh = filter_alterations(bb_scna, "LOH", timing)
        bps.gr = extract_breakpoints(bb_scna, "LOH", timing)
        if (length(bps.gr) == 0) next
        bps.gr = add_gene_names(bps.gr, loh_targets.gr)
        bps.gr = add_replication_time(bps.gr, reptime.gr)

        bps_df = rbind(bps_df, data.frame(
            trajectory = evotype, sample = sid, timing = timing, alteration = "LOH",
            size = width(bps.gr), reptime_score = bps.gr$reptime_score, gene = bps.gr$gene
        ))
    }
    bps_df
}

loh_bps_df = data.table::rbindlist(loh_bps)
p = ggplot(loh_bps_df %>% dplyr::filter(size > 2), aes(x = reptime_score, y = size)) +
    geom_density_2d(alpha = 0.8, color = "black") +
    geom_point(size = 0.5, aes(color = gene, alpha = !is.na(gene))) + 
    geom_vline(xintercept = quantile(reptime.gr$score, c(0.33, 0.66)), linetype = "dashed", color = "black") +
    facet_wrap(~trajectory) + 
    theme(legend.position = "bottom") + 
    scale_y_log10() + 
    scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 0.05), guide = "none") +
    labs(x = "Replication timing score", y = "Size of LOH segment", color = "Gene") + 
    scale_color_manual(values = gene_colours, na.value = "grey") + 
    theme(legend.position = "none")

write_tsv(
  loh_bps_df %>% dplyr::filter(size > 2) %>% dplyr::select(trajectory, size, reptime_score, gene),
  file.path(outdir, "Fig3c_source_data.tsv")
)
ggsave(file.path(figdir, "Fig3c_LOH_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)

# Repeat plot only for clonal
p = ggplot(loh_bps_df %>% dplyr::filter(size > 2, timing == "clonal"), aes(x = reptime_score, y = size)) +
    geom_density_2d(alpha = 0.8, color = "black") +
    geom_point(size = 0.5, color = "grey", alpha = 0.05) + 
    geom_vline(xintercept = quantile(reptime.gr$score, c(0.33, 0.66)), linetype = "dashed", color = "black") +
    facet_wrap(~trajectory) + 
    theme(legend.position = "bottom") + 
    scale_y_log10() + 
    scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 0.05), guide = "none") +
    labs(x = "Replication timing score", y = "Size of clonal LOH segment", color = "Gene") + 
    theme(legend.position = "none")

write_tsv(
  loh_bps_df %>% dplyr::filter(size > 2, timing == "clonal") %>% dplyr::select(trajectory, size, reptime_score, gene),
  file.path(outdir, "EXDF7c_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF7c_LOH_clonal_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)

# Repeat plot only for subclonal
p = ggplot(loh_bps_df %>% dplyr::filter(size > 2, timing == "subclonal"), aes(x = reptime_score, y = size)) +
    geom_density_2d(alpha = 0.8, color = "black") +
    # geom_point_rast(size = 0.5, color = "grey", alpha = 0.05) +
    geom_point(size = 0.5, color = "grey", alpha = 0.05) +
    geom_vline(xintercept = quantile(reptime.gr$score, c(0.33, 0.66)), linetype = "dashed", color = "black") +
    facet_wrap(~trajectory) + 
    theme(legend.position = "bottom") + 
    scale_y_log10() + 
    labs(x = "Replication timing score", y = "Size of subclonal LOH segment", color = "Gene") +
    theme(legend.position = "none")

write_tsv(
  loh_bps_df %>% dplyr::filter(size > 2, timing == "subclonal") %>% dplyr::select(trajectory, size, reptime_score, gene),
  file.path(outdir, "EXDF7f_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF7f_LOH_subclonal_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)

# 2. ANALYSIS OF HD SEGMENTS
# same tumour suppressor genes as above
hd_targets = c("ERG", "TP53", "BRCA1", "PTEN", "FOXP1", "CHD1", "LRP1B", "ZNF292", "RB1", "NKX3-1")
hd_targets.gr = driver_genes.gr[driver_genes.gr$gene %in% hd_targets]

# Gather breakpoints for all samples and annotate with gene names and replication time
hd_bps = foreach(i = 1:nrow(tj_samples)) %dopar% {
    sid = tj_samples$sample[i]
    evotype = tj_samples$trajectory[i]
    bb_fp = list.files(bbdir, pattern = sid, full.names = TRUE)
    bb_scna = read_delim(bb_fp) %>% dplyr::filter(chr %in% 1:22)
    bps_df = data.frame()

    for (timing in c("clonal", "subclonal")){
        hd = filter_alterations(bb_scna, "HD", timing)
        bps.gr = extract_breakpoints(bb_scna, "HD", timing)
        if (length(bps.gr) == 0) next
        bps.gr = add_gene_names(bps.gr, hd_targets.gr)
        bps.gr = add_replication_time(bps.gr, reptime.gr)

        bps_df = rbind(bps_df, data.frame(
            trajectory = evotype, sample = sid, timing = timing, alteration = "HD",
            size = width(bps.gr), reptime_score = bps.gr$reptime_score, gene = bps.gr$gene
        ))
    }
    bps_df
}

hd_bps_df = data.table::rbindlist(hd_bps)
p = ggplot(hd_bps_df %>% dplyr::filter(size > 2), aes(x = reptime_score, y = size)) +
    geom_density_2d(alpha = 0.8, color = "black") +
    # geom_point_rast(size = 0.5, aes(color = gene, alpha = !is.na(gene))) +
    geom_point(size = 0.5, aes(color = gene, alpha = !is.na(gene))) + 
    geom_vline(xintercept = quantile(reptime.gr$score, c(0.33, 0.66)), linetype = "dashed", color = "black") +
    facet_wrap(~trajectory) + 
    theme(legend.position = "bottom") + 
    scale_y_log10() + 
    scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 0.05), guide = "none") +
    labs(x = "Replication timing score", y = "Size of HD segment", color = "Gene") + 
    scale_color_manual(values = gene_colours, na.value = "grey") + 
    theme(legend.position = "none")

write_tsv(
  hd_bps_df %>% dplyr::filter(size > 2) %>% dplyr::select(trajectory, size, reptime_score, gene),
  file.path(outdir, "Fig3d_source_data.tsv")
)
ggsave(file.path(figdir, "Fig3d_HD_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)

# Repeat plot only for clonal
p = ggplot(hd_bps_df %>% dplyr::filter(size > 2, timing == "clonal"), aes(x = reptime_score, y = size)) + 
    geom_density_2d(alpha = 0.8, color = "black") +
    # geom_point_rast(size = 0.5, color = "grey", alpha = 0.05) +
    geom_point(size = 0.5, color = "grey", alpha = 0.05) + 
    geom_vline(xintercept = quantile(reptime.gr$score, c(0.33, 0.66)), linetype = "dashed", color = "black") +
    facet_wrap(~trajectory) + 
    theme(legend.position = "bottom") + 
    scale_y_log10() +
    labs(x = "Replication timing score", y = "Size of clonal HD segment", color = "Gene") + 
    theme(legend.position = "none")

write_tsv(
  hd_bps_df %>% dplyr::filter(size > 2, timing == "clonal") %>% dplyr::select(trajectory, size, reptime_score, gene),
  file.path(outdir, "EXDF7d_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF7d_HD_clonal_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)

# Repeat plot only for subclonal
p = ggplot(hd_bps_df %>% dplyr::filter(size > 2, timing == "subclonal"), aes(x = reptime_score, y = size)) + 
    geom_density_2d(alpha = 0.8, color = "black") +
    # geom_point_rast(size = 0.5, color = "grey", alpha = 0.05) +
    geom_point(size = 0.5, color = "grey", alpha = 0.05) + 
    geom_vline(xintercept = quantile(reptime.gr$score, c(0.33, 0.66)), linetype = "dashed", color = "black") +
    facet_wrap(~trajectory) + 
    theme(legend.position = "bottom") + 
    scale_y_log10() + 
    labs(x = "Replication timing score", y = "Size of subclonal HD segment", color = "Gene") +
    theme(legend.position = "none")

write_tsv(
  hd_bps_df %>% dplyr::filter(size > 2, timing == "subclonal") %>% dplyr::select(trajectory, size, reptime_score, gene),
  file.path(outdir, "EXDF7g_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF7g_HD_subclonal_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)

# 3. ANALYSIS OF GAIN SEGMENTS
event_or_fdr %>% as.data.frame() %>% rownames_to_column("event") %>% dplyr::filter(str_detect(event, "Gain")) %>% dplyr::filter(GE_OR > 1.5) %>% dplyr::arrange(
GE_FDR) %>% head()

# Unclear targets except for CCND1 which is the top hit

gain_targets = c("MYC", "PIK3CA", "MDM4", "CCND1", "FOXA1")

gain_targets.gr = driver_genes.gr[driver_genes.gr$gene %in% gain_targets]

# Gather breakpoints for all samples and annotate with gene names and replication time
gain_bps = foreach(i = 1:nrow(tj_samples)) %dopar% {
    sid = tj_samples$sample[i]
    evotype = tj_samples$trajectory[i]
    bb_fp = list.files(bbdir, pattern = sid, full.names = TRUE)
    bb_scna = read_delim(bb_fp) %>% dplyr::filter(chr %in% 1:22)
    bps_df = data.frame()

    for (timing in c("clonal", "subclonal")){
        gain = filter_alterations(bb_scna, "Gain", timing)
        bps.gr = extract_breakpoints(bb_scna, "Gain", timing)
        if (length(bps.gr) == 0) next
        bps.gr = add_gene_names(bps.gr, gain_targets.gr)
        bps.gr = add_replication_time(bps.gr, reptime.gr)

        bps_df = rbind(bps_df, data.frame(
            trajectory = evotype, sample = sid, timing = timing, alteration = "Gain",
            size = width(bps.gr), reptime_score = bps.gr$reptime_score, gene = bps.gr$gene
        ))
    }
    bps_df
}

gain_bps_df = data.table::rbindlist(gain_bps)
p = ggplot(gain_bps_df %>% dplyr::filter(size > 2), aes(x = reptime_score, y = size)) + 
    geom_density_2d(alpha = 0.8, color = "black") +
    # geom_point_rast(size = 0.5, aes(color = gene, alpha = !is.na(gene))) + 
    geom_point(size = 0.5, aes(color = gene, alpha = !is.na(gene))) +
    geom_vline(xintercept = quantile(reptime.gr$score, c(0.33, 0.66)), linetype = "dashed", color = "black") +
    facet_wrap(~trajectory) + 
    theme(legend.position = "bottom") + 
    scale_y_log10() + 
    scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 0.05), guide = "none") +
    labs(x = "Replication timing score", y = "Size of GAIN segment", color = "Gene") + 
    scale_color_manual(values = gene_colours, na.value = "grey") + 
    theme(legend.position = "none")

write_tsv(
  gain_bps_df %>% dplyr::filter(size > 2) %>% dplyr::select(trajectory, size, reptime_score, gene),
  file.path(outdir, "Fig3e_source_data.tsv")
)
ggsave(file.path(figdir, "Fig3e_GAIN_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)

# Repeat plot only for clonal
p = ggplot(gain_bps_df %>% dplyr::filter(size > 2, timing == "clonal"), aes(x = reptime_score, y = size)) + 
    geom_density_2d(alpha = 0.8, color = "black") +
    # geom_point_rast(size = 0.5, color = "grey", alpha = 0.05) +
    geom_point(size = 0.5, color = "grey", alpha = 0.05) +
    geom_vline(xintercept = quantile(reptime.gr$score, c(0.33, 0.66)), linetype = "dashed", color = "black") +
    facet_wrap(~trajectory) + 
    theme(legend.position = "bottom") + 
    scale_y_log10() + 
    labs(x = "Replication timing score", y = "Size of clonal GAIN segment", color = "Gene") +
    theme(legend.position = "none")

write_tsv(
  gain_bps_df %>% dplyr::filter(size > 2, timing == "clonal") %>% dplyr::select(trajectory, size, reptime_score, gene),
  file.path(outdir, "EXDF7e_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF7e_GAIN_clonal_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)

# Repeat plot only for subclonal
p = ggplot(gain_bps_df %>% dplyr::filter(size > 2, timing == "subclonal"), aes(x = reptime_score, y = size)) + 
    geom_density_2d(alpha = 0.8, color = "black") +
    # geom_point_rast(size = 0.5, color = "grey", alpha = 0.05) +
    geom_point(size = 0.5, color = "grey", alpha = 0.05) +
    geom_vline(xintercept = quantile(reptime.gr$score, c(0.33, 0.66)), linetype = "dashed", color = "black") +
    facet_wrap(~trajectory) + 
    theme(legend.position = "bottom") + 
    scale_y_log10() + 
    labs(x = "Replication timing score", y = "Size of subclonal GAIN segment", color = "Gene") +
    theme(legend.position = "none")

write_tsv(
  gain_bps_df %>% dplyr::filter(size > 2, timing == "subclonal") %>% dplyr::select(trajectory, size, reptime_score, gene),
  file.path(outdir, "EXDF7h_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF7h_GAIN_subclonal_size_vs_reptime.pdf"), p, width = 80/25.4, height = 60/25.3)
