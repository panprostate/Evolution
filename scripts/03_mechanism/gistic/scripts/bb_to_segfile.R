# Transform Battenberg BB files to GISTIC segfile format

rm(list = ls(all = TRUE))

library(argparse)
library(data.table)
library(tidyverse)
library(doMC)
library(foreach)

# Parse command line arguments
parser <- ArgumentParser(description = "Transform Battenberg BB files to GISTIC segfile format")
parser$add_argument("--bb_copynumber_dir", type = "character", required = TRUE, help = "Directory containing Battenberg BB files")
parser$add_argument("--trajectory", type = "character", required = TRUE, help = "Specific trajectory to process")
parser$add_argument("--trajectories_dir", type = "character", required = TRUE, help = "Directory containing trajectory files")
parser$add_argument("--ncores", type = "integer", required = TRUE, help = "Number of cores to use for parallel processing")
args <- parser$parse_args()
bb_copynumber_dir <- args$bb_copynumber_dir
trajectories_dir <- args$trajectories_dir
trajectory <- args$trajectory
ncores <- args$ncores

# FUNCTIONS
source("src/utils.R")

bb_to_gistic_input <- function(bb_fp, nprobes = 10){
    # processed as per BISCUT guidelines which are akin to GISTIC
    # https://github.com/beroukhim-lab/BISCUT-py3/tree/main
    bb_df = read_delim(bb_fp)
    bb_df$sample = extract_ppcg_id(bb_fp, full = F)
    bb_df$width = bb_df$endpos - bb_df$startpos
    bb_df$Avg_CN = ifelse(
        bb_df$frac1_A < 1, 
        bb_df$frac1_A * (bb_df$nMaj1_A + bb_df$nMin1_A) + bb_df$frac2_A * (bb_df$nMaj2_A + bb_df$nMin2_A), 
        bb_df$nMaj1_A + bb_df$nMin1_A
    )
    ploidy = weighted.mean(bb_df$Avg_CN, bb_df$width)
    bb_df$Seg.CN = log2((bb_df$Avg_CN / ploidy) + 1e-5)
    bb_df$num_probes = nprobes
    gistic_segfile = bb_df %>% 
            dplyr::filter(!is.na(Seg.CN), width > 1) %>%
            dplyr::select(sample, chr, startpos, endpos, num_probes, Seg.CN)
    # colnames(gistic_segfile) <- c("Sample", "Chromosome", "Start Position", "End Position", "Num Markers", "Seg.CN")
    return(gistic_segfile)
}


# LOAD DATA
trajectories_fps <- list.files(trajectories_dir, pattern = "PPCG_Feb2026", full = T)
bb_fps <- list.files(bb_copynumber_dir, full = TRUE)

# LOAD DATA
tj_samples = load_trajectories(trajectories_fps)
if (trajectory %in% c("ordering_1", "ordering_2", "ordering_3")) {
    tj <- str_to_title(gsub("_", " ", trajectory))
    samples <- unique(tj_samples$sample[tj_samples$trajectory == tj])
} else if (trajectory == "all") {
    samples <- unique(tj_samples$sample)
} else {
    stop("Invalid trajectory specified. Choose from 'ordering_1', 'ordering_2', 'ordering_3', or 'all'.")
}

bb_fps = bb_fps[str_remove(basename(bb_fps), "_vs.*") %in% samples]

registerDoMC(cores = ncores)
gistic_segfile = foreach(bb_fp = bb_fps, .combine = rbind, .verbose = TRUE) %dopar% {
    bb_to_gistic_input(bb_fp)
}

write_delim(gistic_segfile, 
            paste0(trajectory, "_gistic_input.txt"),
            delim = "\t")