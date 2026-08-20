# Relationship between NRPCC and number of subclones

rm(list = ls(all = TRUE))

# PACKAGES
library(readr)
library(dplyr)
library(stringr)
library(survival)
library(ggthemes)
library(patchwork)

# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")
source("src/dpclust_parsing_functions.R")

# PATHS
outdir <- "outputs/00_preprocessing/evolutionary_metrics"
figdir <- "figures/00_preprocessing/evolutionary_metrics"
dpclust_dir <- "data/processed"
dir.create(outdir, recursive = TRUE); dir.create(figdir, recursive = TRUE)

# GET RELEVANT EVOLUTIONARY METRICS
ccf_subclone <- 0.9
dpc_fps <- list.files(dpclust_dir, pattern = "bestClusterInfo.txt", recursive = T, full = T)

nsubclones <- sapply(dpc_fps, function(dpc_fp){
    return(get_nsubclones(dpc_fp, ccf_subclone))
})

min_ccf <- sapply(dpc_fps, function(dpc_fp) {
    clusters <- read_delim(dpc_fp)
    return(min(clusters$location, na.rm = TRUE))
})

smps <- extract_ppcg_id(dpc_fps, full = F)
nrpcc <- get_nrpcc(smps)

# Association with NRPCC
df <- data.frame(sample = smps, nsubclones = nsubclones, min_ccf = min_ccf, nrpcc = nrpcc)

p1 <- ggplot(df, aes(x = nrpcc, y = as.character(nsubclones))) + 
    geom_boxplot(alpha = .4) +
    geom_vline(xintercept = 10, linetype = "dashed", col = "red") + 
    labs(x = "Number of reads per chromosome copy\n(NRPCC)", y = "Number of subclones detected")

save_ggplot(p1, file.path(figdir, "nsubclones_nrpcc"), w = 50, h = 50)

p2 <- ggplot(df, aes(x = nrpcc, y = min_ccf)) + 
    geom_point(alpha = .4) +
    geom_vline(xintercept = 10, linetype = "dashed", col = "red") + 
    geom_hline(yintercept = 0.4, linetype = "dashed", col = "blue") + 
    ylim(c(0,1)) + 
    labs(x = "NRPCC", y = "Minimum CCF of detected clusters")

save_ggplot(p2, file.path(figdir, "min_ccf_nrpcc"), w = 50, h = 50)


p <- (p1 + p2) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10))

# EXDF1a-b
write_tsv(df %>% dplyr::select(nrpcc, nsubclones, min_ccf), file.path(outdir, "EXDF1ab_source_data.tsv"))
save_ggplot(p, file.path(figdir, "EXDF1ab_nrpcc_power_subclonality"), w = 110, h = 65)
