# Test temporal changes in signature activities per evolutionary subtype


rm(list = ls(all = TRUE))  # Clear environment

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)   
library(qs)
library(data.table)
library(scales)
library(ggh4x)
library(lemon)

# FUNCTIONS ---------------------------------------------------------------
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

# See Gerstung supplement 5.3
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1907-7/MediaObjects/41586_2019_1907_MOESM1_ESM.pdf#page=18.81
pcawg_fold_change <- function(e1, e2) {
  (e2 / (1 - e2)) / (e1 / (1 - e1))
}

rev_fold_change <- function(fc) {
  2^(fc)
} # Log2 is used in dm file

# PATHS ---------------------------------------------------------------
sigtimer_fp <- "outputs/00_preprocessing/sigtimer/sigtimer_nnls.csv"
clinfile <- "data/meta/PPCG_donors_clin_20241217.csv"

outdir <- "outputs/01_landscape/sig_timing"; dir.create(outdir, recursive = TRUE)
figdir <- "figures/01_landscape/sig_timing"; dir.create(figdir, recursive = TRUE)

# LOAD DATA ---------------------------------------------------------------
qc_pass = get_prim_smps()
sigtimer = read_delim(sigtimer_fp, delim = ",") %>% 
    dplyr::rename(sample = id) %>%
    dplyr::filter(epoch == "clonal_subclonal", sigtype == "sbs") %>%
    dplyr::filter(method == "nnls") %>% 
    dplyr::filter(sample %in% qc_pass) %>% 
    dplyr::mutate(grade_group = get_gs_group(sample))

sigtimer$fc <- rev_fold_change(sigtimer$prop)
sigtimer$fc_cln <- dplyr::case_when(sigtimer$fc > 100 ~ 100, sigtimer$fc < 0.01 ~ 0.01, TRUE ~ sigtimer$fc)


# ANALYSIS ---------------------------------------------------------------
# Order signatures by prevalence
sigcount <- dplyr::count(sigtimer, signature)
sigprev <- dplyr::arrange(sigcount, n) %>%
  pull(signature) %>%
  rev() %>%
  unique()

sigtimer$signature <- factor(sigtimer$signature, levels = sigprev)
sigtimer$sigord <- as.numeric(sigtimer$signature)

sigkeep <- dplyr::filter(sigcount, n >= 10) %>%
    dplyr::filter(signature != "error") %>%
    pull(signature) %>%
    unique()

sigtimerk <- dplyr::filter(sigtimer, signature %in% sigkeep)

# Plot fold change of signatures between clonal and subclonal mutations, faceted by signature type and colored by grade group
p = sigtimerk  %>%
    dplyr::filter(epoch == "clonal_subclonal" & !is.na(grade_group)) %>%
    ggplot(aes(x = signature, y = fc_cln, col = grade_group)) + 
    geom_hline(yintercept = 1, linetype = "dashed", alpha = .4) +
    geom_jitter(size = 0.5, alpha = .5, position = position_dodge(width = .8)) +
    geom_boxplot(outlier.shape = NA, alpha = .7) +
    labs(y = "Fold Change", x = "") +
    theme(axis.text.x = element_text(angle = 45, vjust = .9, hjust = .9)) +
    scale_y_log10(n.breaks = 15, labels = label_number(accuracy = 0.01)) +
    scale_color_manual(values = gleason_colours) + # gleason_colours comes from plot_theme.R
    labs(color = "") +
    geom_text(aes(label = n, x = signature, y = 200), size = 2, position = position_dodge(width = .8),
        data = dplyr::count(sigtimerk, signature), inherit.aes = F, check_overlap = F
    ) +
    theme(legend.position='top') + 
    scale_x_discrete(guide = "axis_nested")

write_tsv(
  sigtimerk %>% dplyr::filter(epoch == "clonal_subclonal" & !is.na(grade_group)) %>%
    dplyr::select(grade_group, signature, fc_cln),
  file.path(outdir, "EXDF2b_source_data.tsv")
)
ggsave(file.path(figdir, "EXDF2b_sig_change_clonal_subclonal_by_grade_group.png"), width = 125/25.3, height = 75/25.3, dpi = 300)
ggsave(file.path(figdir, "EXDF2b_sig_change_clonal_subclonal_by_grade_group.pdf"), width = 125/25.3, height = 75/25.3, dpi = 300)

sigtimerk %>% dplyr::group_by(signature) %>%
    dplyr::summarise(median_fc = median(fc, na.rm = TRUE), n = n()) %>%
    arrange(median_fc)
