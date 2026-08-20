# Plot dN / dS estimates
rm(list = ls(all = TRUE))

# LIBRARIES
library(readr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(ggrepel)
library(lemon)
library(data.table)
library(patchwork)

# PATHS
dnds_dir = "outputs/01_landscape/dnds"
figdir = "figures/01_landscape/dnds"
dir.create(figdir)


##### FUNCTIONS ------------
source("src/plot_theme.R")
source("src/plot_functions.R")

read_dnds = function(driver_list, ggroup, timing){
    fp = file.path(dnds_dir, paste0("dnds_", timing), paste0("ordering_all_gs_", ggroup, "_outcome_all"), driver_list, "global_dnds.tsv")
    if (!file.exists(fp)){
        return(data.frame(driver_list = driver_list, group = ggroup, timing = timing, mle = NA, cilow = NA, cihigh = NA))
    }
    dnds = read.table(fp)
    mle = dnds$mle[dnds$name == "wall"]
    cilow = dnds$cilow[dnds$name == "wall"]
    cihigh = dnds$cihigh[dnds$name == "wall"]
    dnds_df = data.frame(
        driver_list = driver_list,
        group = ggroup, 
        timing = timing,
        mle = mle, 
        cilow = cilow, 
        cihigh = cihigh,
        stringsAsFactors = FALSE
        )
    return(dnds_df)
}
  

###### MAIN ---------------

## PPCG LIST OF DRIVERS BY GLEASON
list_drivers = c("prostate_armenia_2018", "prostate_wedge_2018", "prostate_ppcg")
timings = c("clonal", "subclonal", "clonal_and_subclonal")
ggroups = c("1", "2", "3", "4+")

dnds_df = data.frame()
for (driver in list_drivers){
    for (timing in timings){
        for (ggroup in ggroups){
            # read in dnds output
            dnds_df = rbind(dnds_df, read_dnds(driver, ggroup, timing))
        }
    }
}


dodge <- position_dodge2(width = .75)
dnds_df$group = paste("Grade Group", dnds_df$group)
dnds_df$timing = str_to_title(dnds_df$timing)
dnds_df$timing = ifelse(dnds_df$timing == "Clonal_and_subclonal", "All", dnds_df$timing)
dnds_df$driver_list = str_to_title(str_replace_all(dnds_df$driver_list, "_", " "))
dnds_df$driver_list = str_remove(dnds_df$driver_list, "Prostate ")
dnds_df$driver_list = str_replace_all(dnds_df$driver_list, "Ppcg", "PPCG")

p <- dnds_df %>%
    # ensure correct ordering in the plot
    dplyr::mutate(driver_list = factor(driver_list, levels = c("PPCG", "Armenia 2018", "Wedge 2018"))) %>%
    ggplot(aes(y = group, col = timing)) + 
    geom_point(aes(x = mle), size = 3, alpha = .8, position = dodge) + 
    geom_linerange(aes(xmin = cilow, xmax = cihigh), position = dodge) + 
    # geom_segment(aes(x = cilow, xend = cihigh, y = type, yend = type), position = dodge) + 
    geom_vline(xintercept = 1, lty = "dashed") + 
    labs(x = "dN / dS", y = "") + 
    scale_y_discrete(limits = rev(c(
        "Grade Group 4+", 
        "Grade Group 3",
        "Grade Group 2", 
        "Grade Group 1"
    ))) + scale_x_log10() + facet_wrap(~driver_list, nrow = 1) +
    scale_colour_manual(values = c("Subclonal" = "#df2a55", "Clonal" = "#3979bb", "All" = "grey60")) + 
    theme(legend.position = "top")

p <- axes2lemon(p, lt = "both")

write_tsv(dnds_df, file.path("outputs/01_landscape/dnds", "EXDF2a_source_data.tsv"))
save_ggplot(p, file.path(figdir, "EXDF2a_dnds_by_gleason_ppcg"), w = 120, h = 60)
