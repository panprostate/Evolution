# Functions for differentially methylated cytosine post-RnBeads analysis

library(dplyr)
library(stringr)
library(data.table)

#' Filter differentially methylated cytosines
dm_filter <- function(x) x[diffmeth.p.adj.fdr <= 0.05 & abs(mean.diff) >= 0.2]

#' Get name of comparison for a give differential methylation table from RnBeads
#' get_c
get_comparison_name <- function(x) {
  # Get the name of the comparison
  comparison_1 <- stringr::str_remove(names(x)[6], "mean.")
  comparison_2 <- stringr::str_remove(names(x)[7], "mean.")
  paste0(comparison_1, "_vs_", comparison_2)
}

# Create the report
dm_report <- function(x, total_cpgs) {
  # Get the comparison name
  comparison_name <- get_comparison_name(x)

  # Get the number of differentially methylated positions
  number_of_dmcs <- nrow(x)
  number_of_dmcs_hypo <- nrow(x[mean.diff < 0])
  number_of_dmcs_hyper <- nrow(x[mean.diff > 0])

  # Get percentages
  percentage_of_cpgs_dmc <- number_of_dmcs * 100 / total_cpgs
  percentage_of_dmcs_hypo <- number_of_dmcs_hypo * 100 / number_of_dmcs
  percentage_of_dmcs_hyper <- number_of_dmcs_hyper * 100 / number_of_dmcs

  # Report
  data.frame(
    comparison_name = comparison_name,
    number_of_dmcs = number_of_dmcs,
    number_of_dmcs_hypo = number_of_dmcs_hypo,
    number_of_dmcs_hyper = number_of_dmcs_hyper,
    percentage_of_cpgs_dmc = percentage_of_cpgs_dmc,
    percentage_of_dmcs_hypo = percentage_of_dmcs_hypo,
    percentage_of_dmcs_hyper = percentage_of_dmcs_hyper
  )
}

#' Print the report
show_dm_report <- function(report) {
  cat(
    paste0(
      c(
        report$comparison_name, "\n",
        "Total DMCs:", report$number_of_dmcs, "\n",
        "Hypomethylated (in comparision_1):", report$number_of_dmcs_hypo, "\n",
        "Hypermethylated (in comparision_1):", report$number_of_dmcs_hyper, "\n\n"
      ),
      collapse = ""
    )
  )
}

#' Return an object with annotations
dmc_annotation_reporter <- function(dmcs, annotation) {
  dm <- dmcs %>% dplyr::left_join(annotation, by = "cgid")
  dm <- dm %>% dplyr::mutate(mean.diff.direction = ifelse(mean.diff > 0, "hyper", "hypo"))
  comparison_name <- get_comparison_name(dmcs)

  report_cgi <- dm %>%
    dplyr::mutate(cpg_island_label = dplyr::case_when(
      stringr::str_detect(Relation_to_UCSC_CpG_Island, "Shelf") ~ "Shelf",
      stringr::str_detect(Relation_to_UCSC_CpG_Island, "Shore") ~ "Shore",
      stringr::str_detect(Relation_to_UCSC_CpG_Island, "Island") ~ "Island",
      TRUE ~ "Open Sea"
    )) %>%
    dplyr::count(cpg_island_label, mean.diff.direction) %>%
    dplyr::mutate(perc = n / sum(n) * 100) %>%
    dplyr::mutate(comparison_name = comparison_name)

  report_tss <- dm %>%
    dplyr::filter(distance_to_TSS <= 500) %>%
    dplyr::select(cgid, distance_to_TSS, nearestTSS, mean.diff.direction) %>%
    dplyr::count(nearestTSS, mean.diff.direction) %>%
    dplyr::arrange(-n) %>%
    dplyr::filter(n > 1) %>%
    tidyr::pivot_wider(id_cols = nearestTSS, names_from = mean.diff.direction, values_from = n) %>%
    dplyr::mutate(comparison_name = comparison_name)

  report_atac <- dm %>%
    dplyr::count(Corces_ATAC, mean.diff.direction) %>%
    dplyr::mutate(perc = n / sum(n) * 100) %>%
    dplyr::mutate(comparison_name = comparison_name)

  report_pmd <- dm %>%
    dplyr::count(Prostate_PMD_Guo, mean.diff.direction) %>%
    dplyr::mutate(perc = n / sum(n) * 100) %>%
    dplyr::mutate(comparison_name = comparison_name)

  report_chromhmm <- dm %>%
    dplyr::count(ChromHMM_Pomerantz, mean.diff.direction) %>%
    dplyr::mutate(perc = n / sum(n) * 100) %>%
    dplyr::mutate(comparison_name = comparison_name)

  list(
    "comparison_name" = comparison_name,
    "cgi" = report_cgi,
    "tss" = report_tss,
    "atac" = report_atac,
    "pmd" = report_pmd,
    "chromhmm" = report_chromhmm
  )
}

write_annotation_summary <- function(dmc_annotation_report, f_outdir) {
  comparison_name <- dmc_annotation_report$comparison_name
  fields <- setdiff(names(dmc_annotation_report), "comparison_name")
  for (field in fields) {
    outfile <- fs::path(f_outdir, paste0(field, "_", comparison_name, ".csv"))
    data.table::fwrite(
      dmc_annotation_report[[field]],
      outfile
    )
    print(paste0("Written to: ", outfile))
  }
}

#' Read formatted DMC files for each ordering
read_dmc_files <- function(f_outdir, prefix) {
  glob <- paste0("*", prefix, "*")
  data.table::rbindlist(
    lapply(
      fs::dir_ls(f_outdir, glob = glob),
      data.table::fread
    ),
    fill = TRUE
  )
}

clean_annotations <- function(dmc_annotated){
  
  # Add a clean CGI variable
  dmc_annotated <- dmc_annotated %>%
    dplyr::mutate(cpg_island_label = dplyr::case_when(
      stringr::str_detect(Relation_to_UCSC_CpG_Island, "Shelf") ~ "Shelf",
      stringr::str_detect(Relation_to_UCSC_CpG_Island, "Shore") ~ "Shore",
      stringr::str_detect(Relation_to_UCSC_CpG_Island, "Island") ~ "Island",
      TRUE ~ "Open Sea"
    ))

  # Add a clean TSS variable
  valid_tss <- dmc_annotated %>% dplyr::filter(
      distance_to_TSS <= 500
  ) %>% dplyr::count(nearestTSS, mean.diff.direction) %>%
    dplyr::arrange(-n) %>%
    dplyr::filter(n >= 5) %>%
    pull(nearestTSS)
  dmc_annotated <- dmc_annotated %>% dplyr::mutate(
    TSS_500 = dplyr::case_when(
    nearestTSS %in% valid_tss ~ nearestTSS,
    TRUE ~ NA
  ))

}

add_mean_diff_direction <- function(x){
  x %>% dplyr::mutate(mean.diff.direction = ifelse(mean.diff > 0, "hyper", "hypo"))
}

get_average_methylation_at_dm_anno <- function(meth, dm, anno_field, anno_subset, dmr_direction){
  # Get DNAMe at CpGs of interest
  cg_selection = dm$cgid[(dm[[anno_field]] == anno_subset) & (dm[["mean.diff.direction"]]==dmr_direction)]
  n_cpgs = length(cg_selection)
  avmeth = colMeans(meth[cg_selection, ], na.rm=TRUE)
  ncpgs = colSums(!is.na(meth[cg_selection, ]))
  avmeth = melt(avmeth, value.name="meth")
  nmeth = melt(ncpgs, value.name="ncpgs")
  avmeth$PPCG_Meth_Sample_ID = rownames(avmeth)
  nmeth$PPCG_Meth_Sample_ID = rownames(nmeth)
  avmeth = avmeth %>% 
    dplyr::left_join(nmeth, by="PPCG_Meth_Sample_ID") %>%
    dplyr::left_join(
    dplyr::select(rnb_set_f@pheno, PPCG_Meth_Sample_ID, Ordering)
  , by="PPCG_Meth_Sample_ID") %>%
    dplyr::mutate(anno_field=anno_field, anno_subset=anno_subset)

  # Fix orddering names
  avmeth$ordering = avmeth$Ordering
  avmeth = set_ordering_names(avmeth)
  return(avmeth)
}
