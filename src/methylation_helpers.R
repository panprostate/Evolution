ordering_names_recode <- c(
  "ordering_1" = "Canonical",
  "ordering_2" = "Alternative\nAR-driven",
  "ordering_3" = "Alternative\nMYC-driven",
  "normal" = "Normal"
)

set_ordering_names <- function(df) {
  df %>%
    dplyr::mutate(ordering = dplyr::recode(
      ordering,
      "ordering_1" = "Canonical",
      "ordering_2" = "Alternative\nAR-driven",
      "ordering_3" = "Alternative\nMYC-driven",
      "normal" = "Normal"
    ))
}

ggsignif_dunn_comparison <- function(comparison) {
  library(stringr)
  splits <- stringr::str_split(comparison, " - ")
  lapply(
    splits,
    function(x) {
      ordering_names_recode[x]
    }
  )
}

ggsignif_pvalues_to_stars <- function(p_values) {
  sapply(
    p_values,
    function(p) {
      if (p < 0.001) {
        return("***")
      } else if (p < 0.01) {
        return("**")
      } else if (p < 0.05) {
        return("*")
      } else {
        return("N.S.")
      }
    }
  )
}


#' Perform kruskal wallis testing on a set of fields with bonferroni correction
#' @param data A data table in long-format with testing data
#' @param fields A character vector of columns to iteratively test
#' @param response A string literal of the response column to test against fields
kruskal_wallis_test <- function(data, fields, response) {
  # Perform Kruskal-Wallis tests for each field
  kwres <- fields %>%
    lapply(function(field) {
      test_result <- kruskal.test(data[[field]] ~ data[[response]])

      data.frame(
        field = field,
        statistic = unname(test_result$statistic),
        df = unname(test_result$parameter),
        p_value = test_result$p.value
      )
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  kwres$p.adj <- p.adjust(kwres$p_value, method = "bonferroni")
  kwres <- kwres %>%
    dplyr::select(field, statistic, p_value, p.adj) %>%
    dplyr::arrange(p.adj) %>%
    dplyr::mutate(across(where(is.list), ~ map_chr(.x, paste, collapse = ", ")))
  kwres
}

#' Perform post-hoc dunn tests
#' @param data A data table in long-format with testing data
#' @param fields A character vector of columns to iteratively test
#' @param response A string literal of the response column to test against fields
dunn_post_hoc <- function(data, fields, response) {
  library(FSA)
  res <- fields %>%
    lapply(function(field) {
      dunn_res <- FSA::dunnTest(data[[field]] ~ as.factor(data[[response]]))
      dunn_res$res$field <- field
      dunn_res$res
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  res
}

#' Clean RnBeads metadata to select samples for methylation analysis
#' @param df RnBeads metadata from the `@pheno` slot.
find_methylation_qc_fail_samples <- function(df) {
  # Set QC filters
  pass_samples <- df %>%
    dplyr::filter(`sel.1234` == "1") %>%
    dplyr::filter(`Sample type` == "TUM") %>%
    dplyr::filter(!is.na(`Grade group`)) %>%
    dplyr::filter(!is.na(`Age`)) %>%
    dplyr::pull(Matching_WGS_Sample)

  # Return boolean where TRUE if samples fail QC filters
  bool_fail_samples <- !(df$Matching_WGS_Sample %in% pass_samples)

  bool_fail_samples
}

#' Set RnBeads options for differential methylation analysis
#' @param analysis_name String with which to name the analysis
#' @param identifiers_column Field with sample identifiers for analysis
#' @param differential_comparison_columns Fields to use for differential methylation analysis
set_rnbeads_analysis_options <- function(analysis_name, identifiers_column, differential_comparison_columns = c("Ordering")) {
  # Set options
  rnb.options(
    import = TRUE,
    analysis.name = analysis_name,
    identifiers.column = identifiers_column,
    import.table.separator = ",",
    import.sex.prediction = FALSE,
    normalization.method = "none",
    normalization.background.method = "none",
    normalization.plot.shifts = FALSE,
    disk.dump.big.matrices = TRUE,
    disk.dump.bigff = TRUE,
    preprocessing = FALSE,
    qc = FALSE,
    qc.boxplots = FALSE,
    qc.barplots = FALSE,
    qc.negative.boxplot = FALSE,
    qc.sample.batch.size = 50,
    qc.snp.distances = FALSE,
    qc.snp.purity = FALSE,
    qc.snp.boxplot = FALSE,
    qc.snp.heatmap = FALSE,
    filtering.greedycut = FALSE,
    filtering.greedycut.pvalue.threshold = 0.01,
    filtering.missing.value.quantile = 1,
    filtering.snp = "3",
    filtering.sex.chromosomes.removal = FALSE,
    filtering.cross.reactive = T,
    filtering.blacklist = NULL,
    min.group.size = 1,
    exploratory = FALSE,
    exploratory.correlation.qc = FALSE,
    exploratory.beta.distribution = FALSE,
    exploratory.intersample = FALSE,
    exploratory.deviation.plots = FALSE,
    inference = FALSE,
    inference.age.prediction = FALSE,
    differential = TRUE,
    differential.comparison.columns = differential_comparison_columns,
    differential.comparison.columns.all.pairwise = differential_comparison_columns,
    differential.adjustment.celltype = FALSE,
    covariate.adjustment.columns = c("Country", "Age", "Ethnicity", "Grade group", "avg.purity", "EpiCMIT hypo"),
    differential.enrichment.go = FALSE,
    differential.report.sites = TRUE,
    export.to.csv = TRUE,
    export.to.bed = FALSE,
    export.to.trackhub = c()
  )
}
