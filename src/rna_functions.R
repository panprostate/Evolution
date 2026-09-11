
#' Filter genes by expression frequency across samples
#'
#' This function filters out genes with low expression across samples.
#' Genes must have CPM (counts per million) above `cpm_threshold`
#' in at least a fraction `freq` of samples to be kept.
#'
#' @param countData Matrix of raw counts (genes x samples)
#' @param cpm_threshold Numeric. Minimum CPM per sample to consider gene expressed (default 1)
#' @param freq Numeric between 0 and 1. Minimum fraction of samples that must pass cpm_threshold (default 0.2)
#' @return Filtered count matrix with genes passing expression threshold
expression_filter <- function(countData, cpm_threshold = 1, freq = 0.2){
    # Create DGEList object for CPM calculation
    dge <- DGEList(counts = countData)
    # Calculate CPM values
    cpm_values <- cpm(dge)

    # Calculate minimum number of samples required to pass filter
    min_samples <- ceiling(freq * ncol(countData))  

    # Boolean vector of genes to keep
    keep_genes <- rowSums(cpm_values >= cpm_threshold) >= min_samples
    # Subset count matrix to keep only those genes
    filtered_counts <- countData[keep_genes, ]
    return(filtered_counts)
} 


#' Run DESeq2 differential expression analysis
#'
#' This function runs DESeq2 given a biological variable of interest (trajectory),
#' raw count data, sample metadata, and technical covariates to adjust for.
#' It constructs the design formula, creates the DESeqDataSet, runs DESeq,
#' and returns results contrasting TRUE vs FALSE for the trajectory variable.
#'
#' @param trajectory_variable_str Character. Name of the binary trajectory variable in colData (e.g. "is_tj1")
#' @param countData Matrix of raw counts (genes x samples)
#' @param colData Data frame with sample metadata, including trajectory variable and technical covariates
#' @param technical_covariates Character vector of column names in colData representing unwanted variation to adjust for (e.g. RUV factors)
#' @return DESeqResults object with differential expression results for trajectory TRUE vs FALSE
run_deseq2 <- function(trajectory_variable_str, countData, colData, technical_covariates, control_gleason_group = TRUE, return_comparison = TRUE) {

  if (control_gleason_group) {
    # If controlling for Gleason group, add it to the technical covariates to control for it
    technical_covariates <- c(technical_covariates, "grade_group")
    is_na_grade_group <- is.na(colData$grade_group)
    colData <- colData[!is_na_grade_group, ]  # Remove samples with NA Gleason group (error in DESeq2 otherwise)
    countData <- countData[, !is_na_grade_group]  # Filter countData to match colData
  }

  # Construct the design formula string with technical covariates + biological variable
  design = paste0(
    "~", 
    paste0(technical_covariates, collapse = "+"),  # e.g. "~ W1 + W2 + W3"
    paste0("+", trajectory_variable_str)            # e.g. "+ is_tj1"
  )

  # Create DESeqDataSet from rounded count data (DESeq2 expects integers)
  dds <- DESeqDataSetFromMatrix(countData = round(countData),
                                colData = colData, 
                                design = as.formula(design)
                              )

  # Run DESeq2 pipeline (estimate size factors, dispersions, fit model, test)
  dds <- DESeq(dds)

    if (return_comparison) {
        # Extract results for trajectory TRUE vs FALSE comparison
        results <- results(dds, name = paste0(trajectory_variable_str, "TRUE"))
    } else {
        # If not returning specific comparison, return full DESeqDataSet object
        results <- dds
    }
  

  return(results)
}

run_deseq2_by_tj_grade <- function(trajectory_variable_str, grade_group_str, countData, colData, technical_covariates) {

    # Construct the design formula string with technical covariates + biological variable
    design = paste0(
        "~", 
        paste0(technical_covariates, collapse = "+"),  # e.g. "~ W1 + W2 + W3"
        paste0("+", grade_group_str)                   
    )

    # Filter to trajectory of interest
    idx = which(colData[[trajectory_variable_str]] & !is.na(colData[[grade_group_str]]))  # Keep samples with TRUE in trajectory variable and non-NA grade group
    countData = countData[, idx]  # Filter countData to keep only samples with TRUE in trajectory variable
    colData = colData[idx, ]      # Filter colData to match countData
    # Create DESeqDataSet from rounded count data (DESeq2 expects integers)
    dds <- DESeqDataSetFromMatrix(countData = round(countData),
                                  colData = colData, 
                                  design = as.formula(design)
                                  )

    # Run DESeq2 pipeline (estimate size factors, dispersions, fit model, test)
    dds <- DESeq(dds)

    # Extract results
    results <- results(dds)

    return(results)
}

#' Run gene set enrichment analysis with fgsea
#'
#' This function runs fgsea on DESeq2 results given a list of pathway gene sets.
#' It filters DESeq2 results to genes in pathways, prepares ranked gene statistics,
#' runs fgsea with default multilevel method, and returns enrichment results.
#'
#' @param deseq2_results DESeqResults object with gene-level statistics
#' @param pathways Named list of character vectors with genes per pathway
#' @return Data frame with fgsea enrichment results
run_gsea <- function(deseq2_results, pathways){

    deseq2_results = deseq2_results[rownames(deseq2_results) %in% unlist(pathways),]

    # Prepare named vector of statistics (e.g. log2FC or stat)
    ranks <- deseq2_results$log2FoldChange
    names(ranks) <- rownames(deseq2_results)

    # Run fgsea
    fgseaRes <- fgsea(pathways = pathways,
                    stats = ranks)

    # View top enriched pathways
    return(fgseaRes)
}
