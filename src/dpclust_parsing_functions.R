# Given the path to DPClust output files, obtain the number of subclones
# Assumption: Subclone will have CCF < 0.9 (ccf_subclone) 
# and cannot be the clone with largest CCF
# Otherwise, the mutation cluster is clonal
get_nsubclones <- function(dpc_fp, ccf_subclone = .9) {
    df <- read_delim(dpc_fp)
    # Remove the top clone
    df <- dplyr::arrange(df, desc(location))
    df <- df[-1, ]
    # Get all subclones
    return(sum(df$location < ccf_subclone))
}


# Given the path to DPClust output file, obtain the percentage of total
# mutations that are subclonal in a sample
# Assumption: Subclone will have CCF < 0.9 (ccf_subclone) 
# and cannot be the clone with largest CCF
# Otherwise, the mutation is clonal
get_mut_ith <- function(dpc_fp, ccf_subclone = .9) {
    df <- read_delim(dpc_fp)
    # get total number of mutations in DPClust file
    total_muts <- sum(df$no.of.mutations)

    # remove the dominant clone
    # assumption: it is clonal
    df <- dplyr::arrange(df, desc(location))
    df <- df[-1, ]
    # Find subclones
    is_subclone <- df$location < ccf_subclone
    subclonal_muts <- sum(df$no.of.mutations[is_subclone])
    # compute percentage of subclonal mutations
    return(subclonal_muts / total_muts)
}

# Given the path to DPClust output file, obtain the 
# total number of subclonal mutations in a sample
# Assumption: Subclone will have CCF < 0.9 (ccf_subclone) 
# and cannot be the clone with largest CCF
# Otherwise, the mutation is clonal
get_subclone_size <- function(dpc_fp, ccf_subclone) {
    df <- read_delim(dpc_fp)
    # remove dominant clone
    # assumption: it is a clone
    df <- dplyr::arrange(df, desc(location))
    df <- df[-1, ]
    subclonal_muts <- sum(df$no.of.mutations[df$location < ccf_subclone])
    # compute total number of subclonal mutations
    return(subclonal_muts)
}


# Given the path to DPClust output file, obtain the 
# CCF of the largest subclonal expansion in a sample
# Assumption: Subclone will have CCF < 0.9 (ccf_subclone) 
# and cannot be the clone with largest CCF
# Otherwise, the mutation is clonal
get_largest_subclone <- function(dpc_fp, ccf_subclone){
    df <- read_delim(dpc_fp)
    # remove dominant clone
    # assumption: it is a clone
    df <- dplyr::arrange(df, desc(location))
    df <- df[-1, ]
    ccfs_subclones <- df$location[df$location < ccf_subclone]
    if (length(ccfs_subclones) == 0){return(NA)}
    return(max(ccfs_subclones))
}

# Given the path to DPClust output file, obtain the 
# sum of the subclone CCFs in a sample
# Assumption: Subclone will have CCF < 0.9 (ccf_subclone) 
# and cannot be the clone with largest CCF
# Otherwise, the mutation is clonal
get_sum_ccf_subclones <- function(dpc_fp, ccf_subclone){
    df <- read_delim(dpc_fp)
    # remove dominant clone
    # assumption: it is a clone
    df <- dplyr::arrange(df, desc(location))
    df <- df[-1, ]
    ccfs_subclones <- df$location[df$location < ccf_subclone]
    if (length(ccfs_subclones) == 0){return(NA)}
    return(sum(ccfs_subclones, na.rm = T))
}

# Given the path to Battenberg output file, obtain the 
# proportion of genome altered by subclonal SCNAs in a sample
get_subclonal_pga <- function(batt_fp, wgd_aware = TRUE, wgd_samps = NULL) {
    df <- read_delim(batt_fp)
    df <- df[!is.na(df$nMaj1_A),]
    # avoid chromosome X in calculation, which is trickier to estimate
    df <- df[df$chr != "X",]

    df$length <- df$endpos - df$startpos
    total_length <- sum(df$length, na.rm = T)

    # filter to subclonal events
    df <- df[df$frac1_A != 1, ]

    if (wgd_aware){
        smp = extract_ppcg_id(batt_fp, full = F)
        if (!is.null(wgd_samps)){
            is_wgd = smp %in% wgd_samps
        } else {
            is_wgd = smp %in% get_wgd(qc = FALSE)
        }
        if (is_wgd){
            df$is_cna <- !(df$nMaj1_A == 2 & df$nMin1_A == 2)
        } else {
            df$is_cna <- !(df$nMaj1_A == 1 & df$nMin1_A == 1)
        }
    }

    # by definition, if subclonal segment, then we must have at least one subclone
    # with CN different to euploidy / tetraploidy (i.e. no further filtering needed)
    return(sum(df$length) / total_length)
}

# Given the path to Battenberg output file, obtain the number of SCNA events in a sample
get_nscna <- function(batt_fp, wgd_aware = TRUE) {
    df <- read_delim(batt_fp)
    df <- df[!is.na(df$nMaj1_A),]
    # avoid chromosome X in calculation, which is trickier to estimate
    df <- df[df$chr != "X",]
    if (wgd_aware){
        smp = extract_ppcg_id(batt_fp, full = F)
        is_wgd = smp %in% get_wgd()
        if (is_wgd){
            df$is_cna <- !(df$nMaj1_A == 2 & df$nMin1_A == 2)
        } else {
            df$is_cna <- !(df$nMaj1_A == 1 & df$nMin1_A == 1)
        }
    }
    return(sum(df$is_cna, na.rm = T))
}


# Given the path to Battenberg output file, obtain the 
# proportion of genome altered by subclonal SCNAs in a sample
# wgd_aware: if the sample is WGD, we consider as a change anything
# different than 2+2
get_clonal_pga <- function(batt_fp, wgd_aware = TRUE, wgd_samps = NULL) {
    df <- read_delim(batt_fp)
    df <- df[!is.na(df$nMaj1_A),]
    # avoid chromosome X in calculation, which is trickier to estimate
    df <- df[df$chr != "X",]

    df$length <- df$endpos - df$startpos
    total_length <- sum(df$length, na.rm = T)

    # filter to clonal events
    df <- df[df$frac1_A == 1, ]
    # define clonal CNAs as deviation from euploidy / tetraploidy
    if (wgd_aware){
        smp = extract_ppcg_id(batt_fp, full = F)
        if (!is.null(wgd_samps)){
            is_wgd = smp %in% wgd_samps
        } else {
            is_wgd = smp %in% get_wgd(qc = FALSE)
        }
        if (is_wgd){
            df$is_cna <- !(df$nMaj1_A == 2 & df$nMin1_A == 2)
        } else {
            df$is_cna <- !(df$nMaj1_A == 1 & df$nMin1_A == 1)
        }
    }
    df <- df[df$is_cna & !is.na(df$is_cna),]
    return(sum(df$length) / total_length)
}