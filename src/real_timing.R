source("src/spoilsport.R")

# Add column for SNV in pyrimidine context
annotate_tnc <- function(vcf, genome, style="NCBI"){
  # Create GRanges object of SNV trinucleotide context positions
  vcf = vcf[seqnames(vcf) %in% c(1:22, "X", "Y")]
  vcf_ranges = GRanges(as.character(seqnames(vcf)), ranges=GenomicRanges::ranges(vcf))
  vcf_ranges_tri = resize(vcf_ranges, 3, fix="center")
  seqlevelsStyle(vcf_ranges_tri) = "NCBI"
  
  genome = genome[intersect(names(genome),seqnames(vcf_ranges_tri))]
  # Get trinucleotide context from reference genome
  context = data.frame(getSeq(genome, vcf_ranges_tri))[,1]
  rev.context = data.frame(reverseComplement(DNAStringSet(context)))[,1]
  
  # Set the change names
  ref = as.character(ref(vcf))
  alt = as.character(unlist(alt(vcf)))
  change = dplyr::case_when(
    ( ref == "C" & alt == "A") | (ref == "G" & alt == "T") ~ "C>A",
    ( ref == "C" & alt == "G") | (ref == "G" & alt == "C") ~ "C>G",
    ( ref == "C" & alt == "T") | (ref == "G" & alt == "A") ~ "C>T",
    ( ref == "T" & alt == "A") | (ref == "A" & alt == "T") ~ "T>A",
    ( ref == "T" & alt == "C") | (ref == "A" & alt == "G") ~ "T>C",
    ( ref == "T" & alt == "G") | (ref == "A" & alt == "C") ~ "T>G",
    TRUE ~ NA_character_
  )
  
  # Determine the bases either side of the pyrimidine trinucleotide context, using reverse comp for non pyrimidines
  pyrimidine = c("C", "T")
  base_before = ifelse(ref %in% pyrimidine, substring(context, 1, 1), substring(rev.context, 1, 1))
  base_after = ifelse(ref %in% pyrimidine, substring(context, 3), substring(rev.context, 3))
  
  # Annotate trinucleotide context and return
  vcf@info$tnc = paste0(base_before, "_",change,"_",base_after)
  
  return(vcf)
}

# Add column with the SBS signature assignment info
# Returns vcf object but with two new columns in their `info` slot:
# - sbs_sig: SBS signature with highest probability assignment
# - sbs_prob: Probability that mutation caused by sbs_sig
#' @vcf vcf object with mutations where to add SBS assignment
#' @id PPCG Sample ID for the vcf object
#' @sbs_a dataframe with the SBS assignment. Better to be passed than constantly read
#' @r flag to read SBS assignment table from sbs_assign_fp, overriding input @sbs_a
add_sbs_sig <- function(vcf, id, sbs_a = NULL, r = FALSE, target_sigs = c("SBS1", "SBS5", "SBS40a", "SBS40b"), sbs_assign_fp = "data/processed/SNV_signature_annotations__2024-06-21/SNVs_all_signature_20240621.tsv"){
  # process inputs
  if (r){
    sbs_a <- read_delim(sbs_assign_fp)
  }
  if (length(sbs_a) == 0){stop("Please provide SBS assignment table or specify path to be read from")}
  if (!grepl("_DNA", id)){id = paste0(id, "_DNA")}
  if (!id %in% sbs_a$Sample_ID){stop(glue::glue("{id} is not in SBS assignment"))}
  # filter to only sample of interest -> avoid potentially finding repeated muts in diff samples
  sbs_a <- sbs_a[sbs_a$Sample_ID == id, ]
  if (!"m" %in% colnames(sbs_a)){
    sbs_a$m <- paste0(sbs_a$chr, ":", sbs_a$start, "_", sbs_a$ref, "/", sbs_a$alt)
  }

  # match mutations and store SBS signature assignment in vcf info
  m <- rownames(vcf)
  if (!grepl("chr", m[1])){m <- paste0("chr", m)}
  is_asg <- mean(m %in% sbs_a$m)
  if (is_asg < 0.8){stop(glue::glue("{id} has over 20% mutations missing SBS sig assignment"))}
  sbs_a <- sbs_a[sbs_a$m %in% m,]
  idx <- match(sbs_a$m, m)
  idx <- na.omit(idx)
  info(vcf)$sbs_sig <- NA
  info(vcf)$sbs_sig[idx] <- sbs_a$Top_Signature
  info(vcf)$sbs_prob <- NA
  info(vcf)$sbs_prob[idx] <- sbs_a$Top_Signature_Probability
  for (sig in target_sigs){
    if (sig %in% colnames(sbs_a)){
      info(vcf)[[sig]] <- NA
      info(vcf)[[sig]][idx] <- sbs_a[[sig]]
    }
  }
  return(vcf)
}

# Gets all the CpG > T mutations (clock-wise deamination) from vcf file
# Adapted from [Moritz Gerstung et al](https://gerstung-lab.github.io/PCAWG-11/#10_chronological_wgd__mrca)
is_deamination <- function(vcf) grepl("(A_C>T_G)|(C_C>T_G)|(G_C>T_G)|(T_C>T_G)", vcf@info$tnc)

# Get all the mutations where hard assignment to clock-like signature 
is_clock <- function(vcf, clock_sigs = c("SBS1", "SBS5", "SBS40a", "SBS40b")) {vcf@info$sbs_sig %in% clock_sigs}

# parse VCF object and annotate
parse_vcf_timing = function(vcf_fp, genome, add_tnc = TRUE, add_sbs = TRUE, filter_clone = NULL, sbs_a = NULL){
    vcf = readVcf(vcf_fp)
    id = extract_ppcg_id(vcf_fp, full = F)
    if (add_tnc){
        # add into @info the TNC information
        vcf = annotate_tnc(vcf, genome)
        vcf@info$w = is_deamination(vcf)
    }
    if (add_sbs){
        # add into @info the SBS signature assignments
        vcf = add_sbs_sig(vcf, id, sbs_a)
    }
    # filter to clonal period
    if (length(filter_clone) > 0){
        if (filter_clone == "clonal"){
            filter_clone = c("clonal [NA]", "clonal [early]", "clonal [late]")
        }
        vcf = vcf[!is.na(vcf@info$CLS) & vcf@info$CLS %in% filter_clone]
    }
    return(vcf)
}

# get deamination mutations
get_cpg <- function(vcf_fp, genome, filter_clone){
    vcf <- parse_vcf_timing(vcf_fp, genome, add_tnc = TRUE, add_sbs = FALSE, filter_clone = filter_clone, sbs_a = NULL)
    return(sum(vcf@info$w, na.rm = T))
}

# get hard assignment mutations
get_sbs_sig <- function(vcf_fp, genome, filter_clone, sbs_a, target_sigs, simplify = T){
    vcf <- parse_vcf_timing(vcf_fp, genome, add_tnc = FALSE, add_sbs = TRUE, filter_clone = filter_clone, sbs_a = sbs_a)
    o <- sapply(target_sigs, function(sig){sum(vcf@info$sbs_sig == sig, na.rm = T)})
    names(o) <- target_sigs
    if (simplify){
      return(sum(o, na.rm = T))
    } else {
      return(as.data.frame(t(o)))
    }
}

# get total sum probability assignment
# assuming with parse_vcf_timing I can add a column for prob of each sig
get_sbs_sig_sp <- function(vcf_fp, genome, filter_clone, sbs_a, target_sigs, simplify = T){
    vcf <- parse_vcf_timing(vcf_fp, genome, add_tnc = FALSE, add_sbs = TRUE, filter_clone = filter_clone, sbs_a = sbs_a)
    v <- vcf@info
    o <- sapply(target_sigs, function(sig){sum(v[[sig]], na.rm = T)})
    names(o) <- target_sigs
    if (simplify){
      return(sum(o, na.rm = T))
    } else {
      return(as.data.frame(t(o)))
    }
}

power_corr <- function(n.subclonal, id, purdir, cvg){
  pur = get_purity(purdir, id)
  plo = get_ploidy(purdir, id)
  cov = round(cvg[["Tumour.total.depth"]][cvg$PPCG_Sample_ID == id])
  # SpoilSport requires cellular prevelances, so input those
  # 4th argument is the typical number of reads required to call a variant: suggested value 3
  wcc_out = wcc(pur, plo, cov, 3, n.subclonal$Group.1)

  # convert SpoilSport output columns back to corrected CCF and number of mutations
  pcCCF = wcc_out$tp / pur
  pcMut = n.subclonal$x * wcc_out$sf

  return(data.frame(Group.1 = wcc_out$tp, x = pcMut))
}

# Adapted from Gerstung et al, 2020 [https://gerstung-lab.github.io/PCAWG-11/#10_chronological_wgd__mrca]
# Needs to incorporate the correction by power
get_branches_dense = function(vcf_mtr_fp, sbs_a, purdir, cvg, power_correction = TRUE, target_sigs = c("SBS1", "SBS5", "SBS40a", "SBS40b")){
    id = extract_ppcg_id(vcf_mtr_fp, full = F)
    vcf = parse_vcf_timing(vcf_mtr_fp, genome = NULL, add_tnc = FALSE, add_sbs = TRUE, sbs_a = sbs_a)
    w = is_clock(vcf, target_sigs)
    if (sum(w) == 0){return(c(0,0))}
    p <- info(vcf)$pSub[w]
    n.subclonal <- aggregate(p, list(info(vcf)$CNF[w]), sum, na.rm=TRUE)
    n.clonal <- aggregate(1-p, list(info(vcf)$CNF[w]), sum, na.rm=TRUE)

    # Power correction using WCC
    if (power_correction){
      n.subclonal <- power_corr(n.subclonal, id, purdir, cvg)
      n.clonal <- power_corr(n.clonal, id, purdir, cvg)
    }

    b.subclonal = n.subclonal$x %*% (n.subclonal$Group.1) / max(n.subclonal$Group.1)
    b.clonal = sum(n.clonal$x) # Clonal branch (trunk), power adjusted

    # Correction by time-averaged ploidy back to diploid in clonal period
    # And by observed ploidy the subclonal period
    effGenome <- 2/avgWeights(vcf[na.omit(w)])
    b.subclonal <- b.subclonal * 2/get_ploidy(purdir, id)
    b.clonal <- b.clonal * 2/effGenome

    return(c(b.subclonal, b.clonal))
}

get_branches_linear = function(vcf_mtr_fp, sbs_a, purdir, cvg, power_correction = TRUE, target_sigs = c("SBS1", "SBS5", "SBS40a", "SBS40b")){
    id = extract_ppcg_id(vcf_mtr_fp, full = F)
    vcf = parse_vcf_timing(vcf_mtr_fp, genome = NULL, add_tnc = FALSE, add_sbs = TRUE, sbs_a = sbs_a)
    w = is_clock(vcf, target_sigs)
    if (sum(w) == 0){return(c(0,0))}
    p <- info(vcf)$pSub[w]
    n.subclonal <- aggregate(p, list(info(vcf)$CNF[w]), sum, na.rm=TRUE)
    n.clonal <- aggregate(1-p, list(info(vcf)$CNF[w]), sum, na.rm=TRUE)
    # Power correction using WCC
    if (power_correction){
      n.subclonal <- power_corr(n.subclonal, id, purdir, cvg)
      n.clonal <- power_corr(n.clonal, id, purdir, cvg)
    }
    b.subclonal <- sum(n.subclonal$x) # Subclonal branch, power adjusted
    b.clonal <- sum(n.clonal$x) # Clonal branch (trunk), power adjusted
    # Correction by time-averaged ploidy back to diploid in clonal period
    # And by observed ploidy the subclonal period
    effGenome <- 2/avgWeights(vcf[na.omit(w)])
    
    b.subclonal <- b.subclonal * 2/get_ploidy(purdir, id)
    b.clonal <- b.clonal * 2/effGenome
    return(c(b.subclonal, b.clonal))
}

# from Gerstung et al, 2020 [https://gerstung-lab.github.io/PCAWG-11/#10_chronological_wgd__mrca]
cnWeights <- function(vcf){
    t <- if(is.null(info(vcf)$TCN)) (info(vcf)$MajCN + info(vcf)$MinCN) else info(vcf)$TCN
    info(vcf)$MutCN / t
}

# from Gerstung et al, 2020 [https://gerstung-lab.github.io/PCAWG-11/#10_chronological_wgd__mrca]
# It computes the time-averaged effective genome size
avgWeights <- function(vcf, w){
    cnw <- cnWeights(vcf)
    mean(2*cnw[w], na.rm=TRUE)
}

# from computeWgdParamDeam function
# adapted to output only the number of mutations pre and post-WGD in regions 2+0, 2+1 and 2+2
time_clock_wgd <- function(vcf_fs, bb_fs, sbs_a, purdir, ccfdir, id, target_sigs = c("SBS1", "SBS5", "SBS40a", "SBS40b")){
    # load all information from the patient
    vcf_mtr_fp = vcf_fs[str_detect(vcf_fs, id)]
    bb_fp = bb_fs[str_detect(bb_fs, id)]
    bb = readRDS(bb_fp)
    vcf = parse_vcf_timing(vcf_mtr_fp, genome = NULL, add_tnc = FALSE, add_sbs = TRUE, sbs_a = sbs_a)
    purity = get_purity(purdir, id)
    clusters = load_ccf(ccfdir, purity, id)

    # 1. Find segments compatible with WGD
    min.dist <- 0.05
    m <- findMainCluster(bb)
    l <- pmin(bb$time.lo, bb$time - min.dist)
    u <- pmax(bb$time.up, bb$time + min.dist)
    o <- which(l <= m & u >= m)

    # 2. Find clock-like signatures in compatible segments
    # 2+2 segments
    w <- which(info(vcf)$MajCN==2 & info(vcf)$MinCN==2 & sapply(info(vcf)$CNID, length)==1 & is_clock(vcf, target_sigs) & vcf %over% bb[o])
    v <- vcf[w]
    if(nrow(v)<=100) return(c(NA, NA)) # At least 100 SNVs
    
    ###### CHANGE FROM ORIGINAL FUNCTION BELOW
    # CpG > T mutations before and after WGD
    pre_wgd = sum(info(v)$CLS == "clonal [early]", na.rm = TRUE)
    post_wgd = sum(info(v)$CLS == "clonal [late]", na.rm = TRUE) * 2 / get_ploidy(purdir, id) # Correct for ploidy in post-wgd period

    return(c(pre_wgd, post_wgd))
}

mol2time = function(sample_ids, nmuts_before, nmuts_after, mutation_rate, accelerations = c(1,2.5,5,7.5,10,20), event = "MRCA", boot = 100){
    # Get total information from the donors
    total_muts = nmuts_before + nmuts_after
    ages = get_age(sample_ids)

    # Monte-Carlo approach to estimate onset of mutation acceleration in clonal branch
    # Sample initiation ages from uniform distribution between 25% and 5% age diagnosis
    fraction_ages = runif(boot, min = 0.75, max = 0.95)
    
    list_times = foreach(fa=fraction_ages) %dopar% {
      # 1. Calcuate initiation age for each sample
      init_ages = ages * fa     
      # 2. Calculate average mutation rate across all evolutionary period for each rate acceleration
      avg_rate = fa * 1 + (1-fa) * accelerations
      baselines_rates = matrix(1 / avg_rate, nrow = length(total_muts), ncol = length(accelerations), byrow = TRUE) * mutation_rate
      accelerated_rates = baselines_rates * matrix(accelerations, nrow = length(total_muts), ncol = length(accelerations), byrow = TRUE)
      
      
      if (event == "MRCA"){
        # 3. For MRCA, assume full acceleration in subclonal period. Because we only consider mutations AFTER MRCA to estimate the time, 
        # no need to re-scale the mutations occurring in the clonal period
        ts = matrix(nmuts_after, nrow = length(total_muts), ncol = length(accelerations), byrow = FALSE) / accelerated_rates
        colnames(ts) = paste0(accelerations, "x")
        rownames(ts) = sample_ids
        return(ts)
      } else if (event == "WGD"){
        # 4. For WGD, assume acceleration starts at some point in clonal period, and calculate number of mutations before and after acceleration accordingly
        nmuts_before_acceleration = baselines_rates * init_ages
        nmuts_before_acceleration = pmin(nmuts_before_acceleration, matrix(total_muts, nrow = length(total_muts), ncol = length(accelerations), byrow = FALSE), na.rm = TRUE)

        # Divide mutations between the ones that occur before and after acceleration
        nmuts_before = matrix(nmuts_before, nrow = length(total_muts), ncol = length(accelerations), byrow = FALSE)
        nmuts_after = matrix(nmuts_after, nrow = length(total_muts), ncol = length(accelerations), byrow = FALSE)

        pre_to_scale = pmax(matrix(0, nrow = length(total_muts), ncol = length(accelerations)), nmuts_before - nmuts_before_acceleration)
        pre_not_scale = nmuts_before - pre_to_scale
        post_not_scale = nmuts_before_acceleration - nmuts_before + pre_to_scale
        post_to_scale = pmax(matrix(0, nrow = length(total_muts), ncol = length(accelerations)), total_muts - (pre_to_scale + pre_not_scale + post_not_scale))

        stopifnot(all(pre_to_scale + pre_not_scale + post_not_scale + post_to_scale - total_muts < 0.001, na.rm = T))
        ts = post_not_scale / baselines_rates + post_to_scale / accelerated_rates
        colnames(ts) = paste0(accelerations, "x")
        rownames(ts) = sample_ids
        return(ts)
      } else {
        stop("Event must be either MRCA or WGD")
      }
    }
}


mol2time_patient_rate <- function(sample_ids, nmuts_before, nmuts_after, nmuts_subclonal = 0, accelerations = c(1,2.5,5,7.5,10,20), boot = 100, event = "MRCA"){
  
  # nmuts_subclonal is used to enable correct estimation of per-patient mutation rate in WGD case, 
  # where nmuts_before and nmuts_after are the total number of mutations before and after WGD
  # when estimating MRCA, nmuts_subclonal is set to 0, as subclonal muts are passed as nmuts_after

  # Get total information from the donors
  total_muts = nmuts_before + nmuts_after
  ages = get_age(sample_ids)
  # Create patient-specific mutation rates in matrix
  mutation_rate = (total_muts + nmuts_subclonal) / ages
  mutation_rate = matrix(mutation_rate, nrow = length(total_muts), ncol = length(accelerations), byrow = FALSE)

  # Monte-Carlo approach to estimate onset of mutation acceleration in clonal branch
  # Sample initiation ages from uniform distribution between 25% and 5% age diagnosis
  fraction_ages = runif(boot, min = 0.75, max = 0.95)

  foreach(fa=fraction_ages) %dopar%{
    # 1. Calcuate initiation age for each sample
    init_ages = ages * fa     
    # 2. Calculate average mutation rate across all evolutionary period for each rate acceleration
    avg_rate = fa * 1 + (1-fa) * accelerations
    baselines_rates = matrix(1 / avg_rate, nrow = length(total_muts), ncol = length(accelerations), byrow = TRUE) * mutation_rate
    accelerated_rates = baselines_rates * matrix(accelerations, nrow = length(total_muts), ncol = length(accelerations), byrow = TRUE)
    
    
    if (event == "MRCA"){
      # 3. For MRCA, assume full acceleration in subclonal period
      ts = matrix(nmuts_after, nrow = length(total_muts), ncol = length(accelerations), byrow = FALSE) / accelerated_rates
      colnames(ts) = paste0(accelerations, "x")
      rownames(ts) = sample_ids
      return(ts)
    } else if (event == "WGD"){
      # 4. For WGD, assume acceleration starts at some point in clonal period, and calculate number of mutations before and after acceleration accordingly
      nmuts_before_acceleration = baselines_rates * init_ages
      nmuts_before_acceleration = pmin(nmuts_before_acceleration, matrix(total_muts, nrow = length(total_muts), ncol = length(accelerations), byrow = FALSE), na.rm = TRUE)

      # Divide mutations between the ones that occur before and after acceleration
      nmuts_before = matrix(nmuts_before, nrow = length(total_muts), ncol = length(accelerations), byrow = FALSE)
      nmuts_after = matrix(nmuts_after, nrow = length(total_muts), ncol = length(accelerations), byrow = FALSE)

      pre_to_scale = pmax(matrix(0, nrow = length(total_muts), ncol = length(accelerations)), nmuts_before - nmuts_before_acceleration)
      pre_not_scale = nmuts_before - pre_to_scale
      post_not_scale = nmuts_before_acceleration - nmuts_before + pre_to_scale
      post_to_scale = pmax(matrix(0, nrow = length(total_muts), ncol = length(accelerations)), total_muts - (pre_to_scale + pre_not_scale + post_not_scale))

      stopifnot(all(pre_to_scale + pre_not_scale + post_not_scale + post_to_scale - total_muts < 0.001, na.rm = T))
      ts = post_not_scale / baselines_rates + post_to_scale / accelerated_rates
      colnames(ts) = paste0(accelerations, "x")
      rownames(ts) = sample_ids
      return(ts)
    } else {
      stop("Event must be either MRCA or WGD")
    }
  }
}

# Get molecular time of CN alterations
get_mtime <- function(type, clonal_early, clonal_late, clonal_na){
    m2 = clonal_early
    m1 = clonal_late + clonal_na
    # According to methods [here](https://www.nature.com/articles/s41586-019-1907-7#Sec9)
    # NB: Possible to get values over 1 with this code
    if (type == "Mono-allelic Gain"){
      return((3*m2) / (2*m2 + m1))
    } else if (type == "CN-LOH"){
      return(((2*m2) / (2*m2 + m1)))
    } else if (type == "Bi-allelic Gain (WGD)"){
      return(((2*m2) / (2*m2 + m1)))
    } else (
      stop("Please provide valid type of event")
    )
}

# Get mutations clonal early and late in CN segments for CN timing
pull_cn_early_late <- function(smp, bb_fs, vcf_fs, sbs_a, clock_sigs = c("SBS1", "SBS5", "SBS40a", "SBS40b")){    
    # just to pull quickly the mutations clonal early, late...
    muts_by_period <- function(label, mutation_tb, clock_sigs){
        if (!label %in% rownames(mutation_tb)){
            return(0)
        } else {
            return(sum(
                mutation_tb[
                    rownames(mutation_tb) == label, 
                    colnames(mutation_tb) %in% clock_sigs
                    ], na.rm = TRUE 
                )
            )
        }
    }
    return_na <- function(){
      data.frame(
                sample = NA, 
                type = NA,
                chr = NA, 
                start = NA, 
                end = NA, 
                major_cn = NA, 
                minor_cn = NA,
                clock_early = NA, 
                clock_late = NA, 
                clock_na = NA, 
                clock_subclonal = NA, 
                mtime_cn = NA
            )
    }

    # Load data, adding SBS assignments to vcf
    if (!any(grepl(smp, bb_fs)) | !any(grepl(smp, vcf_fs))){return(return_na())}
    bb <- readRDS(bb_fs[grepl(smp, bb_fs)])
    bb <- bb[!is.na(bb$type)]
    vcf_fp <- vcf_fs[grepl(smp, vcf_fs)]
    vcf <- parse_vcf_timing(vcf_fp, genome = NULL, add_tnc = FALSE, sbs_a = sbs_a)
    # Find how many mutations overlap with CN info
    ovs <- findOverlaps(bb, vcf)
    if (length(ovs) == 0){
      return(return_na())
    }

    o <- as.data.frame(rbindlist(
        lapply(1:length(bb), function(i){
            muts <- vcf[ovs@to[ovs@from == i]]  
            chr <- get_chr_from_granges(bb[i])
            start <- start(bb[i]) 
            end <- end(bb[i])
            type <- bb$type[i]
            major_cn <- bb$major_cn[i]
            minor_cn <- bb$minor_cn[i]
            mutation_tb <- table(muts@info$CLS, muts@info$sbs_sig)
            clock_early <- muts_by_period("clonal [early]", mutation_tb, clock_sigs)
            clock_late <- muts_by_period("clonal [late]", mutation_tb, clock_sigs)
            clock_na <- muts_by_period("clonal [NA]", mutation_tb, clock_sigs)
            clock_subclonal <- muts_by_period("subclonal", mutation_tb, clock_sigs)
            mtime_cn <- get_mtime(type, clock_early, clock_late, clock_na)
            return(data.frame(
                sample = smp, 
                type = type,
                chr = chr, 
                start = start, 
                end = end, 
                major_cn = major_cn, 
                minor_cn = minor_cn,
                clock_early = clock_early, 
                clock_late = clock_late, 
                clock_na = clock_na, 
                clock_subclonal = clock_subclonal,
                mtime_cn = mtime_cn
            ))
        })
    ))
    return(o)
}


# Function to count mutations in different adjacent segments
merge_segments <- function(time_segments, distance_threshold = 1e6){
    # quick function to append to output and avoid repetition downstream
    bind_to_output <- function(o, sample, chr, start, end, clock_early, clock_late, clock_na, clock_subclonal, rt_mrca_global_rate, rt_mrca_pt_rate, type, mtime_cn){
        return(rbind(
            o, data.frame(
                sample = sample, chr = chr, start = start, end = end, 
                clock_early = clock_early, clock_late = clock_late, 
                clock_na = clock_na, clock_subclonal = clock_subclonal, rt_mrca_global_rate = rt_mrca_global_rate, 
                rt_mrca_pt_rate = rt_mrca_pt_rate, type = type, mtime_cn = mtime_cn,
                stringsAsFactors = FALSE
            )
        ))
    }

    # function to add mutations and change boundaries of segment after merge
    merge_segment <- function(end, clock_early, clock_late, clock_na, clock_subclonal, type, mtime_cn, o){
        # end is now the end of the merged boundary
        o$end[nrow(o)] <- end

        # we take the weighted mean of the molecular time of segments
        # to merge, weighted by total number of mutations
        if (!is.na(mtime_cn)){
            o$mtime_cn[nrow(o)] <- weighted.mean(
            c(o$mtime_cn[nrow(o)], mtime_cn), 
            # weighted by mutations
            c(
                sum(c(o$clock_early[nrow(o)], o$clock_late[nrow(o)], o$clock_na[nrow(o)])), 
                sum(clock_early, clock_late, clock_na)
                ), 
            na.rm = T
            )
            
        }
        # mutations have to be added to the prior ones
        o$clock_early[nrow(o)] <- o$clock_early[nrow(o)] + clock_early
        o$clock_na[nrow(o)] <- o$clock_na[nrow(o)] + clock_na
        o$clock_late[nrow(o)] <- o$clock_late[nrow(o)] + clock_late
        o$clock_subclonal[nrow(o)] <- o$clock_subclonal[nrow(o)] + clock_subclonal

        # add information of the distinct types of gains found
        if (type != o$type[nrow(o)] & !is.na(type)){
            o$type[nrow(o)] <- paste0(o$type[nrow(o)], "|", type)
        }
        return(o)
    }

    # where we will get merged segments
    o <- data.frame(
        sample = c(), chr = c(), start = c(), end = c(), 
        clock_early = c(), clock_late = c(), clock_na = c(), clock_subclonal = c(), 
        rt_mrca_global_rate = c(), rt_mrca_pt_rate = c(), 
        type = c()
        )

    for (i in 1:nrow(time_segments)){
        print(i)
        # Check whether it should be merged to previous segment
        sample = time_segments$sample[i]; chr = time_segments$chr[i]
        start = time_segments$start[i]; end = time_segments$end[i]
        clock_early = time_segments$clock_early[i]
        clock_late = time_segments$clock_late[i]
        clock_na = time_segments$clock_na[i]
        clock_subclonal = time_segments$clock_subclonal[i]
        rt_mrca_global_rate = time_segments$rt_mrca_global_rate[i]
        rt_mrca_pt_rate = time_segments$mrca_rt_pt_rate[i]
        type = as.character(time_segments$type[i])
        mtime_cn = time_segments$mtime_cn[i]
        if (i == 1){
            o <- bind_to_output(o, sample, chr, start, end, clock_early, clock_late, clock_na, clock_subclonal, rt_mrca_global_rate, rt_mrca_pt_rate, type, mtime_cn)
            next()
        }
        merge <- sample == o$sample[nrow(o)] & chr == o$chr[nrow(o)] & start <= (o$end[nrow(o)] + distance_threshold)
        if (merge){
            o <- merge_segment(end, clock_early, clock_late, clock_na, clock_subclonal, type, mtime_cn, o)
        } else {
            o <- bind_to_output(o, sample, chr, start, end, clock_early, clock_late, clock_na, clock_subclonal, rt_mrca_global_rate, rt_mrca_pt_rate, type, mtime_cn)
        }
    }
    return(o)
}

### HELPER FUNCTIONS FROM GERSTUNG ET AL, 2020
findMainCluster <- function(bb, min.dist=0.05){
    w <- which(bb$n.snv_mnv > 20 & !is.na(bb$time))
#   d <- dist(bb$time[w])
#   ci <- weighted.mean((bb$time.up - bb$time.lo)[w], width(bb)[w])
#   h <- hclust(d, method='average', members=bb$n.snv_mnv[w])
#   c <- cutree(h, h=ci)
#   ww <- c==which.max(table(c))
#   weighted.mean(bb$time[w][ww], 1/((bb$time.up - bb$time.lo + min.dist)[w][ww]), na.rm=TRUE)
    s <- seq(0,1,0.01)
    l2 <- pmin(bb$time.lo, bb$time - min.dist)[w]
    u2 <- pmax(bb$time.up, bb$time + min.dist)[w]
    l1 <- (l2 +  bb$time[w])/2
    u1 <- (u2+  bb$time[w])/2
    wd <- as.numeric(width(bb)[w])
    o <- sapply(s, function(i) sum(wd * ( (l2 <= i & u2 >=i) + (l1 <= i & u1 >= i))))
    s[which.max(o)]
}

piToTime <- function(timing_param, type=c("Mono-allelic Gain","CN-LOH", "Bi-allelic Gain (WGD)")){
    type <- match.arg(type)
    evo <- NA
    w <- timing_param[,"s"]==1 &! timing_param[,"mixFlag"] ## Clonal states, not mixture (CN subclonal)
    M <- sum(w) ## Max multiplicity = Major copy number
    m <- sum(w & timing_param[,"n.m.s"]==2) ## Minor CN
    #evo <- paste0("1:1", paste0(2,":",m), if(M>2) paste0(3,":",) else "", collapse="->")
    t <- timing_param[(M:(M-1)),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE] ## Timing M and M-1
    if(nrow(t)==1) t <- rbind(t, NA)
    if(!any(is.na(t))){
        if(M==4) {
            if(timing_param[M-1,"P.m.sX"] < 0.02){ ## No MCN 3
                if(type=="CN-LOH"){ ## Hotfix for 4+0, treated at 1:1 -> 2:0 -> 4:0
                    t <- timing_param[c(M,M-2),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE]*c(1,0.5) ## Timing M and M-2
                    evo <- "1:1->2:0->4:0"
                }
                else if(type=="Bi-allelic Gain (WGD)"){         
                    if(m==2) {## Hotfix for 4+2 regions, treated at 1:1 -> 2:1 -> 4:2
                        t <- timing_param[c(M,M-2),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE] ## Timing M and M-2
                        t[2,] <- pmax(0,2*t[2,] - t[1,])/3
                        evo <- "1:1->2:1->4:2"
                    } else if(m==4) {## Hotfix for 4+4 regions, treated at 1:1 -> 2:2 -> 4:4
                        t <- timing_param[c(M,M-2),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE]*c(1,0.5) ## Timing M and M-2
                        evo <- "1:1->2:2->4:4"
                    }
                }           
            } else if(type=="Bi-allelic Gain (WGD)"){ ## Can't uniquely time second event
                t[2,] <- NA
            } 
            if(m==3) {
                t[2,] <- NA ## Don'time secondary 4+3 for now, needs more work
            } 
        } else {
            if(M==3 & type=="Bi-allelic Gain (WGD)") {## Hotfix for 3+2 regions, treated at 1:1 -> 2:2 -> 3:2
                t[2,] <- pmax(0,2*t[2,] - t[1,])
                evo <- "1:1->2:2->3:2"
            }
        }
    }
    colnames(t) <- c("","lo","up")
    t[,2] <- pmin(t[,1],t[,2])
    t[,3] <- pmax(t[,1],t[,3])
    t <- pmin(apply(t,2,cumsum),1) ## Times are actually deltas
    if(M < 3) t[2,] <- NA
    t[is.infinite(t)] <- NA
    rownames(t) <- c("",".2nd")
    return(c(t[1,],`.2nd`=t[2,])) ## Linearise
}

bbToTime <- function(bb, timing_param = bb$timing_param, pseudo.count=5){
    sub <- duplicated(bb) 
    covrg <- countQueryHits(findOverlaps(bb, bb)) 
    maj <- sapply(timing_param, function(x) if(length(x) > 0) x[1, "majCNanc"] else NA) #bb$major_cn
    min <- sapply(timing_param, function(x) if(length(x) > 0) x[1, "minCNanc"] else NA) #bb$minor_cn
    type <- sapply(seq_along(bb), function(i){
                if(maj[i] < 2 | is.na(maj[i]) | sub[i] | (maj[i] > 4 & min[i] >= 2)) return(NA)
                type <- if(min[i]==1){ "Mono-allelic Gain" 
                        }else if(min[i]==0){"CN-LOH"}
                        else "Bi-allelic Gain (WGD)"
                return(type)
            })
    time <- t(sapply(seq_along(bb), function(i){
                        if(sub[i] | is.na(type[i])) return(rep(NA,6)) 
                        else piToTime(timing_param[[i]],type[i])
                    }))
    
    res <- data.frame(type=factor(type, levels=c("Mono-allelic Gain","CN-LOH","Bi-allelic Gain (WGD)")), time=time)
    colnames(res) <- c("type","time","time.lo","time.up","time.2nd","time.2nd.lo","time.2nd.up")
    
    # posthoc adjustment of CI's
    res$time.up <- (pseudo.count + bb$n.snv_mnv * res$time.up)/(pseudo.count + bb$n.snv_mnv)
    res$time.lo <- (0 + bb$n.snv_mnv * res$time.lo)/(pseudo.count + bb$n.snv_mnv)
    res$time.2nd.up <- (pseudo.count + bb$n.snv_mnv * res$time.2nd.up)/(pseudo.count + bb$n.snv_mnv)
    res$time.2nd.lo <- (0 + bb$n.snv_mnv * res$time.2nd.lo)/(pseudo.count + bb$n.snv_mnv)
    
    res$time.star <- factor((covrg == 1) + (min < 2 & maj <= 2 | min==2 & maj==2) * (covrg == 1), levels=0:2, labels = c("*","**","***")) ## ***: 2+0, 2+1, 2+2; **: 2<n+1, {3,4}+2; *: subclonal gains
    res$time.star[is.na(res$time)] <- NA
    return(res)
}