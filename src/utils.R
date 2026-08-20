# Extract a PPCG sample ID
extract_ppcg_id <- function(x, full=T){
  if(full){
    return(stringr::str_extract(x, "PPCG[0-9a-zA-Z]+_DNA_vs_PPCG[0-9a-zA-Z]+_DNA"))
  } else {
    return(stringr::str_extract(x, "PPCG[0-9a-zA-Z]+_DNA"))
  }
}

# Extract a PPCG patient ID
extract_ppcg_pt <- function(x){
  return(stringr::str_extract(x, "PPCG[0-9]+"))
}

# Transform a PPCG sample ID to the corresponding patient ID
smp2pt <- function(sample_id){
  return(stringr::str_extract(sample_id, "PPCG[0-9]+"))
}

get_mmrd <- function(sbs_a = read_delim("data/processed/SNV_signature_annotations__2024-06-21/SNVs_all_signature_20240621.tsv"), mmrd_sigs = c("SBS15", "SBS44", "SBS21")){
  is_mmrd = extract_ppcg_pt(sbs_a$Sample_ID[rowSums(sbs_a[mmrd_sigs]) > 0])
  return(unique(is_mmrd))
}

# Get Gleason Grade Group
get_gs_group <- function(ids, meta_fp = "data/meta/PPCG_donors_clin_20241217.csv", asnum = F, collapse_gs4 = T, collapse_gs2 = F){
  meta = read_delim(meta_fp, delim = ",")
  ids = smp2pt(ids)
  m = match(ids, meta$record_id)
  grade_groups = meta$gleason_grade_group[m]
  # Collapse if required
  if (collapse_gs4){
    grade_groups[grade_groups %in% c("4", "5")] = "4+"
  }
  if (collapse_gs2){
    grade_groups[grade_groups %in% c("1", "2")] = "2-"
  }
  if (asnum){
    grade_groups = ifelse(grade_groups == "4+", 4, ifelse(grade_groups == "2-", 2, as.numeric(grade_groups)))
  }
  return(grade_groups)
}

# pull out chromosomes from a GRanges object
get_chr_from_granges <- function(gr){
  rle <- seqnames(gr)
  return(rep(as.character(rle@values), rle@lengths))
}

get_nrpcc <- function(ids, nrpcc_fp = "data/qc/PPCG_Feb2026_NRPCC.txt"){
  qc = read_delim(nrpcc_fp)
  qc$PPCG_Sample_ID = extract_ppcg_id(qc$Barcode, full = FALSE)
  if (!grepl("_DNA", ids[1])){ids = paste0(ids, "_DNA")}
  return(sapply(ids, function(id){qc$NRPCC[qc$PPCG_Sample_ID == id]}))
}

get_wgd <- function(qc = TRUE, wgd_fp = "data/processed/wgd/PPCG_WGD_status_PCAWG_method_20260212.txt", rm_dna = F){
  wgd_log = read.table(wgd_fp, header = T)
  colnames(wgd_log) = c("Barcode", "wgd_status")
  wgd_samps = wgd_log$Barcode[wgd_log$wgd_status == "WGD"]
  wgd_samps = extract_ppcg_id(wgd_samps, full = F)
  if (qc){
    wgd_samps = dplyr::intersect(wgd_samps, qc_pass())
  }
  if (rm_dna){
    wgd_samps = str_remove(wgd_samps, "_DNA")
  }
  return(wgd_samps)
}

qc_pass <- function(qc_fp = "data/qc/qc_pass.csv"){
  qc = read_delim(qc_fp, delim = ",")
  return(qc$PPCG_Sample_ID[qc$ppcg_evo_sel_one])
}

get_age <- function(ids, meta_fp = "data/meta/PPCG_donors_clin_20241217.csv"){
  meta = read_delim(meta_fp, delim = ",")
  ids = smp2pt(ids)
  m = match(ids, meta$record_id)
  return(as.numeric(meta$age_at_tumour_collection[m]))
}

# Get all the primary sample ID in PP-Evo cohort
get_prim_smps <- function(qc = TRUE, rm_dna = F, anatomical_location_fp = "data/meta/Sample_Donor_Tissue_Origin_2024.csv"){
  anatomical_loc = read_delim(anatomical_location_fp, delim = ",")
  prim_smps = anatomical_loc$PPCG_Sample_ID[anatomical_loc$Tissue_Origin == "Primary"]
  if (qc){
    prim_smps = dplyr::intersect(prim_smps, qc_pass())
  }
  if (rm_dna){
    prim_smps = str_remove(prim_smps, "_DNA")
  }
  return(prim_smps)
}

# Get all the metastatic samples ID in PP-Evo cohort
get_met_smps <- function(qc = TRUE, rm_dna = F, anatomical_loc_fp = "data/meta/Sample_Donor_Tissue_Origin_2024.csv"){
  anatomical_loc <- read_delim(anatomical_loc_fp, delim = ",")
  met_smps = anatomical_loc$PPCG_Sample_ID[anatomical_loc$Tissue_Origin == "Metastasis"]
  if (qc){
    met_smps = dplyr::intersect(met_smps, qc_pass())
  }
  if (rm_dna){
    met_smps = str_remove(met_smps, "_DNA")
  }
  return(met_smps)
}

# Get all the primary sample ID with corresponding met in PP-Evo
get_prim_mets <- function(qc = TRUE, rm_dna = F, meta_fp = "data/meta/PPCG_donors_clin_20241217.csv"){
  prim_smps = get_prim_smps(qc, rm_dna)
  clin_meta = read_delim(meta_fp, delim = ",")
  is_met_pt = extract_ppcg_pt(prim_smps) %in% clin_meta$record_id[clin_meta$mets_ind == "mets"]
  return(prim_smps[is_met_pt])
}

# Get all the primary sample ID without corresponding met in PP-Evo 
get_prim_nmets <- function(qc = TRUE, rm_dna = F, followup_filter=0, meta_fp = "data/meta/PPCG_donors_clin_20241217.csv"){
  prim_smps = get_prim_smps(qc, rm_dna)
  clin_meta = read_delim(meta_fp, delim = ",")
  is_met_pt = extract_ppcg_pt(prim_smps) %in% clin_meta$record_id[clin_meta$mets_ind == "mets"]
  if (followup_filter == 0){
    return(prim_smps[!is_met_pt])
  } else {
    is_long_followup = extract_ppcg_pt(prim_smps) %in% clin_meta$record_id[clin_meta$donor_interval_of_last_followup > 365.25*followup_filter]
    return(prim_smps[!is_met_pt & is_long_followup])
  }
}

read_purity <- function(pur_fp){
  purplo = read.table(pur_fp, header=T)
  purity_column = ifelse("purity" %in% colnames(purplo), "purity", "cellularity")
  return(purplo[[purity_column]])
}

read_ploidy <- function(pur_fp){
  purplo = read.table(pur_fp, header=T)
  return(purplo[["ploidy"]])
}

get_purity <- function(purdir, id){
  cc = fs::dir_ls(purdir, glob="*.txt") %>% 
    purrr::keep(~stringr::str_detect(.x, id))
    
  if(length(cc) == 0){
    print("No purity file found")
    return(NA)
  }
  return(read_purity(cc[[1]]))
}


get_ploidy <- function(purdir, id){
  cc = fs::dir_ls(purdir, glob="*.txt") %>% 
    purrr::keep(~stringr::str_detect(.x, id))
    
  if(length(cc) == 0){
    print("No ploidy file found")
    return(NA)
  }

  read_ploidy(cc[[1]])
}

load_trajectories <- function(trajectories_fps, parse_ppcg_id = T) {
    tj1 <- read.table(trajectories_fps[1], header = T) %>%
        {
            if (parse_ppcg_id) {
                extract_ppcg_id(unique(.$Tumour_Name), full = F)
            } else {
                unique(.$Tumour_Name)
            }
        }
    tj2 <- read.table(trajectories_fps[2], header = T) %>%
        {
            if (parse_ppcg_id) {
                extract_ppcg_id(unique(.$Tumour_Name), full = F)
            } else {
                unique(.$Tumour_Name)
            }
        }
    tj3 <- read.table(trajectories_fps[3], header = T) %>%
        {
            if (parse_ppcg_id) {
                extract_ppcg_id(unique(.$Tumour_Name), full = F)
            } else {
                unique(.$Tumour_Name)
            }
        }
    return(
        data.frame(
            sample = c(tj1, tj2, tj3),
            trajectory = c(
                rep("Ordering 1", length(tj1)),
                rep("Ordering 2", length(tj2)),
                rep("Ordering 3", length(tj3))
            )
        )
    )
}

boot_median <- function(x, n, na.rm = TRUE){
  boot_samples <- replicate(n, {
    sample_data <- sample(x, replace = TRUE)  # Resample with replacement
    median(sample_data, na.rm = na.rm)  # Calculate median
  })
  return(boot_samples)  # Return the vector of bootstrapped medians
}

get_relapse_data <- function(clinfile, samples, type = "BCR", only_complete = FALSE) {
    if (type == "BCR"){
      clin_data <- read_delim(clinfile) %>% format_bcr2surv()
    } else if (type == "metastasis") {
      clin_data <- read_delim(clinfile) %>% format_met2surv()
    } else {
      stop("Please provide valid clinical endpoint (i.e. metastasis or BCR)")
    }
    samples <- extract_ppcg_pt(samples)
    samples <- samples[samples %in% clin_data$record_id]
    m <- match(samples, clin_data$record_id)
    if (type == "BCR"){
        relapse <- data.frame(
            record_id = samples,
            relapse_ind = clin_data$relapse_ind[m],
            time2relapse = as.numeric(clin_data$donor_relapse_interval[m])
        )
    } else if (type == "metastasis"){
        relapse <- data.frame(
            record_id = samples, 
            relapse_ind = clin_data$mets_ind[m], 
            time2relapse = as.numeric(clin_data$donor_first_mets_interval_from_m0[m])
        )
    } 
    # Add additional clinical information
    relapse$grade_group <- clin_data$gleason_grade_group[m]
    relapse$psa <- clin_data$psa_at_tumour_collection[m]
    relapse$t_stage <- clin_data$path_t_stage[m]
    relapse$damico_risk <- clin_data$damico_risk[m]
    relapse$age <- clin_data$age_at_tumour_collection[m]
    if (only_complete){
      relapse = relapse[complete.cases(relapse),]
    }
    return(relapse)
}

# Format PPCG data for BCR relapse survival analysis
#' @param clin_data A dataframe containing variables
#' relapse_ind and donor_relapse_interval as formatted in PPCG
#' returns data with relapse_ind and donor_relapse_interval with correct format
format_bcr2surv <- function(clin_data){
  clin_data$relapse_ind = relapse2surv(clin_data$relapse_ind)
  clin_data$donor_relapse_interval = as.numeric(clin_data$donor_relapse_interval)
  # Note: "missing" coerced to NA; desired behaviour
  return(clin_data)
}

# Format PPCG data for metastasis survival analysis
#' @param clin_data A dataframe containing variables
#' mets_ind and donor_first_mets_interval_from_m0 as formatted in PPCG
#' returns data with mets_ind and donor_first_mets_interval_from_m0 with correct format
format_met2surv <- function(clin_data){
  clin_data$mets_ind = met2surv(clin_data$mets_ind)
  clin_data$donor_first_mets_interval_from_m0 = as.numeric(clin_data$donor_first_mets_interval_from_m0)
  # Note: "missing" coerced to NA; desired behaviour
  return(clin_data)
}

# Format PPCG relapse data to binary labels for survival analysis
relapse2surv <- function(relapse_ind){
  # NB: coercion of any other values to NA is intended
  lb = ifelse(relapse_ind == "no relapse", 0, ifelse(relapse_ind == "relapsed", 1, NA))
  return(lb)
}

# Format PPCG metastasis data to binary labels for survival analysis
met2surv <- function(mets_ind){
  # NB: coercion of any other values to NA is intended
  lb = ifelse(mets_ind == "no mets", 0, ifelse(mets_ind == "mets", 1, NA))
  return(lb)
}

# Get country where sequencing was performed
get_country <- function(ids, clinfile = "data/meta/PPCG_donors_clin_20241217.csv"){
  clin_data = read_delim(clinfile, delim = ",")
  ids = smp2pt(ids)
  m = match(ids, clin_data$record_id)
  group = clin_data$group[m]
  country = ifelse(str_detect(group, "^CA"), "CA", group)
  return(country)
}

# Get Cambridge Prognostic Group from grade group, psa and t stage
# https://share.google/q0Hagq9xkzUuOeknS
get_cambridge_prognostic_groups <- function(grade_group, psa, t_stage){
  if (any(is.na(c(grade_group, psa, t_stage)))) {
    return(NA)
  }
  is_T1T2 <- str_detect(t_stage, "^T1|^T2")
  is_T3   <- str_detect(t_stage, "^T3")
  is_T4   <- str_detect(t_stage, "^T4")

  if (grade_group == 5 | is_T4 | (sum(c(grade_group == 4, psa > 20, is_T3)) >= 2)){
    return("CPG5")
  }
  if (grade_group == 4 | psa > 20 | is_T3){
      return("CPG4")
  } 
  if (is_T1T2 & (grade_group == 3 | (grade_group == 2 & psa >= 10 & psa <= 20))){
    return("CPG3")
  }
  if (is_T1T2 & (grade_group == 2 | (psa >= 10 & psa <= 20))){
    return("CPG2")
  }
  if (is_T1T2 & grade_group == 1 & psa < 10){
    return("CPG1")
  } 
  return("Unclassified")  
}
