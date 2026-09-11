library(here)
library(fs)
library(VariantAnnotation)
library(foreach)
library(MutationTimeR)
library(mg14) # Required for plotting
library(RColorBrewer)
library(GenomicRanges)
library(IRanges)
library(logging)
library(glue)
library(dplyr)
library(data.table)
basicConfig()



# Main ----


pipeline_muttime <- function(snvdir, indeldir, ccfdir, purdir, cnadir, wgdfile, outdir, ncores=10){
  # snvdir: directory with a single file for every sample with SNV calls
  # indeldir: directory with a single file for every sample with indel calls
  # ccfdir: directory with a single file for every sample with dpclust calls
  # cnadir: directory with a single file for every sample with CNA calls

  # Note. SNVs and INDELs timed separately as in Nat 2020 paper and difficulty merging VCFs.

  svcfs <- fs::dir_ls(snvdir, glob="*.vcf.gz")
  
  doParallel::registerDoParallel(cores=ncores)
  foreach(sv=svcfs) %dopar% {
    id <- fs::path_file(sv) %>% stringr::str_extract("PPCG[0-9A-Za-d]+.*PPCG[0-9A-Z-a-d]+_DNA")

    # Run Mutation TimeR on SNVs (w/ skip)
    # First, parse SNVs and time
    vcf <- parse_vcf(sv)
    purity <- get_purity(purdir, id)
    ccfs <- load_ccf(ccfdir, purity, id)
    cnas <- load_cnas(cnadir, purity, id)
    cnas_clonal <- load_cnas(cnadir, purity, id, clonal=T) # For WGD
    is.wgd <- call_wgd_gerstung(cnas_clonal)

    missing_data <- any(
      c( 
      !("GRanges" %in% class(cnas)), 
      !("data.frame" %in% class(ccfs)),
      is.na(is.wgd),
      is.na(purity)
      ))
      
    
    if(missing_data){
      loginfo(glue("Skipping {id} due to missing data"))
      return(NA)
    }

    # Run MutationTimeR on SNVs for given sample
    loginfo(glue("Timing SNVs for {id}"))
    snv_outdir <- fs::dir_create(fs::path(outdir, "snvs"))
    mt <- run_ppcg_mt(id, vcf, cnas, ccfs, is.wgd, snv_outdir, gender="male")

    # Run MutatinoTImer on INDELs (w/skip)
    tumid <- stringr::str_extract(id, "PPCG[0-9A-Za-d]+")
    iv <- fs::dir_ls(indeldir, glob=glue("*{tumid}*.vcf.gz"))

    if(length(iv) == 0){
      loginfo(glue("Skipping {id} INDELs due to missing INDEL VCF"))
      return(NA)
    }

    loginfo(glue("Timing INDELs for {id}"))
    ivcf <- parse_indel_vcf(iv)
    indel_outdir <- fs::dir_create(fs::path(outdir, "indels"))
    mt <- run_ppcg_mt(id, ivcf, cnas, ccfs, is.wgd, indel_outdir, gender="male")

    loginfo(glue("Completed MutationTimeR for {id}"))
  }
  doParallel::stopImplicitCluster()

  loginfo("MutationTiming Complete")
}

call_wgd_gerstung <- function(bb){
  # Source: https://github.com/gerstung-lab/PCAWG-11/blob/6fc785dafbb54b00111899b40d93221785aff0a7/code/PCAWG-functions.R#L614C1-L614C70
  bb$total_cn = bb$major_cn + bb$minor_cn
  classWgd <- function(bb) .classWgd(averagePloidy(bb), averageHom(bb))
  .classWgd <- function(ploidy, hom) 2.9 -2*hom <= ploidy
  averageHom <- function(bb){
	sum(width(bb) * (bb$minor_cn == 0) * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
  }
  averagePloidy <- function(bb) {
	c <- if(!is.null(bb$copy_number)) bb$copy_number else bb$total_cn
	sum(width(bb) * c * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
  }
  if(!is.null(bb)){
    return(classWgd(bb))
  } else {
    return(NA)
  }
}

call_wgd <- function(wgdfile, id){

  # Uses old concesus WGD file
  wgd_call <- read.table(
    wgdfile,
    header=T,
  ) %>% dplyr::filter(
    PPCG_Sample_ID == stringr::str_remove(id, "_vs.*")
  ) %>% pull(WGD_2_out_of_3_consensus)

  if(length(wgd_call) == 0){
    return(NA)
  }

  return(wgd_call)
}

load_ccf <- function(ccfdir, purity, glob){
  # Defines the 'clusters' argument for mutation timer.
  # Must be a dataframe with columns: cluster, n_ssms, proportion
  # cluster is just a number identifier, n_ssms is the number of mutations in the cluster,
  # and proportion is purity*CCF. For example, see docstring:
  # https://github.com/gerstung-lab/MutationTimeR/blob/master/R/MutationTime.R#L460
  ccf = fs::dir_ls(ccfdir, glob="*bestClusterInfo.txt", recurse = TRUE) %>% 
    purrr::keep(~stringr::str_detect(.x, glob))

  if(length(ccf) == 0){
    print("No CCF files found")
    # Clonal treatment, as per docstring
    cf = data.frame(cluster=1, proportion=purity, n_ssms=100) 
    return(cf)
  }

  # Format for MutationTimeR
  cf = read.table(ccf, header=T)
  names(cf) = c("cluster", "proportion", "n_ssms")
  cf = dplyr::filter(cf, proportion < 1.5) %>% 
        dplyr::select(cluster, proportion, n_ssms) %>%
        dplyr::mutate(proportion = proportion*purity) %>% # Essential
        dplyr::arrange(-proportion)
  cf$cluster = seq(nrow(cf))

  return(cf)
}

load_cnas <- function(cnadir, purity, glob, clonal=F){
  bbf <- fs::dir_ls(cnadir, glob="*Subclonal_SCNA.txt") %>% 
    purrr::keep(~stringr::str_detect(.x, glob))

  if(length(bbf) == 0){
    print("No CNA files found")
    return(NA)
  }

  if(clonal){
    return(parse_cnas(bbf, purity))
  }

  return(parse_cnas_sc(bbf, purity))
}
# # TEST: Load SNV and check 
# v = "/nemo/project/proj-vanloo-secure/working/nmensah/projects/ppcg/data/raw/data_releases/WGS_Data_Release/Somatic_variants/SNVs/Filtered_SNV_VCFs_20_April_2021/PPCG0041a_DNA_vs_PPCG0041b_DNA.filtered.extended.pon.vcf.gz"
# # Write function to call DBS on SNVs, (temp fix until read-level calls provided).
# x = rowRanges(vcf)
# bychrom = split(x, seqnames(x))
# z = bychrom[[1]]
# DBS <- which(diff(start(z)) == 2)
# ## only on width==1?

# ## test indel
# vcf_file = "/nemo/project/proj-vanloo-secure/working/nmensah/projects/ppcg/data/raw/data_releases/WGS_Data_Release/Somatic_variants/InDels/PPCG_indels_16_Dec_2022/Tier_1/PPCG0004a_DNA_Tier_1.vcf.gz"
# vcf <- parse_indel_vcf(vcf)
# bb_file <- "/nemo/project/proj-vanloo-secure/working/nmensah/projects/ppcg/data/raw/data_releases/WGS_Data_Release/Somatic_variants/SCNA/Subclonal_SCNA_01_June_2020/PPCG0004a_DNA_vs_PPCG0004b_DNA_Subclonal_SCNA.txt"
# cnas <- parse_cnas(bb_file, purity =0.7)
# # Get purity and ploidy for a sample
# # Load local purity/ploidy file
# #get_pp = function(sample_id)p{
# #  ppfile= read.csv(here("data/ppcg_cna_cellularity_ploidy_200601.csv"), header = T, sep=",", colClasses = c("character", "numeric", "numeric"))
# #  ppfile[grepl(sample_id, ppfile$PPCG_Sample_ID),]}

parse_indel_vcf <- function(vcf_file){
    # Setup sample selection
  vcf_samples <- samples(VariantAnnotation::scanVcfHeader(vcf_file))
  tumor_samples <- vcf_samples[grepl("TUMOUR|PPCG.*", vcf_samples)]
  param=ScanVcfParam(samples=tumor_samples)

  vcf = readVcf(as.character(vcf_file), param=param)
  tt_header = DataFrame(Number=c(1,1), Type=c("Integer", "Integer"),
                      Description=c("Tumour alt count", "Tumour ref count"),
                      row.names= c("t_alt_count", "t_ref_count"))
  info(header(vcf)) <- rbind(info(header(vcf)), tt_header)
  info(vcf)$t_alt_count = as.integer(geno(vcf)$NV)
  info(vcf)$t_ref_count = as.integer(geno(vcf)$NR) - as.integer(geno(vcf)$NV)   

  # Check filter is PASS
  vcf <- vcf[rowRanges(vcf)$FILTER=="PASS"]

  return(vcf)
}

#' Parse PPCG VCF file for use with mutation timer
#'
#' @param vcf_file (character) Path to VCF file
#'
#' @return VariantAnnotation::VCF object with variant allele counts
parse_vcf <- function(vcf_file){
  # Setup sample selection
  vcf_samples <- samples(VariantAnnotation::scanVcfHeader(vcf_file))
  tumor_samples <- vcf_samples[grepl("TUMOUR|PPCG.*", vcf_samples)]
  param=ScanVcfParam(samples=tumor_samples)

  # Parse VCF
  param=ScanVcfParam(samples=tumor_samples)
  vcf = readVcf(as.character(vcf_file), param=param)
  tt_header = DataFrame(Number=c(1,1), Type=c("Integer", "Integer"),
                      Description=c("Tumour alt count", "Tumour ref count"),
                      row.names= c("t_alt_count", "t_ref_count"))
  info(header(vcf)) <- rbind(info(header(vcf)), tt_header)
  info(vcf)$t_alt_count = as.integer(geno(vcf)$DPA)
  info(vcf)$t_ref_count = as.integer(geno(vcf)$DPR)

  vcf <- vcf[rowRanges(vcf)$FILTER=="PASS"]

  return(vcf)
}

#' Return clone with highest frequncy in CNA segments from battenberg;
#' a helper for `parse_cnas`
#'
#' @param x (dataframe) A row from the battenberg output
#'
#' @return (integer) Major copy, minor copy and clonal fraction
best_bb_solution <- function(x){
  
  # Assign battenberg row columns to variables
  maj1a=x[["nMaj1_A"]]; min1a=x[["nMin1_A"]]; frac1a= x[["frac1_A"]]
  maj2a=x[["nMaj2_A"]]; min2a=x[["nMin2_A"]]; frac2a=x[["frac2_A"]]
  
  # Return first clone if no subclonal fraction called
  if(is.na(frac2a)){
    return(c(maj1a, min1a, frac1a))
  }
  
  # Else return clonal fraction and copy number
  if(frac1a>=frac2a){
    return(c(maj1a, min1a, frac1a))
  } else {
    return(c(maj2a, min2a, frac2a))
  }
}

#' Parse copy number segments from battenberg
#'
#' @param bb_file (character) Battenberg file with segment frequencies
#' @param purity (numeric) Overall sample purity inferred from battenberg
#'
#' @return (dataframe) Dataframe containing clonal copy number segments and frequency rescaled by purity
parse_cnas <- function(bb_file, purity){
  
  # Read battenberg file as dataframe
  req_cols = c("chr","startpos","endpos","nMaj1_A", "nMin1_A", "frac1_A",
               "nMaj2_A", "nMin2_A", "frac2_A")
  braw = read.table(bb_file, header=TRUE)[, req_cols]
  ## Select the clonal copy number solution for each battenberg row
  bs =  data.frame(t(apply(braw, 1, best_bb_solution)), stringsAsFactors = F)
  names(bs) = c("major_cn", "minor_cn", "clonal_frequency")
  ## Ensure all fields have correct type
  bs$major_cn <- as.integer(bs$major_cn)
  bs$minor_cn <- as.integer(bs$minor_cn)
  bs$clonal_frequency <- as.numeric(bs$clonal_frequency)
  
  # Rescale clonal frequency by tumour purity for mutationTimeR if set
  # Requirement detailed at: 
  #    https://github.com/gerstung-lab/MutationTimeR/issues/5#issuecomment-595717098 
  if (!is.na(purity)){ # Allowing us of NA to skip for testing
  bs$clonal_frequency <- bs$clonal_frequency*purity
  }
  bb = cbind(braw[, c("chr","startpos","endpos")], bs)
  bb = bb[complete.cases(bb),]
  bbgr = makeGRangesFromDataFrame(bb,start.field="startpos", end.field="endpos", keep.extra.columns = T)
  return(bbgr)
  
  # Note: MutationTimeR has a non-public function to load the battenberg file
  # loadBB <- function(file){
  #   tab <- read.table(file, header=TRUE, sep='\t')
  #   GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", tab[-3:-1])
  # }
}


#' Parse copy number segments from battenberg
#' Loads clonal and subclonal profile
#'
#' @param bb_file (character) Battenberg file with segment frequencies
#' @param purity (numeric) Overall sample purity inferred from battenberg
#'
#' @return (dataframe) Dataframe containing clonal copy number segments and frequency rescaled by purity
parse_cnas_sc <- function(bb_file, purity){
  
  # Read battenberg file as dataframe
  req_cols = c("chr","startpos","endpos","nMaj1_A", "nMin1_A", "frac1_A",
               "nMaj2_A", "nMin2_A", "frac2_A")
  braw = read.table(bb_file, header=TRUE)[, req_cols]
  braw = data.table(braw)

  # Create a data.frame with the clonal and subclonal copy number states
  ## Ensure that the same coordinates are given to each clone
  bs = data.table(
    rbind(
      braw[, .(chrom=chr, start=startpos, end=endpos, major_cn=nMaj1_A, minor_cn=nMin1_A, clonal_frequency=frac1_A)],
      braw[, .(chrom=chr, start=startpos, end=endpos, major_cn=nMaj2_A, minor_cn=nMin2_A, clonal_frequency=frac2_A)]
    )
  )
  setkey(bs, chrom, start, end)
  bs = bs[!is.na(clonal_frequency)]

  ## Ensure all fields have correct type
  bs$major_cn <- as.integer(bs$major_cn)
  bs$minor_cn <- as.integer(bs$minor_cn)
  bs$clonal_frequency <- as.numeric(bs$clonal_frequency)
  
  # Rescale clonal frequency by tumour purity for mutationTimeR if set
  # Requirement detailed at: 
  #    https://github.com/gerstung-lab/MutationTimeR/blob/e4e266a494482face03e831322a9df00b1b0ff93/R/MutationTime.R#L468
  #    https://github.com/gerstung-lab/MutationTimeR/issues/5#issuecomment-595717098 
  if (!is.na(purity)){ # Allowing us of NA to skip for testing
  bs$clonal_frequency <- bs$clonal_frequency*purity
  }
  
  # Convert to a granges and return
  bbgr = makeGRangesFromDataFrame(bs,start.field="start", end.field="end", keep.extra.columns = T)
  bbgr = bbgr[order(bbgr)]
  bbgr = bbgr[bbgr$clonal_frequency >0]
  return(bbgr)
  
}

#' Load DPClust results file
#'
#' @param cfile (character) DP clust output
#' @param purity (integer) Overall tumour purity estimated by battenberg
#' @param collapse_superclonal (logical) Whether to collapse superclonal clusters into the clonal cluster. Default TRUE
#' @return (dataframe) SNV clusters, mutation counts and proportions rescaled by purity
parse_cluster <- function(cfile, purity, collapse_superclonal = TRUE, ccf_clonal = 0.9){
  
  # Read DPclst file
  df = read.table(as.character(cfile), header=T, col.names=c("cluster","proportion","n_smss"))
  df = df[,c("cluster","n_smss","proportion")] # Reorder columns as expected
  ## Arrange rows in order of clonality
  df = df[order(df$proportion, decreasing = T), ]
  ## Drop superclonal clusters
  if (collapse_superclonal){
    is_superclonal = df$proportion >= ccf_clonal
    n_smss_superclonal = sum(df$n_smss[is_superclonal], na.rm = T)
    ccf_superclonal = sum(df$proportion[is_superclonal] * df$n_smss[is_superclonal], na.rm = T) / n_smss_superclonal
    df = df[!is_superclonal,]
    df = rbind(data.frame(cluster = 0, n_smss = n_smss_superclonal, proportion = ccf_superclonal), df)
  }
  rownames(df) = seq(nrow(df))
  df$cluster = seq(nrow(df))

  # Rescale CCF to tumour purity
  # Requirement detailed at:
  #   https://github.com/gerstung-lab/MutationTimeR/issues/6#issuecomment-653466479
  if(!is.na(purity)){
    df$proportion = (df$proportion/clonal_prop) * purity
    return(df)
  } else {
    return(df)
  }
}

bb_subclone <- function(x){
  
  # Assign battenberg row columns to variables
  maj1a=x[["nMaj1_A"]]; min1a=x[["nMin1_A"]]; frac1a= x[["frac1_A"]]
  maj2a=x[["nMaj2_A"]]; min2a=x[["nMin2_A"]]; frac2a=x[["frac2_A"]]
  # Return first clone if no subclonal fraction called
  if(is.na(frac2a)){
    return(list(c(x[["chr"]], x[["startpos"]], x[["endpos"]],maj1a, min1a, frac1a)))
  } else {
    return(
      list(
        c(x[["chr"]], x[["startpos"]], x[["endpos"]], maj1a, min1a, frac1a),
        c(x[["chr"]], x[["startpos"]], x[["endpos"]], maj2a, min2a, frac2a)
      )
    )
  }

}


parse_cnas_subclone <- function(bb_file, purity){
  
  # Read battenberg file as dataframe
  req_cols = c("chr","startpos","endpos","nMaj1_A", "nMin1_A", "frac1_A",
               "nMaj2_A", "nMin2_A", "frac2_A")
  braw = read.table(bb_file, header=TRUE)[, req_cols]
  ## Select the clonal copy number solution for each battenberg row
  bt = apply(braw, 1, bb_subclone)
  bt = unlist(bt, recursive=F)
  bs = data.frame(do.call(rbind, bt), stringsAsFactors = F)
  names(bs) = c("chr","startpos","endpos", "major_cn", "minor_cn", "clonal_frequency")
  ## Ensure all fields have correct format
  bs$major_cn <- as.integer(bs$major_cn)
  bs$minor_cn <- as.integer(bs$minor_cn)
  bs$clonal_frequency <- as.numeric(bs$clonal_frequency)
  
  # Rescale clonal frequency by tumour purity for mutationTimeR
  # Requirement detailed at: 
  #    https://github.com/gerstung-lab/MutationTimeR/issues/5#issuecomment-595717098 
  if (!is.na(purity)){ # Allowing us of NA to skip for testing
    bs$clonal_frequency <- bs$clonal_frequency*purity
  }
  bbgr = makeGRangesFromDataFrame(bs,start.field="startpos", end.field="endpos", keep.extra.columns = T)
  return(bbgr)
  
  # Note: MutationTimeR has a non-public function to load the battenberg file
  # loadBB <- function(file){
  #   tab <- read.table(file, header=TRUE, sep='\t')
  #   GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", tab[-3:-1])
  # }
}


# Main function
# Writes to output directory
run_ppcg_mt <- function(sname, vcf, bb, clusters, is.wgd, outdir, gender="male"){

  exp_file <- fs::path(outdir, paste0(sname, "_mutation_timer.vcf")) 
  if(fs::file_exists(exp_file)){
    loginfo("Skipping %s as output exists", sname)
    return()
  }

  # Call mutation TimeR on formatted objects
  ## Note classWgD function is available but we have consensus calls
  mt <- mutationTime(vcf=vcf, cn=bb, clusters=clusters, gender=gender, isWgd=is.wgd, n.boot=100)

  #loginfo("Timing SNVs and gains")
  # Annotate variants with mutation time
  vcf <- addMutTime(vcf, mt$V)
  
  # Add copy gain timings to bb as annotation
  bb_anno = bb
  mcols(bb_anno) <- cbind(mcols(bb_anno),mt$T)
  
  # Write results: annotated vcf file, parsed battenberg file, annotated battenberg object, mutation timer object
  writeVcf(vcf, fs::path(outdir, paste0(sname, "_mutation_timer.vcf")))
  write.csv2(bb, fs::path(outdir, paste0(sname, "_bb.csv")), row.names=F)
  saveRDS(bb_anno, fs::path(outdir, paste0(sname, "_bb.rds")))
  saveRDS(mt, fs::path(outdir, paste0(sname, "_mt.rds")))

  # Save MutationTimeR plot plot
  pdf(file=fs::path(outdir, paste0(sname, "_muttime.pdf")), width=10, height=7)
  plotSample(vcf, bb_anno, title = sname)
  dev.off()
  
}


get_ploidy = function(purdir, id){
  cc = fs::dir_ls(purdir, glob="*.txt") %>%   
      purrr::keep(~stringr::str_detect(.x, id))
      
    if(length(cc) == 0){
      print("No purity/ploidy file found")
      return(NA)
    }

    read.table(cc[[1]], header=T)$ploidy
  }