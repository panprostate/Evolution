# Running dN / dS to identify selective pressures in multiple instances
library(stringr)
library(dndscv)
library(stringr)
library(parallel)
library(doParallel)
library(readr)
library(GenomicRanges)


empty_df = data.frame(
        sampleID = c(), 
        chr = c(),
        pos = c(),
        ref = c(),
        mut = c()
)



# Get all the indel and SNV calls required to run dndscv
#' @snvdir directory containing vcf calls for MutationTimeR SNVs
#' @indeldir directory containing vcf calls for MutationTimeR indels
#' @subclonal flag to filter SNVs by timing (T for "subclonal", F for "clonal", NA for "all")
#' @timing flag to filter SNVs by timing in clonal period (pre for clonal early,post for clonal late)

gather_calls <- function(snvdir, indeldir, ids, subclonal, ncores, timing = NA){    
    snvfs <- fs::dir_ls(snvdir, glob="*.vcf*")
    indelfs <- fs::dir_ls(indeldir, glob="*.vcf*")

    doParallel::registerDoParallel(cores=ncores)
    all_mts = foreach(id=ids) %dopar% {
        snvf = snvfs[grepl(id, snvfs)]
        snvf = snvf[!str_detect(snvf, "tbi$")]
        indelf = indelfs[grepl(id, indelfs)]
        indelf = indelf[!str_detect(indelf, "tbi$")]


        if (!is.na(timing)){
            
            snv_bbf = fs::dir_ls(snvdir, glob=glue("*{id}*bb.rds"))
            if (length(snv_bbf) == 0){
                snvs = empty_df
            } else {
                snvs = format_vcf2dn(snvf, bb_path = snv_bbf, indel = FALSE, timing = timing, subclonal = subclonal)
            }
            indel_bbf = fs::dir_ls(indeldir, glob=glue("*{id}*bb.rds"))
            if (length(indel_bbf) == 0){
                indels = empty_df
            } else {
                 indels = format_vcf2dn(indelf, bb_path = indel_bbf, indel = TRUE, timing = timing, subclonal = subclonal)
            }


        } else {
            if (length(snvf) < 1){
                snvs = empty_df
            } else {
                snvs = format_vcf2dn(snvf, indel = FALSE, subclonal = subclonal)
            }
            if (length(indelf) < 1){
                indels = empty_df
            } else {
                indels = format_vcf2dn(indelf, indel = TRUE, subclonal = subclonal)
            }
        }
        m_df = rbind(snvs,indels)
        return(m_df)
    }
    all_muts = rbindlist(all_mts)
}

# Get all the SNV and indel calls required to run dndscv from HMF vcf files
#' @vcfdir directory containing vcf calls for HMF samples
#' @ncores number of cores to use
gather_calls_hmf <- function(vcf_files,  ncores){
    doParallel::registerDoParallel(cores=ncores)
    all_mts = foreach(vcf_file = vcf_files, .combine = "rbind") %dopar% {
        vcf = readVcf(vcf_file)
        id = gsub(".purple.somatic.vcf.gz", "", basename(vcf_file))
        vcf <- vcf[rowRanges(vcf)$FILTER=="PASS"]
        return(
            data.frame(
                sampleID = id,
                chr = as.character(seqnames(vcf)),
                pos = start(vcf),
                ref = as.character(ref(vcf)),
                mut = as.character(unlist(alt(vcf)))
            )
        )
        
    }
}

# Filter mutations to LOH regions
#' @mts_df dataframe containing all the mutations 
#' @bb_df dataframe containing the bb information for all the samples
#' @loh_filter string indicating whether to filter to "loh" or "non-loh" regions
filter_mutations_by_loh <- function(mts_df, bb_df, loh_filter, ncores, clonal_loh = TRUE){
    if (!loh_filter %in% c("LOH", "non-LOH")){
        stop("Provide adequate value for loh_filter. Possible values: LOH, non-LOH")
    }
    # Remove chromosome X and Y where artificial LOH in male
    bb_df = bb_df[!bb_df$chr %in% c("X", "Y"),] %>% dplyr::filter(!is.na(startpos) &  !is.na(chr) & !is.na(nMin1_A))
    # Define LOH regions
    if (clonal_loh){
        bb_df$is_loh = 
            # Case 1 for clonal LOH: Clonal SCNA status and minor allele = 0
            (bb_df$frac1_A >= 0.9 & bb_df$nMin1_A == 0) | 
            # Case 2 for clonal LOH: Subclonal SCNA status but both subclones with minor = 0
            (bb_df$frac1_A < 0.9 & bb_df$nMin1_A == 0 & bb_df$nMin2_A == 0)
    } else {
        # Subclonal LOH if any of the alleles = 0
        bb_df$is_loh = ifelse(
            # Only one clone
            bb_df$frac1_A > 0.9, bb_df$nMin1_A == 0,
            # Multiple clones
            bb_df$nMin1_A == 0 | bb_df$nMin2_A == 0 
        )
    }

    # Filter sample by sample 
    doParallel::registerDoParallel(cores=ncores)
    ids = dplyr::intersect(unique(mts_df$sampleID), unique(bb_df$sampleID))
    stopifnot(length(ids) > n_distinct(mts_df$sampleID) * 0.9)
    mts_df = mts_df[mts_df$sampleID %in% ids,]
    bb_df = bb_df[bb_df$sampleID %in% ids,]

    filtered_mutations = foreach(id = ids, .combine = "rbind") %dopar% {
        if (loh_filter == "LOH"){
            bb_sub = bb_df[bb_df$sampleID == id & bb_df$is_loh,]
        } else {
            bb_sub = bb_df[bb_df$sampleID == id & !bb_df$is_loh,]
        }
        mts_sub = mts_df[mts_df$sampleID == id,]
        bb.gr = GRanges(
            seqnames = bb_sub$chr,
            ranges = IRanges(start = bb_sub$startpos, end = bb_sub$endpos)
        )
        mts.gr = GRanges(
            seqnames = mts_sub$chr,
            ranges = IRanges(start = mts_sub$pos, end = mts_sub$pos)
        )
        ovs = findOverlaps(mts.gr, bb.gr)
        return(mts_sub[queryHits(ovs),])
    }
    return(filtered_mutations)
}

time_relative_to_wgd <- function(muts_df, wgd_calls, dpclustdir, ncores){
    # 1. Get samples that underwent WGD
    wgd_samples = wgd_calls$sample[wgd_calls$is_wgd == "WGD"]
    muts_df = muts_df[muts_df$sampleID %in% wgd_samples,]
    # 2. Determine multiplicity of mutations
    # a. Load DPClust fit from the sample
    # b. Add multiplicity of mutations to dataframe
    # c. Return updated dataframe
    doParallel::registerDoParallel(cores=ncores)
    muts_df = foreach(sid = unique(muts_df$sampleID), .combine = "rbind") %dopar% {
        muts_sub = muts_df[muts_df$sampleID == sid,]
        dpclust_info = read_delim(list.files(dpclustdir, pattern = glue::glue("{sid}_allDirichletProcessInfo.txt"), recursive = TRUE, full.names = T))
        dpclust_info$mutation_id = paste0(dpclust_info$chr, ":", dpclust_info$end)
        muts_sub$mutation_id = paste0(muts_sub$chr, ":", muts_sub$pos)
        print(glue::glue("{mean(muts_sub$mutation_id %in% dpclust_info$mutation_id)*100}% of mutations found in DPClust info for sample {sid}"))
        muts_sub = merge(muts_sub, dpclust_info[,c("mutation_id", "no.chrs.bearing.mut")])        
        return(muts_sub)
    }
    
    # 3. Split mutations by multiplicity: pre-WGD (m > 1) and post-WGD (m <= 1)
    lmuts = list(
        "pre-WGD" = muts_df %>% dplyr::filter(no.chrs.bearing.mut > 1) %>% dplyr::select(sampleID, chr, pos, ref, mut),
        "post-WGD" = muts_df %>% dplyr::filter(no.chrs.bearing.mut <= 1) %>% dplyr::select(sampleID, chr, pos, ref, mut)
    )
    return(lmuts)
}


# Filter mutations to specific groups of patients
#' @mts_df dataframe containing all the mutations for a given patient
#' @s_df dataframe containing column with sample IDs, trajectory, gs_group and clinical outcome
#' @ordering specific trajectories to which mutations are selected
#' @gs_group specific GS group to which mutations are selected
#' @outcome specific MFS clinical outcome in the patients for which mutations are selected
filter_mutations <- function(mts_df, s_df, ordering = "all", gs_group = "all", outcome = "all", gs_grade_collapsed = "all"){
    # Error handling
    if (!outcome %in% c("no_met", "met", "all")){
        stop("Provide adequate value for clinical outcome")
    } 
    if (!gs_group %in% c("1", "2", "3", "4+", "all")){
        stop("Provide adequate value for gs_group")
    }
    if (!gs_grade_collapsed %in% c("1-2", "3-5", "all")){
        stop("Provide adequate value for gs_grade_collapsed")
    } 
    if (!ordering %in% c("all", paste("Ordering", 1:3))){
        stop("Provide adequate value for ordering")
    }

    # Filter as needed
    if (ordering != "all"){
        s_df = s_df[s_df$trajectory == ordering,]
    }
    if (gs_group != "all"){
        s_df = s_df[s_df$gs_group == gs_group,]
    }
    if (gs_grade_collapsed != "all"){
        s_df = s_df[s_df$gs_group_collapsed == gs_grade_collapsed,]
    } 
    if (outcome != "all"){
        s_df = s_df[s_df$clinical_outcome == outcome,]
    }
    samples = s_df$sample
    return(mts_df[mts_df$sampleID %in% samples,])
}

run_dnds <- function(mts_df, s_df, pathway_genes, pathway, base_outdir, ordering = "all", gs_group = "all", outcome = "all", gs_grade_collapsed = "all", genci = F){
    mts_df = filter_mutations(mts_df, s_df, ordering = ordering, gs_group = gs_group, outcome = outcome, gs_grade_collapsed = gs_grade_collapsed)
    outdir = file.path(base_outdir, glue::glue("ordering_{ordering}_gs_{gs_group}_outcome_{outcome}"), pathway)
    dir.create(outdir, recursive = T)
    tryCatch({
        if (pathway == "all"){
            r = dndscv(mts_df)
        } else {
            target_genes = pathway_genes$gene[pathway_genes$pathway == pathway]
            r = dndscv(mts_df, gene_list = target_genes)
        }
        fs::dir_create(outdir, recurse = TRUE)
        write.table(r$sel_cv, file.path(outdir, ".gene_dnds.tsv"))
        write.table(r$globaldnds, file.path(outdir, "global_dnds.tsv"))

        # Add confidence intervals if desired
        if (genci){
            target_genes = read_delim(gs_path)
            target_genes = target_genes$gene
            r = dndscv(mts_df, gene_list = target_genes, outmats = T)
            confint = geneci(r)
            write_delim(confint, file.path(outdir, ".ci95_genes.tsv"))
        }
    }, error = function(e){
        print(paste0("Error trying to run dnds and output to ", outdir))
    })
}




################################ ARCHIVE - OLD FUNCTIONS

# snvdir <- "data/raw/data_releases/WGS_Data_Release/Somatic_variants/SNVs/Filtered_SNV_VCFs_20_April_2021/"
# indeldir <- "data/raw/data_releases/WGS_Data_Release/Somatic_variants/InDels/PPCG_indels_16_Dec_2022/Tier_1/"
# cnadir <- "data/raw/data_releases/WGS_Data_Release/Somatic_variants/SCNA/Subclonal_SCNA_01_June_2020"
# ccfdir <- "data/raw/pptech_exchange/Working_Groups/PP_Evo/SNV_clonality_DPClust_results_16_02_2023/"
# purdir <- "data/raw/PPCG_Data_Release/210414_dpclust/Single_sample_DPClust/PPCG_898_samples/Cellularity/"
# wgdfile <- "data/raw/data_releases/WGS_Data_Release/Somatic_variants/WGD/WGD_consensus_estimates_25_May_2021.txt"

# Large function to run dnds pipeline for different lists of driver genes
wrapper_dnds_allgene_types = function(snvdir, indeldir, ids, base_outdir, driver_dir, timing = NA, subclonal = NA, ncores, genci = TRUE){
    # all genes
    outdir = file.path(base_outdir, "all")
    dnds_pipeline(
        snvdir = snvdir, 
        indeldir = indeldir, 
        ids = ids, 
        outdir = outdir,
        gr_dir = NULL,  
        gs_path = NULL, 
        subclonal = subclonal,
        timing = timing,
        ncores = ncores, 
        genci = FALSE
    )

    # wedge_2018
    driver_fp = file.path(driver_dir, "prostate_wedge_2018.tsv")
    outdir = file.path(base_outdir, "wedge_2018")
    dnds_pipeline(
        snvdir = snvdir, 
        indeldir = indeldir, 
        ids = ids, 
        outdir = outdir,
        gr_dir = NULL,  
        gs_path = driver_fp, 
        subclonal = subclonal,
        timing = timing,
        ncores = ncores,
        genci = genci
    )

    # armenia 2018
    driver_fp = file.path(driver_dir, "prostate_armenia_2018.tsv")
    outdir = file.path(base_outdir, "armenia_2018")
    dnds_pipeline(
        snvdir = snvdir, 
        indeldir = indeldir, 
        ids = ids, 
        outdir = outdir,
        gr_dir = NULL,  
        gs_path = driver_fp, 
        timing = timing,
        subclonal = subclonal,
        ncores = ncores,
        genci = genci
    )

    # ppcg drivers
    driver_fp = file.path(driver_dir, "prostate_ppcg.tsv")
    outdir = file.path(base_outdir, "ppcg")
    dnds_pipeline(
        snvdir = snvdir, 
        indeldir = indeldir, 
        ids = ids, 
        outdir = outdir,
        gr_dir = NULL,  
        gs_path = driver_fp, 
        timing = timing,
        subclonal = subclonal,
        ncores = ncores, 
        genci = genci
    )

    # cancer genes cosmic
    driver_fp = file.path(driver_dir, "cancergenes_cosmic_cgc81.tsv")
    outdir = file.path(base_outdir, "cosmic_cancergenes")
    dnds_pipeline(
        snvdir = snvdir, 
        indeldir = indeldir, 
        ids = ids, 
        outdir = outdir,
        gr_dir = NULL,  
        gs_path = driver_fp, 
        timing = timing,
        subclonal = subclonal,
        ncores = ncores,
        genci = genci
    )

    # essential genes
    driver_fp = file.path(driver_dir, "essential_blomen.tsv")
    outdir = file.path(base_outdir, "essential_blomen")
    dnds_pipeline(
        snvdir = snvdir, 
        indeldir = indeldir, 
        ids = ids, 
        outdir = outdir,
        gr_dir = NULL,  
        gs_path = driver_fp, 
        timing = timing,
        subclonal = subclonal,
        ncores = ncores,
        genci = genci
    )

    # prostate essential genes cancer lines DepMap
    driver_fp = file.path(driver_dir, "essential_depmap_prostate.tsv")
    outdir = file.path(base_outdir, "essential_depmap_prostate")
    dnds_pipeline(
        snvdir = snvdir, 
        indeldir = indeldir, 
        ids = ids, 
        outdir = outdir,
        gr_dir = NULL,  
        gs_path = driver_fp, 
        timing = timing,
        subclonal = subclonal,
        ncores = ncores,
        genci = genci
    )

    driver_fp = file.path(driver_dir, "hla_genes.tsv")
    outdir = file.path(base_outdir, "hla_genes")
    dnds_pipeline(
        snvdir = snvdir, 
        indeldir = indeldir, 
        ids = ids, 
        outdir = outdir,
        gr_dir = NULL,  
        gs_path = driver_fp, 
        timing = timing,
        subclonal = subclonal,
        ncores = ncores,
        genci = genci
    )
}

dnds_pipeline = function(snvdir, indeldir, ids, outdir, gr_dir = NULL, gs_path = NULL, timing = NA, subclonal = NA, ncores = 10, genci = FALSE){
    # Runs dN/dS in samples 'ids'
    # gr_path is a tsv file containing for each "id" included in the analysis
    # the regions we want to subset dN / dS analysis to
    # gs_path is a list of genes to which to restrict dN / dS analysis
    
    snvfs <- fs::dir_ls(snvdir, glob="*.vcf*")
    indelfs <- fs::dir_ls(indeldir, glob="*.vcf*")

    # ids = path_file(snvfs) %>% stringr::str_extract("^PPCG[0-9A-Za-d]+") # first PPCG id is target sample, second matched normal
    
    doParallel::registerDoParallel(cores=ncores)
    all_mts = foreach(id=ids) %dopar% {
        snvf = snvfs[grepl(id, snvfs)]
        snvf = snvf[!str_detect(snvf, "tbi$")]
        indelf = indelfs[grepl(id, indelfs)]
        indelf = indelf[!str_detect(indelf, "tbi$")]


        if (!is.na(timing)){
            
            snv_bbf = fs::dir_ls(snvdir, glob=glue("*{id}*bb.rds"))
            if (length(snv_bbf) == 0){
                snvs = empty_df
            } else {
                snvs = format_vcf2dn(snvf, bb_path = snv_bbf, indel = FALSE, subclonal = subclonal)
            }
            indel_bbf = fs::dir_ls(indeldir, glob=glue("*{id}*bb.rds"))
            if (length(indel_bbf) == 0){
                indels = empty_df
            } else {
                 indels = format_vcf2dn(indelf, bb_path = indel_bbf, indel = TRUE, subclonal = subclonal)
            }


        } else {
            if (length(snvf) < 1){
                snvs = empty_df
            } else {
                snvs = format_vcf2dn(snvf, indel = FALSE, subclonal = subclonal)
            }
            if (length(indelf) < 1){
                indels = empty_df
            } else {
                indels = format_vcf2dn(indelf, indel = TRUE, subclonal = subclonal)
            }
        }


        m_df = rbind(snvs,indels)
        if (!is.null(gr_dir)){
            gr_path = fs:dir_ls(gr_dir, glob="*.bed")
            m_df = filter_to_gr(m_df, gr_path)
        }
        m_df
    } 

    all_mts = as.data.frame(data.table::rbindlist(all_mts))
    # if (!is.null(gs_path)){
    #     all_mts = filter_to_gr(all_mts, gs_path)
    # }
    tryCatch({
        if (!is.null(gs_path)){
            target_genes = read_delim(gs_path)
            target_genes = target_genes$gene
            r = dndscv(all_mts, gene_list = target_genes)
        } else {
            r = dndscv(all_mts)
        }
        fs::dir_create(outdir, recurse = TRUE)
        write.table(r$sel_cv, file.path(outdir, "gene_dnds.tsv"))
        write.table(r$globaldnds, file.path(outdir, "global_dnds.tsv"))

        # Add confidence intervals if desired
        if (genci){
            if (!is.null(gs_path)){
                target_genes = read_delim(gs_path)
                target_genes = target_genes$gene
                r = dndscv(all_mts, gene_list = target_genes, outmats = T)
                confint = geneci(r)
                write_delim(confint, file.path(outdir, "ci95_genes.tsv"))
            } else {
                r = dndscv(all_mts, outmats = T)
                confint = geneci(r)
                write_delim(confint, file.path(outdir, "ci95_genes.tsv"))
            }
        }
    }, error = function(e){
        print(paste0("Error trying to run dnds and output to ", outdir))
    })
    
    
    # print(head(all_mts))
    # print(nrow(all_mts))
    # r = dndscv(all_mts)
    # fs::dir_create(outdir, recurse = TRUE)
    # write.table(r$sel_cv, file.path(outdir, "gene_dnds.tsv"))
    # write.table(r$globaldnds, file.path(outdir, "global_dnds.tsv"))
}


format_vcf2dn = function(vcf_path, bb_path = NA, timing = NA, indel = FALSE, subclonal = NA){
    # Reads a vcf file and outputs a df with a structure compatible with dN/dS
    id = fs::path_file(vcf_path) %>% stringr::str_extract("^PPCG[0-9A-Za-d]+") # first PPCG id is target sample, second matched normal
    if (indel){vcf = parse_indel_vcf(vcf_path)} else {vcf = parse_vcf(vcf_path)}
    
    # Special processing when we want to run the pipeline for only subclonal or only clonal mutations
    if (!is.na(subclonal)){
        if (subclonal){
            vcf = vcf[na.omit(which(info(vcf)$CLS == "subclonal")),]
        } else {
            vcf = vcf[na.omit(which(str_detect(info(vcf)$CLS, "^clonal"))),]
        }
    }
    
    # Special processing when we want to only focus on preWGD-LOH, postWGD-LOH, preWGD-noLOH and postWGD-noLOH
    if (!is.na(timing)){
        # If we want preWGD, we want clonal [early]; if postWGD we go for clonal [late]
        # if not relative to WGD, just don't filter
        if (str_detect(timing, "^pre")){
            vcf = vcf[na.omit(which(info(vcf)$CLS == "clonal [early]")),]
        } else if (str_detect(timing, "^post")) {
            vcf = vcf[na.omit(which(info(vcf)$CLS == "clonal [late]")),]
        }
        # If we want LOH, we need minor CN @0, otherwise we can take everything that has minor CN > 0
        if (str_detect(timing, "LOH$") & !str_detect(timing, "noLOH")){
            bb = readRDS(bb_path)
            bb = bb[na.omit(which(bb$minor_cn == 0)),]
            ovs = findOverlaps(rowRanges(vcf), bb)
            vcf = vcf[unique(ovs@from),]
        } else if (str_detect(timing, "noLOH")){
            bb = readRDS(bb_path)
            bb = bb[na.omit(which(bb$minor_cn > 0)),]
            ovs = findOverlaps(rowRanges(vcf), bb)
            vcf = vcf[unique(ovs@from),]
        }
    }
    m.gr = rowRanges(vcf)
    # Remove SNPs entries
    # Example rs1421901 in PPCG0001a_DNA_vs_PPCG0001b_DNA.filtered.extended.pon.vcf.gz
    m.gr = m.gr[!grepl("rs[0-9]+", names(m.gr))]

    # return empty data.frame if no mutations
    if (length(m.gr) == 0){return(empty_df)}
    
    # [Format according to dN / dS input]
    # (https://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/dNdScv.html)
    m_df = data.frame(
        sampleID = id, 
        chr = as.vector(seqnames(m.gr)),
        pos = start(m.gr),
        ref = as.character(m.gr$REF),
        mut = as.character(unlist(m.gr$ALT))
    )
    # remove duplicates
    m_df = unique(m_df)
    return(m_df)
}


# Returns a filtered set of mutations located in only certain genomic regions
#' @param m_df dataframe with mutations
#' @param gr_path path to file containing target regions. Should have columns `chr` and `start` (and optionally `end`)
filter_to_gr = function(m_df, gr_path){
  m.gr = GRanges(seqnames = m_df$chr, IRanges(start = m_df$pos, end = m_df$pos))
  # load the path to file with genomic regions / mutations where to subset analysis
  gr_df = read_delim(gr_path)
  if (!all(c("chr", "start") %in% colnames(gr_df))){
    stop("Please check columns `chr` and `start` are included")
  }
  if (! "end" %in% colnames(gr_df)){
    gr_df$end = gr_df$start 
  }
  gr_df$end[is.na(gr_df$end)] = gr_df$start[is.na(gr_df$end)]
  gr.gr = GRanges(seqnames = gr_df$chr, IRanges(start = gr_df$start, end = gr_df$end))
  ovs = findOverlaps(m.gr, gr.gr)
  return(m_df[unique(ovs@from),]) # unique() to avoid duplicated mutations if gr_path file covers same region more than once
}

# FUNCTION TO COMPARE DN/DS OF TWO DATASETS FROM [MARTINCORENA'S GROUP](https://zenodo.org/records/3966023#.YanjS_HP2cZ)
# Function. Comparing dN/dS values between two datasets using the uniform model and removing global differences in dN/dS between two datasets
variable_dNdS_twodatasets = function(dnds1, dnds2, genestotest) {
    
    library("dndscv")
    pvec = rmisvec = rtruvec = rep(NA, length(genestotest)) # Initialising vectors for p-values and for the ratios of wmis and wtru between dataset 1 and 2
    w1 = dnds1$globaldnds$mle; names(w1) = dnds1$globaldnds$name
    w2 = dnds2$globaldnds$mle; names(w2) = dnds2$globaldnds$name
    
    for (g in 1:length(genestotest)) {
        
        # We can implement a simple LRT model based on the uniform dNdS model
        # This is different from the Fisher test in that it uses synonymous mutations (i.e. dN/dS ratios)
        # instead of comparing the contribution of nonsyn muts of a gene *relative* to other genes.
        # Being a uniform model it assumes no considerable changes in the mutation rate variation or coverage
        # across genes in both datasets. But takes into account signature and rate variation between two
        # datasets.
        # H0: wmis1==wmis2 & wtru1==wtru2
        # H1: wmis1!=wmis2 & wtru1!=wtru2
        # This is simply done using obs1, exp1, obs2, exp2 (y1 and y2 vectors below)
        
        y1 = as.numeric(dnds1$genemuts[dnds1$genemuts$gene==genestotest[g],])
        y2 = as.numeric(dnds2$genemuts[dnds2$genemuts$gene==genestotest[g],])
        
        # Global dN/dS ratios from all other genes (to normalise the differences for the gene being tested)        
        ind1 = dnds1$genemuts$gene!=genestotest[g]
        ind2 = dnds2$genemuts$gene!=genestotest[g]
        wmis1_global = sum(dnds1$genemuts$n_mis[ind1])/sum(dnds1$genemuts$exp_mis[ind1])
        wmis2_global = sum(dnds2$genemuts$n_mis[ind2])/sum(dnds2$genemuts$exp_mis[ind2])
        wtru1_global = sum(dnds1$genemuts$n_non[ind1]+dnds1$genemuts$n_spl[ind1])/sum(dnds1$genemuts$exp_non[ind1]+dnds1$genemuts$exp_spl[ind1])
        wtru2_global = sum(dnds2$genemuts$n_non[ind2]+dnds2$genemuts$n_spl[ind2])/sum(dnds2$genemuts$exp_non[ind2]+dnds2$genemuts$exp_spl[ind2])
        
        # MLE dN/dS ratios using the uniform model under H0 and H1
        wmis_mle0 = (y1[3]+y2[3])/(y1[7]*wmis1_global+y2[7]*wmis2_global)
        wtru_mle0 = sum(y1[4:5]+y2[4:5])/sum(y1[8:9]*wtru1_global+y2[8:9]*wtru2_global)
        wmis_mle1 = c(y1[3],y2[3])/c(y1[7]*wmis1_global,y2[7]*wmis2_global)
        wtru_mle1 = c(sum(y1[4:5]),sum(y2[4:5]))/c(sum(y1[8:9]*wtru1_global),sum(y2[8:9]*wtru2_global))
        
        # Observed and predicted counts under H0 and H1
        obs = as.numeric(c(y1[3], sum(y1[4:5]), y2[3], sum(y2[4:5])))
        exp0 = as.numeric(c(y1[7]*wmis1_global*wmis_mle0, sum(y1[8:9])*wtru1_global*wtru_mle0, y2[7]*wmis2_global*wmis_mle0, sum(y2[8:9])*wtru2_global*wtru_mle0))
        exp1 = as.numeric(c(y1[7]*wmis1_global*wmis_mle1[1], sum(y1[8:9])*wtru1_global*wtru_mle1[1], y2[7]*wmis2_global*wmis_mle1[2], sum(y2[8:9])*wtru2_global*wtru_mle1[2])) # Note that exp1 == obs (we only have this line here for confirmation purposes)
        ll0 = c(sum(dpois(x=obs[c(1,3)], lambda=exp0[c(1,3)], log=T)), sum(dpois(x=obs[c(2,4)], lambda=exp0[c(2,4)], log=T)))
        ll1 = c(sum(dpois(x=obs[c(1,3)], lambda=exp1[c(1,3)], log=T)), sum(dpois(x=obs[c(2,4)], lambda=exp1[c(2,4)], log=T)))
        
        # One-sided p-values
        pvals = (1-pchisq(2*(ll1-ll0), df=1))
        if (wmis_mle1[1]<wmis_mle1[2]) { pvals[1] = 1 } else { pvals[1] = pvals[1]/2 }
        if (wtru_mle1[1]<wtru_mle1[2]) { pvals[2] = 1 } else { pvals[2] = pvals[2]/2 }
        
        # Saving the results
        pvec[g] = 1 - pchisq(-2 * sum(log(pvals)), df = 4) # Fisher combined p-value
        rmisvec[g] = wmis_mle1[1]/wmis_mle1[2]
        rtruvec[g] = wtru_mle1[1]/wtru_mle1[2]
        
    }
    out = data.frame(genestotest,pvec,rmisvec,rtruvec)
    return(out)
}