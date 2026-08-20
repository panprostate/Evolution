# A script to format raw PPCG data for analysis with SigTimer

# SigTimer requires clonal and subclonal counts, and active signatures for each sample.
# Active signatures are drawn from SigProfiler whereas counts are from a MutationTimeR VCF.
library(here)
library(readr)
library(glue)
library(fs)
library(qs)
library(stringr)
library(data.table)
library(logger)
library(foreach)
library(doParallel)
library(dplyr)
devtools::load_all("SigTimeR-main/")

cores = 25
# Functions ----
source("src/utils.R")

# Function to load sbs and id data from sigprofiler results
make_sigprofiler_object <- function(sigprofiler_dir){
    sbs_sig = readr::read_tsv(file.path(sigprofiler_dir, "SBS96_12_Signatures", "Decompose_Solution", "Signatures", "Decompose_Solution_Signatures.txt"))
    id_sig = readr::read_tsv(file.path(sigprofiler_dir, "COSMIC_ID83_Decomposed_Solution", "Signatures", "COSMIC_ID83_Signatures.txt"))

    sbs_act = readr::read_tsv(file.path(sigprofiler_dir, "SBS96_12_Signatures", "Decompose_Solution", "Activities", "Decompose_Solution_Activities.txt"))
    id_act = readr::read_tsv(file.path(sigprofiler_dir, "COSMIC_ID83_Decomposed_Solution", "Activities", "COSMIC_ID83_Activities.txt"))

    list(
        "sbs"=list("sig"=sbs_sig, "act"=sbs_act),
        "id"=list("sig"=id_sig, "act"=id_act)
    )
}

# Function to load active signatures for a given sample and signature type
get_active_sigs <- function(sid, sigtype, spo){
    sdata = spo[[sigtype]]
    sdata$act = mutate(sdata$act, Samples = extract_ppcg_id(Samples, full = FALSE))

    if(sigtype=="id"){
        sid_bool = stringr::str_detect(sdata$act[[1]], sid)
        sid_ix = which(sid_bool)[[1]]
    } else {
        sid_bool = (sdata$act[[1]]==sid)
        sid_ix = which(sid_bool)[[1]]
    }
    if(!any(sid_bool)){ # First is sample ids
        return(NULL)
    }
    d = sdata$act[sid_ix, -1] # First is sample ids
    asignames = names(d)[which(d>0)]  
    asigs = sdata$sig[, asignames] # First is signature names
    asigs = as.matrix(asigs)

    if(ncol(asigs)==1){
        return(NULL)
    }

    return(asigs)
}

# Function to laod sigmatrix data across epochs and signature types
make_sigmatrix_object <- function(sigmatgen_dir){

    sbs_f = list(
        "clonal"=fs::dir_ls(sigmatgen_dir, glob="*sbs::clonal.SBS96.all", recurse=T),
        "subclonal"=fs::dir_ls(sigmatgen_dir, glob="*sbs::subclonal.SBS96.all", recurse=T),
        "early"=fs::dir_ls(sigmatgen_dir, glob="*sbs::early.SBS96.all", recurse=T),
        "late"=fs::dir_ls(sigmatgen_dir, glob="*sbs::late.SBS96.all", recurse=T)
    )

    id_f = list(
        "clonal"=fs::dir_ls(sigmatgen_dir, glob="*id::clonal.ID83.all", recurse=T),
        "subclonal"=fs::dir_ls(sigmatgen_dir, glob="*id::subclonal.ID83.all", recurse=T),
        "early"=fs::dir_ls(sigmatgen_dir, glob="*id::early.ID83.all", recurse=T),
        "late"=fs::dir_ls(sigmatgen_dir, glob="*id::late.ID83.all", recurse=T)
    )

    list(
        "sbs"=lapply(sbs_f, fread),
        "id"=lapply(id_f, fread)
    )

}

# Function to load clonal and subclonal counts for a given sample and signature type
get_epoch_counts <- function(sid, sigtype, epochs, smo){
    sdata = smo[[sigtype]][epochs]
    sdata = lapply(sdata, as.data.frame)
    
    search_e1 = stringr::str_detect(names(sdata[[1]]), sid)
    search_e2 = stringr::str_detect(names(sdata[[2]]), sid)

    sample_present = any(search_e1) & any(search_e2)
    
    if(!sample_present){
        return(NULL)
    }

    e1_ix = which(search_e1)
    e2_ix = which(search_e2)

    counts = list(
        sdata[[1]][, which(search_e1)],
        sdata[[2]][, which(search_e2)]
    )

    empty_counts = any( sapply(counts, function(x) sum(x, na.rm=T)==0 )      )
    
    if(empty_counts){
        return(NULL)
    }

    names(counts) = epochs
    counts = lapply(counts, as.matrix)
    return(counts)
}

# Function to run sigtimeR given: ppcg_id, sigtype, epochs and outdir
call_ppcg_sigtimer <- function(sid, sigtype, epochs, spo, smo, outdir){
        
    # Dep res
    st_out = fs::dir_create(fs::path(outdir, sigtype, glue("{epochs[1]}_{epochs[2]}")))
    st_outfile = fs::path(st_out, glue("{sid}.qs"))
    if(fs::file_exists(st_outfile)){
        log_info(glue("SigTimeR output already exists for {sid} : {sigtype} : {epochs[1]} vs {epochs[2]}"))
        return(NULL)
    }

    # Get data for each sample, signature type and epoch pair
    actives = get_active_sigs(sid, sigtype, spo)
    counts = get_epoch_counts(sid, sigtype, epochs, smo)

    if(is.null(actives) | is.null(counts)){
        log_warn(glue("No data for {sid} : {sigtype} : {epochs[1]} vs {epochs[2]}"))
    }

    # Run SigTimeR
    st_data = SigTimeR::STData(counts[[1]], counts[[2]], actives)
    st_nnls = SigTimeR::time(st_data, method="nnls")
    log_info(" NNLS complete for {sid} : {sigtype} : {epochs[1]} vs {epochs[2]}")

    # st_dm = SigTimeR::time(st_data, method="dm")
    # Using NNLS method only
    st_dm = NULL
    log_info(" DM complete for {sid} : {sigtype} : {epochs[1]} vs {epochs[2]}")


    qs::qsave(
        list(
        "sid"=sid,
        "sigtype"=sigtype,
        "epochs"=epochs,
        "nnls"=st_nnls,
        "dm"=st_dm
        ),
        st_outfile)
}

# Main ----

# Setup output directories
outdir = fs::dir_create("outputs/00_preprocessing/sigtimer")

# Setup vars
sigprofiler_dir = "data/processed/ppcg_signature_definitions/"
sigmatgen_dir = "outputs/00_preprocessing/sigmatgen/"

# Get objects for sig activities search
sigpro_object = make_sigprofiler_object(sigprofiler_dir)
sigmat_object = make_sigmatrix_object(sigmatgen_dir)

# Create runtable combining sid, sigtypes, and epochs 
# sids = readr::read_lines("data/files/ppcg_snv_ids.txt")
sids = qc_pass()
sigtypes = c("sbs", "id")
epochs = list(c("clonal", "subclonal"), c("early", "late"))
runtable = expand.grid(sids, sigtypes, epochs, stringsAsFactors=F)
rix = seq(nrow(runtable))

# Run analysis
doParallel::registerDoParallel(cores=cores)
res = foreach(i = rix, .errorhandling="pass") %dopar% {
    sid = runtable[i,1]
    sigtype = runtable[i,2]
    epochs = unlist(runtable[i,3])
    call_ppcg_sigtimer(sid, sigtype, epochs, sigpro_object, sigmat_object, outdir)
    reslog = glue("COMPLETE: {sid} : {sigtype} : {epochs[1]} vs {epochs[2]}")
    return(reslog)
}
 
#print(res)