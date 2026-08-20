# Load qs files for each instance and get the dataset
library(qs)
library(dplyr)
library(data.table)
library(foreach)
doParallel::registerDoParallel(cores=10)

load_sig_qs <- function(x, cls="nnls"){
    d = qs::qread(x)
    res = d[[cls]]@estimate %>% dplyr::filter(info=="log_fc")
    res$id = d$sid
    res$epoch = paste0(d$epochs, collapse="_")
    res$sigtype = d$sigtype
    return(res)
}

load_sig_clonal <- function(x, cls="nnls"){
    d = qs::qread(x)
    res = d[[cls]]@estimate %>% dplyr::filter(info=="mean", epoch == "epoch_1")
    res$id = d$sid
    res$epoch = paste0(d$epochs, collapse="_")
    res$sigtype = d$sigtype
    return(res)
}

sigtimer_root = "outputs/00_preprocessing/sigtimer"
sigfiles = fs::dir_ls(sigtimer_root, glob="*.qs", recurse=T)
sigtimer_nnls = file.path(sigtimer_root, "sigtimer_nnls.csv")

# For each sigtype folder in sigtimer,
# for each subfolder, 
# Call load_sig_qs on all qs files and merge
sf = foreach(x=sigfiles, .combine="rbind") %dopar% {
    load_sig_qs(x, cls="nnls")
}
data.table::fwrite(sf, sigtimer_nnls)
doParallel::stopImplicitCluster()

# Load clonal activity of SBS signatures
sf = foreach(x=sigfiles[grepl("clonal_subclonal", sigfiles)], .combine="rbind") %dopar% {
    load_sig_clonal(x, cls="nnls")
}
data.table::fwrite(sf, file.path(sigtimer_root, "sigtimer_clonal_nnls.csv"))