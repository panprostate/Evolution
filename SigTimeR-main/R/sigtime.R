#'@export
time <- function(data, method, metadata=NA){

  priors = stsigs(data)
  counts = stdata(data)
    
  if(method=="hdp"){

    hmodel = sigtime_hdp(counts, priors)
    estimate = sigtime_hdp_extract(hmodel, priors)
    estimate_fc = calculate_fold_change(estimate)
    estimate = rbind(estimate, estimate_fc)
    
    res = STOut(estimate=estimate, metadata=as.character(metadata), raw=list(hmodel))
    return(res)
  }
  
  if(method=="nnls"){
    # Run NNLS timer on the signature counts
    nnls_raw = sigtime_nnls(counts, priors)
    estimate = sigtime_nnls_extract(nnls_raw)
    estimate_fc = calculate_fold_change(estimate)
    estimate = rbind(estimate, estimate_fc)
    
    res = STOut(estimate=estimate, metadata = as.character(metadata), raw=nnls_raw)
    return(res)
  }
  
  if(method=="dm"){
    # Run NNLS timer on the signature counts
    dm_raw = sigtime_dm(counts, priors)
    estimate = sigtime_dm_extract(dm_raw)
    estimate_fc = calculate_fold_change(estimate)
    estimate = rbind(estimate, estimate_fc)
     
    res = STOut(estimate=estimate, metadata = as.character(metadata), raw=dm_raw)
    return(res)
  }
  
  
  stop("Unrecognised method argument.")
}

calculate_fold_change <- function(st_res){
  headers = c("info", "epoch", "id", "method", "signature", "prop")
  st_wide = tidyr::pivot_wider(
    st_res[st_res$info=="mean",],
    names_from=c("epoch"),
    values_from=c("prop")
  )
  # Add fold-change
  st_wide$prop = log_fold_change(st_wide$epoch_1, st_wide$epoch_2)
  st_wide$info = "log_fc"
  st_wide$epoch = NA
  st_wide$epoch_1=NULL
  st_wide$epoch_2=NULL
  return(st_wide[, headers])
}

log_fold_change <- function(e1, e2){
  log2(
    (e2 /(1-e2)) /
      (e1 /(1-e1)) 
  )
}
