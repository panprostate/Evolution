
sigtime_nnls_extract <- function(nnls_raw, id=NA){
  # Extract from NNLS raw the table of epoch1 and epoch2 results as expected
  e1_means = extract_nnls_means(nnls_raw, 1)
  e2_means = extract_nnls_means(nnls_raw, 2)
  e1_boots = extract_nnls_cis(nnls_raw, 1)
  e2_boots = extract_nnls_cis(nnls_raw, 2)
  
  # Row-bind and return
  extracted = Reduce(rbind, list(e1_means, e2_means, e1_boots, e2_boots))
  extracted = data.frame(extracted)
  extracted$epoch = paste0("epoch_", c(1,2,1,1,2,2))
  extracted$info = c("mean", "mean", "lower", "upper", "lower", "upper")
  extracted$id = as.character(id)
  extracted$method = "nnls"
  
  # Pivot
  extracted = tidyr::pivot_longer(extracted, 1:ncol(e1_means), names_to="signature", values_to="prop")
  extracted = extracted[, c("info", "epoch", "id", "method", "signature", "prop")]
  return(extracted)
}

extract_nnls_means <- function(nnls_raw, ix){
  res = nnls_raw[[ix]]
  weights = res$sol
  sol_scaled = matrix(
    norm_vec(weights),
    ncol=length(weights)
  )
  
  rownames(sol_scaled) = "mean"
  return(sol_scaled)
}

extract_nnls_cis <- function(nnls_raw, ix){
  boots_raw = nnls_raw[[ix]]$boots
  boots = t(apply(boots_raw, 1, function(x) norm_vec(x))) #Normalise
  
  cis = matrix(c(
    sapply(seq(ncol(boots)), function(x) quantile(boots[,x], 0.025, na.rm=T)),
    sapply(seq(ncol(boots)), function(x) quantile(boots[,x], 0.975, na.rm=T))
  ), nrow=ncol(boots))
  rownames(cis) = colnames(boots)
  colnames(cis) = c("lower", "upper")
  cis = t(cis) # Place values on cols as expected
  return(cis)
}
