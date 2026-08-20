
match_hdp_prior_names <- function(hmodel, priors){
  # Create a lookup table that matches prior signatures with dataset
  # Iterate in order of "closest" match
  novel = (comp_categ_distn(hmodel)$mean)
  opts = expand.grid(colnames(priors), rownames(novel))
  names(opts) = c("sig", "novel")
  # Iterate over rows and get cosine similarity
  opts$cos = sapply(seq(nrow(opts)), function(x){
    signame = opts[x,1]
    novname = opts[x,2]
    lsa::cosine(priors[, signame], novel[novname, ])
  })
  # Drop NA
  opts = opts[!is.na(opts$cos),]
  
  # For each signature, extract a row from the dataset
  result = data.frame()
  for( i in seq_along(colnames(priors))){
    x = opts[which.max(opts$cos),]
    result = rbind(result, x)
    opts = opts[!(opts$sig %in% result$sig | opts$novel %in% result$novel), ]
  }
  
  output = as.character(result$novel)
  names(output) = result$sig
  #return(factor(output, levels=colnames(priors)))
  return(output)
}

sigtime_hdp_extract <- function(hmodel, priors, id=NA){
  # Match HDP signatures with prior signature names
  hdp_prior_match = match_hdp_prior_names(hmodel, priors)
  # There is always a component 0 that includes our errors
  hdp_prior_match = c(hdp_prior_match, c("error"="0"))
  
  # Get indices of epochs
  dps = numdp(final_hdpState(chains(hmodel)[[1]]))
  ix = c(dps-1, dps)
  
  # Select mean
  extracted_means =  comp_dp_distn(hmodel)$mean[ix, hdp_prior_match]
  colnames(extracted_means) = names(hdp_prior_match)
  #error = matrix(1-rowSums(extracted_means)); colnames(error) = "error"
  #extracted_means = cbind(extracted_means, error)
  
  # Select lower_ci, upper_ci for each signature
  epoch_1_cis = comp_dp_distn(hmodel)$cred.int[[ix[1]]][, hdp_prior_match]
  epoch_2_cis = comp_dp_distn(hmodel)$cred.int[[ix[2]]][, hdp_prior_match]
  extracted = Reduce(rbind, list(extracted_means, epoch_1_cis, epoch_2_cis))
  extracted = data.frame(extracted)
  
  # Label rows
  extracted$info = c("mean", "mean", "lower", "upper", "lower", "upper")
  extracted$epoch = c("epoch_1", "epoch_2", "epoch_1", "epoch_1", "epoch_2", "epoch_2")
  extracted$id = id
  extracted$method = "hdp"
  
  # Pivot to make one record per signature
  extracted = tidyr::pivot_longer(extracted, names(hdp_prior_match), names_to="signature", values_to="prop")
  extracted = extracted[, c("info", "epoch", "id", "method", "signature", "prop")]
  
  # The final dataframe should have a row for each signature.
  return(extracted)
}

