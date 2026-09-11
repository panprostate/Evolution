
#' Simulate sample clonal-subclonal mutation profile
#' @param sigs Signature matrix
#' @param total N mutations in sample
#' @param clonal_prop clonal proportion of total (numeric). Remainder of mutations will be subclonal
#' @param cs_sig_props proportion of signatures in clonal and subclonal (list)
#' @export
simulate <- function(sigs, total, clonal_prop, cs_sig_props){
  clonal_counts = round(total*clonal_prop, 0)
  subclonal_counts = total - clonal_counts
  clonal_muts <- helper_simulate_mutations(sigs, clonal_counts, cs_sig_props[[1]])
  subclonal_muts <- helper_simulate_mutations(sigs, subclonal_counts, cs_sig_props[[2]])
  return(list(
    params=list(
      sigs=sigs, total=total, clonal_prop=clonal_prop, cs_sig_props=cs_sig_props
    ),
    clonal=clonal_muts,
    subclonal=subclonal_muts,
    total=clonal_muts + subclonal_muts
  ))
}

#' Apply matrix multiplication to get mutation counts from signatures
#' @param sigs Signature matrix
#' @param counts N mutations in (sub)sample
#' @param props proportions of each signature
#' @param mut_classes 96-long character to name mutation classes
helper_simulate_mutations <- function(sigs, counts, props=NULL, mut_classes=NULL){
  
  # Create a vector of mutations per signature
  nsigs = dim(sigs)[2]
  # Set mutations per signature by proportions or equal spread
  if(!is.null(props)){
    stopifnot(nsigs==length(props) | sum(props)<=1)
    multi=counts*props
  }else{
    multi = rep(round(counts/nsigs, 0), nsigs)
  }
  cdif = sum(multi)-counts
  upd_ix = rep_len(seq_along(multi), abs(cdif))
  for(i in upd_ix){
    if(cdif>0){multi[i] = multi[i]-1}else{multi[i]=multi[i]+1}
  }
  
  # Multiply signature activities by mutations and balance
  muts = round(sigs %*% multi, 0)
  mdif = sum(muts)-counts
  upd_m_ix = rep_len(seq_along(muts), abs(mdif))
  for(i in upd_m_ix){
    if(mdif>0){muts[i] = muts[i]-1}else{muts[i]=muts[i]+1}
  }
  
  if(!is.null(mut_classes)){
    rownames(muts) = mut_classes
  }
  
  return(muts)
}
