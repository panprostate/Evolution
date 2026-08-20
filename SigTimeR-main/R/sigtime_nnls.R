
sigtime_nnls <- function(data, sigs, id=NA, reps=1000){
  
  # Get result for epoch_1
  epoch1_nnls = nnls_sol(data[1,], sigs, reps)

  # Get result for ecoch_2
  epoch2_nnls = nnls_sol(data[2,], sigs, reps)
  
  nnls_raw = (list(epoch1=epoch1_nnls, epoch2=epoch2_nnls))
  return(nnls_raw)
}

# Returns the upper and lower CI from bootstrap replicates
nnls_bootstrap <- function(mut_counts, sigs, replic){
  # Function to simulate new counts
  simcounts <- function(x) rmultinom(1, size=sum(x), prob=x)
  
  # Generate bootstrap replicates
  boots = replicate(n=replic,{
    rmuts = simcounts(mut_counts)
    result = nnls_func(rmuts, sigs)
    return(result[[2]])
  })
  
  boots = t(boots)
  colnames(boots) = c(colnames(sigs), "error")

  return(boots)
}

nnls_func <- function(mut_counts, sigs){
  
  obj = (nnls::nnls(sigs, mut_counts))
  sol = coef(obj)
  
  n_unassigned = sum(mut_counts) - sum(sol)
  n_unassigned = ifelse(n_unassigned>0, n_unassigned, 0)
  sol = c(sol, n_unassigned)
  
  names(sol) = c(colnames(sigs), "error")
  return(list(obj, sol))
}

nnls_sol <- function(mut_counts, sigs, reps){
  mut_counts = as.matrix(mut_counts)
  result = nnls_func(mut_counts, sigs)
  nnls_boots = nnls_bootstrap(mut_counts, sigs, reps)
  return(list(obj=result[[1]], sol=result[[2]], boots=nnls_boots))
}


nnls_sol3 <- function(mut_counts, sigs){
  # Calculate signature weights by NNLS
  result = coef(nnls::nnls(sigs, mut_counts))
  
  # Get weights as proportion.
  #   If the weights under-assign counts, we need to include a fraction for unassigned mutations
  n_unassigned = sum(mut_counts) - sum(result)
  all_assigned = n_unassigned <= 0
  if(all_assigned){
    sol = c(norm_vec(result), 0)
  }else{
    sol = norm_vec(c(result, n_unassigned))
  }
  names(sol) = c(colnames(sigs), "error")
  rownames(sol) = "mean"
  
  # Return solution
  return(sol)
}

# Normalise such that sum of vector == 1
norm_vec <- function(x) {x / sum(x)}
