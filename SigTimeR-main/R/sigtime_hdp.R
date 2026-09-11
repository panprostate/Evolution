#' Run heirarchical dirichlet process for one or multiple samples
#' @param data A matrix of mx96 mutation counts, where m is clonal then subclonal for each sample.
#' @param priors A matrix of prior mutational signatures present in the whole dataset.
#' @param infer Attempt to infer novel signatures in the dataset.
#' @return A HDP object
sigtime_hdp <- function(data, priors, infer=FALSE){
  
  # Set function variables
  # Optional: Infer novel signatures with HDP by setting uninformative alpha prior (Gamma dist)
  if(infer){
    prior_alpha_a = c(1,1); prior_alpha_b = c(1,1)
  }else{
    prior_alpha_a = c(1e-5, 1e-5); prior_alpha_b = c(1e5, 1e5)
  }
  n_sigs = dim(priors)[2]
  new_dp_parent_index = 1 + n_sigs + 1
  sample_dp_count = nrow(data)
  wparent_dp_index = new_dp_parent_index : (new_dp_parent_index+sample_dp_count)
  sample_dp_index = wparent_dp_index[-1] # Drop parent
  
  # Initialise HDP with priors
  hmodel <- hdp_prior_init(
    prior_distn = priors,
    prior_pseudoc = rep(1000, n_sigs), # Uniform weight to all signatures
    hh=rep(1,96), # Uniform prior for DP on mutations
    # Set hyperparameter for gamma priors over DP concentration params
    alphaa=prior_alpha_a, # Shape
    alphab=prior_alpha_b # Rate
  )
  
  # Add HDP nodes for each element of our data, along with an additional parent node (required).
  hmodel <- hdp_addconparam(hmodel, alphaa=c(1,1), alphab=c(1,1)) # Add two conparams
  hmodel <- hdp_adddp(
    hdp=hmodel,
    numdp=sample_dp_count+1, # +1 for additional parent node
    # Parent inherits from global DP (1), samples inherit from parent which is after DP of last sig.
    ppindex= c(1, rep(1+n_sigs+1, sample_dp_count)), 
    # Parent inherits from 3rd conparam, remainder from 4th.
    cpindex= c(3, rep(4, sample_dp_count)) 
  )
  
  # Base R apply returns different format if hdp priors are n=2. This breaks Nicola's `hpd_setdata`:
  # https://github.com/nicolaroberts/hdp/blob/c78989b537a3b1c677a10df6c23725dca10bdac5/R/hdp_setdata.R#L46
  # Here, we use a new hpd_setdata function with a fix.
  hmodel <- hdp_setdata(
    hdp=hmodel,
    dpindex=sample_dp_index, # Set indices of DPs for data
    data=data
  )
  
  # Run simulation and extract components. Here, dp index must include the frozen parent.
  hmodel <- run_hdp_multichain(hmodel, dpindex=wparent_dp_index, initcc=n_sigs)
  hmodel <- hdp_extract_components(hmodel)

  return(hmodel)
}

#' Run HDP through parallel chains
run_hdp_multichain <- function(model, dpindex, initcc=2, burnin=2500, n_chains=4){
  chlist <- vector("list", n_chains)
  
  for (i in 1:n_chains){
    # Activating HDP
    activated_hdp <- dp_activate(
      model,
      dpindex=dpindex, # DPs to activate (including parents)
      initcc=initcc, # Number of clusters (signatures) to start with
      seed=i*1e3) 
    
    # Run gibbs sampler over activated nodes
    # Re-assigns cluster (signature) allocation of every data point
    chlist[[i]] <- hdp_posterior(
      activated_hdp,
      burnin=20000,
      n=100,
      space=200,
      cpiter=2,
      seed=i*1e5)
  }
  
  return(hdp_multi_chain(chlist))
}

hdp_plot_chains <- function(hdpsim, type){
  switch(type,
         {lapply(chains(hdpsim), plot_lik, bty="L")},
         {lapply(chains(hdpsim), plot_numcluster, bty="L")},
         {lapply(chains(hdpsim), plot_data_assigned, bty="L")}
  )
}
