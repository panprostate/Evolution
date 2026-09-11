
#' Dirichlet-Multinomial Stan Model
#' Data:
#'  S - number of signatures
#'  N - number of categories in signature (e.g. 96 for SBS)
#'     In stan, the N is implicit given the data vector Y
#'  Y - mutation counts in epoch
#'  alpha dirichlet prior for each signature weight
#'  signatures - reference profiles for active signatures
dm_model <- "
  data {
    int<lower = 0> S;
    int<lower = 0> N;
    int<lower = 0> Y[N];
    vector<lower=0>[S] alpha;
    matrix[N, S] signatures;
  }
  parameters {
    simplex[S] weights;
  }
  model {
    weights ~ dirichlet(alpha);
    Y ~ multinomial(signatures * weights);
  }
"

# Ensure model has been compiled and model RDS is available
model_file = system.file("extdata", "dm_model.stan", package="SigTimeR")
sm = rstan::stan_model(file=model_file, auto_write = T)


#' Time signatures with dirichlet multinomial
sigtime_dm <- function(counts, priors){
  
  # Get result for epoch_1
  epoch1_dm = dm_sol(counts[1,], priors)
  
  # Get result for ecoch_2
  epoch2_dm = dm_sol(counts[2,], priors)
  
  return(list(epoch1=epoch1_dm, epoch2=epoch2_dm))
}

dm_sol <- function(data, sigs){

  n_sigs = dim(sigs)[2]
  n_cats = length(data)
  uniform_alpha_prior = rep_len(1, n_sigs)/n_sigs
  
  dm_data = list(
    S=n_sigs,
    N=n_cats,
    Y=data,
    alpha=uniform_alpha_prior,
    signatures=sigs
  )
  
  dm_sim <- rstan::stan(
    file=model_file,
    data=dm_data,
    chains=4,
    iter=20000*2
  )
  
  #if(return_sim){return(dm_sim)}
  return(list(dm=dm_sim, data=data, sigs=sigs))
}

sigtime_dm_extract <- function(dm_raw, id=NA){
  
  signames = colnames(dm_raw$epoch1$sigs)
  e1_dm = data.frame(rstan::summary(dm_raw$epoch1$dm)[[1]])
  e2_dm = data.frame(rstan::summary(dm_raw$epoch2$dm)[[1]])

  estimate = lapply(
    seq_along(signames),
    function(x){
      e1 = dplyr::select(e1_dm, "mean", lower=X2.5., upper=X97.5.)[x,]
      e1 = dplyr::mutate(e1, epoch="epoch_1", method="dm", id=id, signature=signames[x])
      e2 = dplyr::select(e2_dm, "mean", lower=X2.5., upper=X97.5.)[x,]
      e2 = dplyr::mutate(e2, epoch="epoch_2", method="dm", id=id, signature=signames[x])
      ed = tidyr::pivot_longer(rbind(e1, e2), cols=c("mean","lower","upper"), names_to="info", values_to="prop")
      return(ed[, c("info", "epoch", "id", "method", "signature", "prop")])
    }
  )
  estimate = Reduce(rbind, estimate)
  return(estimate)
  
}