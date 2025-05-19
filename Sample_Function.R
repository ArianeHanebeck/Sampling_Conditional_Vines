#################################
# Functions needed for sampling #
#################################

normalizeRVineMatrix = function(RVM) {
  
  oldOrder = diag(RVM$Matrix)
  Matrix = reorderRVineMatrix(RVM$Matrix)
  
  return(RVineMatrix(Matrix,
                     RVM$family,
                     RVM$par,
                     RVM$par2,
                     names = rev(RVM$names[oldOrder]),
                     check.pars = FALSE))
}

reorderRVineMatrix = function(Matrix, oldOrder = NULL) {
  
  if (length(oldOrder) == 0) {
    oldOrder = diag(Matrix)
  }
  O = apply(t(1:nrow(Matrix)), 2, "==", Matrix)
  
  for (i in 1:nrow(Matrix)) {
    Matrix[O[, oldOrder[i]]] = nrow(Matrix) - i + 1
  }
  return(Matrix)
}

#####################
# Sampling function #
#####################

sample_from_conditional = function(N, STAN, STAN2, RVM, indexcon, ucon, burnin, thin, extreme = 0, ...){
  
  # Default values for burnin and thin
  if(missing(burnin)){
    burnin = 1000
  }
  if(missing(thin)){
    thin = 10
  }
  
  # Extracting the dimensions
  d = length(indexcon)
  d1 = length(which(indexcon == FALSE))
  d2 = d-d1
  
  # Loading Stan
  if (d1 == 1){
    STAN = STAN
  }else{
    STAN = STAN2
  }
  
  # Normalizing RVineMatrix - needed in VineCopula package to compute log-lik
  dataflip = 0
  o = diag(RVM$Matrix)
  if (any(o != length(o):1)) {
    oldRVM = RVM
    RVM = normalizeRVineMatrix(RVM)
    dataflip = 1
  }

  # From the VineCopula package: the different parameters/inputs we need 
  # to compute log-lik in Stan - inserted into Stan as data
  T = 1
  w1 = as.vector(RVM$family)
  w1[is.na(w1)] = 0
  th = as.vector(RVM$par)
  th[is.na(th)] = 0
  th2 = as.vector(RVM$par2)
  th2[is.na(th2)] = 0
  condirect = as.vector(as.numeric(RVM$CondDistr$direct))
  conindirect = as.vector(as.numeric(RVM$CondDistr$indirect))
  maxmat = as.vector(RVM$MaxMat)
  matri = as.vector(RVM$Matrix)
  matri[is.na(matri)] = 0
  maxmat[is.na(maxmat)] = 0
  condirect[is.na(condirect)] = 0
  conindirect[is.na(conindirect)] = 0
  
  # The parameters for the Stan-program
  It = N*thin
  iter_u_2 = It+burnin
  burnin_u_2 = burnin
  chains_u_2 = 4
  adapt_delta_u_2 = 0.8
  max_treedepth_u_2 = 10
  
  # Sample from Stan
  data_stan_u_2 = list(T = T, dataflip = as.integer(dataflip), d = d, 
                       o = as.integer(o), d1 = d1, d2 = d2, family = as.integer(w1),
                       maxmat = as.integer(maxmat),matri = as.integer(matri),
                       condirect = as.integer(condirect),
                       conindirect = as.integer(conindirect),
                       par = as.double(th),par2 = as.double(th2),
                       indexcon = indexcon, ucon = ucon)
  
  init_list_u_2 = list(list(ucalculate = rep(0.5,d1)),list(ucalculate = rep(0.5,d1)),list(ucalculate = rep(0.5,d1)),list(ucalculate = rep(0.5,d1)))
  
  fit_u_2 = sampling(STAN,iter = iter_u_2,warmup = burnin_u_2,chains = chains_u_2,
                   data = data_stan_u_2,init = init_list_u_2, 
                   control = list(adapt_delta = adapt_delta_u_2, 
                                max_treedepth = max_treedepth_u_2),
                   thin = thin, ...)
  return(fit_u_2)
}
