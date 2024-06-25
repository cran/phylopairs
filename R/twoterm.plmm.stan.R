#' twoterm.plmm.stan
#'
#' @description Fits the two-term plmm model of Castillo (2007) in a Bayesian
#'   framework. Bayesian sampling is conducted in the Stan software via the 'rstan'  
#'   package. Users supply vectors of observations and an ultrametric phylogenetic tree.
#'   Users can alter parameters for model-parameter prior distributions and Bayesian sampling
#'   settings. See details.
#'
#' @details The model introduced by Castillo (2007) for analyzing lineage-pair data is a
#'   version of a phylogenetic linear mixed model (plmm) in which there are two 
#'   random-effect terms, one for the 'species 1' and another for the 'species 2' in 
#'   every pair. The model is \code{Y = Xb + Z1u1 + Z1u2 + epsilon}, where u1 and u2 are 
#'   vectors of species-specific random effects that covary according to a phylogenetic 
#'   covariance matrix C. In this model, \code{u1 ~ N(0,physig2*C)}, \code{u2 ~ N(0,physig2*C)}, and
#'   \code{epsilon~N(0, randsig2*I)}, where randsig2 is the scaling parameter for identity 
#'   matrix I and physig2 is the scaling parameter for phylogenetic covariance 
#'   matrix C. See \code{twoterm_plmm_mats} for details on calculating 
#'   the Z1 and Z2 matrices. 
#'
#'   **Prior Distributions for Model Parameters**:
#'   The underlying stan models assume the following prior distributions for model parameters
#'   1. Regression coefficients: Gaussian prior (users can set prior mean and sd)
#'   2. `physig2`: lognormal prior (users can set prior mean and sd)
#'   3. `randsig2`: cauchy prior (users can set location and scale parameters of prior)

#' @usage twoterm.plmm.stan(des, y, sp1s, sp2s, tree, 
#'   iter=2000, chains=4, coef.u=0, coef.sd=10, physig2.u=-1, 
#'   physig2.sd=1, randsig2.loc=0, randsig2.sc=2.5, cores=4, ...)
#'
#' @param des A vector of predictor variable observations OR, in the case of multiple predictors, a matrix in which each column is a vector of observations of a given predictor. Function will add a column of 1s to make this a design matrix whose first column corresponds to the model intercept (unless such a column already exists).
#' @param y A vector of response variable observations. 
#' @param sp1s Character vector naming the "species 1" of every pair. IMPORTANT: name formatting must match that used in the tip.label of the tree
#' @param sp2s Character vector naming the "species 2" of every pair. IMPORTANT: name formatting must match that used in the tip.label of the tree
#' @param tree Phylogenetic tree (a 'phylo' object) for the species included in the analysis (regardless of whether they are species 1 or 2).
#' @param iter Number of iterations to run on each chain; defaults to 2000 (more are often necessary).
#' @param chains Number of Markov chains to run; defaults to 4.
#' @param cores Number of cores to use. 
#' @param ... additional arguments passed to \code{rstan::sampling}, including control parmeters (see rstan::sampling documentation)
#' @param coef.u Mean of the Gaussian prior for each preditor variable coefficient; defaults to 0.
#' @param coef.sd SD of the Gaussian prior for each preditor variable coefficient; defaults to 10.
#' @param physig2.u Mean of the prior distribution (lognormal) for the scale of the phylogenetic component of residual covariance; defaults to -1.
#' @param physig2.sd SD of the prior distribution (lognormal) for the scale of the phylogenetic component of residual covariance; defaults to 1.
#' @param randsig2.loc Location parameter of the prior distribution (cauchy) for the scale of the independent component of residual covariance ; defaults to 0.
#' @param randsig2.sc Scale parameter of the prior distribution (cauchy) for the scale of the independent component of residual covariance ; defaults to 2.5.
#'
#' @return A list containing two elements: (1) the posterior distribution of model parameters, and (2) the log-likelihood of the posteriors for use in downstream analyses (e.g. the calculation of model fitting metrics like loo or waic)
#'
#' @examples 
#' \donttest{
#' # Load a dataset simulated with a non-independent response observations
#' data(data3)
#' # Also load the simulated tree that was used to generate those pairs
#' data(sim.tree1)
#' # Fit an OLS model
#' result1 = linreg.stan(des=data3[,3], y=data3[,4], cores=2)
#' # Fit the twoterm.plmm.stan model
#' result2 = twoterm.plmm.stan(des=data3[,3], y=data3[,4], sp1s=data3[,1], 
#'   sp2s=data3[,2], tree=sim.tree1, cores=2)
#' # Compare posterior parameter estimates
#' result1[[1]]
#' result2[[1]]
#' # Compare the fit of the two models via loo and waic
#' loo1 = loo::loo(result1[[2]])
#' loo2 = loo::loo(result2[[2]])
#' waic1 = loo::waic(result1[[2]])
#' waic2 = loo::waic(result2[[2]])
#' loo1
#' loo2
#' waic1
#' waic2
#' loo::loo_compare(loo1, loo2)
#' loo::loo_compare(waic1, waic2)
#'}

#' @references
#' Castillo, D. M. (2007). Factors contributing to the accumulation of reproductive isolation: A mixed model approach. Ecology and Evolution 7:5808-5820. doi.org/10.1002/ece3.3093

#' @export
twoterm.plmm.stan = function(des, y, sp1s, sp2s, tree, iter=2000, chains=4, coef.u=0, 
  coef.sd=10, physig2.u=-1, physig2.sd=1, randsig2.loc=0, randsig2.sc=2.5, cores=4, ...) {
  # Add an intercept column to design matrix if it's not already present
  if(!all(as.matrix(des)[,1]==1)) des=cbind(rep(1, length(y)), des)
  # Ensure length of species 1s and species 2 vectors match
  if(length(sp1s)!=length(sp2s)) stop("length of 'sp1s' and 'sp2s' vectors must match")
  # Ensure species 1s and species 2s are actually in the tree
  if (any(!sp1s %in% tree$tip.label) | any(!sp2s %in% tree$tip.label)) stop("sps1s and/or sps2s are missing from the tree \n check to ensure the formatting of the names in sp1s and sp2s matches that used in the tree tip.label component")
  # Calculate necessary identity and covariance matrices
  mats = twoterm.plmm.mats(sp.pairs=cbind(sp1s, sp2s), tree=tree)
  # Ensure covariance matrices are valid
  if(suppressWarnings(!all(sapply(mats[3:4], covmat.check)))) stop("at least one invalid covariance matrix calculated; see twoterm_plmm_mats function ")
  # Prep stan data
  stan.dat = list(N=length(y), M=length(unique(sp1s)), P=length(unique(sp2s)), K=ncol(des), Y=y, X=des, C1=mats[[3]], C2=mats[[4]], Z1=mats[[1]], Z2=mats[[2]], 
    coef_mean=coef.u, coef_sd=coef.sd, sig2_mean=physig2.u, sig2_sd=physig2.sd, sig2_loc=randsig2.loc, sig2_sc=randsig2.sc)
  # Fit the model
  fit = sampling(object = stanmodels$twoterm_plmm, data = stan.dat, iter = iter, chains = chains, cores=cores, ...)
  # Get the parameters from the summary
  pars = round(rstan::summary(fit)$summary,2)[c(1:4,dim(fit)[3]),]
  ll = loo::extract_log_lik(fit, parameter_name = "loglik")
  # Package and return the results list
  res = list(pars=pars, ll=ll)
  return(res)
}
