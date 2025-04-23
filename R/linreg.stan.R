#' linreg.stan

#' @description Fits one of three models of linear regression to a dataset in a Bayesian
#'   framework. Bayesian sampling is conducted in the Stan software via the 'rstan'  
#'   package. Users supply vectors of observations and, unless an "ols" model is chosen, 
#'   a covariance matrix. Users can alter parameters for model-parameter prior distributions 
#'   and Bayesian sampling settings. See details. 

#' @details The models that can be chosen are: 
#'   1. "ols": basic ordinary least-squares regression.
#'
#'      \code{Y = XB + epsilon}; where \code{epsilon~N(0,randsig2*I)} and randsig2 is the scaling parameter for identity matrix I.
#'
#'   2. "pgls": phylogenetic generalized least-squares regression in which residuals 
#'      covary according to covariance matrix C. This is a lineage-pair covariance 
#'      matrix (if analyzing lineage-pair data) or a phylogenetic vcv matrix 
#'      (if analyzing species data).
#'
#'       \code{Y = XB + epsilon}; where \code{epsilon~N(0,physig2*C)} and physig2 is the scaling parameter for covariance matrix C.
#'
#'   3. "lp.gls" a mixed-effects version of pgls in which the residuals have both 
#'      uncorrelated and correlated components, where the latter are structured according
#'      to covariance matrix C. This is the 'lp.gls' model from Anderson et al. 2025 (In Revision).
#'
#'      \code{Y = XB + epsilon + u}; where \code{epsilon~N(0,randsig2*I)} and \code{u~N(0,physig2*C)}.
#'
#'
#'   **NOTE**: For models 2 and 3, the user must supply a covariance matrix.
#'
#'   **Prior Distributions for Regression Model Parameters**:
#'   The underlying stan models assume the following prior distributions for model
#'   parameters
#'   1. Regression coefficients: Gaussian prior (users can set prior mean and sd).
#'   2. `physig2`: lognormal prior (users can set prior mean and sd).
#'   3. `randsig2`: cauchy prior (users can set location and scale parameters of prior).

#' @usage linreg.stan(des, y, model="ols", covmat=NULL, iter=2000, 
#'   chains=4, cores=4, coef.u=0, coef.sd=10, physig2.u=-1, physig2.sd=1, 
#'   randsig2.loc=0, randsig2.sc=2.5, ...)
#'   

#' @param des A vector of predictor variable observations OR, in the case of multiple predictors, a matrix in which each column is a vector of observations of a given predictor. Function will add a column of 1s to make this a design matrix whose first column corresponds to the model intercept (unless such a column already exists).
#' @param y A vector of response variable observations. 
#' @param model Choice of "ols", "pgls", or "lp.gls"; defaults to ols.
#' @param covmat Covariance matrix for model residuals (a lineage-pair covariance matrix if analyzing lineage-pair data or a phylogenetic vcv matrix if analyzing species data).
#' @param iter Number of iterations to run on each chain; defaults to 2000 (more are often necessary).
#' @param chains Number of Markov chains to run; defaults to 4.
#' @param cores Number of cores to use. 
#' @param ... additional arguments passed to \code{rstan::sampling()}, including control parmeters (see \code{rstan::sampling()} documentation)
#' @param coef.u Mean of the Gaussian prior for each preditor variable coefficient; defaults to 0.
#' @param coef.sd SD of the Gaussian prior for each preditor variable coefficient; defaults to 10.
#' @param physig2.u Mean of the prior distribution (lognormal) for the scale of the phylogenetic component of residual covariance; defaults to -1.
#' @param physig2.sd SD of the prior distribution (lognormal) for the scale of the phylogenetic component of residual covariance; defaults to 1.
#' @param randsig2.loc Location parameter of the prior distribution (cauchy) for the scale of the independent component of residual covariance ; defaults to 0.
#' @param randsig2.sc Scale parameter of the prior distribution (cauchy) for the scale of the independent component of residual covariance ; defaults to 2.5.

#' @return A list containing two elements: (1) the posterior distribution of model parameters, and (2) the log-likelihood of the posteriors for use in downstream analyses (e.g. the calculation of model fitting metrics like loo or waic). For interpreting model parameters, note that \code{Coef[1]} is the intercept and \code{Coef[2]}, \code{Coef[3]}, ... , \code{Coef[N]} are regression coefficient for the 1st-(N-1)th predictor variables. 

#' @examples 
#' \donttest{
#' ## Fit regression models with and without covariance matrix
#' # Note: data were simulated with Coef[1]=1 (intercept), Coef[2]=0.8 (slope)
#' # Load a data simulated with a non-independent response observations
#' data(data3)
#' # Also load the lineage-pair covariance matrix that arose from those simulations
#' data(sim.cov.pairs)
#' # Fit an OLS model
#' result1 = linreg.stan(des=data3[,3], y=data3[,4], cores=2)
#' # Fit an pgls model
#' result2 = linreg.stan(des=data3[,3], y=data3[,4], model="pgls", covmat=sim.cov.pairs, cores=2)
#'
#' # Compare posterior parameter estimates
#' result1[[1]]
#' result2[[1]]
#'
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
#'
#' # Extend the comparison by fitting a lp.gls model
#' result3 = linreg.stan(des=data3[,3], y=data3[,4], model="lp.gls", 
#'  covmat=sim.cov.pairs, cores=2)
#'
#' # Compare posterior parameter estimates
#' result1[[1]]
#' result2[[1]]
#' result3[[1]]
#'
#' # Compare the fit of the three models via loo and waic
#' loo3 = loo::loo(result3[[2]])
#' waic3 = loo::waic(result3[[2]])
#' loo::loo_compare(loo1, loo2, loo3)
#' loo::loo_compare(waic1, waic2, waic3)
#'}

#' @references
#' Anderson, S. A. S., et al. *In review*. The comparative analysis of lineage-pair data.

#' @export
linreg.stan = function(des, y, model="ols", covmat=NULL, iter=2000, chains=4, cores=4, coef.u=0, 
  coef.sd=10, physig2.u=-1, physig2.sd=1, randsig2.loc=0, randsig2.sc=2.5, ...) {
  # Check if covariance matrix supplied appropriately
  if(is.null(covmat) & model %in% c("pgls", "lp.gls")) stop("User must supply a covariance matrix for this model")
  if(!is.null(covmat) & model == "ols") stop("covariance matrix supplied for 'ols' model. Choose alternative model or don't supply covmat")
  # Add an intercept column to design matrix if it's not already present
  if(!all(as.matrix(des)[,1]==1)) des=cbind(rep(1, length(y)), des)
  # Code the model
  model.type = switch(model, ols = 1, pgls=2, lp.gls=3, stop("'model' given to 'linreg.stan' must be one of 'ols', 'pgls', or 'lp.gls'"))
  # Prep stan data
  stan.dat = list(N=length(y), K=ncol(des), Y=y, X=des, coef_mean=coef.u, coef_sd=coef.sd, sig2_mean=physig2.u, sig2_sd=physig2.sd, sig2_loc=randsig2.loc, sig2_sc=randsig2.sc, model_type=model.type)
  # Add covariance matrix if necessary
  stan.dat$Cp = if(model=="ols") diag(length(y)) else covmat
   # Fit the model
  fit = sampling(object = stanmodels$linreg, data = stan.dat, iter = iter, chains = chains, cores=cores, ...)
  # Get the parameters from the summary
  m = ncol(des)
  pars = if(!is.null(covmat)) rstan::summary(fit)$summary[c(1:(m+2),dim(fit)[3]),] else pars = rstan::summary(fit)$summary[c(1:(m+1),dim(fit)[3]),]
  ll = loo::extract_log_lik(fit, parameter_name = "loglik")
  # Package and return the results list
  res = list(pars=round(pars,2), ll=ll)
  return(res)
}