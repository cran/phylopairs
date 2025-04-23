#' The 'phylopairs' package.
#'
#' @description `phylopairs` provides tools for conducting comparative analyses of lineage-pair traits in a phylogenetically informed context. Regression-style analsyses of pairwise-defined traits like "strength of RI" or "diet niche overlap" collected for numerous pairs of related taxa can yield important insights to biologists, but it is important to recognize that such data are not independent. `phylopairs` provides a function, \code{taxapair.vcv()}, that calculates the expected covariance structure of a pairwise-defined trait given (1) a table of pairs, (2) a phylogeny, and (3) a chosen model (as described in Anderson et al. *in review*). This covariance structure can be used in any number of downstream analyses. `phylopairs` provides tools for a few such analyses, including linear regression, linear mixed models, generalized least squares, and two forms of beta regression for use with bounded response variables. All analyses are conducted in a Bayesian framework using built-in `Stan` software programs interfaced through the `rstan` package. 

#' @name phylopairs-package
#' @aliases phylopairs
#' @useDynLib phylopairs, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @references
#' 		Stan Development Team (NA). RStan: the R interface to Stan. R package version 
#'      2.32.6. https://mc-stan.org
#'
#' 		Anderson et al. *in review*. The comparative analysis of lineage-pair data. 
"_PACKAGE"
NULL
