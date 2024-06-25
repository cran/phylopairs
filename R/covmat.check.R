#' covmat.check
#'
#' @description Tets for validity of a covariance matrix based on four conditions: 
#'   symmetry, diagonal dominance, positive definiteness, and positive variance.
#'
#' @details A valid covariance matrix must be symmetric, diagonally dominant (largest values 
#'   in each row are on the diagonal), positive definite, and have positive variance. 
#'   \code{covmat.check} takes a matrix as input and tests for these four conditions.
#'
#' @usage covmat.check(mat)
#'
#' @param mat A putative covariance matrix. 
#'
#' @return A data.frame containing logical "TRUE" or "FALSE" for each condition. 
#'
#' @examples 
#' # Load sample covariance matrix
#' data(sim.cov.pairs)
#' # Test for validity
#' covmat.check(sim.cov.pairs)

#' @export
covmat.check = function(mat) {
  symmetric=all.equal(mat, t(mat))
  diag.dominant=rep(NA, nrow(mat))
  for(i in 1:nrow(mat)) {
    diag.dominant[i] = i %in% which(mat[i,]==max(mat[i,]))
  }
  diag.dominant=all(diag.dominant)
  pos.definite=all(eigen(mat)$values>0)
  pos.var=all(diag(mat)>=0)
  return(data.frame(symmetric,pos.var,diag.dominant,pos.definite))
}
