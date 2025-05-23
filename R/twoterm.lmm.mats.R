#' twoterm.lmm.mats
#'
#' @description Calculates the four matrices required for fitting the two-term lmm model
#'   of Castillo (2007). See details.
#'
#' @details Castillo (2007) introduced a framework for analyzing lineage-pair data via an
#'   extension of a phylogenetic linear mixed model (plmm) in which there are two 
#'   random-effect terms, one for the 'species 1' and another for the 'species 2' in 
#'   every pair. The model is \code{Y = Xb + Z1u1 + Z1u2 + epsilon}, where b is the vector of 
#'   coefficients, X is a design matrix, u1 and u2 are vectors of species-specific random 
#'   effects that covary according to a phylogenetic covariance matrix C (\code{u1 ~ N(0,C)} 
#'   and \code{u2 ~ N(0,C)}), Z1 and Z2 are design matrices that map the species-specific 
#'   effects to the correct species in each pair, and epsilon is residual error (\code{epsilon~N(0,C)}). 
#'   The C matrix in the model is scaled by parameters that do not come into play here. 
#'
#'   This function takes a table containng columsn for species 1 and species 2 for every
#'   pair and a  phylogenetic tree and returns the matrices Z1 and Z2 as well as pruned 
#'   phylogenetic covariance matrices for  u1 and u2. This pruning is sometimes required  
#'   because not all species found in the  dataset will appear as both 'species 1' and  
#'   'species 2'. Z1 and Z2 will therefore have different sizes and u1 and u2 require  
#'   different covariance matrices. Note that the matrices themselves are pruned, not 
#'   the tree from which they are derived, as the latter could result in the covariance 
#'   between two species being different for the 'species 1' and 'species 2' random effects. 
#'
#' @usage twoterm.lmm.mats(sp.pairs, tree)
#'
#' @param sp.pairs A table (matrix or data.frame) in which the first column contains the names of 'species 1' and the second column contains the names of 'species 2'. Names must be in the same format used in the phylogenetic tree. 
#' @param tree An ultrametric phylogenetic tree ('phylo' object) describing the relationships among the species that appear in the dataset (as either a species 1 or species 2 or both). 
#'
#' @return A list of four matrices: Z1, Z2, cov1, and cov2, where the latter two describe the covariance among the random effects in u1 and u2.
#' @importFrom phytools pbtree

#' @examples 
#' # Simulate a tree
#' lin.tree = phytools::pbtree(n=20)
#' # Generate lineage pairs as the pairwise combinations of species in the tree
#' lin.pairs = data.frame(t(combn(lin.tree$tip.label,2))); colnames(lin.pairs)=c("sp1", "sp2")
#' # Calculate the matrices
#' mats = twoterm.lmm.mats(sp.pairs=lin.pairs, tree=lin.tree)
#' # Check structure of design matrices
#' sapply(mats, dim)
#' head(mats$Z1[,1:5])
#' head(mats$Z2[,1:5])
#' head(mats$cov2[,1:5])
#' head(mats$cov2[,1:5])
#' # Ensure covariance matrices are valid covariance matrices
#' sapply(mats[3:4], covmat.check)

#' @references
#' Castillo, D. M. (2007). Factors contributing to the accumulation of reproductive isolation: A mixed model approach. Ecology and Evolution 7:5808-5820. doi.org/10.1002/ece3.3093

#' @export
twoterm.lmm.mats = function(sp.pairs, tree) {
  # will orient the columns in Z to line up with the phylogenetic covariance matrix
  vc.mat = ape::vcv(tree)
  Z1=Z2=matrix(0, nrow(sp.pairs), ncol(vc.mat))
  colnames(Z1)=colnames(Z2)=colnames(vc.mat)
  for(i in 1:nrow(sp.pairs)) {
    Z1[i,which(colnames(Z1)==sp.pairs[i,1])]=1
    Z2[i,which(colnames(Z2)==sp.pairs[i,2])]=1
  }
  ind1=colnames(vc.mat)%in%unique(sp.pairs[,1])
  ind2=colnames(vc.mat)%in%unique(sp.pairs[,2])
  Z1=Z1[,ind1]
  Z2=Z2[,ind2]
  cov1=vc.mat[ind1,ind1]
  cov2=vc.mat[ind2,ind2]
  return(list("Z1"=Z1, "Z2"=Z2, "cov1"=cov1, "cov2"=cov2))
}
