#' taxapair.vcv
#'
#' @description Calculates an unscaled lineage-pair covariance matrix for use in 
#'   downstream analyses.
#'
#' @details Just as the traits of different species are not independent due to varying 
#'   amounts of shared evolutionary history, so too are pairwise-defined traits like
#'   'strength of RI' and 'range overlap' not independent among related taxonomic pairs.
#'   The function `taxapair.vcv()` calculates the expected covariance structure for 
#'   lineage-pair datasets given the taxa in each pair and a phylogenetic tree
#'   (containing the taxa that appear in the dataset). The exact structure of this 
#'   covariance depends on the underlying model by which phylogenetic signal among
#'   taxa is expected to translate into non-independence among taxonomic pairs 
#'   (Anderson et al. *in review*). 
#'
#'   Let X be some underlying continuous biological character (or set of characters) 
#'   that is defined for each taxon and that has phylogenetic signal. Users can choose
#'   one of four models by which signal in X generates covariance among pairs in a 
#'   lineage-pair trait:
#'
#'   1. sq.diff (default) -- the lineage-pair trait is influenced by the square of the
#'       difference between the two taxa in X.
#'   2. sq.sum -- the lineage-pair trait is linearly by the squared sum of 
#'       the value of X in each taxon.
#'   3. simple.sum -- the lineage-pair trait is influenced by the sum of the 
#'       the value of X in each taxon. IMPORTANT: this model tends to result in 
#'       a singular matrix. 
#'   4. product -- the lineage-pair trait is influenced by the product of the 
#'       the value of X in each taxon.
#'   5. abs.diff -- the lineage-pair trait is influenced by the absolute difference
#'       difference between the two taxa in X.
#'
#'   In each case, it is assumed that there is a linear relationship between 
#'     the value calculated in each model and the response variable.
#'
#'   Note that simple.sum and abs.diff typically result in invalid covariance matrices, so use caution.
#'
#' @usage taxapair.vcv(sp.pairs, tree, spnames=FALSE, dec=2, model="sq.diff", 
#'   regularize=FALSE, regparam=NULL)
#'
#' @param sp.pairs A table (matrix or data.frame) in which the first column contains the names of 'species 1' and the second column contains the names of 'species 2'. Names must be in the same format used in the phylogenetic tree. 
#' @param tree An ultrametric phylogenetic tree ('phylo' object) containing species that appear in the dataset (as either a species 1 or species 2 or both). 
#' @param spnames Logical determining whether to use species names as row and column labels for matrix. If TRUE, then the name of a row or column will be in the form "species1_species2". If FALSE, names are formed from the indices of the species names in the tree "indexOfSpecies1_indexOfSpecies2". Defaults to FALSE.
#' @param dec Number of decimal places to round the values in the matrix; defaults to 2. 
#' @param model One of 'sq.diff', 'sq.sum', 'simple.sum', 'product', or 'abs.diff'. Defaults to 'sq.diff'. See details. 
#' @param regularize Logical indicating whether regularization should be used if resulting matrix is numerically singular; defaults to FALSE. If TRUE, regularization is conducted by adding a small value to the diagonal of the matrix. By default, this value is equal to 1% of the median value on the diagonal, which is continually added until matrix is no longer numerically singular. 
#' @param regparam Custom regularization parameter. Instead of adding the default 1% of the median diagonal value to the diagonal, it will add regparam * that median value. IMPORTANT: it only does this once and does not continue adding the number until the matrix is non-singular. For 2%, write 0.02; for 10%, write 0.10, and so on. 
#' @return A lineage-pair covariance matrix. 

#' @importFrom stats median
#'
#' @examples
#' library(ape)
#' # Load simulated dataset and tree
#' data(data1)
#' data(sim.tree1)
#' # Calculate the lineage-pair covariance matrix
#' linpair.mat = taxapair.vcv(sp.pairs=data1[,1:2], tree=sim.tree1)
#' dim(linpair.mat)
#' # Check the validity of the matrix
#' covmat.check(linpair.mat)

#' @references
#' Anderson, S. A. S., et al. *In review*. The comparative analysis of lineage-pair data.

# requires that sp.pairs has the names of species 1 in the first column and of species 2 in the second column
#' @export
taxapair.vcv = function(sp.pairs, tree, spnames=FALSE, dec=2, model="sq.diff", regularize=FALSE, regparam=NULL) {
  if(is.null(tree$edge.length)) stop("Tree does not have branch lengths; no vc matrix can be returned")
  if(!is.matrix(sp.pairs)) sp.pairs=as.matrix(sp.pairs)
  if(!model%in%c("sq.diff", "sq.sum", "simple.sum","product","abs.diff")) stop("Invalid model chosen for `taxapair.vcv`; see 'details'")
  if(!ape::is.ultrametric(tree, 2)) stop("Tree is not ultrametric; this method requires an ultrametric tree")
  # get the vc matrix for the tree
  vcmat = ape::vcv(tree)
  # make an empty sp.pair vc matrix and label it using edge numbers or taxa names from the tree
  pvcv = matrix(0, nrow=nrow(sp.pairs), ncol=nrow(sp.pairs))
  if(!spnames) colnames(pvcv) = rownames(pvcv) = apply(sp.pairs, 1, 
    function(x) paste(which(tree$tip.label==x[1]), which(tree$tip.label==x[2]),sep="_"))
  if(spnames) colnames(pvcv) = rownames(pvcv) = paste(sp.pairs[,1], sp.pairs[,2],sep="_")
  # for each row, calculate the covariance between pairs
  # note that the species labels in the dataset must match the species labels in the tree
  for(i in 1:nrow(pvcv)) {
    sp1.1 = which(tree$tip.label==sp.pairs[i,][1])
    sp1.2 = which(tree$tip.label==sp.pairs[i,][2])
    for(j in i:ncol(pvcv)) {
      if(model=="sq.diff") {
        if(j==i) {
          pvcv[i,j] = 2*(2*vcmat[1,1] - 2*vcmat[sp1.1,sp1.2])^2
        } else {
          sp2.1 = which(tree$tip.label==sp.pairs[j,][1])
          sp2.2 = which(tree$tip.label==sp.pairs[j,][2])
          pvcv[i,j] = sum(2*vcmat[sp1.1, sp2.1]^2, 2*vcmat[sp1.1, sp2.2]^2, 2*vcmat[sp1.2,sp2.1]^2, 2*vcmat[sp1.2, sp2.2]^2, 
            -4*vcmat[sp1.1, sp2.1]*vcmat[sp1.1, sp2.2], -4*vcmat[sp1.1, sp2.1]*vcmat[sp1.2,sp2.1], -4*vcmat[sp1.1, sp2.2]*vcmat[sp1.2,sp2.2],
            -4*vcmat[sp1.2, sp2.1]*vcmat[sp1.2,sp2.2], 4*vcmat[sp1.1,sp2.1]*vcmat[sp1.2, sp2.2], 4*vcmat[sp1.1, sp2.2]*vcmat[sp1.2, sp2.1])
        }
      }
      if(model=="sq.sum") {
        if(j==i) {
          pvcv[i,j] = 2*(2*vcmat[1,1] + 2*vcmat[sp1.1,sp1.2])^2
        } else {
          sp2.1 = which(tree$tip.label==sp.pairs[j,][1])
          sp2.2 = which(tree$tip.label==sp.pairs[j,][2])
          pvcv[i,j] = sum(2*vcmat[sp1.1, sp2.1]^2, 2*vcmat[sp1.1, sp2.2]^2, 2*vcmat[sp1.2,sp2.1]^2, 2*vcmat[sp1.2, sp2.2]^2, 
            4*vcmat[sp1.1, sp2.1]*vcmat[sp1.1, sp2.2], 4*vcmat[sp1.1, sp2.1]*vcmat[sp1.2,sp2.1], 4*vcmat[sp1.1, sp2.2]*vcmat[sp1.2,sp2.2],
            4*vcmat[sp1.2, sp2.1]*vcmat[sp1.2,sp2.2], 4*vcmat[sp1.1,sp2.1]*vcmat[sp1.2, sp2.2], 4*vcmat[sp1.1, sp2.2]*vcmat[sp1.2, sp2.1])
        }
      }
      if(model=="simple.sum") {
        sp2.1 = which(tree$tip.label==sp.pairs[j,][1])
        sp2.2 = which(tree$tip.label==sp.pairs[j,][2])
        pvcv[i,j] = sum(vcmat[sp1.1, sp2.1], vcmat[sp1.1, sp2.2], vcmat[sp1.2, sp2.1], vcmat[sp1.2, sp2.2])
      }
      if(model=="product") {
        sp2.1 = which(tree$tip.label==sp.pairs[j,][1])
        sp2.2 = which(tree$tip.label==sp.pairs[j,][2])
        pvcv[i,j] = vcmat[sp1.1,sp2.1]*vcmat[sp1.2,sp2.2] + vcmat[sp1.1,sp2.2]*vcmat[sp1.2,sp2.1]
      }
      if(model=="abs.diff") {
        sp2.1 = which(tree$tip.label==sp.pairs[j,][1])
        sp2.2 = which(tree$tip.label==sp.pairs[j,][2])
        pvcv[i,j] = (4/pi)* (vcmat[sp1.1,sp2.1] - vcmat[sp1.1,sp2.2] - vcmat[sp1.2,sp2.1] + vcmat[sp1.2,sp2.2])
      }
    }
  }
  # complete the matrix
  for(i in 2:nrow(pvcv)) pvcv[i, 1:(i-1)] = pvcv[1:(i-1),i]
  # regularize if required
  if(regularize && det(pvcv) <= .Machine$double.eps) {
    if(is.null(regparam)) {
      while(det(pvcv) <= .Machine$double.eps) {
        pvcv = pvcv + diag(median(diag(pvcv))*0.01, nrow(pvcv))
      }
    } else {    
      pvcv = pvcv + diag(median(diag(pvcv))*regparam, nrow(pvcv))
    }
  }
  if(regularize && det(pvcv) <= .Machine$double.eps) {
    stop("You provided a custom 'regparam' but the determinant of matrix is still numerically singular;
      try a new regularization parameter")
  }
  if(!all(eigen(pvcv)$values>0)) {
    message("Important: resulting matrix is not positive semi-definite; not a valid cov. matrix. Try regularization or another model")
  }
  return(pvcv)
}