#' pair.age.in.tree
#'
#' @description Takes a lineage-pair dataset and a phylogenetic tree and returns the 'age' of the node representing the MRCA for each pair. 'Age' is measured from the present and is in whatever units are represented by branch lengths of the tree (which is time in a dated phylogeny).
#'
#' @details Function is a wrapper for \code{ape::node.depth.edgelength()}. That function returns the depths of all the nodes in a tree as measured from the root. \code{phylopairs::node.age()} converts these values to depth as measured from the present and returns the values corresponding to each pair in a user-supplied lineage-pair dataset. 
#'
#' @usage pair.age.in.tree(dataset, tree, taxacolumns)
#'
#' @param dataset A data.frame in which each row corresponds to a pair in the dataset. Must contain two columns of taxa names (one for each taxon in every pair), and the taxa names must be in same format as that used in the tree. 
#' @param tree An ultrametric phylogenetic tree ('phylo' object) containing the species that appear in at least one pair in the dataset. Names must be in the same format as those used in 'dataset'.
#' @param taxacolumns Character vector containing the column names for the two columns containing species names (e.g. c("sp1", "sp2"))
#'
#' @return A numeric vector of node ages ordered to match the rows in 'dataset'.
#'
#' @examples 
#' # Load a dataset and a tree
#' data(data1)
#' data(sim.tree1)
#'
#' # Find the node.ages and add as a column to the dataset
#' data1$age = pair.age.in.tree(dataset=data1, tree=sim.tree1, taxacolumns=c("sp1","sp2"))
#' head(data1)
#'
#' \donttest{
#' # Plot tree with axis in units of branch lengths and perform visual check that dates are correct
#' library(ape)
#' plot(sim.tree1)
#' axisPhylo()
#' }

#' @export
pair.age.in.tree = function(dataset, tree, taxacolumns) {
  # ensure tree is ultrametric
  if(!ape::is.ultrametric(tree, 2)) stop("'tree' must be ultrametric")
  # ensure the dataset has the right column names
  if(all(c(taxacolumns) %in% colnames(dataset))==F) {
    stop("one or more required columns is missing in the dataset, or there is an incorrect column name")
  }
  # Calculate ages for the nodes in the tree
  nodes = round(max(ape::node.depth.edgelength(tree)) - ape::node.depth.edgelength(tree),4)
  # Index the nodes ages corresponding to the mrcas for each pair
  mrcas = apply(dataset[,taxacolumns], 1, function(x) ape::getMRCA(tree, tip=c(x[1], x[2])))
  res = nodes[mrcas]
  return(res)
}