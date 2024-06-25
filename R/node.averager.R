#' node.averager
#'
#' @description Calculate weighted or unweighted node-averaged values of a lineage-pair trait. 
#'
#' @details \code{node.averager()} takes a lineage-pair dataset and a phylogenetic tree and 
#'   returns the average value of the pairwise-defined trait at each node. It calculates, 
#'   at each node in the tree, the average value of a pairwise-defined trait for all pairs 
#'   whose species span the node. 
#'
#'   The simple or 'unweighted' average is calculated as introduced by Coyne 
#'   and Orr (1989). The 'weighted' node averaging procedure was introduced by  
#'   Fitzpatrick (2002) and discussed in Fitzpatrick and Turelli (2006). In weighted 
#'   averaging, the trait values for the pairs spanning a node are first halved K-1 
#'   times, where K is the number of nodes between the species in a pair. These halved 
#'   values are then summed to get arrive at a weighted average for a node.
#'
#'   **Note:** For datasets containing many individual species, the best available
#'   tree might be very large. By default, \code{node.averager()} prunes the tree to 
#'   contain just the species represented in the dataset. Beware that pruning can 
#'   affect both the number of nodes and the node averaged values (for example by 
#'   altering the number of nodes between pairs when calculating weighted node averages).
#'
#' @usage node.averager(dataset, tree, taxacolumns, varb, 
#'   weighted=FALSE, prune=TRUE, av=TRUE)
#'
#' @param dataset A data.frame in which each row corresponds to a pair in the dataset. Must contain two columns of taxa names (one for each taxon in every pair) and at least one data column with the pairwise-defined trait that is to be averaged. Taxa names must be in same format as that used in the tree. 
#' @param tree An ultrametric phylogenetic tree ('phylo' object) containing the species that appear in at least one pair in the dataset. Names must be in the same format as those used in 'dataset'.
#' @param taxacolumns Character vector containing the column names for the two columns containing species names (e.g. c("sp1", "sp2"))
#' @param varb The variable to be averaged (e.g. "RI", "range_overlap", etc.)
#' @param weighted Logical indicating whether weighted node averages are to be calculated; defaults to FALSE. 
#' @param prune Logical indicating whether tree should be pruned to contain just the species represented in the dataset; defaults to TRUE. 
#' @param av Logical indicating whether to average the values of multiple entries for the same pair, should they appear in the dataset; defaults to TRUE. If set to FALSE, function will stop in the case of more than one entry in the dataset corresponding to the same pair.
#'
#' @return Numeric vector of node averages, named according to node indices.
#'
#' @examples
#' # Load simulated dataset and tree
#' data(data1)
#' data(sim.tree1)
#' 
#' # Perform node averaging
#' unwtd <- node.averager(dataset = data1, tree = sim.tree1, varb = "pred", 
#'   taxacolumns = c("sp1", "sp2"))
#' wtd <- node.averager(dataset = data1, tree = sim.tree1, varb = "pred", 
#'   taxacolumns = c("sp1", "sp2"), weighted = TRUE)
#' 
#' # Compare outcomes of weighted and unweighted node averaging
#' unwtd
#' wtd
#' summary(unwtd)
#' summary(wtd)
#' 
#' # Calculate data loss
#' nrow(data1) - length(wtd)
#' nrow(data1) - length(unwtd)
#' 
#' \donttest{
#' # Plot tree and node labels
#' library(ape)
#' plot(sim.tree1)
#' nodelabels()
#' }

#' @references
#'    Coyne, J. A., Orr, H. A. 1989. Patterns of speciation in Drosophila. Evolution 43:362-381.
#' 
#'    Fitzpatrick, B. M. 2002. Molecular correlates of reproductive isolation. Evolution 56:191-198.
#' 
#'    Fitzpatrick, B. M., Turelli. 2006. The geography of mammalian speciation: mixed signals from phylogenies and range maps. Evolution 60:601-615.

#' @export
node.averager = function(dataset, tree, taxacolumns, varb, weighted=FALSE, prune=TRUE, av=TRUE) {
  # Ensure the dataset has the right column names
  if(!all(c(taxacolumns, varb) %in% colnames(dataset))) {
    stop("one or more required columns is missing in the dataset, or there is an incorrect column name given as an argument")
  }
  sp1=dataset[,taxacolumns[1]]
  sp2=dataset[,taxacolumns[2]]
  # Prune the tree to have just the taxa that are in the dataset (required due to extreme size and number of internal nodes of some trees)
  if(prune) tree = ape::drop.tip(phy=tree, tip=tree$tip.label[which(!tree$tip.label %in% c(dataset[,taxacolumns[1]],dataset[,taxacolumns[2]]))])
  # Find the nodelabels (i.e. tip label numbers) in the tree for each species in the dataset
  # warning: I think the lines below will always work, but its worth doing a quick visual check via tree plotting to make sure that the order of tip.labels matches the order of nodelabels (i.e. tip numbers)
  sp1.tip = sapply(sp1, function(x) which(tree$tip.label==x))
  sp2.tip = sapply(sp2, function(x) which(tree$tip.label==x))
  # Get the index of each internal node in the tree
  nodes = unique(tree$edge[,1])
  # For each node, calculate the average of the trait value for pairs spanning the node
  res = rep(NA, length(nodes)); names(res)=nodes
  for(i in 1:length(nodes)) {
    # Get the index in the tree for the tips that descend from just this node
    sps = sapply(ape::extract.clade(tree, nodes[i])$tip.label, function(x) which(tree$tip.label==x))
    # What are the possible pairwise combinations for these tips?
    pairwise.sp=utils::combn(sps,2)
    ## Not all pairwise combinations will be used going forward; we only need the combinations of taxa that occur on OPPOSITE SIDES OF THIS NODE
    # Such pairwise combinations will contain taxa that share this node as their MRCA
    inds = which(apply(pairwise.sp, 2, function(x) ape::getMRCA(tree, x)==nodes[i]))
    # Trim the pairwise combos matrix to just these relevant pairs (unless this is a sister node)
    if(length(inds)>1) pairwise.sp = pairwise.sp[,inds]
    # For each relevant pair, find the value of the variable for the two taxa in the dataset and, if weighted=TRUE, the number of nodes between them in the tree
    node.res=rep(NA, ncol(pairwise.sp))
    for(k in 1:ncol(pairwise.sp)) {
      # Find variable (varb) between these taxa
      var.indexed = dataset[,varb][sp1.tip %in% pairwise.sp[,k] & sp2.tip %in% pairwise.sp[,k]]
      if(length(var.indexed)==0) var.indexed=NA
      if(length(var.indexed)>1 & av) var.indexed=mean(var.indexed)
      if(length(var.indexed)>1 & !av) stop(paste("multiple values in dataset for the pair:", dataset[,taxacolumns][sp1.tip %in% pairwise.sp[,k] & sp2.tip %in% pairwise.sp[,k]], sep=" "))
      # Find the number of nodes between these taxa in the tree (subtract 2 from the length of the "nodepath", which in 'ape' includes the two tips)
      if(weighted) {
        num.nodes=length(ape::nodepath(tree, from=pairwise.sp[1,k], to=pairwise.sp[2,k]))-2
        var.indexed = var.indexed*(0.5)^(num.nodes-1)
      }
      node.res[k]=var.indexed
    }
    # Average the values in the node.res vector
    if(!all(is.na(node.res))) {
      res[i] = if(weighted) sum(node.res, na.rm=TRUE) else mean(node.res, na.rm=TRUE)
    }
  }
  return(res)
}
