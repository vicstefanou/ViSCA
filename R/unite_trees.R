#' Unite multiple trees
#'
#' Unites multiple lineage or division trees from the same movie in one.
#' This function is useful for combining multiple selection criteria with logical OR operator(s).
#'
#' Nodes and edges of the resulting tree are the union of nodes and edges of the trees, respectively.
#' Therefore, the resulting tree is disconnected.
#' Any missing attributes among the trees are filled with \code{NA} values in the resulting tree.
#'
#' @param trees A list with the trees to be united.
#' Each element of the list is a lineage or division tree (an object of class \code{"igraph"}).
#'
#' @return The united tree, an object of class \code{"igraph"}.
#'
#' @seealso \code{\link{select_subtree}} for combining multiple selection criteria with logical AND operator(s).
#' @export
#' @import igraph

unite_trees <- function(trees) {

  ## subtree of tree myDT with cells of generation <= 5 or generation 9
  # myCriteria <- list()
  # myCriteria[[1]] <- list(attr = "generation", val = 5, op = "<=")
  # myDT1 <- select_subtree(tree = myDT, criteria = myCriteria)
  # myCriteria[[1]] <- list(attr = "generation", val = 9, op = "==")
  # myDT2 <- select_subtree(tree = myDT, criteria = myCriteria)
  # myTrees <- list()
  # myTrees[[1]] <- myDT1
  # myTrees[[2]] <- myDT2
  # myDT1.2 <- unite_trees(trees = myTrees)

  # find all possible attributes between all trees
  attrs <- NULL
  for (i in 1:length(trees)) {
    nodes <- as_data_frame(x = trees[[i]], what = "vertices")
    attrs <- unique(c(attrs, colnames(nodes)))
  }

  nodesNew <- NULL
  edgesNew <- NULL

  for (i in 1:length(trees)) {

    nodes <- as_data_frame(x = trees[[i]], what = "vertices")
    edges <- as_data_frame(x = trees[[i]], what = "edges")

    missing_attrs <- setdiff(attrs, colnames(nodes))  # Find names of missing attributes
    nodes[missing_attrs] <- NA   # Add them, filled with NA
    nodes <- nodes[attrs]        # change order of attributes to be the same between all trees

    nodesNew <- as.data.frame(unique(rbind(nodesNew, nodes)))
    edgesNew <- as.data.frame(unique(rbind(edgesNew, edges)))

  }

  treeNew <- graph_from_data_frame(d = as.data.frame(edgesNew), directed = TRUE, vertices = nodesNew)

  return(treeNew)

}
