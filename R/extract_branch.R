#' Extract a branch from a tree
#'
#' Extracts a branch from a lineage or division tree. The branch for extraction is defined by its root cell.
#'
#' @param tree The lineage or division tree, an object of class \code{"igraph"}.
#'
#' @param cell The label of the cell in the \code{tree} which is the root of the branch to be extracted, a character string.
#' It can be any cell, as returned from \code{\link{get_cells}}.
#'
#' @return A named list with the following components:
#' \item{treeNew}{The new tree with the branch extracted, an object of class \code{"igraph"}.}
#' \item{branch}{The extracted connected motherless branch, an object of class \code{"igraph"}.}
#'
#' @seealso \code{\link{add_branch}} for connecting a motherless branch to a lineage tree.
#' @export
#' @import igraph

extract_branch <- function(tree, cell) {

  ############# arguments check ############

  possible_cells <- get_cells(tree = tree, type = "all")

  if (!(cell %in% possible_cells)) {
    stop(paste("Selected cell", paste("\"", cell, "\"", sep = ""), "does not exist\n"))
  }

  #####################################

  branch_cells <- subcomponent(graph = tree, v = cell, mode = "out")
  treeNew_cells <- V(tree)[!(V(tree)$name %in% V(tree)[branch_cells]$name)]

  branch <- induced_subgraph(graph = tree, vids = branch_cells)
  treeNew <- induced_subgraph(graph = tree, vids = treeNew_cells)

  return(list(treeNew = treeNew, branch = branch))

}
