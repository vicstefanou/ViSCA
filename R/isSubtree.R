#' Check if a tree is a subtree of another tree
#'
#' Checks if a lineage or division tree is a subtree of another lineage or division tree,
#' regarding only the nodes of the trees.
#'
#' @param subtree The lineage or division tree which is supposed to be the subtree of \code{tree},
#' an object of class \code{"igraph"}.
#'
#' @param tree The lineage or division tree which is supposed to be the supertree of \code{subtree},
#' an object of class \code{"igraph"}.
#'
#' @param treeT A character string naming the type of \code{subtree} and \code{tree}:
#' \itemize{
#' \item \code{"LT"} if \code{subtree} and \code{tree} are lineage trees
#' \item \code{"DT"} if \code{subtree} and \code{tree} are division trees
#' }
#'
#' @param type A character string naming the type of cells regarding which the statement will be checked:
#' \itemize{
#' \item \code{"all"} for all cells (including any existing imaginary \emph{root} cells)
#' \item \code{"nr"} for all non-root cells (excluding any existing imaginary \emph{root} cells)
#' \item \code{"inc"} for all cells included in the analysis
#' }
#'
#' @return A logical value (\code{TRUE} or \code{FALSE})
#' indicating whether the \code{subtree} is a subtree of the \code{tree} or not.
#'
#' @export
#' @import igraph

isSubtree <- function(subtree, tree, treeT = c("LT", "DT"), type = c("all", "nr", "inc")) {

  ################## arguments check ####################

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
  }

  if (!(type %in% c("all", "nr", "inc"))) {
    stop("type must be \"all\" / \"nr\" / \"inc\"\n")
  }

  #######################################################

  subtree_cells <- V(subtree)[get_cells(tree = subtree, treeT = treeT, type = type)]$cellName
  tree_cells  <- V(tree)[get_cells(tree = tree, treeT = treeT, type = type)]$cellName

  if (length(setdiff(subtree_cells, tree_cells)) == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }

}
