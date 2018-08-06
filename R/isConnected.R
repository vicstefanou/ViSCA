#' Check if a tree is connected
#'
#' Checks if a lineage or division tree is connected.
#' A tree is connected if only one motherless cell exists (i.e. the root of the tree).
#'
#' @param tree The lineage or division tree, an object of class \code{"igraph"}.
#'
#' @return A logical value (\code{TRUE} or \code{FALSE})
#' indicating whether the \code{tree} is connected or not.
#'
#' @export
#' @import igraph

isConnected <- function(tree) {

  V(tree)$degree <- degree(graph = tree, v = V(tree), mode = "in")
  roots <- V(tree)[V(tree)$degree == 0]$name

  if (length(roots) == 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }

}
