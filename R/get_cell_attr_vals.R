#' Get attributes' values of a cell
#'
#' Returns the values of the attributes of a cell in a lineage or division tree.
#'
#' @param tree The lineage or division tree where the cell specified in \code{cell} belongs,
#' an object of class \code{"igraph"}.
#'
#' @param treeT A character string naming the type of \code{tree}:
#' \itemize{
#' \item \code{"LT"} if \code{tree} is a lineage tree
#' \item \code{"DT"} if \code{tree} is a division tree
#' }
#'
#' @param cell The label of the cell in the \code{tree} whose attributes' values will be returned, a character string.
#' It can be any non-root cell, as returned from \code{\link{get_cells}}.
#'
#' @return The attributes' values of the \code{cell}, a named list.
#'
#' @export
#' @import igraph

get_cell_attr_vals <- function(tree, treeT = c("LT", "DT"), cell) {

  ################## arguments check #######################

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
  }

  possible_cells <- get_cells(tree = tree, treeT = treeT, type = "nr")

  if (!(cell %in% possible_cells)) {
    stop(paste("Selected cell ", cell, " does not exist\n", sep = ""))
  }

  ######################

  return(vertex_attr(graph = tree, index = cell))

}
