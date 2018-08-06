#' Add an attribute
#'
#' Adds an attribute to a lineage or division tree.
#'
#' @param tree The lineage or division tree, an object of class \code{"igraph"}.
#'
#' @param attr The name of the new attribute, a character string.
#' It should not contain the underscore character \code{"_"}.
#'
#' @param vals The values of the new attribute,
#' a vector of numeric values, logical values or character strings.
#' Length of the vector must be equal to the number of all cells in the \code{tree},
#' as returned from \code{\link{get_cells}}.
#'
#' @return The updated tree with the new attribute \code{attr} added, an object of class \code{"igraph"}.
#'
#' @export
#' @import igraph

add_attr <- function(tree, attr, vals) {

  ############ arguments check ###########################

  if (length(unlist(strsplit(attr, "_"))) > 1) { # "_" exists
    stop("attr should not contain \"_\"\n")
  }

  if (length(vals) != length(get_cells(tree = tree, type = "all"))) {
    stop("Length of vals must be equal to the number of all cells in tree")
  }

  #######################################
  
  tree <- set_vertex_attr(graph = tree, name = attr, index = V(tree), value = vals)

  return(tree)

}
