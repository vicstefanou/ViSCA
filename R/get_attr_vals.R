#' Get attribute's values
#'
#' Returns the values of an attribute for all cells of a lineage or division tree,
#' as returned from \code{\link{get_cells}}.
#'
#' @param tree The lineage or division tree, an object of class \code{"igraph"}.
#'
#' @param attr The name of the attribute in the \code{tree}, a character string.
#' It can be any numeric, boolean or character string attribute, as returned from \code{\link{get_attr_names}}.
#'
#' @return A vector with the values of attribute \code{attr}.
#'
#' @export
#' @import igraph

get_attr_vals <- function(tree, attr) {

  ############# arguments check ############

  if (!(attr %in% vertex_attr_names(tree))) {
    stop(paste("attr \"", attr, "\" does not exist", sep = ""))
  }

  ##################

  cells <- get_cells(tree = tree, type = "all")

  return(vertex_attr(graph = tree, name = attr, index = cells))

}
