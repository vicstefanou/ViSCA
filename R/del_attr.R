#' Delete an attribute
#'
#' Deletes an attribute from a lineage or division tree.
#'
#' @param tree The lineage or division tree, an object of class \code{"igraph"}.
#'
#' @param attr The name of the attribute in the \code{tree}, a character string.
#' It can be any numeric, boolean or character string attribute, as returned from \code{\link{get_attr_names}}.
#'
#' @return The updated tree with the attribute \code{attr} deleted, an object of class \code{"igraph"}.
#'
#' @export
#' @import igraph

del_attr <- function(tree, attr) {

  ############ arguments check ###########################

  if (!(attr %in% vertex_attr_names(tree))) {
    stop(paste("attr \"", attr, "\" does not exist", sep = ""))
  }

  #######################################

  tree <- delete_vertex_attr(graph = tree, name = attr)

  return(tree)

}
