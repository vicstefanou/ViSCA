#' Get attributes' names
#'
#' Returns the names of the attributes of a lineage or division tree.
#'
#' @param tree The lineage or division tree, an object of class \code{"igraph"}.
#'
#' @param type A character string naming the type of attributes to be returned:
#' \itemize{
#' \item \code{"n"} for numeric attributes
#' \item \code{"b"} for boolean attributes
#' \item \code{"c"} for character string attributes
#' }
#'
#' @return The correponding attributes' names, a vector of character strings.
#'
#' @export
#' @import igraph

get_attr_names <- function(tree, type = c("n", "b", "c")) {

  attrs <- NULL

  if (type == "n") {
    for (v_attr in vertex_attr_names(tree)) {
      if (is.numeric(vertex_attr(graph = tree, name = v_attr, index = V(tree)))) {
        attrs <- c(attrs, v_attr)
      }
    }
  } else if (type == "b") {
    for (v_attr in vertex_attr_names(tree)) {
      if (is.logical(vertex_attr(graph = tree, name = v_attr, index = V(tree)))) {
        attrs <- c(attrs, v_attr)
      }
    }
  } else { # type == "c"
    for (v_attr in vertex_attr_names(tree)) {
      if (is.character(vertex_attr(graph = tree, name = v_attr, index = V(tree)))) {
        attrs <- c(attrs, v_attr)
      }
    }
  }

  return(attrs)

}
