#' Get the label of an instant of a cell
#'
#' Returns the label of an instant of a cell in a lineage tree,
#' given the label of the cell in the corresponding division tree.
#'
#' @param DT The division tree where the cell specified in \code{cell} belongs,
#' an object of class \code{"igraph"}.
#' @param cell The label of the cell in the \code{DT}, a character string.
#' It can be any non-root cell, as returned from \code{\link{get_cells}}.
#' @param instant The instant of the cell specified in \code{cell},
#' a non-zero positive integer value.
#'
#' @return The label of the \code{cell} \code{instant} in the corresponding lineage tree of the \code{DT},
#' a character string.
#' In case \code{cell} has less instants than the specified \code{instant}, \code{NULL} is returned.
#'
#' @seealso \code{\link{get_DT_cell}} for the reverse.
#' @export
#' @import igraph

get_LT_cell <- function(DT, cell, instant) {

  ############# arguments check ############

  possible_cells <- get_cells(tree = DT, treeT = "DT", type = "nr")

  if (!(cell %in% possible_cells)) {
    stop(paste("Selected cell", paste("\"", cell, "\"", sep = ""), "does not exist\n"))
  }

  #####################################

  cell_instants <- unlist(strsplit(V(DT)[cell]$cellInstants, ", "))

  if (abs(instant) > length(cell_instants)) {
    warning(paste("cell", paste("\"", cell, "\"", sep = ""), "has", length(cell_instants), "cell instants\n"))
	return(NULL)
  } else {
    LT_cellLabel <- cell_instants[instant]
	return(LT_cellLabel)
  }

}
