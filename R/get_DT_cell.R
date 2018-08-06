#' Get the label of the cell of a cell instant
#'
#' Returns the label of a cell in a division tree,
#' given the label of an instant of the cell in the corresponding lineage tree.
#'
#' @param LT The lineage tree where the cell specified in \code{cell} belongs,
#' an object of class \code{"igraph"}.
#' @param cell The label of the cell in the \code{LT}, a character string.
#' It can be any non-root cell, as returned from \code{\link{get_cells}}.
#' @param DT The corresponding division tree of the \code{LT},
#' an object of class \code{"igraph"}.
#'
#' @return A named list with the following components:
#' \item{cell}{The label of the \code{cell} in the \code{DT}, a character string.}
#' \item{instant}{The instant of the cell, a non-zero positive integer value.}
#' In case \code{cell} is not found in the \code{DT}, \code{NULL} is returned.
#'
#' @seealso \code{\link{get_LT_cell}} for the reverse.
#' @export
#' @import igraph

get_DT_cell <- function(LT, cell, DT) {

  ############# arguments check ############

  possible_cells <- get_cells(tree = LT, treeT = "LT", type = "nr")

  if (!(cell %in% possible_cells)) {
    stop(paste("Selected cell", paste("\"", cell, "\"", sep = ""), "does not exist\n"))
  }

  #####################################

  findCell <- function(x) {
    cell_instants <- unlist(strsplit(x, ", "))
    instant <- which(cell_instants == cell, arr.ind = TRUE)
    if (length(instant) == 0) {
      return(0)
    } else {
      return(instant)
    }
  }

  a <- unlist(lapply(V(DT)$cellInstants, findCell))
  DTind <- which(a != 0, arr.ind = TRUE)

  if (length(DTind) == 0) {
    return(NULL)
  } else {
    return(list(cell = V(DT)$name[DTind], instant = a[DTind]))
  }

}
