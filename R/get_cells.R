#' Get cells of a tree
#'
#' Returns the labels of the cells in a lineage or division tree.
#'
#' @param tree The lineage or division tree, an object of class \code{"igraph"}.
#'
#' @param treeT A character string naming the type of \code{tree}:
#' \itemize{
#' \item \code{"LT"} if \code{tree} is a lineage tree
#' \item \code{"DT"} if \code{tree} is a division tree
#' }
#' This argument is ignored in case \code{type = "all"}
#'
#' @param type A character string naming the type of cells to be returned:
#' \itemize{
#' \item \code{"all"} for all cells (including any existing imaginary \emph{root} cells)
#' \item \code{"nr"} for all non-root cells (excluding any existing imaginary \emph{root} cells)
#' \item \code{"inc"} for all cells included in the analysis
#' }
#'
#' @return The labels of the corresponding cells, a vector of character strings.
#'
#' @export
#' @import igraph

get_cells <- function(tree, treeT = c("LT", "DT"), type = c("all", "nr", "inc")) {

  ################## arguments check ####################

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
  }

  if (!(type %in% c("all", "nr", "inc"))) {
    stop("type must be \"all\" / \"nr\" / \"inc\"\n")
  }

  #######################################################

  if (type == "all") {
    cells <- V(tree)$name
  } else if (type == "nr") {
    if (treeT == "LT") {
      cells <- V(tree)[V(tree)$frame != -1]$name
    } else { # treeT == "DT"
      cells <- V(tree)[V(tree)$generation != -1]$name
    }
  } else { # type == "inc"
    if (treeT == "LT") {
      cells <- V(tree)[V(tree)$frame != -1]$name
    } else { # treeT == "DT"
      cells <- V(tree)[V(tree)$isConsidered == TRUE]$name
    }
  }

  return(cells)

}
