#' Get family cells of a cell
#'
#' Returns the labels of the family cells of a cell in a lineage or division tree.
#'
#' @param tree The connected lineage or division tree where the cell specified in \code{cell} belongs, an object of class \code{"igraph"}.
#'
#' @param treeT A character string naming the type of \code{tree}:
#' \itemize{
#' \item \code{"LT"} if \code{tree} is a lineage tree
#' \item \code{"DT"} if \code{tree} is a division tree
#' }
#'
#' @param cell The label of the cell in the \code{tree} whose family cells will be returned, a character string.
#' It can be any non-root cell, as returned from \code{\link{get_cells}}.
#'
#' @param type A character string naming the type of the family cells to be returned:
#' \itemize{
#' \item \code{"m"} for mother
#' \item \code{"gm"} for grandmother
#' \item \code{"d"} for daughters
#' \item \code{"gd"} for granddaughters
#' \item \code{"s"} for sibling
#' \item \code{"c"} for cousins
#' }
#'
#' @return The labels of the correponding family cells of the \code{cell}, a vector of character strings.
#' In case no family cells are found, \code{NULL} is returned.
#'
#' @seealso \code{\link{isConnected}} for checking if a tree is connected.
#' @export
#' @import igraph

get_cell_fam <- function(tree, treeT = c("LT", "DT"), cell, type = c("m", "gm", "d", "gd", "s", "c")) {

  ################## arguments check #######################

  if (!isConnected(tree = tree)) {
    stop("tree is disconnected\n")
  }

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
  }

  if (!(type %in% c("m", "gm", "d", "gd", "s", "c"))) {
    stop("type must be \"m\" / \"gm\" / \"d\" / \"gd\" / \"s\" / \"c\"\n")
  }

  possible_cells <- get_cells(tree = tree, treeT = treeT, type = "nr")

  if (!(cell %in% possible_cells)) {
    stop(paste("Selected cell \"", cell, "\" does not exist\n", sep = ""))
  }

  ###############

  if (type == "m") { # "mother"
    cells <- V(tree)[unlist(ego(graph = tree, order = 1, nodes = cell, mode = "in", mindist = 1))]$name
  } else if (type == "gm") { # "grandmother"
    cells <- V(tree)[unlist(ego(graph = tree, order = 2, nodes = cell, mode = "in", mindist = 2))]$name
  } else if (type == "d") { # "daughters"
    cells <- V(tree)[unlist(ego(graph = tree, order = 1, nodes = cell, mode = "out", mindist = 1))]$name
  } else if (type == "gd") { # "granddaughters"
    cells <- V(tree)[unlist(ego(graph = tree, order = 2, nodes = cell, mode = "out", mindist = 2))]$name
  } else if (type == "s") { # "sibling"
    mother <- V(tree)[unlist(ego(graph = tree, order = 1, nodes = cell, mode = "in", mindist = 1))]$name
    if (mother %in% possible_cells) {
      children <-  V(tree)[unlist(ego(graph = tree, order = 1, nodes = mother, mode = "out", mindist = 1))]$name
      cells <- children[children != cell]
    } else {
      cells <- character()
    }
  } else { # type == "c" --> "cousins"
    grandmother <- V(tree)[unlist(ego(graph = tree, order = 2, nodes = cell, mode = "in", mindist = 2))]$name
    if (grandmother %in% possible_cells) {
      mother <- V(tree)[unlist(ego(graph = tree, order = 1, nodes = cell, mode = "in", mindist = 1))]$name
      mothers <-  V(tree)[unlist(ego(graph = tree, order = 1, nodes = grandmother, mode = "out", mindist = 1))]$name
      aunt <- mothers[mothers != mother]
      if (aunt %in% possible_cells) {
        cells <-  V(tree)[unlist(ego(graph = tree, order = 1, nodes = aunt, mode = "out", mindist = 1))]$name
      } else {
        cells <- character()
      }
    } else {
      cells <- character()
    }
  }

  if (length(cells <- cells[cells %in% possible_cells]) == 0) {
    return(NULL)
  } else {
    return(cells)
  }

}
