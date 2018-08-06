#' Add a branch to an LT
#'
#' Adds a lineage tree as a branch to another lineage tree at a specific position,
#' defined by the cell which will be the mother of the root of the branch.
#'
#' @section Prerequisites:
#' See the \emph{Prerequisites} section of \code{\link{get_cand_mother_cells}}.
#'
#' @param LT The lineage tree to which the \code{branch} will be added, an object of class \code{"igraph"}.
#'
#' @param branch The connected lineage tree (motherless branch) which will be added as a branch to the \code{LT}, an object of class \code{"igraph"}.
#'
#' @param cell The label of the cell in the \code{LT} where the \code{branch} will be added, a character string.
#' It can be any valid candidate mother cell, as returned from \code{\link{get_cand_mother_cells}}.
#'
#' @param cell_list A list containing all the cell instants of the movie.
#'
#' @param col_list A list containing all the colony instants of the movie.
#'
#' @param Ncols Number of colonies in the movie, a non-zero positive integer value.
#'
#' @return The new LT with the branch added, an object of class \code{"igraph"}.
#'
#' @seealso \code{\link{extract_branch}} for extracting a branch from a tree,
#' \code{\link{isConnected}} for checking if a tree is connected.
#' @export
#' @import igraph

add_branch <- function(LT, branch, cell,
                       cell_list, col_list, Ncols) {

  ############# arguments check ############

  if (!isConnected(tree = branch)) {
    stop("branch is disconnected\n")
  }

  possible_cells <- get_cells(tree = LT, treeT = "LT", type = "all")
  if (!(cell %in% possible_cells)) {
    stop(paste("Selected cell", paste("\"", cell, "\"", sep = ""), "does not exist\n"))
  }

  if (!is.null(m <- get_cand_mother_cells(branch = branch, LT = LT,
                                          cell_list = cell_list, col_list = col_list, Ncols = Ncols,
                                          Nd = 2, Ncands = NULL, show = FALSE))) {
    if (!(cell %in% m)) {
      stop(paste("Selected cell", paste("\"", cell, "\"", sep = ""), "is not a candidate mother cell\n"))
    }

  } else {
    stop("Selected branch has no candidate mother cell\n")
  }

  if (!is.null(d <- get_cell_fam(tree = LT, treeT = "LT", cell = cell, type = "d"))) {
    if (length(d) == 2) {
      stop(paste("Selected cell", paste("\"", cell, "\"", sep = ""), "has already 2 daughters\n"))
    }
  }

  ##########################################

  # track/correct colony (for branch)
  V(branch)$colony <- V(LT)[cell]$colony

  # concatenate edges and nodes of LT and branch
  nodes_LT <- as_data_frame(x = LT, what = "vertices")
  edges_LT <- as_data_frame(x = LT, what = "edges")

  nodes_branch <- as_data_frame(x = branch, what = "vertices")
  edges_branch <- as_data_frame(x = branch, what = "edges")

  nodes <- as.data.frame(rbind(nodes_LT, nodes_branch), stringsAsFactor = FALSE)
  edges <- as.data.frame(rbind(edges_LT, edges_branch), stringsAsFactor = FALSE)

  # add edge between LT and branch
  branch_root <- V(branch)[degree(graph = branch, v = V(branch), mode = "in") == 0]$name
  edgeNew <- c(as.numeric(cell), as.numeric(branch_root))
  edges <- as.data.frame(rbind(edges, edgeNew), stringsAsFactor = FALSE)

  LTnew <- graph_from_data_frame(d = edges, directed = TRUE, vertices = nodes)

  return(LTnew)

}
