#' Find candidate mothers
#'
#' Returns the candidate mother cells for the root of a motherless branch in a lineage tree.
#'
#' @section Prerequisites:
#' This function can be used by \emph{BaSCA} users, importing the data with \code{\link{import_basca}}.
#' \cr\cr
#' Users of \emph{Oufti} or \emph{SuperSegger} 
#' who imported the data with \code{\link{import_oufti}} or \code{\link{import_ss}}, respectively,
#' are \bold{excluded} from using this function, as no colony list was returned.
#' \cr\cr
#' If \code{\link{import_json}} was used for importing the data,
#' it is necessary that cell list elements have the \code{centroid} and \code{colId} components
#' and colony list elements have the \code{ULcorner} component.
#' See \code{\link{import_json}} for more details.
#' In other case, this function cannot be used (throws an error).
#'
#' Candidate mothers are cells of the same colony in the previous frame.
#'
#' @param branch The connected lineage tree (motherless branch) for which to find candidate mothers, an object of class \code{"igraph"}.
#'
#' @param LT The lineage tree where the candidate mothers of the \code{branch} will belong, an object of class \code{"igraph"}.
#'
#' @param Nd Maximum number of daughters for a cell in the \code{LT} to be considered as a candidate mother, a positive integer value.
#' The default value is \code{1}.
#'
#' @param Ncands Maximum number of candidate mothers to be returned, a non-zero positive integer value.
#' The default value is \code{5}. Use value \code{NULL} to return all candidate mothers.
#'
#' @param cell_list A list containing all the cell instants of the movie.
#' @param col_list A list containing all the colony instants of the movie.
#' @param Ncols Number of colonies in the movie, a non-zero positive integer value.
#'
#' @param show A logical value (\code{TRUE} or \code{FALSE}) indicating whether \code{\link{view_cell}}
#' will be called for the root of the \code{branch} and the candidate mothers.
#' This capability is useful in order to choose the desired mother cell and call \code{\link{add_branch}}.
#' The default value is \code{FALSE}.
#' See the \emph{Prerequisites} section of \code{\link{view_cell}} in order to be able to use value \code{TRUE}.
#'
#' @return The labels of the candidate mothers, a vector of character strings.
#' Candidate mothers are sorted by ascending distance between their centroid and
#' the centroid of the root of the \code{branch}, which implies best matching.
#' In case no candidate mothers are found, \code{NULL} is returned.
#'
#' @seealso \code{\link{isConnected}} for checking if a tree is connected.
#' @export
#' @import igraph
#' @importFrom stats dist

get_cand_mother_cells <- function(branch, LT,
                                  Nd = 1, Ncands = 5,
                                  cell_list, col_list, Ncols,
                                  show = FALSE) {

  ################## arguments check #######################

  if (!isConnected(tree = branch)) {
    stop("branch is disconnected\n")
  }

  ### prerequisites
  if ((m <- checkField(myList = cell_list, fieldName = "centroid")) != "") {
    stop(paste("cell list:", m))
  }
  if ((m <- checkField(myList = cell_list, fieldName = "colId")) != "") {
    stop(paste("cell list:", m))
  }
  if ((m <- checkField(myList = col_list, fieldName = "ULcorner")) != "") {
    stop(paste("colony list:", m))
  }

  #########################################

  V(branch)$degree <- degree(graph = branch, v = V(branch), mode = "in")
  branchRoot <- V(branch)[V(branch)$degree == 0]$name

  branchRoot_cell_list_ID <- as.numeric(branchRoot) - (Ncols + 1)

  branchRoot_prev_colId <- which(sapply(col_list, function(x) x$colName) %in% col_list[[cell_list[[branchRoot_cell_list_ID]]$colId]]$prev_colName)

  V(LT)$degree <- degree(LT, v = V(LT), mode = "out")

  LTcands <- V(LT)[as.numeric(V(LT)$colId) %in% branchRoot_prev_colId &
                     V(LT)$degree <= Nd]$name

  ##### find top N candidate cells based on distances between the centroids of the cells ########

  if (length(LTcands) == 0) {
    return(NULL)
  }

  cands_centroid <- NULL
  for (cand in LTcands) {
    cand_cell_list_ID <- as.numeric(cand) - (Ncols + 1)
    cand_ULcorner <- col_list[[cell_list[[cand_cell_list_ID]]$colId]]$ULcorner
    cand_centroid <- cell_list[[cand_cell_list_ID]]$centroid + cand_ULcorner
    cands_centroid <- rbind(cands_centroid, cand_centroid)
  }

  branchRoot_ULcorner <- col_list[[cell_list[[branchRoot_cell_list_ID]]$colId]]$ULcorner
  branchRoot_centroid <- cell_list[[branchRoot_cell_list_ID]]$centroid + branchRoot_ULcorner

  distances <-  NULL
  for (i_cand in 1:dim(cands_centroid)[1]) {
    distances <- c(distances, dist(rbind(branchRoot_centroid, cands_centroid[i_cand, ])))
  }

  LTcands <- LTcands[order(distances)]

  if (!is.null(Ncands)) {
    if (length(LTcands) > Ncands) {
      LTcands <- LTcands[1:Ncands]
    }
  }

  if (show) {
    view_cell(LT = branch,
              cells = branchRoot,
              cell_list = cell_list, col_list = col_list, Ncols = Ncols)
    view_cell(LT = LT,
              cells = LTcands,
              cell_list = cell_list, col_list = col_list, Ncols = Ncols)
  }

  return(LTcands)

}
