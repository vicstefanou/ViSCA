#' Create cells' FLT
#'
#' Creates the cells' Forest of Lineage Trees (FLT) given a list containing all the cell instants of the movie.
#'
#' Obligatory attributes are:
#' \itemize{
#' \item \code{"name"} is the ID of the cell, a non-zero positive integer number stored as a character string.
#' Value \code{"1"} corresponds to the main \emph{root} cell.
#' Values \code{"1+<i>"} correspond to the colonies' \emph{root} cells,
#' where \code{"<i>"} is the colony ID.
#' \item \code{"cellName"} is the name of the cell, a character string.
#' The main \emph{root} cell is named as \code{"root"}.
#' The colonies' \emph{root} cells are named as \code{"colony<i>"},
#' where \code{"<i>"} is the colony ID.
#' \item \code{"frame"} is the ID of the frame of the cell, a non-zero positive integer number
#' \item \code{"colony"} is the ID of the colony of the cell, a non-zero positive integer number
#' }
#'
#' Any other numeric or boolean attribute of a cell is also stored as node attribute to the FLT.
#' The imaginary \emph{root} cells have value \code{-1} in all numeric attributes and value \code{FALSE} in all boolean attributes.
#'
#' Apart from the cell instants of the movie,
#' an imaginary main \emph{root} cell as well as imaginary \emph{root} cells for each colony are also added to the FLT.
#' Colonies' \emph{root} cells are daughters of the main \emph{root} cell.
#' Cell instants of the first frame of the movie are daughters of the corresponding colony's \emph{root} cell.
#' The imaginary \emph{root} cells are used to facilitate the tree representation of the movie and the colony tracking
#' and are automatically excluded from the analysis.
#'
#' @param cell_list A list containing all the cell instants of the movie.
#'
#' @param Ncols Number of colonies in the movie, a non-zero positive integer value.
#'
#' @return A named list with the following components:
#' \item{LTmain}{The main part of the overall FLT,
#' a connected lineage tree containing the imaginary \emph{root} cells (object of class \code{"igraph"}).}
#' \item{branches}{A list with the motherless branches of the overall FLT.
#' Each branch (element of the list) is a connected lineage tree (object of class \code{"igraph"}).
#' Motherless branches arise from tracking errors, in case a cell instant (root of the branch) fails to be
#' connected to any cell instant of the previous frame, or when a cell instant (root of the branch) just entered the field of view.}
#'
#' @seealso \code{\link{save_tree}} for saving a tree on disc,
#' \code{\link{add_branch}} for connecting a motherless branch to a lineage tree.
#' @export
#' @import igraph
#' @importFrom utils setTxtProgressBar txtProgressBar

createFLT <- function(cell_list, Ncols) {

  ###################################### create edgesLT ################################

  edgesLT <- NULL

  # tree root
  for (i_colony in 1:Ncols) {
    edgesLT <- rbind(edgesLT, c(1, i_colony + 1))
  }

  # colonies roots
  cells <- which(sapply(cell_list, function(x) x$frame) == 1)
  for (cell in cells) {
    edgesLT <- rbind(edgesLT, c(cell_list[[cell]]$colony + 1, cell + Ncols + 1))
  }

  pb <- txtProgressBar(min = 0, max = length(cell_list), style = 3) ### set progress bar
  ipb <- 0

  # cells
  for (i_cell in 1:length(cell_list)) {

    ipb <- ipb + 1
    setTxtProgressBar(pb, ipb) ### update progress bar

    daughters <- cell_list[[i_cell]]$daughterIds
    if (!is.null(daughters)) {
      for (daughter_name in daughters) {
        daughter_id <- which(sapply(cell_list, function(x) x$cellName == daughter_name))
        if (length(daughter_id) != 0) {  # daughter found
          edgesLT <- rbind(edgesLT, c(i_cell + Ncols + 1, daughter_id + Ncols + 1))
        }
      }
    }

  }

  close(pb) ### close progress bar
  cat("\n")

  ###################################### create nodesLT ################################

  nodesLT <- as.matrix(c(1:(length(cell_list) + Ncols + 1))) # all cells
  nodesLT <- as.data.frame(nodesLT[order(nodesLT), ], stringsAsFactors = FALSE)
  colnames(nodesLT) <- "cellID"

  cells <- nodesLT[-(1:(Ncols + 1)), "cellID"] - (Ncols + 1) # all except roots

  # cellName attribute
  roots_cellName <- "root"
  for (i_colony in 1:Ncols) {
    roots_cellName <- c(roots_cellName, paste("colony", i_colony, sep = ""))
  }
  nodesLT[, "cellName"] <- c(roots_cellName, unlist(unname(sapply(cell_list, function(x) x$cellName))[cells]))

  # numeric attributes
  numeric_attrs <- names(cell_list[[1]])[sapply(cell_list[[1]], function(x) class(x) == "numeric" || class(x) == "integer")]
  for (attr in numeric_attrs) {
    if (attr == "colId") {
      nodesLT[, "colId"] <- as.character(c(rep(-1, Ncols + 1), unlist(unname(sapply(cell_list, function(x) x$colId))[cells])))
    } else {
      nodesLT[, attr] <- c(rep(-1, Ncols + 1), unlist(unname(sapply(cell_list, function(x) x[attr]))[cells]))
    }
  }

  # boolean attributes
  boolean_attrs <- names(cell_list[[1]])[sapply(cell_list[[1]], function(x) is.logical(x))]
  for (attr in boolean_attrs) {
    nodesLT[, attr] <- c(rep(FALSE, Ncols + 1), unlist(unname(sapply(cell_list, function(x) x[attr]))[cells]))
  }

  ######################################### create LT ##############################

  overallLT <- graph_from_data_frame(d = as.data.frame(edgesLT), directed = TRUE, vertices = nodesLT)
  branches <- decompose(graph = overallLT) # returns a list of the clusters graphs
  i_wanted_gr <- which.max(clusters(overallLT)$csize)
  LTmain <- branches[[i_wanted_gr]]  # keep the biggest graph, which also happens to contain the root cells
  branches[[i_wanted_gr]] <- NULL

  # track/correct colonies (for LTmain only)
  for (i_colony in 1:Ncols) {
    V(LTmain)[subcomponent(graph = LTmain, v = i_colony + 1, mode = "out")]$colony <- i_colony
  }

  return(list(LTmain = LTmain, branches = branches))

}
