#' Merge two cells
#'
#' Merges two cells from the same or different lineage trees.
#'
#' After the merge operation, \code{cell1} is replaced in both \code{LT1} and \code{cell_list} by the resulting cell.
#' Daughter branches of \code{cell2} are extracted from \code{LT2} (or \code{LT1} in case \code{LT2 = NULL}),
#' with successive calls of \code{\link{extract_branch}}.
#' These motherless branches (lineage trees) are added as daughter branches to \code{cell1},
#' with successive calls of \code{\link{add_branch}}, until \code{cell1} has two daughters.
#' Finally, \code{cell2} is deleted from \code{LT2} (or \code{LT1} in case \code{LT2 = NULL})
#' with \code{\link{extract_branch}}.
#'
#' @section Prerequisites:
#' This function can be used by \emph{BaSCA} users \bold{only},
#' importing the data with \code{\link{import_basca}}.
#'
#' @param LT1 The lineage tree where the cell specified in \code{cell1} belongs, an object of class \code{"igraph"}.
#' @param LT2 The lineage tree where the cell specified in \code{cell12} belongs, an object of class \code{"igraph"}.
#' When the default value \code{NULL} is used, \code{cell2} belongs to the \code{LT}.
#' @param cell1 The label of the first cell in the \code{LT1} to be merged, a character string.
#' It can be any non-root cell, as returned from \code{\link{get_cells}}.
#' @param cell2 The label of the second cell in the \code{LT2} (or \code{LT1} in case \code{LT2 = NULL}) to be merged, a character string.
#' It can be any valid candidate merge cell, as returned from \code{\link{get_cand_merge_cells}}.
#' @param cell_list A list containing all the cell instants of the movie.
#' @param col_list A list containing all the colony instants of the movie.
#' @param Ncols Number of colonies in the movie, a non-zero positive integer value.
#' @param pixelR The pixel ratio in units of length, a non-zero positive numeric value.
#'
#' @param matFolder A character string naming the absolute path of the directory where
#' the \code{.mat} file generated by \emph{BaSCA} is saved (excluding the last \code{"/"}).
#' The default value is the current working directory \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{"/"} on Windows.
#' @param matFileName A character string naming the \code{.mat} file generated by \emph{BaSCA}
#' (including the suffix \code{".mat"}).
#' The filename is relative to the \code{matFolder}.
#' @param exeFolder A character string naming the absolute path of the installation folder of the MATLAB executable
#' (excluding the last \code{"/"}).
#' \cr\cr
#' NOTE: The components should be separated by \code{"/"} on Windows.
#' @param mcrFolder A character string naming the absolute path of the installation folder of the Matlab Compiler Runtime (MCR)
#' (excluding the last \code{"/"}).
#' \cr\cr
#' NOTE: The components should be separated by \code{"/"} on Windows.
#'
#' @param show A logical value (\code{TRUE} or \code{FALSE}) indicating whether \code{\link{view_cell}}
#' will be called for \code{cell1} and \code{cell2} before the merge operation
#' and the resulting cell after the merge operation.
#' This capability is useful in order to see the result of the function.
#' The default value is \code{TRUE}.
#'
#' @return A named list with the following components:
#' \item{LT1}{The updated LT1 with \code{cell1} replaced and the possible daughter branches added/replaced,
#' an object of class \code{"igraph"}.}
#' \item{LT2}{The updated LT2 with \code{cell2} and its daughter branches deleted,
#' an object of class \code{"igraph"}.
#' \code{NULL} is returned in case LT2 was anywise \code{NULL} or if it ended up with no cells.}
#' \item{cell_list}{The updated cell_list with \code{cell1} replaced.}
#' \item{branches}{A list with the remaining motherless branches.
#' Each branch (element of the list) is an object of class \code{"igraph"}.}
#'
#' @seealso \code{\link{split_cell}} for the reverse.
#' @export
#' @import igraph

merge_cells <- function(LT1, LT2 = NULL,
                        cell1, cell2,
                        cell_list, col_list, Ncols,
                        pixelR,
                        matFolder, matFileName,
                        exeFolder, mcrFolder,
                        show = TRUE) {

  #exeFolder = "/usr/BaSCA/mat_correct"
  #mcrFolder = "/usr/local/MATLAB/MATLAB_Compiler_Runtime"

  ############ arguments check ###########################

  if (is.null(LT2)) {

    possible_cells <- get_cells(tree = LT1, treeT = "LT", type = "nr")
    if (length(m <- setdiff(c(cell1, cell2), possible_cells)) != 0) {
      stop(paste("Selected cell(s)", toString(paste("\"", m, "\"", sep = "")), "do not exist\n"))
    }

  } else {

    possible_cells  <- get_cells(tree = LT1, treeT = "LT", type = "nr")
    if (!(cell1 %in% possible_cells)) {
      stop(paste("Selected cell", paste("\"", cell1, "\"", sep = ""), "does not exist\n"))
    }

    possible_cells  <- get_cells(tree = LT2, treeT = "LT", type = "nr")
    if (!(cell2 %in% possible_cells)) {
      stop(paste("Selected cell", paste("\"", cell2, "\"", sep = ""), "does not exist\n"))
    }

  }

  if (!is.null(m <- get_cand_merge_cells(LT = LT1, cell = cell1, LTcand = LT2,
                                         cell_list = cell_list, col_list = col_list, Ncols = Ncols, show = FALSE))) {
    if (!(cell2 %in% m)) {
      stop(paste("Selected cell", paste("\"", cell1, "\"", sep = ""), "cannot be merged with cell", paste("\"", cell2, "\"", sep = ""), "\n"))
    }
  } else {
    stop(paste("Selected cell", paste("\"", cell1, "\"", sep = ""), "cannot be merged with any cell\n"))
  }

  ###################### 2 to 1
  branches <- list()

  if (show) {
    if (is.null(LT2)) {
      view_cell(LT = LT1, cells = c(cell1, cell2),
                cell_list = cell_list, col_list = col_list, Ncols = Ncols)
    } else {
      view_cell(LT = LT1, cells = cell1,
                cell_list = cell_list, col_list = col_list, Ncols = Ncols)
      view_cell(LT = LT2, cells = cell2,
                cell_list = cell_list, col_list = col_list, Ncols = Ncols)
    }
  }

  cell_list_ID1 <- as.numeric(cell1) - (Ncols + 1)
  cell_list_ID2 <- as.numeric(cell2) - (Ncols + 1)

  myFrame <- cell_list[[cell_list_ID1]]$frame
  myColony <- cell_list[[cell_list_ID1]]$colony
  myColId <- cell_list[[cell_list_ID1]]$colId
  correctColony <- V(LT1)[cell1]$colony # if LT1 is mainLT, it is already correct

  mergeCellProps <- callBascaMatlabExe(Nsplit = 0,
                                       cell_list_IDs = c(cell_list_ID1, cell_list_ID2),
                                       cell_list = cell_list,
                                       col_list = col_list,
                                       mat_folder = matFolder,
                                       mat_fileName = matFileName,
                                       install_folder_exe = exeFolder,
                                       install_folder_MCR = mcrFolder)

  if (dim(mergeCellProps$cellProps)[3] != 1) {
    stop("Cells could not be merged\n")
  }

  # update cell_list, LT1 for cell1
  cell_list[[cell_list_ID1]] <- mergeCellProps$cellProps[, , 1]
  cell_list[[cell_list_ID1]] <- updateCellInBascaList(x = cell_list[[cell_list_ID1]],
                                                      frame = myFrame,
                                                      colony = myColony,
                                                      colId = myColId,
                                                      col_list = col_list,
                                                      pixelRatio = pixelR)

  LT1 <- updateCellInBascaLT(LT = LT1,
                             cell = cell1,
                             cell_list = cell_list,
                             Ncols = Ncols,
                             colony = correctColony)


  # add children of cell2 as children to cell1
  if (is.null(LT2)) {
    children <- get_cell_fam(tree = LT1, treeT = "LT", cell = cell2, type = "d")
  } else {
    children <- get_cell_fam(tree = LT2, treeT = "LT", cell = cell2, type = "d")
  }

  if (!is.null(children)) {
    for (child in children) {

      if (is.null(LT2)) {
        br <- extract_branch(tree = LT1, cell = child)$branch
      } else {
        br <- extract_branch(tree = LT2, cell = child)$branch
      }

      d <- get_cell_fam(tree = LT1, treeT = "LT", cell = cell1, type = "d")
      if (is.null(d) || length(d) == 1) { # cell1 has 0 or 1 children -> add as child
        LT1 <- add_branch(LT = LT1,
                          branch = br,
                          cell = cell1,
                          cell_list = cell_list,
                          col_list = col_list,
                          Ncols = Ncols)
      } else { # cell1 has 2 children -> add to list with motherless branches
        branches[[length(branches) + 1]] <- br
      }

    }
  }

  # delete cell2
  if (is.null(LT2)) {
    LT1 <- extract_branch(tree = LT1, cell = cell2)$treeNew
  } else {
    LT2 <- extract_branch(tree = LT2, cell = cell2)$treeNew
    if (gorder(LT2) == 0) { # LT2 is now empty
      LT2 <- NULL
    }
  }

  if (show) {
    view_cell(LT = LT1, cells = cell1,
              cell_list = cell_list, col_list = col_list, Ncols = Ncols)
  }

  cat("Checking updated cell list...\n")
  cell_list <- checkCellList(cell_list = cell_list, col_list = col_list)

  return(list(LT1 = LT1, LT2 = LT2, cell_list = cell_list, branches = branches))

}
