#' View a cell
#'
#' Displays the images of a set of given cells in the Plots Pane of RStudio.
#'
#' A separate image for each cell specified in \code{cells} is generated.
#' The cell is viewed in its colony and is marked as red.
#' The rest cells of the colony are marked as white.
#'
#' @section Prerequisites:
#' This function can be used by \emph{BaSCA} users, importing the data with \code{\link{import_basca}}.
#' \cr\cr
#' Users of \emph{Oufti} or \emph{SuperSegger} 
#' who imported the data with \code{\link{import_oufti}} or \code{\link{import_ss}}, respectively,
#' are \bold{excluded} from using this function, as no colony list was returned.
#' \cr\cr
#' If \code{\link{import_json}} was used for importing the data,
#' it is necessary that cell list elements have the \code{pixelList} and \code{colId} components
#' and colony list elements have the \code{colImage} component.
#' See \code{\link{import_json}} for more details.
#' In other case, this function cannot be used (throws an error).
#'
#' @param LT The lineage tree where the cells specified in \code{cells} belong, an object of class \code{"igraph"}.
#' @param cells The labels of the cells in the \code{LT} to be viewed, a vector of character strings.
#' They can be any non-root cells, as returned from \code{\link{get_cells}}.
#' @param cell_list A list containing all the cell instants of the movie.
#' @param col_list A list containing all the colony instants of the movie.
#' @param Ncols Number of colonies in the movie, a non-zero positive integer value.
#'
#' @export
#' @import igraph
#' @importFrom graphics par plot rasterImage

view_cell <- function(LT, cells,
                      cell_list, col_list, Ncols) {

  ############ arguments check ###########################

  possible_cells <- get_cells(tree = LT, treeT = "LT", type = "nr")

  if (length(m <- setdiff(cells, possible_cells)) != 0) {
    stop(paste("Selected cell(s)", toString(paste("\"", m, "\"", sep = "")), "do not exist\n"))
  }

  ### prerequisites
  if ((m <- checkField(myList = cell_list, fieldName = "pixelList")) != "") {
    stop(paste("cell list:", m))
  }
  if ((m <- checkField(myList = cell_list, fieldName = "colId")) != "") {
    stop(paste("cell list:", m))
  }
  if ((m <- checkField(myList = col_list, fieldName = "colImage")) != "") {
    stop(paste("colony list:", m))
  }

  ######################

  for (cell in cells) {

    cell_list_ID <- as.numeric(cell) - (Ncols + 1)
    cell_pixels <- cell_list[[cell_list_ID]]$pixelList
    im <- col_list[[cell_list[[cell_list_ID]]$colId]]$colImage

    im[im == 0] <- "black"
    im[im == 1] <- "snow3"
    im[cell_pixels] <- "red"

    oldPar <- par()
    par(mar = c(0, 0, 2, 0))

    # blank plot
    plot(NA,
         main = paste("cell \"", cell, "\"", sep = ""),
         xlab = "", ylab = "",
         xlim = c(0, 1), ylim = c(0, 1),
         bty = "n", axes = 0, xaxs = 'i', yaxs = 'i')

    rasterImage(EBImage::Image(data = im, colormode = 'Color'),
                xleft = 0, ybottom = 0, xright = 1, ytop = 1)

    par(mar = oldPar$mar)

  }

}
