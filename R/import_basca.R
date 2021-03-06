#' Import \emph{BaSCA} data
#'
#' Imports a \code{.mat} file generated by \emph{BaSCA} and converts it into
#' a cell list and colony list containing all the cell instants and colony instants of the movie, respectively.
#'
#' @param file A character string naming the \code{.mat} file generated by \emph{BaSCA}
#' (including the suffix \code{".mat"}) from which the data is to be imported.
#' If it does not contain an absolute path, the file name is relative to the current working directory, \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{"/"} on Windows.
#'
#' @param pixelR The pixel ratio in units of length, a non-zero positive numeric value.
#'
#' @param ringW The width of the ring (in \emph{pixels}) which marks the boundary pixels of a colony, a non-zero positive integer value.
#' This value is used to decide whether a cell is on the boundary of its colony or not.
#' The default value is \code{10}.
#'
#' @return A named list with the following components:
#' \item{col_list}{A list containing all the colony instants of the movie.
#' Each element of the list is a named list with the following components:
#' \itemize{
#' \item \code{colName} is the name of the colony instant, a character string in the format \code{"f<frame>_c<colony>"},
#' where \code{"<frame>"} and \code{"<colony>"} is the ID of the frame and colony (in the frame) of the colony instant, respectively.
#' \item \code{prev_colName} is a vector of character strings containing the \code{colName}
#' of the corresponding colony instant(s) in the previous frame.
#' For colony instants that do not have a corresponding colony instant in the previous frame, this is equal to \code{"f0_c0"}.
#' \item \code{next_colName} is a character string containing the \code{colName}
#' of the corresponding colony instant in the next frame.
#' For colony instants that do not have a corresponding colony instant in the next frame, this is equal to \code{"f00_c0"}.
#' \item \code{colImage} is the mask of the box surrounding the colony instant, a matrix of \code{0} and \code{1}.
#' \code{1}s denote the pixels of cells and \code{0}s the background pixels.
#' \item \code{ULcorner} is an 1x2 matrix of non-zero integer values
#' denoting the upper-left pixel of the box surrounding the colony instant in global (frame) coordinates.
#' The first integer represents the row and the second the column of the pixel.
#' \item \code{colBoundaryPixels} is a Nx2 matrix of non-zero integer values
#' denoting the boundary pixels of the colony instant in colony coordinates
#' (i.e. relative to the \code{colImage}).
#' Each one of the N rows indicates a boundary pixel of the colony instant detected based on the \code{ringW} argument.
#' The first column represents the row and the second the column of the boundary pixel.
#' \item \code{colCentroid} is an 1x2 matrix of non-zero numeric values
#' denoting the centroid (geometric center) of the colony instant in colony coordinates
#' (i.e. relative to the \code{colImage}).
#' }
#' }
#' \item{cell_list}{A list containing all the cell instants of the movie.
#' Each element of the list is a named list with the following components:
#' \itemize{
#' \item \code{cellName} is the name of the cell, a character string in the format \code{"x<------->_y<------->_f<frame>"}.
#' \item \code{frame} is the ID of the frame of the cell, a non-zero positive integer number.
#' \item \code{colony} is the ID of the colony of the cell in the \code{frame}, a non-zero positive integer number.
#' \item \code{daughterIds} is a vector of character strings containing the \code{cellName}
#' of the linked cell(s) in the next frame,
#' or \code{NULL} in case no such cells exist.
#' \item \code{colId} is a \emph{pointer} to the corresponding colony instant of the cell in the \code{col_list},
#' a non-zero positive integer value.
#' \cr\cr
#' Colonies that entered the field of view at a time point and did not exist from the beginning of the movie (i.e. from the first frame)
#' should not have tracked cells, until they merge (if this is the case) with another existing colony.
#' This means that no element should \emph{point} to such colony instants.
#' \item \code{pixelList} is a Nx2 matrix of non-zero integer values
#' denoting the pixels of the cell in colony coordinates
#' (i.e. relative to the \code{colImage} of the \code{\ifelse{html}{\out{colId<sup>th</sup>}}{\eqn{colId^{th}}}}
#' element in the \code{col_list}).
#' Each one of the N rows indicates a pixel of the cell.
#' The first column represents the row and the second the column of the pixel.
#' \item \code{boundaryPixelList} is a Nx2 matrix of non-zero integer values
#' denoting the boundary pixels of the cell in colony coordinates
#' (i.e. relative to the \code{colImage} of the \code{\ifelse{html}{\out{colId<sup>th</sup>}}{\eqn{colId^{th}}}}
#' element in the \code{col_list}).
#' Each one of the N rows indicates a boundary pixel of the cell.
#' The first column represents the row and the second the column of the boundary pixel.
#' \item \code{centroid} is an 1x2 matrix of non-zero numeric values
#' denoting the centroid (geometric center) of the cell in colony coordinates
#' (i.e. relative to the \code{colImage} of the \code{\ifelse{html}{\out{colId<sup>th</sup>}}{\eqn{colId^{th}}}}
#' element in the \code{col_list}).
#' It is the \emph{mean} of the \code{pixelList} by column.
#' \item \code{length} is the length of the cell in units of length, a non-zero positive numeric value.
#' \item \code{width} is the width of the cell in units of length, a non-zero positive numeric value.
#' \item \code{LW} is the length-to-width ratio, a non-zero positive numeric value.
#' \item \code{area} is the area of the cell in squared units of length, a non-zero positive numeric value.
#' \item \code{perimeter} is the perimeter of the cell in units of length, a non-zero positive numeric value.
#' \item \code{minorAxis} is the short axis of the ellipse surrounding the cell in units of length, a non-zero positive numeric value.
#' \item \code{majorAxis} is the long axis of the ellipse surrounding the cell in units of length, a non-zero positive numeric value.
#' \item \code{eccentricity} is a numeric value in the range \code{[0, 1]}
#' defining the eccentricity of the cell.
#' \item \code{orientation} is a numeric value in the range \code{[0, 360)}
#' defining the orientation of the cell in \emph{degrees}.
#' \item \code{solidity} is a numeric value in the range \code{[0, 1]}
#' defining the solidity of the cell.
#' \item \code{distFromCentroid} is the distance of the centroid of the cell from the centroid of its colony
#' (i.e. the euclidean distance between the \code{centroid} and
#' the \code{colCentroid} of the \code{\ifelse{html}{\out{colId<sup>th</sup>}}{\eqn{colId^{th}}}}
#' element in the \code{col_list})
#' in units of length, a non-zero positive numeric value.
#' \item \code{isOnBoundary} is a logical value (\code{TRUE} or \code{FALSE})
#' indicating whether the cell is on the boundary of its colony or not.
#' It is \code{TRUE} in case there is more than 66% overlapping between the \code{pixelList} and
#' the \code{colBoundaryPixels} of the \code{\ifelse{html}{\out{colId<sup>th</sup>}}{\eqn{colId^{th}}}}
#' element in the \code{col_list},
#' \code{FALSE} otherwise.
#' \item \code{fluorescenceInt<i>} is a Nx1 matrix of positive numeric values indicating the
#' cell fluorescence intensity of channel \code{<i>}.
#' The values correspond to each one of the pixels in \code{pixelList} (1-1 correspondence).
#' \item \code{fluorescenceInt<i>.mean} is the \emph{mean} cell fluorescence intensity of channel \code{<i>},
#' a positive numeric value.
#' \item \code{fluorescenceInt<i>.std} is the \emph{standard deviation} of the cell fluorescence intensity of channel \code{<i>},
#' a positive numeric value.
#' \item \code{fluorescenceInt<i>.coverage} is the percentage of the cell area covered by fluorescence of channel \code{<i>},
#' a numeric value in the range \code{[0, 1]}.
#' }
#' The \code{fluorescenceInt<i>*} components are included in case they were exctracted by \emph{BaSCA}.
#' }
#' \item{Nframes}{Number of frames in the movie, a non-zero positive integer value.
#' IDs of frames are in the range \code{[1, Nframes]}.}
#' \item{Ncols}{Number of colonies in the movie, a non-zero positive integer value.
#' IDs of colonies are in the range \code{[1, Ncols]}.
#' This value corresponds to the number of colonies at the start of the movie.}
#' \item{frameH}{Frame image height in \emph{pixels}, a non-zero positive integer value.}
#' \item{frameW}{Frame image width in \emph{pixels}, a non-zero positive integer value.}
#'
#' @references
#' A. Balomenos, P. Tsakanikas, Z. Aspridou, A. Tampakaki, K. Koutsoumanis and E. Manolakos,
#' \emph{“Image analysis driven single-cell analytics for systems microbiology”},
#' BMC Systems Biology, vol. 11, no. 1, 2017.
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar

import_basca <- function(file, pixelR, ringW = 10) {

  cat("Loading file...\n")
  mat_file <- R.matlab::readMat(file)
  data <- mat_file$frame.for.R[, , 1]

  pb <- txtProgressBar(min = 0, max = ncol(data), style = 3) ### set progress bar
  ipb <- 0

  # unlist daughterIds
  for (i_frame in 1:ncol(data)) {

    ipb <- ipb + 1
    setTxtProgressBar(pb, ipb) ### update progress bar

    for (i_colony in 1:dim(data[ , i_frame]$colonyProps)[3]) {

      if (is.null(ncol(data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1]))) {  # one cell in the colony
        if (!is.null(data[,i_frame]$colonyProps[, , i_colony]$cellProps[, , 1]$daughterIds)) {   # daughterIds field exists
          if ((length(data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1]$daughterIds) != 0 )) {  # not empty list
            data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1]$daughterIds <- as.vector(unlist(data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1]$daughterIds))
          }
        }
      } else {
        for (i_cell in 1:ncol(data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1])) {
          if (!is.null(data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1][, i_cell]$daughterIds)) {   # daughterIds field exists
            if ((length(data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1][, i_cell]$daughterIds) != 0)) {  # not empty list
              data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1][, i_cell]$daughterIds <- as.vector(unlist(data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1][, i_cell]$daughterIds))
            }
          }
        }
      }
    }
  }

  close(pb) ### close progress bar
  cat("\n")

  #####################

  cell_list <- col_list <- list()
  curr_cell <- curr_colony <- 0

  cat("Creating cell and colony lists...\n")
  pb <- txtProgressBar(min = 0, max = ncol(data), style = 3) ### set progress bar
  ipb <- 0

  for (i_frame in 1:ncol(data))  {

    ipb <- ipb + 1
    setTxtProgressBar(pb, ipb) ### update progress bar

    if (i_frame == 1){
      frameH <- as.vector(data[, i_frame]$x) #nrows
      frameW <- as.vector(data[, i_frame]$y) #ncols
    }

    for (i_colony in 1:dim(data[, i_frame]$colonyProps)[3]) {

      curr_colony <- curr_colony + 1

      ULcorner <- data[, i_frame]$colonyProps[, , i_colony]$bBoxULCorner # (x,y) -> (col,row) in frame
      ULcorner <- cbind(ULcorner[, 2], ULcorner[, 1]) # (row,col) in frame

      colCentroid <- data[, i_frame]$colonyProps[, , i_colony]$centroid # (x,y) -> (col,row) in frame
      colCentroid <- cbind(colCentroid[, 2], colCentroid[, 1]) # (row,col) in frame
      colCentroid <- colCentroid - ULcorner # (row,col) in colImage

      distanceMatrix <- EBImage::distmap(data[, i_frame]$colonyProps[, , i_colony]$bwColonyMask, metric = 'euclidean')
      outterRing <- distanceMatrix < ringW & distanceMatrix != 0
      colBoundaryPixels <- which(outterRing == TRUE, arr.ind = TRUE) # (row,col) in colImage

      if (length(data[, i_frame]$colonyProps[, , i_colony]$correspondingColPrevFrameInd) != 0 ) {
        prev_colName <- paste("f", i_frame - 1, "_c", data[, i_frame]$colonyProps[, , i_colony]$correspondingColPrevFrameInd, sep = "")
      } else {
        prev_colName <- "f0_c0"
      }

      if (data[, i_frame]$colonyProps[, , i_colony]$correspondingColNextFrameInd != 0 ) {
        next_colName <- paste("f", i_frame + 1, "_c", data[, i_frame]$colonyProps[, , i_colony]$correspondingColNextFrameInd, sep = "")
      } else {
        next_colName <- "f00_c0"
      }

      col_list[[curr_colony]] <- list(colName = paste("f", i_frame, "_c", i_colony, sep = ""),
                                      prev_colName = prev_colName,
                                      next_colName = next_colName,
                                      colImage = data[, i_frame]$colonyProps[, , i_colony]$bwColony,
                                      ULcorner = ULcorner,
                                      colBoundaryPixels = colBoundaryPixels,
                                      colCentroid = colCentroid)

      if (is.null(ncol(data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1]))) {  # one cell in the colony
        curr_cell <- curr_cell + 1
        cell_list[[curr_cell]] <- data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1]
        cell_list[[curr_cell]] <- updateCellInBascaList(x = cell_list[[curr_cell]],
                                                        frame = i_frame,
                                                        colony = i_colony,
                                                        colId = curr_colony,
                                                        col_list = col_list,
                                                        pixelRatio = pixelR,
                                                        updateDaughterIds = TRUE)
      } else {
        for (i_cell in 1:ncol(data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1])) {
          curr_cell <- curr_cell + 1
          cell_list[[curr_cell]] <- data[, i_frame]$colonyProps[, , i_colony]$cellProps[, , 1][, i_cell]
          cell_list[[curr_cell]] <- updateCellInBascaList(x = cell_list[[curr_cell]],
                                                          frame = i_frame,
                                                          colony = i_colony,
                                                          colId = curr_colony,
                                                          col_list = col_list,
                                                          pixelRatio = pixelR,
                                                          updateDaughterIds = TRUE)
        }
      }
    }
  }

  close(pb) ### close progress bar
  cat("\n")

  cat("Checking colony list...\n")
  checkColList(col_list = col_list, frameH = frameH, frameW = frameW)
  cat("Checking cell list...\n")
  cell_list <- checkCellList(cell_list = cell_list, col_list = col_list)

  Nframes <- length(unique(unlist(sapply(cell_list, function(x) x$frame))))
  Ncols <- length(unique(unlist(sapply(cell_list, function(x) x$colony)[unlist(sapply(cell_list, function(x) x$frame == 1))])))

  return(list(col_list = col_list, cell_list = cell_list,
              Nframes = Nframes, Ncols = Ncols,
              frameH = frameH, frameW = frameW))

}
