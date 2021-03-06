#' Import custom-made data
#'
#' Imports a \code{.json} file containing all the cell instants of the movie
#' and another \code{.json} file containing all the colony instants of the movie.
#'
#' The \code{file_cols} must be a JSON array of objects.
#' Each object denotes a colony instant of the movie and should contain the following key-value pairs
#' (otherwise an error is produced):
#' \itemize{
#' \item \code{colName} is the name of the colony instant, a character string in the format \code{"f<frame>_c<colony>"},
#' where \code{"<frame>"} and \code{"<colony>"} is the ID of the frame and colony (in the frame) of the colony instant, respectively.
#' \item \code{prev_colName} is a JSON array of character strings containing the \code{colName}
#' of the corresponding colony instant(s) in the previous frame.
#' In case one such colony instant exists, it can also be a character string.
#' For colony instants that do not have a corresponding colony instant in the previous frame, it should be equal to \code{"f0_c0"}.
#' \item \code{next_colName} is a character strings containing the \code{colName}
#' of the corresponding colony instant in the next frame.
#' For colony instants that do not have a corresponding colony instant in the next frame, it should be equal to \code{"f00_c0"}.
#' \item \code{colImage} is the mask of the box surrounding the colony instant, a JSON array of H arrays,
#' where H is the height of the box.
#' Each element of the array is a JSON array of W \code{0}s and/or \code{1}s,
#' where W is the width of the box.
#' \code{1}s denote the pixels of cells and \code{0}s the background pixels.
#' \cr\cr
#' NOTE: This key-value pair is not necessary.
#' \item \code{ULcorner} is a JSON array of 2 non-zero integer values
#' denoting the upper-left pixel of the box surrounding the colony instant in global (frame) coordinates.
#' The first integer represents the row and the second the column of the pixel.
#' \cr\cr
#' NOTE: This key-value pair is not necessary unless the key \code{colImage} is contained (an error is produced).
#' }
#' The \code{file_cells} must be a JSON array of objects.
#' Each object denotes a cell instant of the movie.
#' Key-value pair(s) denoting numeric or boolean attribute(s) should be contained in every object.
#' The following key-value pairs are also required to be contained in every object
#' (otherwise an error is produced):
#' \itemize{
#' \item \code{cellName} is the name of the cell, a character string in the format \code{"<cell>_f<frame>"}.
#' \item \code{frame} is the ID of the frame of the cell, a non-zero positive integer number.
#' \item \code{colony} is the ID of the colony of the cell in the \code{frame},
#' a non-zero positive integer number.
#' If the whole frame is treated as a single colony, value \code{1} must be used for all cells.
#' \item \code{daughterIds} is a JSON array of character strings containing the \code{cellName}
#' of the linked cell(s) in the next frame.
#' In case one such cell exists, it can also be a character string.
#' In case no such cells exist, it can either be \code{NULL} or omitted.
#' \item \code{colId} is a \emph{pointer} to the corresponding colony instant of the cell in the \code{file_cols},
#' a non-zero positive integer value.
#' \cr\cr
#' Colonies that entered the field of view at a time point and did not exist from the beginning of the movie (i.e. from the first frame)
#' should not have tracked cells, until they merge (if this is the case) with another existing colony.
#' This means that no object should \emph{point} to such colony instants.
#' \cr\cr
#' NOTE: This key-value pair is not necessary unless \code{file_cols != NULL} (an error is produced) and
#' should be omitted if \code{file_cols = NULL} (a warning is produced).
#' \item \code{pixelList} is a JSON array of arrays.
#' Each element of the array is a JSON array of 2 non-zero integer values,
#' indicating a pixel of the cell in colony coordinates
#' (i.e. relative to the \code{colImage} key of the \code{\ifelse{html}{\out{colId<sup>th</sup>}}{\eqn{colId^{th}}}}
#' object in \code{file_cols}).
#' The first value represents the row and the second the column of the pixel.
#' \cr\cr
#' NOTE: This key-value pair is not necessary.
#' It should be omitted if \code{file_cols = NULL} or if objects in \code{file_cols} do not contain the key \code{colImage}
#' (a warning is produced).
#' }
#'
#' @param file_cells,file_cols Character strings naming the \code{.json} files
#' (including the suffix \code{".json"}) containing all the cell and colony instants of the movie, respectively,
#' from which the data is to be imported.
#' If a string does not contain an absolute path, the file name is relative to the current working directory, \code{getwd()}.
#' See the \emph{Details} field for information about the format of these files.
#' \cr\cr
#' When the default value \code{file_cols = NULL} is used, no colony list is imported.
#' \cr\cr
#' NOTE: The components should be separated by \code{"/"} on Windows.
#'
#' @param frameH Frame image height in \emph{pixels}, a non-zero positive integer value.
#' This argument is ignored in case \code{file_cols = NULL}.
#' @param frameW Frame image width in \emph{pixels}, a non-zero positive integer value.
#' This argument is ignored in case \code{file_cols = NULL}.
#'
#' @return A named list with the following components:
#' \item{col_list}{A list containing all the colony instants of the movie,
#' or \code{NULL} if \code{file_cols = NULL}.
#' Each element of the list is a named list having as components the corresponding key-value pairs in the \code{file_cols}.}
#' \item{cell_list}{A list containing all the cell instants of the movie.
#' Each element of the list is a named list having as components the corresponding key-value pairs in the \code{file_cells}.
#' \cr\cr
#' In case the key \code{pixelList} is contained, \code{centroid} is also computed
#' as the \emph{mean} of the \code{pixelList} by column.
#' \code{centroid} is an 1x2 matrix of non-zero numeric values
#' denoting the centroid (geometric center) of the cell in colony coordinates
#' (i.e. relative to the \code{colImage} of the \code{\ifelse{html}{\out{colId<sup>th</sup>}}{\eqn{colId^{th}}}}
#' element in the \code{col_list}).}
#' \item{Nframes}{Number of frames in the movie, a non-zero positive integer value.
#' IDs of frames are in the range \code{[1, Nframes]}.}
#' \item{Ncols}{Number of colonies in the movie, a non-zero positive integer value.
#' IDs of colonies are in the range \code{[1, Ncols]}.
#' This value corresponds to the number of colonies at the start of the movie.}
#'
#' @export


import_json <- function(file_cells, file_cols = NULL, frameH, frameW) {

  if (!is.null(file_cols)) {
    col_list <- jsonlite::fromJSON(txt = file_cols, simplifyDataFrame = FALSE) # as nested list, not as data frame
	cat("Checking colony list...\n")
    checkColList(col_list = col_list, frameH = frameH, frameW = frameW)
  } else {
    col_list <- NULL
  }

  cell_list <- jsonlite::fromJSON(txt = file_cells, simplifyDataFrame = FALSE) # as nested list, not as data frame
  cat("Checking cell list...\n")
  cell_list <- checkCellList(cell_list = cell_list, col_list = col_list)

  Nframes <- length(unique(unlist(sapply(cell_list, function(x) x$frame))))
  Ncols <- length(unique(unlist(sapply(cell_list, function(x) x$colony)[unlist(sapply(cell_list, function(x) x$frame == 1))])))

  return(list(col_list = col_list, cell_list = cell_list,
              Nframes = Nframes, Ncols = Ncols))

}
