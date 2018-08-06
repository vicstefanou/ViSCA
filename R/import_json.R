#' Import custom-made cell and colony list
#'
#' Imports a \code{.json} file containing all the cell instants of the movie
#' and another \code{.json} file containing all the colony instants of the movie.
#'

#'
#' The \code{file_cells} should be a JSON array of objects.
#' Each object denotes a cell instant of the movie and should contain the following key-value pairs:
#' \itemize{
#' \item \code{cellName} is the name of the cell, a character string in the format \code{"<cell>_f<frame>"}.
#' \item \code{frame} is the ID of the frame of the cell, a non-zero positive integer number.
#' \item \code{colony} is the ID of the colony of the cell in the \code{frame},
#' a non-zero positive integer number.
#' If the whole frame is treated as a single colony, value \code{1} should be used for all cells.
#' \item \code{daughterIds} is a JSON array of character strings containing the \code{cellName}
#' of the linked cell(s) in the next frame.
#' In case one such cell exists, it can also be a character string.
#' In case no such cells exist, it can either be \code{NULL} or omitted.
#' \item \code{colId} is a pointer to the corresponding colony instant of the cell in the \code{file_cols},
#' a non-zero positive integer value.
#' \cr\cr
#' NOTE: This key-value pair should be omitted if \code{file_cols = NULL}. A warning message is produced.
#'
#' \item \code{pixelList} is a JSON array of arrays.
#' Each element of the array is a JSON array of 2 non-zero integer values,
#' indicating a pixel of the cell in colony coordinates
#' (i.e. relative to the \code{colImage} of the \code{colId} element
#' in the colony list returned from \code{\link{import_col_list}}).
#' The first value represents the row and the second the column of the pixel.
#' \cr\cr
#' NOTE: This key-value pair should be omitted if \code{file_cols = NULL} or
#' if objects in \code{file_cols} do not contain the key \code{colImage}.
#' A warning message is produced.
#' }
#'

#'
#' \itemize{
#' \item \code{colName} is the name of the colony instant, a character string in the format \code{"f<frame>_c<colony>"},
#' where \code{"<frame>"} and \code{"<colony>"} is the ID of the frame and colony (in the frame) of the colony instant, respectively.
#'
#' \item \code{prev_colName} is a vector of character strings containing the \code{colName}
#' of the corresponding colony instant(s) in the previous frame,
#' or \code{NA} for colony instants of the first frame.
#'
#' \item \code{next_colName} is a vector of character strings containing the \code{colName}
#' of the corresponding colony instant(s) in the next frame.
#' For colony instants of the last frame this is equal to \code{"f<Nframes+1>_c0"}.
#'
#' \item \code{colImage} is the mask of the box surrounding the colony instant, a matrix of \code{0} and \code{1}.
#' \code{1}s denote the pixels of cells and \code{0}s the background pixels.
#'
#' \item \code{ULcorner} is a vector of 2 non-zero integer values
#' denoting the upper-left pixel of the box surrounding the colony instant in global (frame) coordinates.
#' The first integer represents the row and the second the column of the pixel.
#' }
#'
#'
#' @param file_cells,file_cols Character strings naming the \code{.json} files
#' (including the suffix \code{".json"}) containing all the cell and colony instants of the movie, respectively,
#' from which the data is to be imported.
#' If a string does not contain an absolute path, the file name is relative to the current working directory, \code{getwd()}.
#' \cr\cr
#' When the default value \code{file_cols = NULL} is used, no colony list is imported.
#' \cr\cr
#' NOTE: The components should be separated by \code{/} (not \code{\\}) on Windows.
#'
#' @return A named list with the following components:
#' \item{cell_list}{A list containing all the cell instants of the movie.
#' Each element of the list is a named list having as components the corresponding key-value pairs in the \code{file_cells}.}
#' \item{cell_list}{A list containing all the colony instants of the movie.
#' Each element of the list is a named list having as components the corresponding key-value pairs in the \code{file_cols}.}
#' \item{Nframes}{Number of frames in the movie, a non-zero positive integer value.
#' IDs of frames are in the range \code{[1, Nframes]}.}
#' \item{Ncols}{Number of colonies in the movie, a non-zero positive integer value.
#' IDs of colonies are in the range \code{[1, Ncols]}.
#' This value corresponds to the number of colonies at the start of the movie.}
#'
#' @export


import_json <- function(file_cells, file_cols = NULL) {

  ### col_list ###

  if (!is.null(file_cols)) {

    col_list <- jsonlite::fromJSON(txt = file_cols,
                                    simplifyDataFrame = FALSE) # as nested list, not as data frame

    if ((m <- checkField(myList = col_list, fieldName = "colName")) != "") {
      stop(paste("colony list:", m))
    } else {
      ###
    }

    if ((m <- checkField(myList = col_list, fieldName = "prev_colName")) != "") {
      stop(paste("colony list:", m))
    } else {
      ###
    }

    if ((m <- checkField(myList = col_list, fieldName = "next_colName")) != "") {
      stop(paste("colony list:", m))
    } else {
      ###
    }

    # if (checkField(myList = col_list, fieldName = "ULcorner") == "" &
    #     checkField(myList = col_list, fieldName = "colImage") == "") {
    #
    # }

  } else {
    col_list <- NULL
  }

  ### cell_list ###

  cell_list <- jsonlite::fromJSON(txt = file_cells,
                                  simplifyDataFrame = FALSE) # as nested list, not as data frame

  if ((m <- checkField(myList = cell_list, fieldName = "cellName")) != "") {
    stop(paste("cell list:", m))
  } else {
    ###
  }

  if ((m <- checkField(myList = cell_list, fieldName = "frame")) != "") {
    stop(paste("cell list:", m))
  } else {
    ###
  }

  if ((m <- checkField(myList = cell_list, fieldName = "colony")) != "") {
    stop(paste("cell list:", m))
  } else {
    ###
  }

  if ((m <- checkField(myList = cell_list, fieldName = "daughterIds")) != "") {
    stop(paste("cell list:", m))
  } else {
    ###
  }


  if (!is.null(col_list)) {
    if ((m <- checkField(myList = cell_list, fieldName = "colId")) != "") {
      stop(paste("colony list is imported\ncell list:", m))
    } else {
      ###
    }
  } else {
    if (checkField(myList = cell_list, fieldName = "colId") == "") {
      warning("colony list is not imported\n
              cell list: Component \"colId\" is useless and should be ommited\n")
    }
  }


  if (!is.null(col_list)) {
    if ((m <- checkField(myList = col_list, fieldName = "colImage")) == "") {
      if (checkField(myList = cell_list, fieldName = "pixelList") == "") {
        ###
        # centroid
      }
    } else {
      if (checkField(myList = cell_list, fieldName = "pixelList") == "") {
        warning(paste("colony list: ", m,
                      "cell list: Component \"pixelList\" is useless and should be ommited\n", sep = ""))
      }
    }
  } else {
    if (checkField(myList = cell_list, fieldName = "pixelList") == "") {
      warning("colony list is not imported\n
              cell list: Component \"pixelList\" is useless and should be ommited\n")
    }
  }


  Nframes <- length(unique(unlist(sapply(cell_list, function(x) x$frame))))
  Ncols <- length(unique(unlist(sapply(cell_list, function(x) x$colony)[unlist(sapply(cell_list, function(x) x$frame == 1))])))

  return(list(cell_list = cell_list, col_list = col_list,
              Nframes = Nframes, Ncols = Ncols))

}
