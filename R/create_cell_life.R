#' Visualize cell life
#'
#' Creates the life images of a set of given cells.
#'
#' Separate life image(s) for each cell specified in \code{cells} are generated.
#' The cell is viewed in its colony for every frame of its lifespan and is marked as red.
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
#' @param DT The division tree where the cells specified in \code{cells} belong, an object of class \code{"igraph"}.
#'
#' @param cells The labels of the cells in the \code{DT} whose life images will be created,
#' a vector of character strings.
#' They can be any non-root cells, as returned from \code{\link{get_cells}}.
#' The default value \code{"all"} stands for all non-root cells.
#'
#' @param cell_list A list containing all the cell instants of the movie.
#' @param col_list A list containing all the colony instants of the movie.
#' @param Ncols Number of colonies in the movie, a non-zero positive integer value.
#'
#' @param saveT A character string naming the options of image saving:
#' \itemize{
#' \item \code{"sep"} saves a separate \code{.png} file for each frame of the cell's lifespan,
#' named as \code{"<cell>_f<first_frame>_<i>.png"}.
#' Image files for each cell will be saved in a directory named after the cell.
#' \item \code{"comb"} saves a single \code{.png} file with all frames of the cell's lifespan combined
#' in \emph{left-to-right, top-to-bottom} order, named as \code{"<cell>_f<first_frame>.png"}.
#' }
#' \code{"<cell>"} is the label of the cell, \code{"<first_frame>"} is its birth frame and
#' \code{"<i>"} is the instant of the cell starting from \code{1}.
#'
#' @param savePath  A character string naming the absolute path of the directory
#' where the cell life images will be saved (excluding the last \code{"/"}).
#' The default value is the current working directory \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{"/"} on Windows.
#'
#' @export
#' @import igraph
#' @importFrom grDevices dev.off png
#' @importFrom utils setTxtProgressBar txtProgressBar


create_cell_life <- function(DT, cells = "all",
                             cell_list, col_list, Ncols,
                             saveT = c("sep", "comb"), savePath = getwd()) {

  ############ arguments check ###########################

  if (!(saveT %in% c("sep", "comb"))) {
    stop("saveT must be \"sep\" / \"comb\"\n")
  }

  possible_cells <- get_cells(tree = DT, treeT = "DT", type = "nr")

  if (length(cells) != 1) {
    if (length(m <- setdiff(cells, possible_cells)) != 0) {
      stop(paste("Selected cell(s)", toString(paste("\"", m, "\"", sep = "")), "do not exist\n"))
    }
  } else {
    if (!(cells %in% c(possible_cells, "all"))) {
      stop(paste("Selected cell", cells, "does not exist\n"))
    }
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

  centerImage <- function(image, height, width) {

    extra_value <- "black"

    if ((width - ncol(image)) %% 2 == 0) {
      extra_columns_left <- extra_columns_right <- matrix(extra_value, nrow = nrow(image), ncol = (width - ncol(image)) / 2)
    } else {
      extra_columns_left <- matrix(extra_value, nrow = nrow(image), ncol = (width - ncol(image)) %/% 2)
      extra_columns_right <- matrix(extra_value, nrow = nrow(image), ncol = (width - ncol(image)) %/% 2 + 1)
    }

    centered_im <- cbind(extra_columns_left, image, extra_columns_right)

    if ((height - nrow(image)) %% 2 == 0) {
      extra_rows_up <- extra_rows_down <- matrix(extra_value, nrow = (height - nrow(image)) / 2, ncol = width)
    } else {
      extra_rows_up <- matrix(extra_value, nrow = (height - nrow(image)) %/% 2, ncol = width)
      extra_rows_down <- matrix(extra_value, nrow = (height - nrow(image)) %/% 2 + 1, ncol = width)
    }

    centered_im <- rbind(extra_rows_up, centered_im, extra_rows_down)

    return(centered_im)

  }

  ######################

  if (length(cells) == 1){
    if (cells == "all") {
      cells <- possible_cells
    }
  }

  oldDir <- getwd()
  if (saveT == "sep") {
    dir.create(file.path(savePath, cell), showWarnings = FALSE)
    setwd(file.path(savePath, cell))
  } else {
    setwd(savePath)
  }

  pb <- txtProgressBar(min = 0, max = length(cells), style = 3) ### set progress bar
  ipb <- 0

  for (cell in cells) {

    ipb <- ipb + 1
    setTxtProgressBar(pb, ipb) ### update progress bar

    first_frame <- unlist(strsplit(V(DT)[cell]$frame, ", "))[1]

    LT_IDs <- as.numeric(unlist(strsplit(V(DT)[cell]$cellInstants, ", ")))

    cell_list_IDs <- LT_IDs - (Ncols + 1)

    ######################### find max x,y

    colony_list_IDs <- unname(unlist(sapply(cell_list, function(x) x$colId)))[cell_list_IDs]
    N_height_pixels <- max(unname(unlist(sapply(col_list, function(x) nrow(x$colImage))))[colony_list_IDs])
    N_width_pixels <- max(unname(unlist(sapply(col_list, function(x) ncol(x$colImage))))[colony_list_IDs])

    #########################

    life <- list()
    i <- 0

    for (cell_instant in cell_list_IDs) {

      i <- i + 1

      cell_pixels <- cell_list[[cell_instant]]$pixelList
      im <- col_list[[cell_list[[cell_instant]]$colId]]$colImage

      im[im == 0] <- "black"
      im[im == 1] <- "snow3"
      im[cell_pixels] <- "red"

      centered_im <- centerImage(image = im, height = N_height_pixels, width = N_width_pixels)

      life[[i]] <- EBImage::Image(data = centered_im, colormode = 'Color')

      if (saveT == "sep") {
        png(filename = paste(cell, "_f", first_frame, "_", i, ".png", sep = ""),
            width = N_width_pixels * 2,
            height = N_height_pixels * 2)
        EBImage::display(x = life[[i]], method = "raster", all = TRUE)
        dev.off()
      }

    }

    if (saveT == "comb") {
      png(filename = paste(cell, "_f", first_frame, ".png", sep = ""),
          width = 1500,
          height = 1500)
      EBImage::display(x = EBImage::combine(life), method = "raster", all = TRUE)
      dev.off()
    }

  }

  close(pb) ### close progress bar
  cat("\n")

  if (saveT == "sep") {
    cat(paste("Image files were saved in ", getwd(), "/<cell>/<cell>_f<first_frame>_<i>.png", "\n", sep = ""))
  } else { #saveT == "comb"
    cat(paste("Image files were saved in ", getwd(), "/<cell>_f<first_frame>.png", "\n", sep = ""))
  }

  setwd(oldDir)

}
