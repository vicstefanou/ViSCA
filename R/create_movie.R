#' Create movie
#'
#' Creates the movie of the segmented cells.
#'
#' A separate image for each frame existing in the \code{overallLT} is generated.
#'
#' @section Prerequisites:
#' This function can be used by \emph{BaSCA} or \emph{SuperSegger} users,
#' importing the data with \code{\link{import_basca}} or \code{\link{import_ss}}, respectively.
#' \cr\cr
#' If \code{\link{import_json}} was used for importing the data,
#' it is necessary that cell list elements have the \code{pixelList} and \code{colId} components
#' and colony list elements have the \code{ULcorner} component.
#' See \code{\link{import_json}} for more details.
#' In other case, this function cannot be used (throws an error).
#' \cr\cr
#' Users of \emph{Oufti} who imported the data with \code{\link{import_oufti}}
#' are also \bold{excluded} from using this function, as no colony list was returned.
#'
#' \cr\cr
#' \emph{FFmpeg} is \bold{required} in the system in order to automatically
#' convert the series of the generated \code{.png} files to an \code{.mp4} file,
#' at 5 fps (frames per second).
#' Image files for more than two consequtive frames starting from frame \code{1}
#' must have been generated.
#'
#' @param overallLT The lineage tree the cells of which will be visualized in the created movie,
#' an object of class \code{"igraph"}.
#' Cells that do not belong to the \code{LT} are colored white.
#' \cr\cr
#' For visualizing all cells extracted from the analysis of the movie,
#' use the tree returned from \code{\link{unite_trees}}
#' when called for all parts of the overall FLT (main part and possible motherless branches).
#'
#' @param LT The lineage tree the cells of which will be colored in the created movie,
#' an object of class \code{"igraph"}.
#' This tree must be a subtree of the \code{overallLT} regarding the included cells in the analysis.
#'
#' @param cell_list A list containing all the cell instants of the movie.
#' @param col_list A list containing all the colony instants of the movie.
#' @param frameH Frame image height in \emph{pixels}, a non-zero positive integer value.
#' @param frameW Frame image width in \emph{pixels}, a non-zero positive integer value.
#' @param Ncols Number of colonies in the movie, a non-zero positive integer value.
#'
#' @param attrC The name of the attribute in the \code{LT} by which the cells will be colored, a character string.
#' It can be any numeric or boolean attribute, as returned from \code{\link{get_attr_names}}.
#' Coloring is applied to all non-root cells of the \code{LT}, as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in this attribute which are colored gray.
#' When the default value \code{""} (the empty character) is used, all cells of the \code{LT} are colored gray.
#' Cells that belong to the \code{overallLT} but not to the \code{LT} are always colored white.
#'
#' @param unitC The unit of \code{attrC}, a character string.
#' It should be in the format \code{"<string>,<number>"},
#' where \code{",<number>"} represents the power and is optional
#' (e.g. \code{"m"} for meters and \code{"cm,3"} for cubic centimeters).
#' The default value is the empty character \code{""}, which implies that \code{attrC} is in arbitrary units.
#' This argument is ignored in case \code{attrC = ""}.
#'
#' @param NC Number of colonies in the movie (if \code{attrC = "colony"}) or
#' number of generations in the movie (if \code{attrC = "generation"}), a non-zero positive integer value.
#' This argument is ignored in case \code{attrC != "colony"}
#' and \code{attrC != "generation"}.
#'
#' @param savePath  A character string naming the absolute path of the directory
#' where the generated \code{.png} files will be saved.
#' The default value is the current working directory \code{getwd()}.
#' Image files are named as \code{"frame_<i>.png"}, where \code{"<i>"} is the frame ID.
#' The movie file is named as \code{"movie.mp4"}.
#' All files will be saved in folder \code{"<attrC>"}
#' (or \code{"noColor"} in case \code{attrC = ""}), created under the specified directory.
#' \cr\cr
#' NOTE: The components should be separated by \code{/} (not \code{\\}) on Windows.
#'
#' @seealso \code{\link{isSubtree}} for checking if a tree is a subtree of another tree.
#' @export
#' @import igraph
#' @importFrom graphics layout par plot rasterImage
#' @importFrom grDevices dev.off png
#' @importFrom utils setTxtProgressBar txtProgressBar

create_movie <- function(overallLT, LT,
                         cell_list, col_list,
                         frameH, frameW,
                         Ncols,
                         attrC = "", unitC = "", NC = NULL,
                         savePath = getwd()) {

  ############ arguments check ###########################

  if (length(get_cells(tree = overallLT, treeT = "LT", type = "nr")) == 0) {
    stop("overallLT is empty\n")
  }

  if (length(get_cells(tree = LT, treeT = "LT", type = "nr")) == 0) {
    stop("LT is empty\n")
  }

  if (!isSubtree(subtree = LT, tree = overallLT, treeT = "LT", type = "inc")) {
    stop("LT is not a subtree of overallLT")
  }

  numeric_attrs <- get_attr_names(tree = overallLT, type = "n")
  boolean_attrs <- get_attr_names(tree = overallLT, type = "b")

  if (!(attrC %in% c(numeric_attrs, boolean_attrs, ""))) {
    stop(paste("Wrong attrC \"", attrC, "\"\n", sep = ""))
  }

  if (attrC %in% c("colony", "generation") & is.null(NC)) {
    stop("Specify number of colors NC\n")
  }

  ### prerequisites
  if ((m <- checkField(myList = cell_list, fieldName = "pixelList")) != "") {
    stop(paste("cell list:", m))
  }
  if ((m <- checkField(myList = cell_list, fieldName = "colId")) != "") {
    stop(paste("cell list:", m))
  }
  if ((m <- checkField(myList = col_list, fieldName = "ULcorner")) != "") {
    stop(paste("colony list:", m))
  }

  ######################
  oldDir <- getwd()

  if (attrC != "") {
    dir.create(file.path(savePath, attrC), showWarnings = FALSE)
    setwd(file.path(savePath, attrC))
  } else {
    dir.create(file.path(savePath, "noColor"), showWarnings = FALSE)
    setwd(file.path(savePath, "noColor"))
  }


  ######################### Change cell colors by a value

  V(overallLT)$color <- "snow3"

  LT_cells  <- get_cells(tree = LT, treeT = "LT", type = "nr")
  V(overallLT)[LT_cells]$color <- "gray30"

  if (attrC != "") {
    colored_cells  <- V(LT)[V(LT)$name %in% LT_cells &
                              !is.na(vertex_attr(graph = LT, name = attrC, index = V(LT)))]$name

    if (is.null(cellColors <- getCellColors(tree = LT, cells = colored_cells, attr = attrC, N_colors = NC))) {
      attrC <- ""
    } else {
      V(overallLT)[colored_cells]$color <- cellColors
    }
  }

  ############################

  possible_frames <- sort(unique(V(overallLT)$frame))
  possible_frames <- possible_frames[possible_frames != -1]

  pb <- txtProgressBar(min = 0, max = length(possible_frames), style = 3) ### set progress bar
  ipb <- 0

  for (i_frame in possible_frames) {

    ipb <- ipb + 1
    setTxtProgressBar(pb, ipb) ### update progress bar

    im <- matrix("black", nrow = frameW, ncol = frameH)

    frame_cells <- V(overallLT)[V(overallLT)$frame == i_frame]$name

    for (cell in frame_cells) {

      cell_list_ID <- as.numeric(cell) - (Ncols + 1)

      cell_pixels <- cell_list[[cell_list_ID]]$pixelList
      ULcorner <- col_list[[cell_list[[cell_list_ID]]$colId]]$ULcorner

      pixels <- cbind(cell_pixels[,1] + ULcorner[1], cell_pixels[,2] + ULcorner[2])
      im[pixels] <- V(overallLT)[cell]$color

    }

    png(filename = paste("frame_", i_frame, ".png", sep = ""),
        width = 750 * 2,
        height = round(750 * frameH / frameW) * 2,
        res = 200)

    if (attrC != "") {
      if (attrC %in% c("colony", "generation", boolean_attrs)) {

        par(mar = c(0, 0, 0, 4.2), bg = "black", fg = "white")
        plot(NA, xlim = 0:1, ylim = 0:1, bty = "n", axes = 0, xaxs = 'i', yaxs = 'i')
        rasterImage(EBImage::Image(data = im, colormode = 'Color'), 0, 0, 1, 1)

        whichColors <- sort(unique(vertex_attr(graph = tree, name = attrC, index = colored_cells)))
        plotColorLegend(attr = attrC, whichColors = whichColors, N_colors = NC, size_lab = 0.6)

      } else {

        layout(matrix(c(1, 2,
                        1, 2), nrow = 2, byrow = TRUE), widths = c(10, 1))

        par(mar = rep(0, 4), bg = "black", fg = "white")
        plot(NA, xlim = 0:1, ylim = 0:1, bty = "n", axes = 0, xaxs = 'i', yaxs = 'i')
        rasterImage(EBImage::Image(data = im, colormode = 'Color'), 0, 0, 1, 1)

        plotColormap(values = vertex_attr(graph = LT, name = attrC, index = colored_cells),
                     attr = attrC, unit = unitC,
                     size_lab = 0.8, size_values = 0.6)

      }
    } else {
      par(mar = rep(0, 4), bg = "black", fg = "white")
      plot(NA, xlim = 0:1, ylim = 0:1, bty = "n", axes = 0, xaxs = 'i', yaxs = 'i')
      rasterImage(EBImage::Image(data = im, colormode = 'Color'), 0, 0, 1, 1)
    }


    dev.off()

  }

  close(pb) ### close progress bar
  cat("\n")

  cat(paste("Image files were saved in ", getwd(), "/frame_<i>.png", "\n", sep = ""))

  ############################ create movie  ####################################

  if (length(possible_frames) >= 2) { # multiple frames

    if (possible_frames[1] == 1) { # start from frame 1

      diff_f <- diff(possible_frames)

      if (length(unique(diff_f)) == 1) { # consequtive frames by a value

        if (unique(diff_f) == 1) { # consequtive frames by 1

          system(paste("ffmpeg", "-r", 5, "-i", "frame_%d.png", "-vcodec libx264 -crf 15 -pix_fmt yuv420p", paste(getwd(), "/movie.mp4", sep = "")))

          cat(paste("Movie file was saved in ", getwd(), "/movie.mp4", "\n", sep = ""))

        } else {
          warning("Unable to create movie\nFrames in overallLT are not consequtive\n")
        }
      } else {
        warning("Unable to create movie\nFrames in overallLT are not consequtive\n")
      }

    } else {
      warning("Unable to create movie\nFrame 1 does not exist in overallLT\n")
    }
  } else {
    warning("Unable to create movie\nSingle frame exists in overallLT\n")
  }

  ################################################################################

  setwd(oldDir)

}
