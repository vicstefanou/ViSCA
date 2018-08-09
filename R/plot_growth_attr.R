#' Plot raw growth curves of an attribute
#'
#' Plots the raw single-cell growth curves (unfitted time-series data) of a numeric attribute
#' for a set of given cells in a division tree.
#'
#' A common plot for all cells specified in \code{cells} is generated.
#' \cr\cr
#' x-axis represents the time in \emph{frames},
#' from the start (value \code{1}) to end of the movie (value \code{Nframes}).
#' \cr\cr
#' Each single-cell growth curve is randomly colored.
#'
#' @param DT The division tree where the cells specified in \code{cells} belong, an object of class \code{"igraph"}.
#'
#' @param LT The corresponding lineage tree of the \code{DT}, an object of class \code{"igraph"}.
#'
#' @param attr The name of the attribute in the \code{LT}, a character string.
#' It can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' except for \code{"colony"}, \code{"generation"}, \code{"frame"} and \code{"age"}.
#'
#' @param unit The unit of \code{attr}, a character string.
#' It should be in the format \code{"<string>,<number>"},
#' where \code{",<number>"} represents the power and is optional
#' (e.g. \code{"m"} for meters and \code{"cm,3"} for cubic centimeters).
#' The default value is the empty character \code{""}, which implies that \code{attr} is in arbitrary units.
#'
#' @param cells The labels of the cells in the \code{DT} whose growth curve will be plotted,
#' a vector of character strings.
#' They can be any cells that are included in the analysis, as returned from \code{\link{get_cells}}.
#' The default value \code{"all"} stands for all cells in the \code{DT}.
#'
#' @param Nframes Number of frames in the movie, a non-zero positive integer value.
#'
#' @param save A logical value (\code{TRUE} or \code{FALSE}) indicating whether the generated plot
#' will be saved in a \code{.png} file or displayed in the Plots Pane of RStudio, respectively.
#' The default value is \code{FALSE}.
#'
#' @param savePars A named list specifying the parameters of the generated image file.
#' This argument is ignored in case \code{save = FALSE}.
#' Elements of the list are the following:
#' \describe{
#' \item{\code{w}}{The width of the image file in \emph{pixels}, a non-zero positive integer value.
#' The default value is \code{2500}.}
#' \item{\code{h}}{The height of the image file in \emph{pixels}, a non-zero positive integer value.
#' The default value is \code{2000}.}
#' \item{\code{res}}{The resolution of the image file in \emph{pixels} per \emph{inch} (ppi), a non-zero positive integer value.
#' The smaller this value, the larger the plot area in inches, and the smaller the text relative to the graph itself.
#' The default value is \code{300}.}
#' \item{\code{path}}{A character string naming the directory where the image file will be saved (excluding the last \code{"/"}).
#' If it does not contain an absolute path, the image file will be saved relative to the current working directory \code{getwd()}.
#' The default value is the current working directory \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{"/"} on Windows.}
#' \item{\code{name}}{The image file name, a character string.
#' The suffix \code{".png"} is added automatically.
#' The default value is \code{"my_growth_attr"}.}
#' }
#'
#' @export
#' @import igraph
#' @importFrom graphics lines par plot
#' @importFrom grDevices colors dev.off png

plot_growth_attr <- function(DT, LT,
                             attr, unit = "",
                             cells = "all",
                             Nframes,
                             save = FALSE, savePars = list(w = 2500, h = 2000, res = 300, path = getwd(), name = "my_growth_attr")) {

  ############ arguments check ###########################

  numeric_attrs <- get_attr_names(tree = LT, type = "n")
  if (!(attr %in% numeric_attrs) || attr %in% c("colony", "generation", "frame", "age")) {
    stop(paste("Wrong attr \"", attr, "\"\n", sep = ""))
  }

  possible_cells <- get_cells(tree = DT, treeT = "DT", type = "inc")
  if (length(cells) != 1) {
    if (length(m <- setdiff(cells, possible_cells)) != 0) {
      stop(paste("Selected cell(s) ", toString(paste("\"", m, "\"", sep = "")), " do not exist\n", sep = ""))
    }
  } else {
    if (!(cells %in% c(possible_cells, "all"))) {
      stop(paste("Selected cell ", cells, " does not exist\n", sep = ""))
    }
  }

  ######################

  if (length(cells) == 1){
    if (cells == "all") {
      cells <- possible_cells
    }
  }

  if (length(cells) == 0) {
    stop("No cells in DT\n")
  }

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, ".png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  #### find x, y limits #######

  xLimits <- c(1, Nframes)
  yLimits <- range(as.numeric(unlist(strsplit(vertex_attr(graph = DT, name = attr, index = cells), ", "))))

  #############################

  # removed shades of gryy -> max 433 colors
  V(DT)[cells]$color <- sample(colors()[grep('gr(a|e)y', colors(), invert = T)], length(cells), replace = TRUE)

  oldPar <- par()
  par(mar = c(4.5, 4.5, 1, 1))

  cells_plotted <- NULL

  for (cell in cells) {

    values <- as.numeric(unlist(strsplit(vertex_attr(graph = DT, name = attr, index = cell), ", ")))
    time <- c(V(DT)[cell]$birthTime:V(DT)[cell]$divisionTime)

    if (is.null(cells_plotted)) { # first cell to plot

      plot(x = time, y = values,
           main = "", xlab = "Time (frames)", ylab = createAxisLab(attr = attr, unit = unit),
           xlim = xLimits, ylim = yLimits,
           cex.main = 1.5, cex.lab = 1, font = 2, font.lab = 2, las = 1,
           type = "l", lwd = 2, col = V(DT)[cell]$color,
           xaxt = "n")

      axis(1, at = c(1:Nframes), labels = c(1:Nframes), font = 2)

    } else {

      lines(x = time, y = values,
            xlim = xLimits, ylim = yLimits,
            lwd = 2, col = V(DT)[cell]$color)

    }
    cells_plotted <- c(cells_plotted, cell)
  }

  par(mar = oldPar$mar)

  if (save) {
    dev.off()
  }

}
