#' Plot single-cell growth curve of an attribute
#'
#' Plots the raw (and fitted) single-cell growth curve (time-series data) of a numeric attribute
#' for a set of given cells in a division tree.
#'
#' A separate plot for each cell specified in \code{cells} is generated.
#' \cr\cr
#' x-axis represents the life of the cell in \emph{frames},
#' starting from its birth (value \code{1}).
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
#' They can be any non-root cells, as returned from \code{\link{get_cells}}.
#' The default value \code{"all"} stands for all cells in the \code{DT}.
#'
#' @param model A character string naming the optional fitted growth curve(s) to be plotted,
#' if the corresponding parameters have been already successfully estimated
#' using \code{\link{add_attr_growth_fit_pars}}:
#' \itemize{
#' \item \code{""} (the empty string) for plotting just the raw growth curve
#' \item \code{"lin"} for plotting the fitted \emph{linear} model
#' \item \code{"exp"} for plotting the fitted \emph{exponential} model
#' \item \code{"both"} for plotting both fitted \emph{linear} and \emph{exponential} model
#' }
#'
#' @param frameR Frame rate of the movie in \emph{frames} per \emph{minute}, a non-zero positive numeric value.
#' This argument is ignored in case \code{model = ""}.
#'
#' @param save A logical value (\code{TRUE} or \code{FALSE}) indicating whether the generated plot(s)
#' will be saved in \code{.png} file(s) or displayed in the Plots Pane of RStudio, respectively.
#' The default value is \code{FALSE}.
#'
#' @param savePars A named list specifying the parameters of each generated image file.
#' This argument is ignored in case \code{save = FALSE}.
#' Elements of the list are the following:
#' \describe{
#' \item{\code{w}}{The width of the image file in \emph{pixels}, a non-zero positive integer value.
#' The default value is \code{1500}.}
#' \item{\code{h}}{The height of the image file in \emph{pixels}, a non-zero positive integer value.
#' The default value is \code{1000}.}
#' \item{\code{res}}{The resolution of the image file in \emph{pixels} per \emph{inch} (ppi), a non-zero positive integer value.
#' The smaller this value, the larger the plot area in inches, and the smaller the text relative to the graph itself.
#' The default value is \code{150}.}
#' \item{\code{path}}{A character string naming the directory where the image file will be saved (excluding the last \code{"/"}).
#' If it does not contain an absolute path, the image file will be saved relative to the current working directory \code{getwd()}.
#' The default value is the current working directory \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{"/"} on Windows.}
#' \item{\code{name}}{The image file name, a character string.
#' The suffix \code{".png"} is added automatically.
#' The default value is \code{"my_growth_attr_cell"}.}
#' }
#'
#' @export
#' @import igraph
#' @importFrom graphics axis legend lines par plot
#' @importFrom grDevices dev.off png
#' @importFrom stats as.formula

plot_growth_attr_cell <- function(DT, LT,
                                  attr, unit = "",
                                  cells = "all",
                                  model = c("", "lin", "exp", "both"),
                                  frameR,
                                  save = FALSE, savePars = list(w = 1500, h = 1000, res = 150, path = getwd(), name = "my_growth_attr_cell")) {

  ############ arguments check ###########################

  if (!(model %in% c("", "lin", "exp", "both"))) {
    stop("model must be \"\" / \"lin\" / \"exp\ / \"both\"\n")
  }

  numeric_attrs <- get_attr_names(tree = LT, type = "n")
  if (!(attr %in% numeric_attrs) || attr %in% c("colony", "generation", "frame", "age")) {
    stop(paste("Wrong attr \"", attr, "\"\n", sep = ""))
  }

  possible_cells <- get_cells(tree = DT, treeT = "DT", type = "nr")
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

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, "_%03d.png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  if (length(cells) == 1){
    if (cells == "all") {
      cells <- possible_cells
    }
  }

  oldPar <- par()
  par(mar = c(4.5, 4.5, 4, 1))

  linFormula <- as.formula("values ~ a * time + b")
  expFormula <- as.formula("values ~ value0 * exp(k * time)")

  for (cell in cells) {

    values <- as.numeric(unlist(strsplit(vertex_attr(graph = DT, name = attr, index = cell), ", ")))
    time_raw <- seq(1, V(DT)[cell]$lifeFrames, 1) # plot in frames
    time <- seq(1, V(DT)[cell]$lifeFrames, length.out = 1000) # plot in frames

    linear_plotted <- exponential_plotted <- FALSE

    if (model %in% c("lin", "both")) {

      linear_parameter <- paste(attr, "a", sep = "_")

      if (linear_parameter %in% vertex_attr_names(DT)) {
        if (!is.na(vertex_attr(graph = DT, name = linear_parameter, index = cell))) {

          a <- vertex_attr(graph = DT, name = paste(attr, "a", sep = "_"), index = cell) / (frameR * 60)
          b <- vertex_attr(graph = DT, name = paste(attr, "b", sep = "_"), index = cell)

          linear_plotted <- TRUE

        } else {
          warning(paste(attr, " of cell \"", cell, "\" does not fit a linear model\n", sep = ""))
        }
      } else {
        warning("Parameter(s) of linear model have not been estimated\nCall function add_attr_growth_fit_pars()\n")
      }

    }

    if (model %in% c("exp", "both")) {

      exp_parameter <- paste(attr, "k", sep = "_")

      if (exp_parameter %in% vertex_attr_names(DT)) {
        if (!is.na(vertex_attr(graph = DT, name = exp_parameter, index = cell))) {

          k <- vertex_attr(graph = DT, name = paste(attr, "k", sep = "_"), index = cell) / (frameR * 60)
          value0 <- vertex_attr(graph = DT, name = paste(attr, "0", sep = "_"), index = cell)

          exponential_plotted <- TRUE

        } else {
          warning(paste(attr, " of cell \"", cell, "\" does not fit an exponential model\n", sep = ""))
        }
      } else {
        warning("Parameter(s) of exponential model have not been estimated\nCall function add.attr.growth.fit.pars()\n")
      }

    }

    if (linear_plotted & exponential_plotted) {
      ylim <- range(c(values, eval(linFormula[[3]]), eval(expFormula[[3]])))
    } else if (linear_plotted) {
      ylim <- range(c(values, eval(linFormula[[3]])))
    } else if (exponential_plotted) {
      ylim <- range(c(values, eval(expFormula[[3]])))
    } else {
      ylim <- range(values)
    }

    plot(x = time_raw, y = values,
         main = paste("cell \"", cell, "\"", sep = ""), xlab = "Cell Life (frames)", ylab = createAxisLab(attr = attr, unit = unit),
         ylim = ylim,
         cex.main = 1.5, cex.lab = 1, font = 2, font.lab = 2, las = 1,
         type = "b", lty = 1, lwd = 2,
         pch = 16, col = "slateblue4", cex = 1.2)

    if (linear_plotted) {
      lines(time, eval(linFormula[[3]]), lwd = 3, col = "coral2")
    }

    if (exponential_plotted) {
      lines(time, eval(expFormula[[3]]), lwd = 3, col = "chartreuse3")
    }

    axis(1, at = time_raw, labels = time_raw, font = 2)

    if (linear_plotted & exponential_plotted) {
      legend("topleft", cex = 1,
             legend = c("fitted linear model", "fitted exponential model", "raw"),
             pch = c(NA, NA, 16),
             col = c("coral2", "chartreuse3", "slateblue4"),
             lwd = 3,
             xpd = TRUE, inset = 0.01)
    } else if (linear_plotted) {
      legend("topleft", cex = 1,
             legend = c("fitted linear model", "raw"),
             pch = c(NA, 16),
             col = c("coral2", "slateblue4"),
             lwd = 3,
             xpd = TRUE, inset = 0.01)
    } else if (exponential_plotted) {
      legend("topleft", cex = 1,
             legend = c("fitted exponential model", "raw"),
             pch = c(NA, 16),
             col = c("chartreuse3", "slateblue4"),
             lwd = 3,
             xpd = TRUE, inset = 0.01)
    } else {
      legend("topleft", cex = 1,
             legend = "raw",
             pch = 16,
             col = "slateblue4",
             lwd = 3,
             xpd = TRUE, inset = 0.01)
    }

  }

  par(mar = oldPar$mar)

  if (save) {
    dev.off()
  }

}
