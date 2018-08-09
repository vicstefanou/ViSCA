#' Plot a tree
#'
#' Plots a lineage or division tree.
#'
#' Nodes represent all cells, as returned from \code{\link{get_cells}}.
#' By default, all nodes are white dot points and edges are colored darkgray.
#'
#' @param tree The connected lineage or division tree, an object of class \code{"igraph"}.
#'
#' @param treeT A character string naming the type of \code{tree}:
#' \itemize{
#' \item \code{"LT"} if \code{tree} is a lineage tree
#' \item \code{"DT"} if \code{tree} is a division tree
#' }
#'
#' @param attrC The name of the attribute in the \code{tree} by which the cells will be colored, a character string.
#' It can be any numeric or boolean attribute, as returned from \code{\link{get_attr_names}}.
#' Coloring is applied to all non-root cells, as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in this attribute which are colored gray.
#' Any existing imaginary \emph{root} cells are always colored white.
#' When the default value \code{""} (the empty character) is used, nodes' color is the default.
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
#' @param attrS The name of the attribute in the \code{tree} by which the cells will be shaped, a character string.
#' It can be any boolean attribute, as returned from \code{\link{get_attr_names}}.
#' Cells are represented by dots or squares if their value in this attribute is \code{TRUE} or \code{FALSE}, respectively.
#' When the default value \code{""} (the empty character) is used, nodes' shape is the default.
#'
#' @param cellsC The labels of the cells in the \code{tree} which will be colored red,
#' a vector of character strings.
#' They can be any cells, as returned from \code{\link{get_cells}}.
#' The default value is the empty character \code{""}, which stands for no cell.
#'
#' @param colorCol A logical value (\code{TRUE} or \code{FALSE}) indicating whether
#' the edges of the tree will be colored based on the colony or not, respectively.
#' The value \code{TRUE} can be used only if \code{tree} contains the imaginary main \emph{root} cell.
#' When the default value \code{FALSE} is used, edges' color is the default.
#'
#' @param Ncols Number of colonies in the movie, a non-zero positive integer value.
#' This argument is ignored in case \code{colorCol = FALSE}.
#'
#' @param showLabels A logical value (\code{TRUE} or \code{FALSE}) indicating whether
#' the labels of the cells will be shown on the plot or not, respectively.
#' The default value is \code{FALSE}.
#'
#' @param sizeLabel Size of cell labels, a non-zero positive numeric value.
#' The default value is \code{0.1}.
#' This argument is ignored in case \code{showLabels = FALSE}.
#'
#' @param circular A logical value (\code{TRUE} or \code{FALSE}) indicating whether
#' the tree will be plotted in circular or classical tree layout, respectively.
#' The default value is \code{FALSE}.
#'
#' @param showLegends A logical value (\code{TRUE} or \code{FALSE}) indicating whether
#' the explanatory legends and title will be shown on the plot or not, respectively.
#' The default value is \code{FALSE}.
#'
#' @param sizeL Size of explanatory legends and title, a non-zero positive numeric value.
#' The default value is \code{1}.
#' This argument is ignored in case \code{showLegends = FALSE}.
#'
#' @param sizeV Size of vertices, a non-zero positive numeric value.
#' The default value is \code{0.5}.
#'
#' @param sizeE Width of edges, a non-zero positive numeric value.
#' The default value is \code{0.1}.
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
#' The default value is \code{6000}.}
#' \item{\code{h}}{The height of the image file in \emph{pixels}, a non-zero positive integer value.
#' The default value is \code{6000}.}
#' \item{\code{res}}{The resolution of the image file in \emph{pixels} per \emph{inch} (ppi), a non-zero positive integer value.
#' The smaller this value, the larger the plot area in inches, and the smaller the text relative to the graph itself.
#' The default value is \code{600}.}
#' \item{\code{path}}{A character string naming the directory where the image file will be saved (excluding the last \code{"/"}).
#' If it does not contain an absolute path, the image file will be saved relative to the current working directory \code{getwd()}.
#' The default value is the current working directory \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{"/"} on Windows.}
#' \item{\code{name}}{The image file name, a character string.
#' The suffix \code{".png"} is added automatically.
#' The default value is \code{"my_tree"}.}
#' }
#'
#' @seealso \code{\link{isConnected}} for checking if a tree is connected.
#' @export
#' @import igraph
#' @importFrom graphics layout legend par text
#' @importFrom grDevices dev.off png

plot_tree <- function(tree, treeT = c("LT", "DT"),
                      attrC = "", unitC = "", NC = NULL,
                      attrS = "",
                      cellsC = "",
                      colorCol = FALSE, Ncols,
                      showLabels = FALSE, sizeLabel = 0.1,
                      circular = TRUE,
                      showLegends = TRUE, sizeL = 1,
                      sizeV = 0.5, sizeE = 0.1,
                      save = FALSE, savePars = list(w = 6000, h = 6000, res = 600, path = getwd(), name = "my_tree")) {

  ################## arguments check #######################

  if (length(get_cells(tree = tree, treeT = treeT, type = "nr")) == 0) {
    stop("tree is empty\n")
  }

  if (!isConnected(tree = tree)) {
    stop("tree is disconnected\n")
  }

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
  }

  numeric_attrs <- get_attr_names(tree = tree, type = "n")
  boolean_attrs <- get_attr_names(tree = tree, type = "b")

  if (!(attrC %in% c(numeric_attrs, boolean_attrs, ""))) {
    stop(paste("Wrong attrC \"", attrC, "\"\n", sep = ""))
  }

  if (attrC %in% c("colony", "generation") & is.null(NC)) {
    stop("Specify number of colors NC\n")
  }

  if (colorCol == TRUE & !("1" %in% V(tree)$name)) {
    stop("Root is missing from tree\nSet colorCol to FALSE\n")
  }

  if (!(attrS %in% c(boolean_attrs, ""))) {
    stop(paste("Wrong attrS \"", attrS, "\"\n", sep = ""))
  }

  possible_cells <- get_cells(tree = tree, type = "all")
  if (length(cellsC) != 1) {
    if (length(m <- setdiff(cellsC, possible_cells)) != 0) {
      stop(paste("Selected cell(s) ", toString(paste("\"", m, "\"", sep = "")), " do not exist\n", sep = ""))
    }
  } else {
    if (!(cellsC %in% c(possible_cells, ""))) {
      stop(paste("Selected cell ", cellsC, " does not exist\n", sep = ""))
    }
  }

  #######################

  myPlotTree <- function(mar) {

    oldPar <- par()
    par(mar = mar)

    plot.igraph(x = tree,
                layout = layout_as_tree(tree, circular = circular),
                edge.width = sizeE,
                edge.color = E(tree)$color,
                edge.arrow.size = 0,
                vertex.shape = V(tree)$shape,
                vertex.size = sizeV,
                vertex.color = V(tree)$color,
                vertex.frame.color = NA,
                vertex.label = ifelse(showLabels, V(tree)$name, NA),
                vertex.label.cex = ifelse(showLabels, sizeLabel, NA),
                vertex.label.font = ifelse(showLabels, 2, NA),
                vertex.label.color = ifelse(showLabels, "white", NA))

    if (showLegends == TRUE) {
      corners <- par("usr")
      if (treeT == "LT") {
        text(x = mean(corners[1:2]), y = corners[4], labels = expression(bold("Lineage Tree")), cex = 1.5 * sizeL)
      } else { # treeT == "DT"
        text(x = mean(corners[1:2]), y = corners[4], labels = expression(bold("Division Tree")), cex = 1.5 * sizeL)
      }
    }

    par(mar = oldPar$mar)

  }

  #######################

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, ".png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  if (attrS == "") {
    V(tree)$shape <- "circle"
  } else {
    shapes <- c("square", "circle")
    V(tree)$shape <- shapes[as.numeric(vertex_attr(graph = tree, name = attrS, index = V(tree))) + 1] # FALSE -> 1, TRUE -> 2
  }

  if (attrC == "") {
    V(tree)$color <- "snow3"
  } else {
    V(tree)$color <- "gray30"

    root_cells <- V(tree)[!(V(tree)$name %in% get_cells(tree = tree, treeT = treeT, type = "nr"))]$name
    V(tree)[root_cells]$color <- "snow3"

    colored_cells <- V(tree)[V(tree)$name %in% get_cells(tree = tree, treeT = treeT, type = "inc") &
                               !is.na(vertex_attr(graph = tree, name = attrC, index = V(tree)))]$name

    if (is.null(cellColors <- getCellColors(tree = tree, cells = colored_cells, attr = attrC, N_colors = NC))) {
      V(tree)$color <- "snow3"
      attrC <- ""
    } else {
      V(tree)[colored_cells]$color <-  cellColors
    }

  }

  if (length(cellsC) == 1) {
    if (cellsC != "") {
      V(tree)[cellsC]$color <- "red"
    }
  } else {
    V(tree)[cellsC]$color <- "red"
  }


  E(tree)$color <- "darkgray"

  if (colorCol == TRUE) {

    edgeColors <- getGroupColors(type = "colony", N_colors = Ncols)

    possible_cols <- sort(unique(V(tree)$colony))
    possible_cols <- possible_cols[possible_cols != -1]

    for (i_colony in possible_cols) {

      edges_to <- subcomponent(graph = tree, v = as.character(i_colony + 1), mode = "out")$name
      edges_from <- V(tree)[unlist(ego(graph = tree, order = 1, nodes = edges_to, mode = "in", mindist = 1))]$name

      ei <- get.edge.ids(graph = tree, vp = c(t(cbind(edges_from, edges_to))))
      E(tree)[ei]$color <- edgeColors[i_colony]

    }
  }


  oldPar <- par()
  par(bg = "black", fg = "white")

  if (showLegends == TRUE) {

    if (attrC != "" || attrS != "") { # color or shape

      if (attrC %in% c("colony", "generation", boolean_attrs, "")) {

        myPlotTree(mar = c(0, 0, 2 * sizeL, 7 * sizeL))

        if (attrS == "") { # only color

          whichColors <- sort(unique(vertex_attr(graph = tree, name = attrC, index = colored_cells)))
          plotColorLegend(attr = attrC, whichColors = whichColors, N_colors = NC, size_lab = sizeL)

        } else if (attrC == "") { # only shape

          corners <- par("usr")
          legend(x = corners[2], y = mean(corners[3:4]),
                 title =  as.expression(bquote(bold(.(attrS)))),
                 legend = c("FALSE", "TRUE"),
                 pch = c(15, 16),
                 yjust = 0.5,
                 cex = sizeL)

        } else { # color and shape

          whichColors <- sort(unique(vertex_attr(graph = tree, name = attrC, index = colored_cells)))
          lp <- plotColorLegend(attr = attrC, N_colors = NC, size_lab = sizeL)

          legend(x = mean(c(lp$rect$left, lp$rect$left + lp$rect$w)), y = lp$rect$top - lp$rect$h,
                 title = as.expression(bquote(bold(.(attrS)))),
                 legend = c("FALSE", "TRUE"),
                 pch = c(15, 16),
                 yjust = 1, xjust = 0.5,
                 cex = sizeL)

        }

      } else {

        def.par <- par(no.readonly = TRUE)

        if (attrS == "") { # only color

          layout(matrix(c(1, 2,
                          1, 2), nrow = 2, byrow = TRUE), widths = c(10, 1.5 * sizeL))

          myPlotTree(mar = c(0, 0, 2 * sizeL, 0))

          plotColormap(values = vertex_attr(graph = tree, name = attrC, index = colored_cells),
                       attr = attrC, unit = unitC,
                       size_lab = sizeL, size_values = 0.7 * sizeL)

        } else { # color and shape

          layout(matrix(c(1, 2,
                          1, 2), nrow = 2, byrow = TRUE), widths = c(10, 3 * sizeL))

          myPlotTree(mar = c(0, 0, 2 * sizeL, 0))

          plotColormap(values = vertex_attr(graph = tree, name = attrC, index = colored_cells),
                       attr = attrC, unit = unitC,
                       size_lab = sizeL, size_values = 0.7 * sizeL)

          corners <- par("usr")
          legend(x = mean(corners[1:2]), y = corners[3] - 2 ,
                 title = as.expression(bquote(bold(.(attrS)))),
                 legend = c("FALSE", "TRUE"), pch = c(15, 16),
                 yjust = 1, xjust = 0.5,
                 cex = sizeL)

        }

        par(def.par)

      }

    } else { # no color or shape (only title)

      myPlotTree(mar = c(0, 0, 2 * sizeL, 0))

    }

  } else { # no legends and title

    myPlotTree(mar = rep(0, 4))

  }

  par(bg = oldPar$bg, fg = oldPar$fg)

  if (save) {
    dev.off()
  }

}
