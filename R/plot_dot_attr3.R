#' Create scatter plot of three attributes
#'
#' Creates the X-Y-Z scatter plot of three numeric attributes of a lineage or division tree.
#'
#' Each data point (x,y,z) represents the value of \code{attr1}, \code{attr2} and \code{attr3} of a cell, respectively.
#' The horizontal axis of the plot corresponds to the x-axis, the vertical to the z-axis and the diagonal to the y-axis.
#' The scatter plot is created for all cells that are included in the analysis,
#' as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in \code{attr1} or \code{attr2} or \code{attr3}.
#' \cr\cr
#' By default, all data points are purple dot points.
#' \cr\cr
#' The linear regression plane is also drawn on the plot.
#' X, Y are the predictor variables and Z is the response.
#' The parameters of the regression plane are found using the \emph{linear least squares} method
#' provided by \code{\link[stats]{lm}}.
#'
#' @param tree The lineage or division tree, an object of class \code{"igraph"}.
#'
#' @param treeT A character string naming the type of \code{tree}:
#' \itemize{
#' \item \code{"LT"} if \code{tree} is a lineage tree
#' \item \code{"DT"} if \code{tree} is a division tree
#' }
#'
#' @param attr1,attr2,attr3 The names of the attributes in the \code{tree}, character strings.
#' Each one can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' except for \code{"colony"}, \code{"generation"} and \code{"frame"}.
#'
#' @param unit1,unit2,unit3 The units of the corresponding \code{attr}, character strings.
#' Each one should be in the format \code{"<string>,<number>"},
#' where \code{",<number>"} represents the power and is optional
#' (e.g. \code{"m"} for meters and \code{"cm,3"} for cubic centimeters).
#' The default value of each is the empty character \code{""},
#' which implies that the corresponding \code{attr} is in arbitrary units.
#'
#' @param attrC The name of the attribute in the \code{tree} by which the cells will be colored, a character string.
#' It can be any numeric or boolean attribute, as returned from \code{\link{get_attr_names}}.
#' Coloring is applied to all depicted cells,
#' except for cells with \code{NA} value in this attribute which are colored gray.
#' When the default value \code{""} (the empty character) is used, data points' color is the default.
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
#' @param sizeL Size of explanatory legends, a non-zero positive numeric value.
#' The default value is \code{0.8}.
#' This argument is ignored in case \code{attrC = ""} or
#' \code{attrC} is a boolean attribute or
#' \code{attrC = "colony"} or \code{attrC = "generation"}.
#'
#' @param attrS The name of the attribute in the \code{tree} by which the cells will be shaped, a character string.
#' It can be any boolean attribute, as returned from \code{\link{get_attr_names}}.
#' Cells are represented by dots or squares if their value in this attribute is \code{TRUE} or \code{FALSE}, respectively.
#' When the default value \code{""} (the empty character) is used, data points' shape is the default.
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
#' The default value is \code{350}.}
#' \item{\code{path}}{A character string naming the directory where the image file will be saved (excluding the last \code{"/"}).
#' If it does not contain an absolute path, the image file will be saved relative to the current working directory \code{getwd()}.
#' The default value is the current working directory \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{"/"} on Windows.}
#' \item{\code{name}}{The image file name, a character string.
#' The suffix \code{".png"} is added automatically.
#' The default value is \code{"my_dot_attr3"}.}
#' }
#'
#' @return A named list with the following components:
#' \item{Ncells}{Number of cells, a non-zero positive integer value.}
#' \item{regression}{A named list with the following components:
#' \itemize{
#' \item \code{a1} is the slope of the regression plane in the direction of x, a numeric value
#' \item \code{a2} is the slope of the regression plane in the direction of y, a numeric value
#' \item \code{b} is the y-intercept of the regression plane, a numeric value
#' \item \code{r2} is the R-squared coefficient of the regression plane
#' as returned from \code{\link[stats]{summary.lm}}, a numeric value in the range \code{[0, 1]}
#' }
#' In case less than 3 unique data points exist, \code{NULL} is returned, instead.}
#' In case no cells exist, no plot is generated and \code{NULL} is returned.
#'
#' @export
#' @import igraph
#' @importFrom graphics layout legend lines mtext par
#' @importFrom grDevices dev.off png
#' @importFrom stats na.omit

plot_dot_attr3 <- function(tree, treeT = c("LT", "DT"),
                           attr1, unit1 = "", attr2, unit2 = "", attr3, unit3 = "",
                           attrC = "", unitC = "", NC = NULL,
                           attrS = "",
                           sizeL = 0.8,
                           save = FALSE, savePars = list(w = 2500, h = 2000, res = 350, path = getwd(), name = "my_dot_attr3")) {

  ################## arguments check #######################

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
  }

  numeric_attrs <- get_attr_names(tree = tree, type = "n")
  boolean_attrs <- get_attr_names(tree = tree, type = "b")

  if (!(attr1 %in% numeric_attrs & attr2 %in% numeric_attrs & attr3 %in% numeric_attrs & length(unique(c(attr1, attr2, attr3))) == 3) ||
      attr1 %in% c("colony", "generation", "frame") ||
      attr2 %in% c("colony", "generation", "frame") ||
      attr3 %in% c("colony", "generation", "frame")) {
    stop(paste("Wrong cell attributes \"", attr1, "\", \"", attr2, "\", \"", attr3, "\"\n", sep = ""))
  }

  if (!(attrC %in% c(numeric_attrs, boolean_attrs, ""))) {
    stop(paste("Wrong attrC \"", attrC, "\"\n", sep = ""))
  }

  if (attrC %in% c("colony", "generation") & is.null(NC)) {
    stop("Specify number of colors NC\n")
  }

  if (!(attrS %in% c(boolean_attrs, ""))) {
    stop(paste("Wrong attrS \"", attrS, "\"\n", sep = ""))
  }

  #################################################

  addRegressionPlane <- function(myPlot) {

    myModel <- lm(df$var3 ~ df$var1 + df$var2)
    a1 <- unname(coef(myModel)[2])
    a2 <- unname(coef(myModel)[3])
    if (!is.na(a1) & !is.na(a2)) { # y = a1*x1 +a2*x2 + b
      b <- unname(coef(myModel)[1])
      r2 <- summary(myModel)$r.squared
      myPlot$plane3d(myModel, draw_line = TRUE, draw_polygon = TRUE)
      regression <- list(a1 = a1, a2 = a2, b = b, r2 = r2)
    } else { # y = a2*x2 + b or y = a1*x1 + b or y = b
      regression <- NULL # less than 3 unique data points
    }

    return(regression = regression)

  }


  plot_scatter3d <- function(mar) {

    oldPar <- par()
    par(las = 1)

    myXlab <- createAxisLab(attr = attr1, unit = unit1)
    myYlab <- createAxisLab(attr = attr2, unit = unit2)
    myZlab <- createAxisLab(attr = attr3, unit = unit3)

    myPlot <- scatterplot3d::scatterplot3d(x = df$var1, y = df$var2, z = df$var3,
                                    main = "", xlab = myXlab, ylab = "", zlab = myZlab,
                                    cex.axis = 0.8, cex.lab = 0.8, font = 2, font.lab = 2,
                                    pch = df$shape, color = df$color,
                                    grid = TRUE, col.grid = "gray", box = FALSE,
                                    y.margin.add = 0,
                                    angle = 50,
                                    mar = mar)

    # Add ylab
    plot_ylabel <- myPlot$xyz.convert
    body(plot_ylabel) <- substitute({
      mem.par <- par(mar = mar, usr = usr)
      mtext(x, side = 4, line = -3, at = z.min, adj  = 0, las = 3, cex = 0.8)
    })
    plot_ylabel(myYlab)

    # Add box lines
    plot_boxlines <- myPlot$box3d
    body(plot_boxlines) <- substitute({
      mem.par <- par(mar = mar, usr = usr)
      lines(c(x.min + yx.f * y.max, x.min + yx.f * y.max), c(z.min + yz.f * y.max, z.max + yz.f * y.max), ...)
      lines(c(x.max + yx.f * y.max, x.max + yx.f * y.max), c(z.min + yz.f * y.max, z.max + yz.f * y.max), ...)
      lines(c(x.min, x.min + y.max * yx.f), c(z.max, z.max + y.max * yz.f), ...)
      lines(c(x.min + yx.f * y.max, x.max + yx.f * y.max), c(z.max + yz.f * y.max, z.max + yz.f * y.max), ...)
    })
    plot_boxlines()

    regression <- addRegressionPlane(myPlot = myPlot)

    par(mar = oldPar$mar)

    return(regression)

  }

  #################################################

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, ".png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  cells <- get_cells(tree = tree, treeT = treeT, type = "inc")

  var1 <- vertex_attr(graph = tree, name = attr1, index = cells)
  var2 <- vertex_attr(graph = tree, name = attr2, index = cells)
  var3 <- vertex_attr(graph = tree, name = attr3, index = cells)

  df <- data.frame(cell = cells, var1 = var1, var2 = var2, var3 = var3, stringsAsFactors = FALSE)
  df <- na.omit(df)

  if (nrow(df) == 0) {
    if (save) {
      dev.off()
    }
    return(NULL)
  }

  if (attrC == "") {
    df$color <- "slateblue4"
  } else {

    colored_cells <- V(tree)[V(tree)$name %in% df$cell &
                               !is.na(vertex_attr(graph = tree, name = attrC, index = V(tree)))]$name

    if (is.null(cellColors <- getCellColors(tree = tree, cells = colored_cells, attr = attrC, N_colors = NC))) {
      df$color <- "slateblue4"
      attrC <- ""
    } else {
      df2 <- data.frame(cell = colored_cells,
                        color = cellColors,
                        stringsAsFactors = FALSE)
      df <- merge(x = df, y = df2, by = "cell", all = TRUE)
      df[which(is.na(df$color), arr.ind = TRUE), ]$color <- "gray30"
    }

  }

  if (attrS == "") {
    df$shape <- 16
  } else {
    shapes <- c(15, 16)
    df$shape <- shapes[as.numeric(vertex_attr(graph = tree, name = attrS, index = df$cell)) + 1] # FALSE -> 1, TRUE -> 2
  }

  df <- df[sample(nrow(df)), ] # suffle rows

  if (attrC != "") {

    if (attrC %in% c("colony", "generation", boolean_attrs)) {

      regression <- plot_scatter3d(mar = c(3, 3.5, 1, 6))

      if (attrS != "") {
        legend("topleft", cex = 0.8, title = as.expression(bquote(bold(.(attrS)))),
               legend = c("FALSE", "TRUE"), pch = c(15, 16), xpd = TRUE, inset = 0.01)
      }

      whichColors <- sort(unique(vertex_attr(graph = tree, name = attrC, index = df$cell)))
      plotColorLegend(attr = attrC, whichColors = whichColors, N_colors = NC, size_lab = sizeL)

    } else {

      def.par <- par(no.readonly = TRUE)
      layout(matrix(c(1, 2,
                      1, 2), nrow = 2, byrow = TRUE), widths = c(10, 1))

      regression <- plot_scatter3d(mar = c(3, 3.5, 1, 0))

      if (attrS != "") {
        legend("topleft", cex = 0.8, title = as.expression(bquote(bold(.(attrS)))),
               legend = c("FALSE", "TRUE"), pch = c(15, 16), xpd = TRUE, inset = 0.01)
      }

      plotColormap(values = vertex_attr(graph = tree, name = attrC, index = df$cell),
                   attr = attrC, unit = unitC,
                   size_lab = 1, size_values = 0.7)

      par(def.par)

    }

  } else {

    regression <- plot_scatter3d(mar = c(3, 3.5, 1, 2))

    if (attrS != "") {
      legend("topleft", cex = 0.8, title = as.expression(bquote(bold(.(attrS)))),
             legend = c("FALSE", "TRUE"), pch = c(15, 16), xpd = TRUE, inset = 0.01)
    }

  }

  if (save) {
    dev.off()
  }

  return(list(Ncells = nrow(df), regression = regression))

}
