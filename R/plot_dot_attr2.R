#' Create scatter plot of two attributes
#'
#' Creates the X-Y scatter plot of two numeric attributes of a lineage or division tree.
#'
#' Each data point (x,y) represents the value of \code{attr1} and \code{attr2} of a cell, respectively.
#' The scatter plot is created for all cells that are included in the analysis,
#' as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in \code{attr1} or \code{attr2}.
#' \cr\cr
#' By default, all data points are dot points colored based on the density of the (X,Y) variable value.
#' 2D binning is applied on the range of all \code{attr1} and \code{attr2} values.
#' Density represents the counts (number of data points) in each bin.
#' \cr\cr
#' The linear regression line is also drawn on the plot.
#' X is the predictor variable and Y is the response.
#' The parameters of the regression line are found using the \emph{linear least squares} method
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
#' @param attr1,attr2 The names of the attributes in the \code{tree}, character strings.
#' Each one can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' except for \code{"colony"}, \code{"generation"} and \code{"frame"}.
#'
#' @param unit1,unit2 The units of the corresponding \code{attr}, character strings.
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
#' When the default value \code{""} (the empty character) is used, data points' colors are the default.
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
#' When the default value \code{""} (the empty character) is used, data points' shape is the default.
#'
#' @param Nbins Number of equally spaced bins to be used for \code{attr1} and \code{attr2},
#' a vector of two non-zero positive integer values \code{>=2}.
#' The first value is used for \code{attr1} and the second for \code{attr2}.
#' These values are indicative and may change automatically depending on the values of \code{attr1} and \code{attr2},
#' producing a warning message.
#' The default value is \code{c(20, 20)}.
#' This argument is ignored in case \code{attrC != ""}.
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
#' \item{\code{path}}{A character string naming the directory where the image file will be saved.
#' If it does not contain an absolute path, the image file will be saved relative to the current working directory \code{getwd()}.
#' The default value is the current working directory \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{/} (not \code{\\}) on Windows.}
#' \item{\code{name}}{The image file name, a character string.
#' The suffix \code{".png"} is added automatically.
#' The default value is \code{"my_dot_attr2"}.}
#' }
#'
#' @return A named list with the following components:
#' \item{Ncells}{Number of cells, a non-zero positive integer value.}
#' \item{r}{The Pearson correlation coefficient (a numeric value in the range \code{[-1, 1]})
#' or \code{NA} in case less than 2 unique data points exist.}
#' \item{regression}{A named list with the following components:
#' \itemize{
#' \item \code{a} is the slope of the regression line, a numeric value
#' \item \code{b} is the y-intercept of the regression line, a numeric value
#' \item \code{r2} is the R-squared coefficient of the regression line
#' as returned from \code{\link[stats]{summary.lm}}, a numeric value in the range \code{[0, 1]}
#' }
#' In case less than 2 unique data points exist, \code{NULL} is returned, instead.}
#' In case no cells exist, no plot is generated and \code{NULL} is returned.
#'
#' @export
#' @import igraph
#' @import ggplot2
#' @importFrom graphics plot plot.new
#' @importFrom grDevices dev.off png
#' @importFrom stats cor na.omit

plot_dot_attr2 <- function(tree, treeT = c("LT", "DT"),
                           attr1, unit1 = "", attr2, unit2 = "",
                           attrC = "", unitC = "", NC = NULL,
                           attrS = "",
                           Nbins = c(20, 20),
                           save = FALSE, savePars = list(w = 2500, h = 2000, res = 350, path = getwd(), name = "my_dot_attr2")) {

  ################## arguments check #######################

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
  }

  numeric_attrs <- get_attr_names(tree = tree, type = "n")
  boolean_attrs <- get_attr_names(tree = tree, type = "b")

  if (!(attr1 %in% numeric_attrs & attr2 %in% numeric_attrs & attr1 != attr2) ||
      attr1 %in% c("colony", "generation", "frame") ||
      attr2 %in% c("colony", "generation", "frame")) {
    stop(paste("Wrong cell attributes \"", attr1, "\", \"", attr2, "\"\n", sep = ""))
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

  if ((m <- length(Nbins)) != 2) {
    stop("Nbins must be a vector of two values\n")
  }
  if (Nbins[1] < 2) {
    stop("Nbins[1] must be >=2\n")
  }
  if (Nbins[2] < 2) {
    stop("Nbins[2] must be >=2\n")
  }

  ##################################################

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, ".png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  cells <- get_cells(tree = tree, treeT = treeT, type = "inc")

  var1 <- vertex_attr(graph = tree, name = attr1, index = cells)
  var2 <- vertex_attr(graph = tree, name = attr2, index = cells)

  df <- data.frame(cell = cells, var1 = var1, var2 = var2, stringsAsFactors = FALSE)
  df <- na.omit(df)

  if (nrow(df) == 0) {
    if (save) {
      dev.off()
    }
    return(NULL)
  }

  if (attrC != "") {
    df$color <- vertex_attr(graph = tree, name = attrC, index = df$cell)
    if (attrC %in% c("colony", "generation") || attrC %in% boolean_attrs) {
      df$color <- as.factor(df$color)
    }
  } else {
    Nbins[1] <- findNbins(values = df$var1, Nbins = Nbins[1])
    Nbins[2] <- findNbins(values = df$var2, Nbins = Nbins[2])
    if (Nbins[1] >= 2 & Nbins[2] >= 2) {
      dens <- gplots::hist2d(x = df$var1, y = df$var2, nbins = Nbins, show = FALSE)
      ix <- findInterval(x = df$var1, vec = dens$x.breaks, rightmost.closed = TRUE)
      iy <- findInterval(x = df$var2, vec = dens$y.breaks, rightmost.closed = TRUE)
      df$color <- dens$counts[cbind(ix, iy)]  # counts --> integer values
    } else {
	    df$color <- 0
    }
  }

  if (attrS != "") {
    df$shape <- as.factor(vertex_attr(graph = tree, name = attrS, index = df$cell))
  }

  df <- df[sample(nrow(df)), ] # suffle rows

  myPlot <- ggplot(data = df, aes_string(x = "var1", y = "var2")) +
    labs(x = createAxisLab(attr = attr1, unit = unit1),
         y =  createAxisLab(attr = attr2, unit = unit2)) +
    theme(axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10),
          axis.title = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 10))

  resultColor <- getCellColorsDot(df = df, myPlot = myPlot, tree = tree, attr = attrC, unit = unitC, N_colors = NC)
  myPlot <- resultColor$myPlot
  colored <- resultColor$colored

  if (colored == TRUE) {
    if (attrS != "") {
      myPlot <- myPlot +
        labs(shape = createAxisLab(attr = attrS))  +
        scale_shape_manual(values = c(15, 16)) +
        geom_point(aes_string(color = "color", shape = "shape"))
    } else {
      myPlot <- myPlot + geom_point(aes_string(color = "color"))
    }
  } else {
    myPlot <- myPlot + geom_point(color = "slateblue4")
  }

  myRegression <- addRegressionLine(myPlot = myPlot, df = df)

  myPlot <- myRegression$myPlot
  regression <- myRegression$regression

  if (colored == FALSE) {
    plot(myPlot)
  } else {

    myPlot <- myPlot +
      theme(axis.text.y = element_text(size = 10, angle = 90, hjust = 0.5),
            legend.text = element_text(size = 7),
            legend.position = c(0,1), legend.justification = c(0,1),
            legend.direction = "horizontal", legend.box = "horizontal")

    blankPlot <- createBlankPlot()

    plot1 <- ggplot(data = df, aes_string(x = "var1")) +
      geom_histogram(breaks = dens$x.breaks, fill = "lightgray", color = "gray14") +
      labs(y = "Counts", x = "") +
      theme_minimal()+
      theme(axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10, angle = 90, hjust = 0.5),
            axis.title = element_text(face = "bold", size = 12))

    plot2 <- ggplot(data = df, aes_string(x = "var2")) +
      geom_histogram(breaks = dens$y.breaks, fill = "lightgray", color = "gray14") +
      coord_flip() +
      labs(y = "Counts", x = "") +
      theme_minimal()+
      theme(axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10, angle = 90, hjust = 0.5),
            axis.title = element_text(face = "bold", size = 12))

    lay <- rbind(c(1, 1, 1, 2),
                 c(3, 3, 3, 4),
                 c(3, 3, 3, 4),
                 c(3, 3, 3, 4))

    nP <- gridExtra::arrangeGrob(plot1, blankPlot, myPlot, plot2, layout_matrix = lay)

    plot.new()
    grid::grid.draw(nP)

  }

  if (nrow(df) == 1) {
    r <- NA # single cell
  } else {
    r <- tryCatch(cor(df[, c("var1", "var2")], use = "all.obs", method = "pearson")["var1", "var2"],
                  warning = function(w) {
                    # less than 2 unique data points
                    return(NA)
                  })
  }

  if (save) {
    dev.off()
  }

  return(list(Ncells = nrow(df), r = r, regression = regression))

}

