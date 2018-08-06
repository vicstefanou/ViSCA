#' Create dot plot of an attribute by time
#'
#' Creates the dot plot of a numeric attribute of a lineage or division tree by time.
#'
#' Each data point (x,y) represents the time in \emph{frames} and the value of \code{attr} of a cell, respectively.
#' The dot plot is created for all cells that are included in the analysis,
#' as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in \code{attr}.
#' \cr\cr
#' Data points are dot points colored based on the density of the Y variable value at every time point x.
#' Binning is applied on the range of all \code{attr} values.
#' Density represents the counts (number of data points) in each bin at every time point x.
#' Ultimately, the plot is a color representation of the histogram of \code{attr} in time.
#'
#' @param tree The lineage or division tree, an object of class \code{"igraph"}.
#'
#' @param treeT A character string naming the type of \code{tree}:
#' \itemize{
#' \item \code{"LT"} if \code{tree} is a lineage tree
#' \item \code{"DT"} if \code{tree} is a division tree
#' }
#'
#' @param attr The name of the attribute in the \code{tree}, a character string.
#' In case \code{treeT = "LT"}, it can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' except for \code{"colony"}, \code{"generation"} and \code{"frame"}.
#' In case \code{treeT = "DT"}, it can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' having the suffix \code{"_birth"} or \code{"_division"}.
#'
#' @param unit The unit of \code{attr}, a character string.
#' It should be in the format \code{"<string>,<number>"},
#' where \code{",<number>"} represents the power and is optional
#' (e.g. \code{"m"} for meters and \code{"cm,3"} for cubic centimeters).
#' The default value is the empty character \code{""}, which implies that \code{attr} is in arbitrary units.
#'
#' @param Nbins Number of equally spaced bins to be used for \code{attr}, a non-zero positive integer value \code{>=2}.
#' This value is indicative and may change automatically depending on the values of \code{attr},
#' producing a warning message.
#' The default value is \code{20}.
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
#' The default value is \code{"my_dot_time_attr"}.}
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
#' @importFrom graphics hist plot.new
#' @importFrom grDevices dev.off png
#' @importFrom stats na.omit

plot_dot_time_attr <- function(tree, treeT = c("LT", "DT"),
                               attr, unit = "",
                               Nbins = 20,
                               save = FALSE, savePars = list(w = 2500, h = 2000, res = 350, path = getwd(), name = "my_dot_time_attr")) {

  ################## arguments check #######################

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
  }

  numeric_attrs <- get_attr_names(tree = tree, type = "n")
  if (treeT == "LT") {
    if (!(attr %in% numeric_attrs) || attr %in% c("colony", "generation", "frame")) {
      stop(paste("Wrong attr \"", attr, "\"\n", sep = ""))
    }
  } else { # treeT == "DT"
    if (!(attr %in% numeric_attrs & (grepl("_birth", attr) || grepl("_division", attr)))) {
      stop(paste("Wrong attr \"", attr, "\"\n", sep = ""))
    }
  }

  if (Nbins < 2) {
    stop("Nbins must be >=2\n")
  }

  ##################################################

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, ".png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  if (treeT == "LT") {
    time_attr <- "frame"
  } else { # treeT == "DT"
    if (grepl("_birth", attr)) {
      time_attr <- "birthTime"
    } else { # grepl("_division", attr)
      time_attr <- "divisionTime"
    }
  }

  cells <- get_cells(tree = tree, treeT = treeT, type = "inc")

  var1 <- vertex_attr(graph = tree, name = time_attr, index = cells)
  var2 <- vertex_attr(graph = tree, name = attr, index = cells)

  df <- data.frame(var1 = var1, var2 = var2, stringsAsFactors = FALSE)
  df <- na.omit(df)

  if (nrow(df) == 0) {
    if (save) {
      dev.off()
    }
    return(NULL)
  }

  myPlot <- ggplot(data = df, aes_string(x = "var1", y = "var2")) +
    labs(x = createAxisLab(attr = "Time", unit = "frames"),
         y =  createAxisLab(attr = attr, unit = unit)) +
    theme(axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10),
          axis.title = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 10))

  Nbins <- findNbins(values = df$var2, Nbins = Nbins) # bin overall range of var2

  if (Nbins >= 2) {
    for (fr in sort(unique(df$var1))) { # find hist (counts) separately for each frame
      dens <- hist(x = df$var2[df$var1 == fr],
                   breaks = seq(min(df$var2), max(df$var2), length.out = Nbins + 1),
                   plot = FALSE)
      i <- findInterval(x = df$var2[df$var1 == fr], vec = dens$breaks, rightmost.closed = TRUE)
      df$color[df$var1 == fr] <- dens$counts[i]  # counts --> integer values
    }
  } else {
    df$color <- 0
  }

  resultColor <- getCellColorsDot(df = df, myPlot = myPlot, tree = tree, attr = "")
  myPlot <- resultColor$myPlot
  colored <- resultColor$colored

  if (colored == TRUE) {
    myPlot <- myPlot + geom_point(aes_string(color = "color"))
  } else {
    myPlot <- myPlot + geom_point(color = "slateblue4")
  }

  myRegression <- addRegressionLine(myPlot = myPlot, df = df)

  myPlot <- myRegression$myPlot
  regression <- myRegression$regression

  myPlot <- myPlot +
    theme(axis.text.y = element_text(size = 10, angle = 90, hjust = 0.5),
          legend.text = element_text(size = 7),
          legend.position = c(0,1), legend.justification = c(0,1),
          legend.direction = "horizontal", legend.box = "horizontal")

  plot1 <- ggplot(data = df, aes_string(x = "var1")) +
    geom_bar() +
    labs(y = "Counts", x = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10, angle = 90, hjust = 0.5),
          axis.title = element_text(face = "bold", size = 12))

  if (colored == FALSE) {

    lay <- rbind(c(1, 1, 1),
                 c(2, 2, 2),
                 c(2, 2, 2),
                 c(2, 2, 2))

    nP <- gridExtra::arrangeGrob(plot1, myPlot, layout_matrix = lay)

  } else {

    blankPlot <- createBlankPlot()

    plot2 <- ggplot(data = df, aes_string(x = "var2")) +
      geom_histogram(breaks = seq(min(df$var2), max(df$var2), length.out = Nbins + 1),
                     fill = "lightgray", color = "gray14") +
      coord_flip() +
      labs(y = "Counts", x = "") +
      theme_minimal() +
      theme(axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10, angle = 90, hjust = 0.5),
            axis.title = element_text(face = "bold", size = 12))

    lay <- rbind(c(1, 1, 1, 2),
                 c(3, 3, 3, 4),
                 c(3, 3, 3, 4),
                 c(3, 3, 3, 4))

    nP <- gridExtra::arrangeGrob(plot1, blankPlot, myPlot, plot2, layout_matrix = lay)

  }

  plot.new()
  grid::grid.draw(nP)

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
