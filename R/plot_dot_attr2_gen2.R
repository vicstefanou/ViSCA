#' Create scatter plot of attribute(s) between two generations
#'
#' Creates the X-Y scatter plot of the same or different numeric attribute
#' between cells of two specific generations of a division tree.
#'
#' Each data point (x,y) represents the corresponding attribute value of
#' the older and younger cell, respectively.
#' The scatter plot is created for all specified cells that are included in the analysis,
#' as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in \code{attr1} (or \code{attr2} in case \code{attr2 != ""}).
#' \cr\cr
#' Data points are purple dot points.
#' \cr\cr
#' The linear regression line is also drawn on the plot.
#' X is the predictor variable and Y is the response.
#' The parameters of the regression line are found using the \emph{linear least squares} method
#' provided by \code{\link[stats]{lm}}.
#'
#' @param DT The connected division tree, an object of class \code{"igraph"}.
#'
#' @param gen1 The ID of the first generation in the \code{DT}, a positive integer value.
#' @param gen2 The ID of the second generation in the \code{DT}, a positive integer value.

#' @param attr1 The name of the first attribute in the \code{DT}, a character string.
#' It can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' except for \code{"colony"} and \code{"generation"}.
#' When \code{attr2 = ""}, this attribute is used for cells of both generations \code{gen1} and \code{gen2}.
#' In other case, this attribute is used for cells of generation \code{gen1}.
#'
#' @param unit1 The unit of \code{attr1}, a character string.
#' It should be in the format \code{"<string>,<number>"},
#' where \code{",<number>"} represents the power and is optional
#' (e.g. \code{"m"} for meters and \code{"cm,3"} for cubic centimeters).
#' The default value is the empty character \code{""}, which implies that \code{attr1} is in arbitrary units.
#'
#' @param attr2 The name of the second attribute in the \code{DT}, a character string.
#' When the default value \code{""} (the empty string) is used,
#' the same attribute \code{attr1} is used for cells of both generations \code{gen1} and \code{gen2}.
#' In other case, it can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' except for \code{"colony"} and \code{"generation"}, and it is used for cells of generation \code{gen2}.
#'
#' @param unit2 The unit of \code{attr2}, a character string.
#' It should be in the format \code{"<string>,<number>"},
#' where \code{",<number>"} represents the power and is optional
#' (e.g. \code{"m"} for meters and \code{"cm,3"} for cubic centimeters).
#' The default value is the empty character \code{""}, which implies that \code{attr2} is in arbitrary units.
#' This argument is ignored in case \code{attr2 = ""}.
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
#' The default value is \code{"my_dot_attr2_gen2"}.}
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
#' @seealso \code{\link{isConnected}} for checking if a tree is connected.
#' @export
#' @import igraph
#' @import ggplot2
#' @importFrom graphics abline plot
#' @importFrom grDevices dev.off png
#' @importFrom stats cor na.omit

plot_dot_attr2_gen2 <- function(DT,
                                gen1, gen2,
                                attr1, unit1 = "", attr2 = "", unit2 = "",
                                save = FALSE, savePars = list(w = 2500, h = 2000, res = 350, path = getwd(), name = "my_dot_attr2_gen2")) {

  ################## arguments check #######################

  if (!isConnected(tree = DT)) {
    stop("DT is disconnected\n")
  }

  numeric_attrs <- get_attr_names(tree = DT, type = "n")

  if (!(attr1 %in% numeric_attrs) || attr1 %in% c("colony", "generation")) {
    stop(paste("Wrong attr1 \"", attr1, "\"\n", sep = ""))
  }

  if (!(attr2 %in% c(numeric_attrs, "")) || attr2 %in% c("colony", "generation")) {
    stop(paste("Wrong attr2 \"", attr2, "\"\n", sep = ""))
  }

  possible_gens <- unique(V(DT)$generation)
  possible_gens <- possible_gens[possible_gens != -1]
  if (!(gen1 %in% possible_gens & gen2 %in% possible_gens & gen1 != gen2)) {
    stop("Wrong generations\n")
  }

  ####################################################################

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, ".png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  if (gen1 < gen2) {
    g1 <- gen1
    g2 <- gen2
  } else {
    g1 <- gen2
    g2 <- gen1
  }

  cells2 <- V(DT)[V(DT)$name %in% get_cells(tree = DT, treeT = "DT", type = "inc") &
                    V(DT)$generation == g2]$name
  cells1 <- V(DT)[unlist(ego(DT, g2 - g1, nodes = cells2, mode = "in", mindist = g2 - g1))]$name

  var1 <- vertex_attr(graph = DT, name = attr1, index = cells1)

  if (attr2 == "") {
    var2 <- vertex_attr(graph = DT, name = attr1, index = cells2)
  } else {
    var2 <- vertex_attr(graph = DT, name = attr2, index = cells2)
  }

  df <- data.frame(var1 = var1, var2 = var2, stringsAsFactors = FALSE)
  df <- na.omit(df)

  if (nrow(df) == 0) {
    if (save) {
      dev.off()
    }
    return(NULL)
  }

  df <- df[sample(nrow(df)), ] # suffle rows

  if (attr2 == "") {
    myTitle <- createAxisLab(attr = attr1, unit = unit1)
    myXlab <- paste("generation", g1)
    myYlab <- paste("generation", g2)
  } else {
    myTitle <- ""
    myXlab <- createAxisLab(attr = attr1, unit = unit1, string = paste("- generation", g1))
    myYlab <- createAxisLab(attr = attr2, unit = unit2, string = paste("- generation", g2))
  }

  myPlot <- ggplot(data = df, aes_string(x = "var1", y = "var2")) +
    labs(title = myTitle,
         x = myXlab,
         y =  myYlab) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10),
          axis.title = element_text(face = "bold", size = 12))

  myPlot <- myPlot + geom_point(color = "slateblue4")

  myRegression <- addRegressionLine(myPlot = myPlot, df = df)

  myPlot <- myRegression$myPlot
  regression <- myRegression$regression

  plot(myPlot)

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
