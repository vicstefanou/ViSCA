#' Create scatter plot of an attribute between family cells
#'
#' Creates the X-Y scatter plot of a numeric attribute between family cells of a division tree.
#' The scatter plot can be created for specific generation(s) or for the whole population.
#'
#' A common scatter plot for all generations specified in \code{gens} is created.
#' In the general case, each data point (x,y) represents the attribute value of
#' the older and younger cell of the family, respectively.
#' The scatter plot is created for all specified family cells that are included in the analysis,
#' as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in \code{attr}.
#' \cr\cr
#' Data points are dot points colored based on the generation of the younger cell.
#' Note that when \code{type = "s"} or \code{type = "c"}, there is no distinction between older and younger cell,
#' since both cells are of the same generation.
#' \cr\cr
#' The linear regression line is also drawn on the plot.
#' X is the predictor variable and Y is the response.
#' The parameters of the regression line are found using the \emph{linear least squares} method
#' provided by \code{\link[stats]{lm}}.
#'
#' @param DT The connected division tree, an object of class \code{"igraph"}.
#'
#' @param attr The name of the attribute in the \code{DT}, a character string.
#' It can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' except for \code{"colony"} and \code{"generation"}.
#'
#' @param unit The unit of \code{attr}, a character string.
#' It should be in the format \code{"<string>,<number>"},
#' where \code{",<number>"} represents the power and is optional
#' (e.g. \code{"m"} for meters and \code{"cm,3"} for cubic centimeters).
#' The default value is the empty character \code{""}, which implies that \code{attr} is in arbitrary units.
#'
#' @param type A character string naming the type of the family cells for which to create the scatter plot:
#' \itemize{
#' \item \code{"s"} for siblings
#' \item \code{"c"} for cousins
#' \item \code{"md"} for mother and daughter
#' \item \code{"gmgd"} for grand-mother and grand-daughter
#' }
#'
#' @param gens The IDs of the generations which will be included in the scatter plot,
#' a vector of non-zero positive integer values.
#' When \code{type = "md"} or \code{type = "gmgd"},
#' each value denotes the generation of the younger cell of the family (i.e. daughter or grand-daughter, respectively).
#' Acceptable values are in the range \code{[1, Ngens-1]} in case \code{type = "s"} or \code{type = "md"}
#' and in the range \code{[2, Ngens-1]} in case \code{type = "c"} or \code{type = "gmgd"}.
#' The default value \code{-1} stands for all existing generations in the \code{DT} (whole population).
#'
#' @param Ngens Number of generations in the movie, a non-zero positive integer value.
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
#' The default value is \code{"my_dot_attr_fam"}.}
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
#' @import data.table
#' @importFrom graphics plot
#' @importFrom grDevices dev.off png
#' @importFrom stats cor na.omit
#' @importFrom utils combn

plot_dot_attr_fam <- function(DT,
                              attr, unit = "",
                              type = c("s", "c", "md", "gmgd"),
                              gens = -1, Ngens,
                              save = FALSE, savePars = list(w = 2500, h = 2000, res = 350, path = getwd(), name = "my_dot_attr_fam")) {

  ################## arguments check #######################

  if (!isConnected(tree = DT)) {
    stop("DT is disconnected\n")
  }

  if (!(type %in% c("s", "c", "md", "gmgd"))) {
    stop("type must be \"s\" / \"c\" / \"md\" / \"gmgd\"\n")
  }

  numeric_attrs <- get_attr_names(tree = DT, type = "n")
  if (!(attr %in% numeric_attrs) || attr %in% c("colony", "generation")) {
    stop(paste("Wrong attr \"", attr, "\"\n", sep = ""))
  }

  possible_gens <- unique(V(DT)$generation)
  possible_gens <- possible_gens[possible_gens != -1]
  if (type %in% c("s", "md")) {
    possible_gens <- possible_gens[possible_gens >= 1]
  } else { # type %in% c("c","gmgd")
    possible_gens <- possible_gens[possible_gens >= 2]
  }
  if (length(gens) != 1) {
    if (length(m <- setdiff(gens, possible_gens)) != 0) {
      stop(paste("Selected generation(s) ", toString(m), " are wrong or do not exist\n", sep = ""))
    }
  } else {
    if (!(gens %in% c(possible_gens, -1))) {
      stop(paste("Selected generation ", gens, " does not exist\n", sep = ""))
    }
  }

  ####################################################################

  # to avoid error "no visible binding for global variable ..." in devtools::check()
  child <- mother <- grandmother <- NULL

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, ".png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  if (length(gens) == 1){
    if (gens == -1) {
      gens <- possible_gens
    }
  } else {
    gens <- sort(gens)
  }

  children <- V(DT)[V(DT)$name %in% get_cells(tree = DT, treeT = "DT", type = "inc") &
                      V(DT)$generation %in% gens]$name

  if (length(children) == 0) {
    if (save) {
      dev.off()
    }
    return(NULL)
  }

  if (type == "s") { # "Siblings"

    mothers <- V(DT)[unlist(ego(graph = DT, order = 1, nodes = children, mode = "in", mindist = 1))]$name

    dt <- as.data.table(cbind(children, mothers))
    colnames(dt) <- c("child", "mother")

    myDT <- dt[, toString(child), by = mother]

    cells <- NULL
    for (i in 1:dim(myDT)[1]) {
      b <- unlist(strsplit(myDT$V1[i], ", "))
      if (length(b) > 1) {
        cells <- cbind(cells, combn(b, 2))
      }
    }

    cells1 <- cells[1, ]
    cells2 <- cells[2, ]

    xtitle <- "first sibling"
    ytitle <- "second sibling"
    myText <- "siblings"

  } else if (type == "c") { # "Cousins"

    mothers <- V(DT)[unlist(ego(graph = DT, order = 1, nodes = children, mode = "in", mindist = 1))]$name
    grandmothers <- V(DT)[unlist(ego(graph = DT, order = 2, nodes = children, mode = "in", mindist = 2))]$name

    dt <- as.data.table(cbind(children, mothers, grandmothers))
    colnames(dt) <- c("child", "mother", "grandmother")

    gm_c <- dt[, toString(child), by = grandmother]
    colnames(gm_c) <- c("grandmother", "children")
    gm_m <- dt[, toString(mother), by = grandmother]
    colnames(gm_m) <- c("grandmother", "mothers")

    myDT <- merge(gm_c, gm_m)

    cells <- NULL

    for (i in 1:dim(myDT)[1]) { # for every grandmother

      grandchildren <- NULL
      d <- NULL
      unique_mothers <- unique(unlist(strsplit(myDT$mothers[i], ", ")))
      N_mothers <- length(unique_mothers)

      if (N_mothers > 1) { # grandchildren from more than one mother => there are cousins

        for (j in 1:N_mothers) { # find grandchildren from different mothers
          grandchildren[j] <- toString(unlist(strsplit(myDT$children[i], ", "))[which(unlist(strsplit(myDT$mothers[i], ", ")) == unique_mothers[j])])
        }

        comb_mothers <- combn(c(1:N_mothers), 2)

        for (j in 1:dim(comb_mothers)[2]) {   #for all combinations of "mothers"
          d <- c(d, as.vector(outer(unlist(strsplit(grandchildren[comb_mothers[1, j]], ", ")), unlist(strsplit(grandchildren[comb_mothers[2, j]], ", ")), paste, sep = ", ")))
        }

        for (j in 1:length(d)) {
          cells <- rbind(cells, unlist(strsplit(d[j], ", ")))
        }
      }
    }

    cells1 <- cells[, 1]
    cells2 <- cells[, 2]

    xtitle <- "first cousin"
    ytitle <- "second cousin"
    myText <- "cousins"

  } else if (type == "md") { # "MotherDaughter"

    cells2 <- children
    cells1 <- V(DT)[unlist(ego(graph = DT, order = 1, nodes = children, mode = "in", mindist = 1))]$name #mothers

    xtitle <- "mother"
    ytitle <- "daughter"
    myText <- "mother and daughter"

  } else { # type == "gmgd" --> "GrandmotherGranddaughter"

    cells2 <- children
    cells1 <- V(DT)[unlist(ego(graph = DT, order = 2, nodes = children, mode = "in", mindist = 2))]$name #grandmothers

    xtitle <- "grand-mother"
    ytitle <- "grand-daughter"
    myText <- "grand-mother and grand-daughter"

  }

  ##################

  var1 <- vertex_attr(graph = DT, name = attr, index = cells1)
  var2 <- vertex_attr(graph = DT, name = attr, index = cells2)

  df <- data.frame(var1 = var1, var2 = var2, stringsAsFactors = FALSE)
  df <- na.omit(df)

  if (nrow(df) == 0) {
    if (save) {
      dev.off()
    }
    return(NULL)
  }

  df$color <- as.factor(V(DT)[cells2]$generation) # color based on generation of daughter/grand-daughter -> always exists

  df <- df[sample(nrow(df)), ] # suffle rows

  if (length(gens) == 1) { # one generation
    myText2 <- paste("generation", gens)
  } else { # multiple generations
    diff_g <- diff(gens)
    if (length(gens) > 2 & length(unique(diff_g)) == 1) { # consequtive generations by a value
      if (unique(diff_g) == 1) { # consequtive generations by 1
        myText2 <- paste("generations ", min(gens), "-", max(gens), sep = "")
      } else {
        myText2 <- paste("generations", toString(gens))
      }
    } else {
      myText2 <- paste("generations", toString(gens))
    }
  }

  myPlot <- ggplot(data = df, aes_string(x = "var1", y = "var2")) +
    labs(title = createAxisLab(attr = attr, unit = unit),
         subtitle = myText2,
         x = xtitle,
         y = ytitle) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10),
          axis.title = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 10))

  myPlot <- getCellColorsDot(df = df, myPlot = myPlot, tree = DT, attr = "generation", N_colors = Ngens)$myPlot

  myPlot <- myPlot +
    geom_point(aes_string(color = "color"))

  if (length(levels(df$color)) == 1) {
    myPlot <- myPlot + theme(legend.position="none")
  }

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
