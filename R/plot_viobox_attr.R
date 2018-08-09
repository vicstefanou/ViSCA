#' Create violin plot or boxplot for an attribute
#'
#' Creates the violin plot or boxplot of a numeric attribute of a lineage or division tree.
#' The plot can be created for each colony or generation or for the whole population.
#'
#' A common plot for all groups specified in \code{groups} is generated.
#' Each violin plot or boxplot is created considering all cells that are included in the analysis,
#' as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in \code{attr}.
#' \cr\cr
#' Color denotes the corresponding group.
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
#' It can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' except for \code{"colony"}, \code{"generation"} and \code{"frame"}.
#'
#' @param unit The unit of \code{attr}, a character string.
#' It should be in the format \code{"<string>,<number>"},
#' where \code{",<number>"} represents the power and is optional
#' (e.g. \code{"m"} for meters and \code{"cm,3"} for cubic centimeters).
#' The default value is the empty character \code{""}, which implies that \code{attr} is in arbitrary units.
#'
#' @param grouped A character string naming the grouping method:
#' \itemize{
#' \item \code{"col"} for grouping by colony
#' \item \code{"gen"} for grouping by generation
#' \item \code{"pop"} for no grouping (whole population)
#' }
#'
#' @param groups The IDs of the groups for which to create the violin plot or boxplot, a vector of positive integer values.
#' This argument is ignored in case \code{grouped = "pop"}.
#' The default value \code{-1} stands for all existing groups in the \code{tree}.
#'
#' @param Ngroups Number of colonies in the movie (if \code{grouped = "col"}) or
#' number of generations in the movie (if \code{grouped = "gen"}), a non-zero positive integer value.
#' This argument is ignored in case \code{grouped = "pop"}.
#'
#' @param type A character string naming the type of the plot to be created:
#' \itemize{
#' \item \code{"vio"} for violin plot. The black dot represents the \emph{median}.
#' \item \code{"box"} for boxplot. The black dots represent the outliers (i.e. datapoints out of the 1st/3rd quantile).
#' }
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
#' The default value is \code{"my_viobox_attr"}.}
#' }
#'
#' @return A dataframe with the following columns:
#' \enumerate{
#' \item \code{group} is the ID of the group (a positive integer value)
#' or \code{-2} in case \code{grouped = "pop"}.
#' \item \code{Ncells} is the number of cells, a positive integer value.
#' }
#' For groups with \code{Ncells < 3}, no plot is generated.
#' In case no cells exist, no plot is generated and \code{NULL} is returned.
#'
#' @export
#' @import igraph
#' @import ggplot2
#' @import data.table
#' @importFrom graphics plot
#' @importFrom grDevices dev.off png
#' @importFrom stats median

plot_viobox_attr <- function(tree, treeT = c("LT", "DT"),
                             attr, unit = "",
                             grouped = c("col", "gen", "pop"), groups = -1, Ngroups,
                             type = c("vio", "box"),
                             save = FALSE, savePars = list(w = 2500, h = 2000, res = 300, path = getwd(), name = "my_viobox_attr")) {

  ################## arguments check ####################

  if (!(type %in% c("vio", "box"))) {
    stop("type must be \"vio\" / \"box\"\n")
  }

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
  }

  numeric_attrs <- get_attr_names(tree = tree, type = "n")
  if (!(attr %in% numeric_attrs) || attr %in% c("colony", "generation", "frame")) {
    stop(paste("Wrong attr \"", attr, "\"\n", sep = ""))
  }

  if (!(grouped %in% c("col", "gen", "pop"))) {
    stop("grouped must be \"col\" / \"gen\" / \"pop\"\n")
  }
  if (grouped == "col") {
    grouped <- "colony"
  } else if (grouped == "gen"){
    grouped <- "generation"
  } else { # grouped == "pop"
    grouped <- "population"
  }

  if (grouped != "population") {
    possible_groups <- unique(vertex_attr(graph = tree, name = grouped, index = V(tree)))
    possible_groups <- possible_groups[possible_groups != -1]
    if (length(groups) != 1) {
      if (length(m <- setdiff(groups, possible_groups)) != 0) {
        stop(paste("Selected ", grouped, "(s)", toString(m), " do not exist\n", sep = ""))
      }
    } else {
      if (!(groups %in% c(possible_groups, -1))) {
        stop(paste("Selected ", grouped, " ", groups, " does not exist\n", sep = ""))
      }
    }
  }

  ##########################################

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, ".png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  if (grouped == "population") { # all population

    cells <- V(tree)[V(tree)$name %in% get_cells(tree = tree, treeT = treeT, type = "inc") &
                       !is.na(vertex_attr(graph = tree, name = attr, index = V(tree)))]$name

    if (length(cells) == 0) {
      if (save) {
        dev.off()
      }
      return(NULL)
    }

    myGroups <- as.factor("Population")
    values <- vertex_attr(graph = tree, name = attr, index = cells)
    df <- data.frame(var1 = values, var2 = myGroups, stringsAsFactors = FALSE)

    myData <- data.frame(group = -2,
                         Ncells = length(cells),
                         stringsAsFactors = FALSE)

    if (length(cells) >= 3) { # if Ncells < 3 -> error only in violin plot

      isPlot <- TRUE

      myColors <- getGroupColors()

      myPlot <- ggplot(data = df, aes_string(x = "var2", y = "var1", fill = "var2")) +
        labs(x = createAxisLab(attr = grouped),
             y = createAxisLab(attr = attr, unit = unit)) +
        coord_flip() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(face = "bold", color = "black", size = 10),
              axis.title = element_text(color = "black", size = 12))

    } else {
      isPlot <- FALSE
    }

  } else {

    if (length(groups) == 1) {
      if (groups == -1) {
        groups <- possible_groups
      }
    }

    groups <- sort(groups)

    cells <- V(tree)[V(tree)$name %in% get_cells(tree = tree, treeT = treeT, type = "inc") &
                       !is.na(vertex_attr(graph = tree, name = attr, index = V(tree))) &
                       vertex_attr(graph = tree, name = grouped, index = V(tree)) %in% groups]$name

    if (length(cells) == 0) {
      if (save) {
        dev.off()
      }
      return(NULL)
    }

    myGroups <- as.factor(vertex_attr(graph = tree, name = grouped, index = cells))
    values <- vertex_attr(graph = tree, name = attr, index = cells)
    df <- data.frame(var1 = values, var2 = myGroups, stringsAsFactors = FALSE)


    myData <- as.data.frame(as.data.table(df)[, .N, by = "var2"], stringsAsFactors = FALSE)
    colnames(myData) <- c("group", "Ncells")
    myData$group <- as.numeric(levels(myData$group))
    # add groups with 0 cells
    groups_noCells <- setdiff(myData$group, groups)
    for (i_group in groups_noCells) {
      myData <- rbind(myData, data.frame(group = i_group,
                                         Ncells = 0,
                                         stringsAsFactors = FALSE))
    }
    myData <- myData[order(myData$group), ]

    plotted_groups <- myData[myData$Ncells >= 3, ]$group # if Ncells < 3 -> error only in violin plot

    if (length(plotted_groups) != 0) {

      isPlot <- TRUE

      df <- df[as.numeric(levels(df$var2)) %in% plotted_groups, ]

      graphColors <- getGroupColors(type = grouped, N_colors = Ngroups)
      if (grouped == "colony") {
        myColors <- graphColors[plotted_groups]
      } else {
        myColors <- graphColors[plotted_groups + 1]
      }

      myPlot <- ggplot(data = df, aes_string(x = "var2", y = "var1", fill = "var2")) +
        labs(x = createAxisLab(attr = grouped),
             y = createAxisLab(attr = attr, unit = unit))  +
        coord_flip() +
        theme(axis.text.x = element_text(face = "bold", size = 10),
              axis.text.y = element_text(face = "bold", size = 10),
              axis.title = element_text(face = "bold", size = 12))

    } else {
      isPlot <- FALSE
    }

  }

  if (isPlot) {
    if (type == "vio") {
      myPlot <- myPlot +
        geom_violin(color = "gray14", trim = FALSE) +
        stat_summary(fun.y = median, geom = "point", shape = 16, size = 2, color = "black") # add median
    } else { # type == "box"
      myPlot <- myPlot +
        geom_boxplot(color = "gray14")
    }

    myPlot <- myPlot +
      scale_fill_manual(values = myColors, guide = FALSE)

    plot(myPlot)
  }

  if (save) {
    dev.off()
  }

  return(myData)

}
