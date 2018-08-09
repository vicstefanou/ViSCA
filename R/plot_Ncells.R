#' Plot cell counts
#'
#' Calculates and plots the number of cells (cell counts) in a lineage or division tree.
#' Cell counts can be plotted for each colony, generation or frame.
#'
#' A common plot for all groups specified in \code{groups} is generated.
#' Cell counts are calculated considering all cells that are included in the analysis,
#' as returned from \code{\link{get_cells}}.
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
#' @param grouped A character string naming the grouping method:
#' \itemize{
#' \item \code{"col"} for grouping by colony
#' \item \code{"gen"} for grouping by generation
#' \item \code{"frame"} for grouping by frame, which is acceptable in case \code{treeT = "LT"}.
#' }
#'
#' @param groups The IDs of the groups for which to plot the cell counts, a vector of positive integer values.
#' The default value \code{-1} stands for all existing groups in the \code{tree}.

#' @param Ngroups Number of colonies in the movie (if \code{grouped = "col"}) or
#' number of generations in the movie (if \code{grouped = "gen"}) or
#' number of frames in the movie (if \code{grouped = "frame"}), a non-zero positive integer value.
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
#' The default value is \code{4000}.}
#' \item{\code{h}}{The height of the image file in \emph{pixels}, a non-zero positive integer value.
#' The default value is \code{2000}.}
#' \item{\code{res}}{The resolution of the image file in \emph{pixels} per \emph{inch} (ppi), a non-zero positive integer value.
#' The smaller this value, the larger the plot area in inches, and the smaller the text relative to the graph itself.
#' The default value is \code{250}.}
#' \item{\code{path}}{A character string naming the directory where the image file will be saved (excluding the last \code{"/"}).
#' If it does not contain an absolute path, the image file will be saved relative to the current working directory \code{getwd()}.
#' The default value is the current working directory \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{"/"} on Windows.}
#' \item{\code{name}}{The image file name, a character string.
#' The suffix \code{".png"} is added automatically.
#' The default value is \code{"my_Ncells"}.}
#' }
#'
#' @return A dataframe with the following columns:
#' \enumerate{
#' \item \code{group} is the ID of the group, a positive integer value.
#' \item \code{Ncells} is the number of cells, a positive integer value.
#' }
#' In case no cells exist, no plot is generated and \code{NULL} is returned.
#'
#' @export
#' @import igraph
#' @import ggplot2
#' @import data.table
#' @importFrom graphics plot
#' @importFrom grDevices dev.off png

plot_Ncells <- function(tree, treeT = c("LT", "DT"),
                        grouped = c("col", "gen", "frame"), groups = -1, Ngroups,
                        save = FALSE, savePars = list(w = 4000, h = 2000, res = 250, path = getwd(), name = "my_Ncells")) {

  ################## arguments check ####################

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
  }

  if (!(grouped %in% c("col", "gen", "frame"))) {
    stop("grouped must be \"col\" / \"gen\" / \"frame\"\n")
  }
  if (grouped == "frame" & treeT == "DT") {
    stop("grouped can be \"frame\" only when treeT is \"LT\"\nFor treeT = \"DT\" set grouped to \"col\" / \"gen\"\n")
  }
  if (grouped == "col") {
    grouped <- "colony"
  } else if (grouped == "gen"){
    grouped <- "generation"
  } #else { # grouped == "frame"
    #grouped <- "frame"
  #}

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

  ###################################################

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, ".png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  if (length(groups) == 1)  {
    if (groups == -1) { # all groups
      groups <- possible_groups
    }
  }

  groups <- sort(groups)

  cells <- V(tree)[V(tree)$name %in% get_cells(tree = tree, treeT = treeT, type = "inc") &
                     vertex_attr(graph = tree, name = grouped, index = V(tree)) %in% groups]$name

  if (length(cells) == 0) {
    if (save) {
      dev.off()
    }
    return(NULL)
  }

  myGroups <- as.factor(vertex_attr(graph = tree, name = grouped, index = cells))
  df <- data.frame(var1 = myGroups, stringsAsFactors = FALSE)

  if (grouped == "colony") {
    graphColors <- getGroupColors(type = grouped, N_colors = Ngroups)
    myColors <- graphColors[as.numeric(levels(df$var1))]
    myTitle <- paste("Distribution of", ifelse(treeT == "LT", "cell instants", "cells"), "across colonies")
  } else if (grouped == "generation") {
    graphColors <- getGroupColors(type = grouped, N_colors = Ngroups)
    myColors <- graphColors[as.numeric(levels(df$var1)) + 1]
    myTitle <- paste("Distribution of", ifelse(treeT == "LT", "cell instants", "cells"), "across generations")
  } else {
    myColors <- getGroupColors()
    myTitle <- "Distribution of cell instants across frames"
  }

  myPlot <- ggplot(data = df, aes_string(x = "var1")) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10),
          axis.title = element_text(face = "bold", size = 12)) +
    labs(title = myTitle,
         x = createAxisLab(attr = grouped),
         y = "Number of Cells") +
    geom_bar(fill = myColors, stat = "count")

  plot(myPlot)

  if (save) {
    dev.off()
  }

  myData <- as.data.frame(as.data.table(df)[, .N, by = "var1"], stringsAsFactors = FALSE)
  colnames(myData) <- c("group", "Ncells")
  myData$group <- as.numeric(myData$group)
  # add groups with 0 cells
  groups_noCells <- setdiff(myData$group, groups)
  for (i_group in groups_noCells) {
    myData <- rbind(myData, data.frame(group = i_group,
                                       Ncells = 0,
                                       stringsAsFactors = FALSE))
  }
  myData <- myData[order(myData$group), ]

  return(myData)

}
