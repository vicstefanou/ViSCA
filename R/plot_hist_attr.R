#' Create histogram of an attribute
#'
#' Creates the histogram of a numeric attribute of a lineage or division tree.
#' The histogram can be created for each colony or generation or for the whole population.
#'
#' A separate plot for each group specified in \code{groups} is generated.
#' Each histogram is created considering all cells that are included in the analysis,
#' as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in \code{attr}.
#' \cr\cr
#' y-axis represents the cell counts.
#' The range of x-axis values depicted in each plot is common and
#' is calculated as the range of \code{attr} values of all specified groups.
#' Binning is applied on the same common range.
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
#' @param groups The IDs of the groups for which to create the histogram, a vector of positive integer values.
#' This argument is ignored in case \code{grouped = "pop"}.
#' The default value \code{-1} stands for all existing groups in the \code{tree}.
#'
#' @param Ngroups Number of colonies in the movie (if \code{grouped = "col"}) or
#' number of generations in the movie (if \code{grouped = "gen"}), a non-zero positive integer value.
#' This argument is ignored in case \code{grouped = "pop"}.
#'
#' @param Nbins Number of equally spaced bins to be used for \code{attr}, a non-zero positive integer value \code{>=2}.
#' This value is indicative and may change automatically depending on the values of \code{attr},
#' producing a warning message.
#' The default value is \code{20}.
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
#' The default value is \code{2000}.}
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
#' The default value is \code{"my_hist_attr"}.}
#' }
#'
#' @return A dataframe with the following columns:
#' \enumerate{
#' \item \code{group} is the ID of the group (a positive integer value)
#' or \code{-2} in case \code{grouped = "pop"}.
#' \item \code{Ncells} is the number of cells, a positive integer value.
#' }
#' For groups with \code{Ncells = 0}, no plot is generated.
#' In case less than 2 unique cells exist, no plot is generated.
#' In case no cells exist, \code{NULL} is returned.
#'
#' @export
#' @import igraph
#' @import ggplot2
#' @importFrom graphics plot
#' @importFrom grDevices dev.off png

plot_hist_attr <- function(tree, treeT = c("LT", "DT"),
                           attr, unit = "",
                           grouped = c("col", "gen", "pop"), groups = -1, Ngroups,
                           Nbins = 20,
                           save = FALSE, savePars = list(w = 2000, h = 2000, res = 250, path = getwd(), name = "my_hist_attr")) {

  ################## arguments check ####################

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

  if (Nbins < 2) {
    stop("Nbins must be >=2\n")
  }

  ###################################################

  plot_hist <- function(values.range, Nbins, group, groupColor) {

    if (Nbins >= 2) {

      myPlot <- ggplot(data = df, aes_string(x = "var1")) +
        labs(title = paste("Histogram of", attrTitle),
             subtitle = ifelse(grouped == "population", "Population", paste(grouped, group)),
             x = createAxisLab(attr = attr, unit = unit),
             y = "Number of Cells") +
        theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
              plot.subtitle = element_text(size = 12, hjust = 0.5),
              axis.text.x = element_text(face = "bold", size = 10),
              axis.text.y = element_text(face = "bold", size = 10),
              axis.title = element_text(face = "bold", size = 12)) +
        xlim(values.range) +
        geom_histogram(breaks = seq(values.range[1], values.range[2], length.out = Nbins + 1),
                       fill = groupColor, color = "gray14")

      plot(myPlot)

    }

  }

  ###################################################

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, "_%03d.png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  attrTitle <- getAttrForTitle(attr = attr)

  if (grouped == "population") { # all population

    cells <- V(tree)[V(tree)$name %in% get_cells(tree = tree, treeT = treeT, type = "inc") &
                       !is.na(vertex_attr(graph = tree, name = attr, index = V(tree)))]$name

    if (length(cells) == 0) {
      if (save) {
        dev.off()
      }
      return(NULL)
    }

    values <- vertex_attr(graph = tree, name = attr, index = cells)

    df <- data.frame(var1 = values, stringsAsFactors = FALSE)

    plot_hist(values.range = range(values),
              Nbins = findNbins(values = values, Nbins = Nbins),
              group = -2,
              groupColor = getGroupColors())

    myData <- data.frame(group = -2,
                         Ncells = length(cells),
                         stringsAsFactors = FALSE)

  } else {

    if (length(groups) == 1)  {
      if (groups == -1) { # all groups
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

    values <- vertex_attr(graph = tree, name = attr, index = cells)
    values.range <- range(values)
    Nbins <- findNbins(values = values, Nbins = Nbins)
    graphColors <- getGroupColors(type = grouped, N_colors = Ngroups)

    myData <- data.frame(group = integer(),
                         Ncells = integer(),
                         stringsAsFactors = FALSE)

    for (i_group in groups) {

      if (grouped == "colony") {
        groupColor <- graphColors[i_group]
      } else { # group_by == "generation"
        groupColor <- graphColors[i_group + 1]
      }

      cells <- V(tree)[V(tree)$name %in% get_cells(tree = tree, treeT = treeT, type = "inc") &
                         !is.na(vertex_attr(graph = tree, name = attr, index = V(tree))) &
                         vertex_attr(graph = tree, name = grouped, index = V(tree)) == i_group]$name

      if (length(cells) != 0) {

        values <- vertex_attr(graph = tree, name = attr, index = cells)
        df <- data.frame(var1 = values, stringsAsFactors = FALSE)

        plot_hist(values.range = values.range,
                  Nbins = Nbins,
                  group = i_group,
                  groupColor = groupColor)
      }

      myData <- rbind(myData, data.frame(group = i_group,
                                         Ncells = length(cells),
                                         stringsAsFactors = FALSE))

    }

  }

  if (save) {
    dev.off()
  }

  return(myData)

}
