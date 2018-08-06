#' Plot mean and sd of an attribute per frame
#'
#' Computes and plots the \emph{mean} and \emph{standard deviation} of a numeric attribute of a lineage or division tree
#' per frame. The plot can be created for each colony or generation or for the whole population.
#'
#' A separate plot for each group specified in \code{groups} is generated.
#' Each plot is created considering all cells that are included in the analysis,
#' as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in \code{attr}.
#' \cr\cr
#' x-axis represents the time in \emph{frames}.
#' The range of x and y-axis values depicted in each plot is common and
#' is calculated as the range of the corresponding values of all specified groups.
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
#' @param grouped A character string naming the grouping method:
#' \itemize{
#' \item \code{"col"} for grouping by colony
#' \item \code{"gen"} for grouping by generation
#' \item \code{"pop"} for no grouping (whole population)
#' }
#'
#' @param groups The IDs of the groups for which to create the plot, a vector of positive integer values.
#' This argument is ignored in case \code{grouped = "pop"}.
#' The default value \code{-1} stands for all existing groups in the \code{tree}.
#'
#' @param Ngroups Number of colonies in the movie (if \code{grouped = "col"}) or
#' number of generations in the movie (if \code{grouped = "gen"}), a non-zero positive integer value.
#' This argument is ignored in case \code{grouped = "pop"}.
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
#' The default value is \code{"my_time_attr"}.}
#' }
#'
#' @return A named list with the following components:
#' \item{group}{The ID of the group (a positive integer value)
#' or \code{-2} in case \code{grouped = "pop"}.}
#' \item{data}{A dataframe with the following columns:
#' \enumerate{
#' \item \code{frame} is the frame ID, a non-zero positive integer value.
#' \item \code{Ncells} is the number of cells, a positive integer value.
#' \item \code{mean} is the \emph{mean} of \code{attr}, a numeric value.
#' \item \code{sd} is the \emph{standard deviation} of \code{attr}, a positive numeric value
#' or \code{NA} in case \code{Ncells = 1}.
#' }
#' For groups with no cells, no plot is generated and \code{NULL} is returned.
#' }
#'
#' @export
#' @import igraph
#' @import ggplot2
#' @importFrom graphics plot
#' @importFrom grDevices dev.off png
#' @importFrom stats na.omit sd

plot_time_attr <- function(tree, treeT = c("LT", "DT"),
                           attr, unit = "",
                           grouped = c("col", "gen", "pop"), groups = -1, Ngroups,
                           save = FALSE, savePars = list(w = 2500, h = 2000, res = 350, path = getwd(), name = "my_time_attr")) {

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

  ###################################################

  find_df2 <- function() {

    df2 <- data.frame(frame = integer(),
                      Ncells = integer(),
                      mean = double(),
                      sd = double(),
                      min = double(),
                      max = double(),
                      stringsAsFactors = FALSE)

    for (fr in sort(unique(df$var1))) { # no run if df$var1 = integer(0) (no cells) -> returns empty df2

      fr_mean <- mean(df$var2[df$var1 == fr])
      fr_sd <- sd(df$var2[df$var1 == fr])

      df2 <- rbind(df2, data.frame(frame = fr,
                                   Ncells = length(df$var2[df$var1 == fr]),
                                   mean = fr_mean,
                                   sd = fr_sd,
                                   min = ifelse(is.na(fr_sd), fr_mean, fr_mean - fr_sd),
                                   max = ifelse(is.na(fr_sd), fr_mean, fr_mean + fr_sd),
                                   stringsAsFactors = FALSE))

    }

    return(df2)

  }


  plot_group_data <- function(data) {

    df <- data$df
    group <- data$group
    groupColor <- data$color

    myPlot <- ggplot(data = df, aes_string(x = "frame", y = "mean")) +
      labs(title = paste("mean \u00B1 sd of", attrTitle, "per frame"),
           subtitle = ifelse(grouped == "population", "Population", paste(grouped, group)),
           x = createAxisLab(attr = "Time", unit = "frames"),
           y =  createAxisLab(attr = attr, unit = unit)) +
      theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face = "bold", size = 10),
            axis.title = element_text(face = "bold", size = 12),
            legend.title = element_text(face = "bold", size = 10)) +
      xlim(time.range) +
      ylim(values.range)


    if (nrow(df) > 2) {
      myPlot <- myPlot +
        geom_ribbon(aes_string(ymin = "min", ymax = "max"), fill = groupColor) +
        geom_line(color = "gray14", size = 0.8)
    }

    myPlot <- myPlot + geom_point()

    plot(myPlot)

  }

  ###################################################

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, "_%03d.png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  attrTitle <- getAttrForTitle(attr = attr)

  if (treeT == "LT") {
    attrT <- "frame"
  } else { # treeT == "DT"
    if (grepl("_birth", attr)) {
      attrT <- "birthTime"
    } else { # grepl("_division", attr)
      attrT <- "divisionTime"
    }
  }

  plotData <- list()
  myData <- list()

  if (grouped == "population") {

    cells <- get_cells(tree = tree, treeT = treeT, type = "inc")

    var1 <- vertex_attr(graph = tree, name = attrT, index = cells)
    var2 <- vertex_attr(graph = tree, name = attr, index = cells)

    df <- data.frame(var1 = var1, var2 = var2, stringsAsFactors = FALSE)
    df <- na.omit(df)

    df2 <- find_df2()

    if (length(df2) != 0) {
      plotData[[length(plotData) + 1]] <- list(df = df2, group = -2, color = getGroupColors())
      myData[[length(myData) + 1]] <- list(group = -2, data = df2[, c("frame", "Ncells", "mean", "sd")])
    } else {
      myData[[length(myData) + 1]] <- list(group = -2, data = NULL)
    }

  } else {

    if (length(groups) == 1)  {
      if (groups == -1) { # all groups
        groups <- possible_groups
      }
    }

    groups <- sort(groups)

    graphColors <- getGroupColors(type = grouped, N_colors = Ngroups)

    for (i_group in groups) {

      if (grouped == "colony") {
        groupColor <- graphColors[i_group]
      } else { # group_by == "generation"
        groupColor <- graphColors[i_group + 1]
      }

      cells <- V(tree)[V(tree)$name %in% get_cells(tree = tree, treeT = treeT, type = "inc") &
                         vertex_attr(graph = tree, name = grouped, index = V(tree)) == i_group]$name

      var1 <- vertex_attr(graph = tree, name = attrT, index = cells)
      var2 <- vertex_attr(graph = tree, name = attr, index = cells)

      df <- data.frame(var1 = var1, var2 = var2, stringsAsFactors = FALSE)
      df <- na.omit(df)

      df2 <- find_df2()

      if (length(df2) != 0) {
        plotData[[length(plotData) + 1]] <- list(df = df2, group = i_group, color = groupColor)
        myData[[length(myData) + 1]] <- list(group = i_group, data = df2[, c("frame", "Ncells", "mean", "sd")])
      } else {
        myData[[length(myData) + 1]] <- list(group = i_group, data = NULL)
      }

    }
  }

  if (length(plotData) != 0) {

    dfs <- lapply(plotData, "[[", "df")
    time.range <- range(unlist(sapply(dfs, "[[", "frame")))
    values.range <- c(min(unlist(sapply(dfs, "[[", "min"))), max(unlist(sapply(dfs, "[[", "max"))))

    for (i in 1:length(plotData)) {
      plot_group_data(data = plotData[[i]])
    }
  }

  if (save) {
    dev.off()
  }

  return(myData)

}
