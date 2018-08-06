#' Plot fitted growth curves of an attribute
#'
#' Plots the fitted single-cell growth curves (fitted time-series data) of a numeric attribute
#' for all cells in a division tree.
#' The growth curves can be plotted per colony or generation or for the whole population.
#' The average growth curve of each group is also computed and plotted separately.
#'
#' For each group specified in \code{groups}, two types of plots are generated.
#' The first plot depicts the single-cell growth curves of all cells that are included in the analysis,
#' as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} values in the attributes/parameters of the \code{model},
#' as estimated by \code{\link{add_attr_growth_fit_pars}}.
#' By default, each single-cell growth curve is randomly colored.
#' \cr\cr
#' The second plot depicts the average growth curve of the group, +/- one standard deviation.
#' The following curves are drawn:
#' \itemize{
#' \item When \code{model = "lin"}:
#' \itemize{
#' \item \code{y = a_mean * t + b_mean}
#' \item \code{y = (a_mean + a_sd) * t + (b_mean + b_sd)}
#' \item \code{y = (a_mean - a_sd) * t + (b_mean - b_sd)}
#' }
#' \item When \code{model = "exp"}:
#' \itemize{
#' \item \code{y = y0_mean * e^(k_mean * t)}
#' \item \code{y = (y0_mean + y0_sd) * e^((k_mean + k_sd) * t)}
#' \item \code{y = (y0_mean - y0_sd) * e^((k_mean - k_sd) * t)}
#' }
#' }
#' The parameters of these curves are computed based on the corresponding parameters of the single-cell growth curves.
#' See the \emph{Value} field for more details.
#' Color of the area between the curves denotes the corresponding group.
#' \cr\cr
#' In both types of plots, x-axis represents the time in \emph{hours} in the range \code{[0, dur]}.
#' Value \code{0} is considered to be the birth time of each cell.
#' The range of y-axis values is common among the plots of the same type and
#' is calculated as the range of the corresponding values of all specified groups (excluding the outliers).
#'
#' @param DT The division tree, an object of class \code{"igraph"}.
#'
#' @param LT The corresponding lineage tree of the \code{DT}, an object of class \code{"igraph"}.
#'
#' @param attr The name of the attribute in the \code{LT}, a character string.
#' It can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' except for \code{"colony"}, \code{"generation"}, \code{"frame"} and \code{"age"}.
#'
#' @param unit The unit of \code{attr}, a character string.
#' It should be in the format \code{"<string>,<number>"},
#' where \code{",<number>"} represents the power and is optional
#' (e.g. \code{"m"} for meters and \code{"cm,3"} for cubic centimeters).
#' The default value is the empty character \code{""}, which implies that \code{attr} is in arbitrary units..
#'
#' @param model A character string naming the type of the fitted growth models to be plotted,
#' \itemize{
#' \item \code{"lin"} for plotting the fitted \emph{linear} models
#' \item \code{"exp"} for plotting the fitted \emph{exponential} models
#' }
#' Parameters of the model must have already been estimated using \code{\link{add_attr_growth_fit_pars}}.
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
#' The default value \code{-1} stands for all existing groups in the \code{DT}.
#'
#' @param Ngroups Number of colonies in the movie (if \code{grouped = "col"}) or
#' number of generations in the movie (if \code{grouped = "gen"}), a non-zero positive integer value.
#' This argument is ignored in case \code{grouped = "pop"}.
#'
#' @param attrC The name of the attribute in the \code{DT} by which the cells' growth curves will be colored,
#' a character string.
#' It can be any numeric or boolean attribute, as returned from \code{\link{get_attr_names}}.
#' Coloring is applied to all depicted cells of all groups specified in \code{groups},
#' except for cells with \code{NA} value in this attribute which are colored gray.
#' When the default value \code{""} (the empty character) is used, cells' colors are the default.
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
#' @param dur The time span in \code{hours}, a non-zero positive numeric value.
#' The default value is \code{1}.
#'
#' @param sizeL Size of explanatory legends, a non-zero positive numeric value.
#' The default value is \code{0.7}.
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
#' The default value is \code{1500}.}
#' \item{\code{h}}{The height of the image file in \emph{pixels}, a non-zero positive integer value.
#' The default value is \code{1000}.}
#' \item{\code{res}}{The resolution of the image file in \emph{pixels} per \emph{inch} (ppi), a non-zero positive integer value.
#' The smaller this value, the larger the plot area in inches, and the smaller the text relative to the graph itself.
#' The default value is \code{150}.}
#' \item{\code{path}}{A character string naming the directory where the image file will be saved.
#' If it does not contain an absolute path, the image file will be saved relative to the current working directory \code{getwd()}.
#' The default value is the current working directory \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{/} (not \code{\\}) on Windows.}
#' \item{\code{name}}{The image file name, a character string.
#' The suffix \code{".png"} is added automatically.
#' The default value is \code{"my_growth_attr_fit"}.}
#' }

#' @return A dataframe with the following columns:
#' \itemize{
#' \item When \code{model = "lin"}:
#' \enumerate{
#' \item \code{group} is the ID of the group (a positive integer value)
#' or \code{-2} in case \code{grouped = "pop"}
#' \item \code{Ncells} is the number of cells, a positive integer value
#' \item \code{a_mean} is the \emph{mean} of \code{"<attr>_a"}
#' (a non-zero positive numeric value in units of \code{attr} per \emph{hour}),
#' or \code{NA} in case \code{Ncells = 0}
#' \item \code{a_sd} is the \emph{standard deviation} of \code{"<attr>_a"}
#' (a non-zero positive numeric value in units of \code{attr} per \emph{hour}),
#' or \code{NA} in case \code{Ncells = 0} or \code{Ncells = 1}
#' \item \code{b_mean} is the \emph{mean} of \code{"<attr>_b"}
#' (a positive numeric value in units of \code{attr})
#' or \code{NA} in case \code{Ncells = 0}
#' \item \code{b_sd} is the \emph{standard deviation} of \code{"<attr>_b"}
#' (a positive numeric value in units of \code{attr})
#' or \code{NA} in case \code{Ncells = 0} or \code{Ncells = 1}
#' }
#' \item When \code{model = "exp"}:
#' \enumerate{
#' \item \code{group} is the ID of the group (a positive integer value)
#' or \code{-2} in case \code{grouped = "pop"}
#' \item \code{Ncells} is the number of cells, a positive integer value
#' \item \code{k_mean} is the \emph{mean} of \code{"<attr>_k"}
#' (a non-zero positive numeric value in units of \code{attr} per \emph{hour}),
#' or \code{NA} in case \code{Ncells = 0}
#' \item \code{k_sd} is the \emph{standard deviation} of \code{"<attr>_k"}
#' (a non-zero positive numeric value in units of \code{attr} per \emph{hour}),
#' or \code{NA} in case \code{Ncells = 0} or \code{Ncells = 1}
#' \item \code{y0_mean} is the \emph{mean} of \code{"<attr>_0"}
#' (a non-zero positive numeric value in units of \code{attr})
#' or \code{NA} in case \code{Ncells = 0}
#' \item \code{y0_sd} is the \emph{standard deviation} of \code{"<attr>_0"}
#' (a non-zero positive numeric value in units of \code{attr})
#' or \code{NA} in case \code{Ncells = 0} or \code{Ncells = 1}
#' }
#' }
#' \cr\cr
#' For groups with \code{Ncells = 0}, no plot of first or second type is generated.
#' For groups with \code{Ncells = 1}, no plot of second type is generated.
#' In case no cells exist, no plot is generated and \code{NULL} is returned.
#'
#' @export
#' @import igraph
#' @importFrom graphics layout lines mtext par plot polygon
#' @importFrom grDevices colors dev.off png
#' @importFrom stats as.formula sd


plot_growth_attr_fit <- function(DT, LT,
                                 attr, unit = "",
                                 model = c("lin", "exp"),
                                 grouped = c("col", "gen", "pop"), groups = -1, Ngroups,
                                 attrC = "", unitC = "", NC = NULL,
                                 dur = 1,
                                 sizeL = 0.7,
                                 save = FALSE, savePars = list(w = 1500, h = 1000, res = 150, path = getwd(), name = "my_growth_attr_fit")) {

  ################## arguments check #######################

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

  numeric_attrs <- get_attr_names(tree = LT, type = "n")
  if (!(attr %in% numeric_attrs) || attr %in% c("colony", "generation", "frame", "age")) {
    stop(paste("Wrong attr \"", attr, "\"\n", sep = ""))
  }

  if (model == "lin") {
    model <- "linear"
    parameter <- paste(attr, "a", sep = "_")
  } else if (model == "exp") {
    model <- "exponential"
    parameter <- paste(attr, "k", sep="_")
  } else {
    stop("model must be \"lin\" / \"exp\"\n")
  }

  if (!(parameter %in% vertex_attr_names(DT))) {
    stop("Parameter(s) of the model have not been estimated\nCall function add_attr_growth_fit_pars()\n")
  }

  if (grouped != "population") {
    possible_groups <- unique(vertex_attr(graph = DT, name = grouped, index = V(DT)))
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

  numeric_attrs <- get_attr_names(tree = DT, type = "n")
  boolean_attrs <- get_attr_names(tree = DT, type = "b")

  if (!(attrC %in% c(numeric_attrs, boolean_attrs, ""))) {
    stop(paste("Wrong attrC \"", attrC, "\"\n", sep = ""))
  }

  if (attrC %in% c("colony", "generation") & is.null(NC)) {
    stop("Specify number of colors NC\n")
  }

  #####################################

  find_y_upper_lim <- function(groups) { # from all defined groups

    time <- dur

    ########################### cells plot ############################

    v1 <- vertex_attr(graph = DT, name = paste(attr, v1_name, sep = "_"), index = cells)
    v2 <- vertex_attr(graph = DT, name = paste(attr, v2_name, sep = "_"), index = cells)
    max_values_cells <- eval(myFormula[[3]])

    if (model == "linear") {
      y_upper_lim_cells <- ceiling(max(max_values_cells))
    } else {
      if (length(unique(v1)) >= 2) { # sd not 0
        z_scores_cells <- abs(v1 - mean(v1)) / sd(v1)
        y_upper_lim_cells <- ceiling(max(max_values_cells[z_scores_cells < 1]))
      } else {
        y_upper_lim_cells <- ceiling(max(max_values_cells))
      }
    }

    ########################### group plot ############################

    if (length(groups) == 1) { # one group or population

      groups_mean_v1 <- mean(vertex_attr(graph = DT, name = paste(attr, v1_name, sep = "_"), index = cells))
      groups_mean_v2 <- mean(vertex_attr(graph = DT, name = paste(attr, v2_name, sep = "_"), index = cells))

      if (length(cells) >= 1) { # sd not NA
        groups_sd_v1 <- sd(vertex_attr(graph = DT, name = paste(attr, v1_name, sep = "_"), index = cells))
        groups_sd_v2 <- sd(vertex_attr(graph = DT, name = paste(attr, v2_name, sep = "_"), index = cells))
      } else {
        groups_sd_v1 <- 0
        groups_sd_v2 <- 0
      }

    } else {

      groups_mean_v1 <- groups_mean_v2 <- groups_sd_v1 <- groups_sd_v2 <- NULL

      for (i_group in groups) {

        group_cells <- V(DT)[V(DT)$name %in% cells & vertex_attr(graph = DT, name = grouped, index = V(DT)) == i_group]$name

        if (length(group_cells) != 0) {

          groups_mean_v1 <- c(groups_mean_v1, mean(vertex_attr(graph = DT, name = paste(attr, v1_name, sep = "_"), index = group_cells)))
          groups_mean_v2 <- c(groups_mean_v2, mean(vertex_attr(graph = DT, name = paste(attr, v2_name, sep = "_"), index = group_cells)))

          if (length(group_cells) >= 1) { # sd not NA
            groups_sd_v1 <- c(groups_sd_v1, sd(vertex_attr(graph = DT, name = paste(attr, v1_name, sep = "_"), index = group_cells)))
            groups_sd_v2 <- c(groups_sd_v2, sd(vertex_attr(graph = DT, name = paste(attr, v2_name, sep = "_"), index = group_cells)))
          } else {
            groups_sd_v1 <- c(groups_sd_v1, 0)
            groups_sd_v2 <- c(groups_sd_v2, 0)
          }

        }

      }
    }

    v1 <- groups_mean_v1 + groups_sd_v1
    v2 <- groups_mean_v2 + groups_sd_v2
    max_values_groups <- eval(myFormula[[3]])

    if (model == "linear") {
      y_upper_lim_groups <- ceiling(max(max_values_groups))
    } else {
      if (length(unique(v1)) >= 2) { # sd not 0
        z_scores_groups <- abs(v1 - mean(v1)) / sd(v1)
        y_upper_lim_groups <- ceiling(max(max_values_groups[z_scores_groups < 1]))
      } else {
        y_upper_lim_groups <- ceiling(max(max_values_groups))
      }
    }

    ####################################################################

    return(list(y_upper_lim_groups = y_upper_lim_groups, y_upper_lim_cells = y_upper_lim_cells))

  }



  plot_group_data <- function(cells, group, groupColor) {

    if (length(cells) != 0) {

      oldPar <- par()
      par(mar = c(4.5, 4.5, 4, 1))

      mean_v1 <- mean(vertex_attr(graph = DT, name = paste(attr, v1_name, sep = "_"), index = cells))
      mean_v2 <- mean(vertex_attr(graph = DT, name = paste(attr, v2_name, sep = "_"), index = cells))

      v1 <- mean_v1
      v2 <- mean_v2
      y_middle <- eval(myFormula[[3]])

      plot(x = time, y = y_middle,
           ylim = y_lim_groups,
           main = "Average cell Growth Curve", xlab = bquote(bold("Time (h)")), ylab = createAxisLab(attr = attr, unit = unit),
           cex.main = 1.5, cex.lab = 1, font = 2, font.lab = 2, las = 1,
           type = "n")

      if (grouped == "population") {
        mtext(paste("Population (fitted ", model, " model)", sep = ""), side = 3, cex = 1)
      } else {
        mtext(paste(grouped, " ", group, " (fitted ", model, " model)", sep = ""), side = 3, cex = 1)
      }

      if (length(cells) != 1) { # sd not NA

        sd_v1 <- sd(vertex_attr(graph = DT, name = paste(attr, v1_name, sep = "_"), index = cells))
        sd_v2 <- sd(vertex_attr(graph = DT, name = paste(attr, v2_name, sep = "_"), index = cells))

        myData <- rbind(myData, data.frame(group = group,
                                           Ncells = length(cells),
                                           v1_mean = mean_v1,
                                           v1_sd = sd_v1,
                                           v2_mean = mean_v2,
                                           v2_sd = sd_v2,
                                           stringsAsFactors = FALSE))

        if (sd_v1 != 0 || sd_v2 != 0) {

          v1 <- mean_v1 + sd_v1
          v2 <- mean_v2 + sd_v2
          y_up <- eval(myFormula[[3]])

          lines(time, y_up,
                ylim = y_lim_groups,
                type = "l", lwd = 3, lty = 3, col = groupColor)

          v1 <- mean_v1 - sd_v1
          v2 <- mean_v2 - sd_v2
          y_down <- eval(myFormula[[3]])

          lines(x = time, y = y_down,
                ylim = y_lim_groups,
                type = "l", lwd = 3, lty = 3, col = groupColor)

          polygon(x = c(time, rev(time)), y = c(y_up, rev(y_down)), col = "lightgray", border = NA)

        }

      } else {

        myData <- rbind(myData, data.frame(group = group,
                                           Ncells = 1,
                                           v1_mean = mean_v1,
                                           v1_sd = NA,
                                           v2_mean = mean_v2,
                                           v2_sd = NA,
                                           stringsAsFactors = FALSE))
      }

      lines(x = time, y = y_middle,
            ylim = y_lim_groups,
            type = "l", lwd = 3, col = groupColor)

      par(mar = oldPar$mar)

    } else {

      myData <- rbind(myData, data.frame(group = group,
                                         Ncells = 0,
                                         v1_mean = NA,
                                         v1_sd = NA,
                                         v2_mean = NA,
                                         v2_sd = NA,
                                         stringsAsFactors = FALSE))
    }

    return(myData)

  }


  plot_cells_data <- function(cells, pop_cells, group) {

    if (length(cells) != 0) {

      if (attrC != "") {
        if (attrC %in% c("colony", "generation", boolean_attrs)) {
          oldPar <- par()
          par(mar = c(4.5, 4.5, 4, 6))
        } else {
          def.par <- par(no.readonly = TRUE)
          layout(matrix(c(1, 2,
                          1, 2), nrow = 2, byrow = TRUE), widths = c(10, 1))
          oldPar <- par()
          par(mar=c(4.5, 4.5, 4, 1))
        }
      } else {
        # Cell colors random for each group (removed shades of gray -> max 433 colors)
        V(DT)[cells]$color <- sample(colors()[grep('gr(a|e)y', colors(), invert = T)], length(cells), replace = TRUE)
        oldPar <- par()
        par(mar = c(4.5, 4.5, 4, 1))
      }

      #################################################################

      cells_plotted <- NULL

      cells <- sample(cells) # suffle cells

      for (cell in cells) {

        v1 <- vertex_attr(graph = DT, name = paste(attr, v1_name, sep = "_"), index = cell)
        v2 <- vertex_attr(graph = DT, name = paste(attr, v2_name, sep = "_"), index = cell)

        if (is.null(cells_plotted)) { # first cell to plot

          plot(x = time, y = eval(myFormula[[3]]),
               ylim = y_lim_cells,
               main = "Single-cell Growth Curves", xlab = bquote(bold("Time (h)")), ylab = createAxisLab(attr = attr, unit = unit),
               cex.main = 1.5, cex.lab = 1, font = 2, font.lab = 2, las = 1,
               type = "l", lwd = 2, col = V(DT)[cell]$color)

          if (grouped == "population") {
            mtext(paste("Population (fitted ", model, " models)", sep = ""), side = 3, cex = 1)
          } else {
            mtext(paste(grouped, " ", group, " (fitted ", model, " models)", sep = ""), side = 3, cex = 1)
          }

        } else {

          lines(x = time, y = eval(myFormula[[3]]),
                ylim = y_lim_cells,
                type = "l", lwd = 2, col = V(DT)[cell]$color)

        }

        cells_plotted <- c(cells_plotted, cell)
      }

      #################################################################

      if (attrC != "") {
        if (attrC %in% c("colony", "generation", boolean_attrs)) {
          whichColors <- sort(unique(vertex_attr(graph = DT, name = attrC, index = pop_cells)))
          plotColorLegend(attr = attrC, whichColors = whichColors, N_colors = NC, size_lab = sizeL)
          par(mar = oldPar$mar)
        } else {
          par(mar = oldPar$mar)
          plotColormap(values = vertex_attr(graph = DT, name = attrC, index = pop_cells),
                       attr = attrC, unit = unitC,
                       size_lab = 1, size_values = 0.7)
          par(def.par)
        }
      } else {
        par(mar = oldPar$mar)
      }

    }

  }

  #####################################

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, "_%03d.png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  if (model == "linear") {
    v1_name <- "a"
    v2_name <- "b"
    myFormula <- as.formula("values ~ v1 * time + v2") # "values ~ a * time + b"
  } else {
    v1_name <- "k"
    v2_name <- "0"
    myFormula <- as.formula("values ~ v2 * exp(v1 * time)") # "values ~ value0 * exp(k * time)"
  }

  time <- seq(0, dur, length.out = 1000)

  myData <- data.frame(group = integer(),
                       Ncells = integer(),
                       v1_mean = double(),
                       v1_sd = double(),
                       v2_mean = double(),
                       v2_sd = double(),
                       stringsAsFactors = FALSE)

  if (grouped == "population") {

    cells <- V(DT)[V(DT)$name %in% get_cells(tree = DT, treeT = "DT", type = "inc") &
                     !is.na(vertex_attr(graph = DT, name = paste(attr, v1_name, sep = "_"), index = V(DT)))]$name

    if (length(cells) == 0) {
      return(NULL)
    }

    if (attrC != "") {
      # Cell colors based on all cells in the specified groups
      V(DT)[cells]$color <- "gray"

      colored_cells <- V(DT)[V(DT)$name %in% cells &
                               !is.na(vertex_attr(graph = DT, name = attrC, index = V(DT)))]$name

      if (is.null(cellColors <- getCellColors(tree = DT, cells = colored_cells, attr = attrC, N_colors = NC))) {
        attrC <- ""
      } else {
        V(DT)[colored_cells]$color <-  cellColors
      }
    }

    y_upper_lim <- find_y_upper_lim(groups = -2)
    y_lim_groups <- c(0, y_upper_lim$y_upper_lim_groups)
    y_lim_cells <- c(0, y_upper_lim$y_upper_lim_cells)

    myData <- plot_group_data(cells = cells,
                              group = -2,
                              groupColor = getGroupColors())

    plot_cells_data(cells = cells,
                    pop_cells = cells,
                    group = -2)

  } else {

    if (length(groups) == 1) {
      if (groups == -1) {
        groups <- possible_groups
      }
    }

    groups <- sort(groups)

    cells <- V(DT)[V(DT)$isConsidered == TRUE &
                     !is.na(vertex_attr(graph = DT, name = paste(attr, v1_name, sep = "_"), index = V(DT))) &
                     vertex_attr(graph = DT, name = grouped, index = V(DT)) %in% groups]$name

    if (length(cells) == 0) {
      return(NULL)
    }

    if (attrC != "") {
      # Cell colors based on all cells in the specified groups
      V(DT)[cells]$color <- "gray"

      colored_cells <- V(DT)[V(DT)$name %in% cells &
                               !is.na(vertex_attr(graph = DT, name = attrC, index = V(DT)))]$name

      if (is.null(cellColors <- getCellColors(tree = DT, cells = colored_cells, attr = attrC, N_colors = NC))) {
        attrC <- ""
      } else {
        V(DT)[colored_cells]$color <-  cellColors
      }
    }

    y_upper_lim <- find_y_upper_lim(groups = groups)
    y_lim_groups <- c(0, y_upper_lim$y_upper_lim_groups)
    y_lim_cells <- c(0, y_upper_lim$y_upper_lim_cells)

    graphColors <- getGroupColors(type = grouped, N_colors = Ngroups)

    for (i_group in groups) {

      if (grouped == "colony") {
        groupColor <- graphColors[i_group]
      } else { # group_by == "generation"
        groupColor <- graphColors[i_group + 1]
      }

      group_cells <- V(DT)[V(DT)$name %in% cells & vertex_attr(graph = DT, name = grouped, index = V(DT)) == i_group]$name

      myData <- plot_group_data(cells = group_cells,
                                group = i_group,
                                groupColor = groupColor)

      plot_cells_data(cells = group_cells,
                      pop_cells = cells,
                      group = i_group)

    }

  }

  if (save) {
    dev.off()
  }

  if (model == "linear") {
    colnames(myData) <- c("group", "Ncells", "a_mean", "a_sd", "b_mean", "b_sd")
  } else { # model == "exponential"
    colnames(myData) <- c("group", "Ncells", "k_mean", "k_sd", "y0_mean", "y0_sd")
  }

  return(myData)

}
