#' Plot the PDF of an attribute
#'
#' Computes and plots the Probability Density Function (PDF) of a numeric attribute of a lineage or division tree.
#' The PDF can be plotted for each colony or generation or for the whole population.
#' Each PDF is computed by fitting a distribution model (\emph{Normal}, \emph{Gamma} or \emph{Lognormal})
#' to the corresponding data.
#'
#' Each PDF is computed considering all cells that are included in the analysis,
#' as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in \code{attr}.
#' \cr\cr
#' The range of x-axis (attribute) values depicted in each 2D plot is common and
#' is calculated as the range of values of all groups specified in \code{groups}.
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
#' @param groups The IDs of the groups for which to plot the PDF, a vector of positive integer values.
#' This argument is ignored in case \code{grouped = "pop"}.
#' The default value \code{-1} stands for all existing groups in the \code{tree}.
#'
#' @param Ngroups Number of colonies in the movie (if \code{grouped = "col"}) or
#' number of generations in the movie (if \code{grouped = "gen"}), a non-zero positive integer value.
#' This argument is ignored in case \code{grouped = "pop"}.
#'
#' @param model A character string naming the distribution model to be fitted:
#' \itemize{
#' \item \code{"norm"} is for fitting the \emph{Normal} distribution.
#' \item \code{"gamma"} is for fitting the \emph{Gamma} distribution.
#' \item \code{"lnorm"} is for fitting the \emph{Lognormal} distribution.
#' \item \code{"auto"} is for finding the best-fit distribution.
#' This is accomplished by fitting separately the \emph{Normal}, \emph{Gamma} and \emph{Lognormal} distribution.
#' The best-fit distribution is then chosen using the Bayesian Inference Criterion (BIC),
#' according to which the best model is the one with the lowest numeric BIC value.
#' }
#' Each model is fitted using the \emph{maximum likelyhood estimation (MLE)} method
#' provided by \code{\link[fitdistrplus]{fitdist}}.
#' \cr
#' Note that the \emph{Gamma} and \emph{Lognormal} distributions can be fitted to attributes with
#' non-zero positive numeric values.
#' Zero values are automatically replaced by value \code{1e-6}.
#' For negative values, an error is produced.
#'
#' @param plot3D A logical value (\code{TRUE} or \code{FALSE}) indicating whether
#' a 3D or 2D plot will be generated, respectively.
#' When the default value \code{TRUE} is used, a common 3D plot for all groups specified in \code{groups} is generated.
#' When the value \code{FALSE} is used, a separate 2D plot for each group specified in \code{groups} is generated.
#' This argument is ignored (regarded as \code{FALSE}) in case \code{grouped = "pop"}
#' or if only one group is specified in \code{groups}.
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
#' The default value is \code{"my_pdf_attr"}.}
#' }
#'
#' @return A dataframe with the following columns:
#' \enumerate{
#' \item \code{group} is the ID of the group (a positive integer value)
#' or \code{-2} in case \code{grouped = "pop"}.
#' \item \code{Ncells} is the number of cells, a positive integer value.
#' \item \code{distr} is a character string naming the distribution model that was fitted:
#' \code{"norm"} for \emph{Normal},
#' \code{"gamma"} for \emph{Gamma} and
#' \code{"lnorm"} for \emph{Lognormal} distribution or
#' \code{NA} if no distribution was fitted (less than 2 unique values of \code{attr} exist).
#' \item \code{mean} is the \emph{\mu} parameter (\emph{mean}) of the \emph{Normal} distribution (a numeric value),
#' or \code{NA} in case \code{distr != "norm"}.
#' \item \code{sd} is the \emph{\sigma} parameter (\emph{standard deviation}) of the \emph{Normal} distribution (a non-zero positive numeric value),
#' or \code{NA} in case \code{distr != "norm"}.
#' \item \code{shape} is the \emph{\alpha} parameter (\emph{shape}) of the \emph{Gamma} distribution (a non-zero positive numeric value),
#' or \code{NA} in case \code{distr != "gamma"}.
#' \item \code{rate} is the \emph{\beta} parameter (\emph{rate}) of the \emph{Gamma} distribution (a non-zero positive numeric value),
#' or \code{NA} in case \code{distr != "gamma"}.
#' \item \code{meanlog} is the \emph{\mu} parameter of the \emph{Lognormal} distribution (a numeric value),
#' or \code{NA} in case \code{distr != "lnorm"}.
#' \item \code{sdlog} is the \emph{\sigma} parameter of the \emph{Lognormal} distribution (a non-zero positive numeric value),
#' or \code{NA} in case \code{distr != "lnorm"}.
#' \item \code{BIC} is the \emph{BIC} value of the fitted distribution (a numeric value),
#' or \code{NA} in case \code{distr = NA}.
#' \item \code{dBIC} is a character string summarizing the strength of the chosen distribution model
#' specified in \code{distr} against the other models with higher BIC values.
#' Value is \code{NA} in case \code{model != "auto"} or if \code{distr = NA}.
#' \cr\cr
#' The format of the string is \code{"<dBIC_norm>, <dBIC_gamma>, <dBIC_lnorm>"}.
#' Each \code{<dBIC_model>} value is rounded.
#' The larger a \code{<dBIC_model>} value, the stronger the evidence that attribute \code{attr} of the \code{group}
#' follows the chosen \code{distr} distribution against the \code{<model>} distribution.
#' Values \code{>10} typically indicate strong preference to the chosen distribution.
#' }
#' For groups with \code{distr = NA}, no plot is generated.
#' In case no cells exist, no plot is generated and \code{NULL} is returned.
#'
#' @export
#' @import igraph
#' @import ggplot2
#' @importFrom graphics lines mtext persp plot polygon segments text
#' @importFrom grDevices adjustcolor dev.off png trans3d
#' @importFrom stats dgamma dlnorm dnorm

plot_pdf_attr <- function(tree, treeT = c("LT", "DT"),
                          attr, unit = "",
                          grouped = c("col", "gen", "pop"), groups = -1, Ngroups,
                          model = c("norm", "gamma", "lnorm", "auto"),
                          plot3D = TRUE,
                          save = FALSE, savePars = list(w = 2000, h = 2000, res = 250, path = getwd(), name = "my_pdf_attr")) {

  ################## arguments check ####################

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
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

  numeric_attrs <- get_attr_names(tree = tree, type = "n")
  if (!(attr %in% numeric_attrs) || attr %in% c("colony", "generation", "frame")) {
    stop(paste("Wrong attr \"", attr, "\"\n", sep = ""))
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

  possible_models <- c("norm", "gamma", "lnorm")
  if (!(model %in% c(possible_models, "auto"))) {
    m <- paste(possible_models, collapse = "\" / \"")
    stop(paste("model must be \"", m ,"\"/ \"auto\"\n", sep = ""))
  }

  ###################################################

  find_pdf <- function(plotData, myData, values, group, groupColor) {

    if (length(unique(x.values)) >= 2 & length(unique(values)) >= 2) {

      if (model != "auto") {

        if (model == "norm") {
          lower_parameters <- c(-Inf, 0)
        } else if (model == "gamma") {
          lower_parameters <- c(0, 0)
          values[values == 0] <- 1e-6
        } else if (model == "lnorm") {
          lower_parameters <- c(-Inf, 0)
          values[values == 0] <- 1e-6
        }

        myModel  <- fitdistrplus::fitdist(data = values,
                                          distr = model,
                                          method = "mle",
                                          lower = lower_parameters)

        if (model == "norm") {
          z.values <- dnorm(x.values,
                            mean = unname(myModel$estimate["mean"]),
                            sd = unname(myModel$estimate["sd"]),
                            log = FALSE)
        } else if (model == "gamma") {
          z.values <- dgamma(x.values,
                             shape = unname(myModel$estimate["shape"]),
                             rate = unname(myModel$estimate["rate"]),
                             log = FALSE)
        } else if (model == "lnorm") {
          z.values <- dlnorm(x.values,
                             meanlog = unname(myModel$estimate["meanlog"]),
                             sdlog = unname(myModel$estimate["sdlog"]),
                             log = FALSE)
        }

        plotData[[length(plotData) + 1]] <- list(z = z.values, col = groupColor, groupName = group)

        myData <- rbind(myData, data.frame(group = group,
                                           Ncells = length(values),
                                           distr = model,
                                           mean = ifelse(model == "norm",
                                                         unname(myModel$estimate)[1],
                                                         NA),
                                           sd = ifelse(model == "norm",
                                                       unname(myModel$estimate)[2],
                                                       NA),
                                           shape = ifelse(model == "gamma",
                                                          unname(myModel$estimate)[1],
                                                          NA),
                                           rate = ifelse(model == "gamma",
                                                         unname(myModel$estimate)[2],
                                                         NA),
                                           meanlog = ifelse(model == "lnorm",
                                                            unname(myModel$estimate)[1],
                                                            NA),
                                           sdlog = ifelse(model == "lnorm",
                                                          unname(myModel$estimate)[2],
                                                          NA),
                                           BIC = myModel$bic,
                                           dBIC = NA,
                                           stringsAsFactors = FALSE))

      } else {

        myModels <- list()

        for (i_d in 1:length(possible_models)) {

          if (possible_models[i_d] == "norm") {
            lower_parameters <- c(-Inf, 0)
          } else if (possible_models[i_d] == "gamma") {
            lower_parameters <- c(0, 0)
            values[values == 0] <- 1e-6
          } else if (possible_models[i_d] == "lnorm") {
            lower_parameters <- c(-Inf, 0)
            values[values == 0] <- 1e-6
          }

          myModels[[i_d]]  <- fitdistrplus::fitdist(data = values,
                                                    distr = possible_models[i_d],
                                                    method = "mle",
                                                    lower = lower_parameters)

        }

        bics <- sapply(myModels, "[[", "bic")
        # BIC may be negative! -> minimum numeric BIC! not by absolute value!
        # eg for the same number of data points (n) and number of parameters (k) for all models
        # BIC = ln(n)*k - 2*ln(L) -> min BIC <-> min  log likelyhood L
        best_fit <- which(bics == min(bics))

        if (possible_models[best_fit] == "norm") {
          z.values <- dnorm(x.values,
                            mean = unname(myModels[[best_fit]]$estimate["mean"]),
                            sd = unname(myModels[[best_fit]]$estimate["sd"]),
                            log = FALSE)
        } else if (possible_models[best_fit] == "gamma") {
          z.values <- dgamma(x.values,
                             shape = unname(myModels[[best_fit]]$estimate["shape"]),
                             rate = unname(myModels[[best_fit]]$estimate["rate"]),
                             log = FALSE)
        } else if (possible_models[best_fit] == "lnorm") {
          z.values <- dlnorm(x.values,
                             meanlog = unname(myModels[[best_fit]]$estimate["meanlog"]),
                             sdlog = unname(myModels[[best_fit]]$estimate["sdlog"]),
                             log = FALSE)
        }

        plotData[[length(plotData) + 1]] <- list(z = z.values, col = groupColor, groupName = group)

        Dbics <- abs(bics - min(bics)) # abs because BICs may be negative

        myData <- rbind(myData, data.frame(group = group,
                                           Ncells = length(values),
                                           distr = possible_models[best_fit],
                                           mean = ifelse(possible_models[best_fit] == "norm",
                                                         unname(myModels[[best_fit]]$estimate)[1],
                                                         NA),
                                           sd = ifelse(possible_models[best_fit] == "norm",
                                                       unname(myModels[[best_fit]]$estimate)[2],
                                                       NA),
                                           shape = ifelse(possible_models[best_fit] == "gamma",
                                                          unname(myModels[[best_fit]]$estimate)[1],
                                                          NA),
                                           rate = ifelse(possible_models[best_fit] == "gamma",
                                                         unname(myModels[[best_fit]]$estimate)[2],
                                                         NA),
                                           meanlog = ifelse(possible_models[best_fit] == "lnorm",
                                                            unname(myModels[[best_fit]]$estimate)[1],
                                                            NA),
                                           sdlog = ifelse(possible_models[best_fit] == "lnorm",
                                                          unname(myModels[[best_fit]]$estimate)[2],
                                                          NA),
                                           BIC = myModels[[best_fit]]$bic,
                                           dBIC = toString(round(Dbics)),
                                           stringsAsFactors = FALSE))

      }
    } else {

      myData <- rbind(myData, data.frame(group = group,
                                         Ncells = length(values),
                                         distr = NA,
                                         mean = NA,
                                         sd = NA,
                                         shape = NA,
                                         rate = NA,
                                         meanlog = NA,
                                         sdlog = NA,
                                         BIC = NA,
                                         dBIC = NA,
                                         stringsAsFactors = FALSE))
    }

    return(list(plotData = plotData, myData = myData))

  }


  find_pdfs <- function(plotData, myData, graphColors) {

    for (i_group in groups) {

      if (grouped == "colony") {
        groupColor <- graphColors[i_group]
      } else { # group_by == "generation"
        groupColor <- graphColors[i_group + 1]
      }

      cells <- V(tree)[V(tree)$name %in% get_cells(tree = tree, treeT = treeT, type = "inc") &
                         !is.na(vertex_attr(graph = tree, name = attr, index = V(tree))) &
                         vertex_attr(graph = tree, name = grouped, index = V(tree)) == i_group]$name

      values <- vertex_attr(graph = tree, name = attr, index = cells)

      results <- find_pdf(plotData = plotData,
                          myData = myData,
                          values = values,
                          group = i_group,
                          groupColor = groupColor)

      plotData <- results$plotData
      myData <- results$myData

    }

    return(list(plotData = plotData, myData = myData))

  }


  plot3D_pdf <- function() {

    if (length(plotData) != 0) {

      xlim <- c(floor(min(x.values) * 2) / 2, ceiling(max(x.values) * 2) / 2) # floor/ceiling to half
      ylim <- c(0, 101)
      zlim <- c(0, ceiling(max(unlist(sapply(plotData, "[[", "z"))) * 10 / 2) / 10 * 2) # ceiling to the first even decimal

      n <- 2

      # plot 3D box
      mat <- persp(x = xlim,
                   y = ylim,
                   z = matrix(0, n, n), # empty
                   zlim = zlim,
                   main = paste("Probability Density Functions of", attrTitle), xlab = "", ylab = "", zlab = "",
                   theta = 45, phi = 15,
                   d = 5, expand = 1,
                   shade = 0.3,
                   box = FALSE, axes = FALSE)

      mtext(paste(distr_name, " fitted distributions per ",
                  toupper(substr(grouped, 1, 1)), substr(grouped, 2, nchar(grouped)), sep = ""),
            side = 3, cex = 1)

      # x label
      text(trans3d((xlim[2] - xlim[1]) / 2, ylim[1] - 20, zlim[1], mat),
           labels = createAxisLab(attr = attr, unit = unit),
           srt = 337, cex = 0.9, font = 2)
      # y label
      text(trans3d(xlim[2], (ylim[2]-ylim[1]) / 2, zlim[1], mat),
           labels = paste("\n\n\n\n", toupper(substr(grouped, 1, 1)), substr(grouped, 2, nchar(grouped)), sep = ""),
           srt = 22, cex = 0.9, font = 2)
      # z label
      text(trans3d(xlim[1], ylim[1], (zlim[2] - zlim[1]) / 2, mat),
           labels = "Probability Density\n\n\n",
           srt = 92, cex = 0.9, font = 2)

      # plot lines for box
      C <- trans3d(xlim[1], ylim[1], seq(zlim[1], zlim[2], by = 0.1), mat)
      lines(C, lwd = 1)

      C <- trans3d(xlim[1], ylim[2], seq(zlim[1], zlim[2], by = 0.1), mat)
      lines(C, lwd = 1)

      C <- trans3d(xlim[2], ylim[2], seq(zlim[1], zlim[2], by = 0.1), mat)
      lines(C, lwd = 1)

      C <- trans3d(xlim[1], seq(ylim[1], ylim[2]), zlim[2], mat)
      lines(C, lwd = 1)

      C <- trans3d(seq(xlim[1], xlim[2], by = 0.1), ylim[2], zlim[2], mat)
      lines(C, lwd = 1)

      # plot lines for groups
      if (length(plotData) == 1) {
        y.positions <- 50
      } else if (length(plotData) == 2) {
        y.positions <- c(30, 70)
      } else {
        y.positions <- seq(10, 90, length = length(plotData))
      }

      for(i in 1:length(plotData)) {
        C <- trans3d(xlim, y.positions[i], zlim[1], mat)
        lines(C, lwd = 1)
      }

      # x ticks
      if (length(x.values[x.values == round(x.values)]) == length(x.values)) { # integer
        x.axis <- round(seq(xlim[1], xlim[2], length.out = 5))
      } else {
        if (xlim[2] - xlim[1] <= 1) {
          x.axis <- round(seq(xlim[1], xlim[2], length.out = 5), 2)
        } else if (xlim[2] - xlim[1] <= 5) {
          x.axis <- round(seq(xlim[1], xlim[2], length.out = 5), 1)
        } else {
          x.axis <- round(seq(xlim[1], xlim[2], length.out = 5))
        }
      }
      labels <- as.character(x.axis)

      tick.start <- trans3d(x.axis, ylim[1], zlim[1], mat)
      tick.end <- trans3d(x.axis, ylim[1] - 2, zlim[1], mat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

      label.pos <- trans3d(x.axis, ylim[1] - 5, zlim[1], mat)
      text(label.pos$x, label.pos$y, labels = labels, adj = c(0, NA), srt = 340, cex = 0.9, font = 2)

      # y ticks
      labels <- paste("\n", as.character(unlist(sapply(plotData, "[[", "groupName"))), sep = "")
      label.pos <- trans3d(xlim[2], y.positions[1:length(plotData)], zlim[1], mat)
      text(label.pos$x, label.pos$y, labels = labels, adj = c(0, NA), cex = 0.7, font = 2)

      # z ticks
      if (zlim[2] - zlim[1] <= 1) {
        z.axis <- round(seq(zlim[1], zlim[2], length.out = 5), 2)
      } else {
        z.axis <- round(seq(zlim[1], zlim[2], length.out = 5), 1)
      }

      if (xlim[1] != 0) {
        tick.start <- trans3d(xlim[1], ylim[1], z.axis, mat)
        tick.end <- trans3d(xlim[1], ylim[1] + 1.5, z.axis, mat)
        segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

        labels <- as.character(z.axis)
        label.pos <- trans3d(xlim[1], ylim[1] + 3, z.axis, mat)
        text(label.pos$x, label.pos$y, labels = labels, adj = c(0, NA), srt = 10, cex = 0.9, font = 2)
      } else {
        tick.start <- trans3d(xlim[1], ylim[1], z.axis[-1], mat)
        tick.end <- trans3d(xlim[1], ylim[1] + 1.5, z.axis[-1], mat)
        segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

        labels <- as.character(z.axis)[-1]
        label.pos <- trans3d(xlim[1], ylim[1] + 3, z.axis[-1], mat)
        text(label.pos$x, label.pos$y, labels = labels, adj = c(0, NA), srt = 10, cex = 0.9, font = 2)
      }

      # plot pdf of groups
      for(i in length(plotData):1) {

        xp <- x.values
        yp <- rep(y.positions[i], length(x.values))
        z0 <- rep(0, length(x.values)) # "upsos" pou na ksekinaei to plot
        zp <- plotData[[i]]$z

        C <- trans3d(x = c(xp, rev(xp)), y = c(yp, yp), z = c(zp, z0), pmat =  mat)
        polygon(C, border = NA, col = adjustcolor(plotData[[i]]$col, alpha.f = 0.7), density = 100)

        C <- trans3d(x = xp, y = yp, z = zp, pmat = mat)
        lines(C, col = adjustcolor(plotData[[i]]$col, alpha.f = 0.5), lwd = 2)

      }

    }

  }


  plot2D_pdf <- function() {

    if (length(plotData) != 0) {

      for(i in 1:length(plotData)) {

        df <- data.frame(var1 = x.values, var2 = plotData[[i]]$z, stringsAsFactors = FALSE)

        myPlot <- ggplot(data = df, aes_string(x = "var1", y = "var2")) +
          labs(title = paste("Probability Density Function of", attrTitle),
               subtitle = ifelse(grouped == "population",
                                 paste("Population (", distr_name, " fitted distribution)", sep = ""),
                                 paste(grouped, " ", plotData[[i]]$groupName, " (", distr_name, " fitted distribution)", sep = "")),
               x = createAxisLab(attr = attr, unit = unit),
               y = "Probability Density") +
          theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                plot.subtitle = element_text(size = 12, hjust = 0.5),
                axis.text.x = element_text(face = "bold", size = 10),
                axis.text.y = element_text(face = "bold", size = 10),
                axis.title = element_text(face = "bold", size = 12)) +
          geom_ribbon(aes(ymin = 0, ymax = df$var2), fill = plotData[[i]]$col) +
          geom_line(color = "gray14", size = 0.8)

        plot(myPlot)

      }

    }

  }

  ###################################################

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, "_%03d.png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  if (model == "norm") {
    distr_name <- "normal"
  } else if (model == "gamma") {
    distr_name <- "gamma"
  } else if (model == "lnorm") {
    distr_name <- "lognormal"
  } else if (model == "auto") {
    distr_name <- "best"
  }

  attrTitle <- getAttrForTitle(attr = attr)

  myData <- data.frame(group = integer(),
                       Ncells = integer(),
                       distr = character(),
                       mean = double(),
                       sd = double(),
                       shape = double(),
                       rate = double(),
                       meanlog = double(),
                       sdlog = double(),
                       BIC = double(),
                       dBIC = character(),
                       stringsAsFactors = FALSE)

  stp <- 500

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

    if (length(values[values < 0]) >= 1 & (model != "norm")) {
      if (save) {
        dev.off()
      }
      stop(paste("Unable to fit", model, "distribution to negative values\n"))
    }

    if (length(values[values == round(values)]) == length(values)) { # integer
      x.values <- seq(min(values), max(values))
    } else {
      x.values <- seq(min(values), max(values), length.out = stp)
    }

    results <- find_pdf(plotData = list(),
                        myData = myData,
                        values = values,
                        group = -2,
                        groupColor = getGroupColors())

    plotData <- results$plotData
    plot2D_pdf()

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

    if (length(values[values < 0]) >= 1 & (model != "norm")) {
      if (save) {
        dev.off()
      }
      stop(paste("Unable to fit", model, "distribution to negative values\n"))
    }

    if (length(values[values == round(values)]) == length(values)) { # integer
      x.values <- seq(min(values), max(values))
    } else {
      x.values <- seq(min(values), max(values), length.out = stp)
    }

    results <- find_pdfs(plotData = list(),
                         myData = myData,
                         graphColors = getGroupColors(type = grouped, N_colors = Ngroups))

    plotData <- results$plotData
    if (length(groups) == 1) {
      plot2D_pdf()
    } else {
      if (plot3D) {
        plot3D_pdf()
      } else {
        plot2D_pdf()
      }
    }

  }

  if (save) {
    dev.off()
  }

  return(results$myData)

}
