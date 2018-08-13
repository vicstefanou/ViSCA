#' Plot growth curves of cell counts
#'
#' Computes and plots the growth curves of cell counts (number of cells) in a lineage tree.
#' The growth curve can be plotted for each colony or for the whole population.
#' Each growth curve is computed by fitting the \emph{Baranyi and Roberts} model to the corresponding data.
#'
#' The parameters of the \emph{Baranyi and Roberts} model are found using the
#' \emph{non-linear least squares} method provided by \code{\link[stats]{nls}},
#' considering all non-root cells, as returned from \code{\link{get_cells}}.
#' \cr\cr
#' When \code{cols = -1}, a common plot for all colonies is generated.
#' In other cases, a separate plot for each colony specified in \code{cols} is generated.
#' \cr\cr
#' x-axis represents the time in \emph{hours} from the start (value \code{0}) to end of the movie.
#' y-axis represents the cell counts in logarithmic scale.
#' The range of y-axis values depicted in each plot is common and
#' is calculated as the range of values of all specified colonies,
#' with minimum upper limit \code{2} (i.e. 20 cells).
#' \cr\cr
#' Color denotes the corresponding colony.
#'
#' @param LT The lineage tree, an object of class \code{"igraph"}.
#'
#' @param DT The corresponding division tree of the \code{LT},
#' an object of class \code{"igraph"}.
#'
#' @param cols The IDs of the colonies for which to fit the \emph{Baranyi and Roberts} model,
#' a vector of non-zero positive integer values.
#' The default value \code{-1} stands for all existing colonies in the \code{LT}.
#' Use value \code{-2} for fitting the model to the whole population.
#'
#' @param Ncols Number of colonies in the movie, a non-zero positive integer value.
#' This argument is ignored in case \code{cols = -2}.
#' @param Nframes Number of frames in the movie, a non-zero positive integer value.
#' @param frameR Frame rate of the movie in \emph{frames} per \emph{minute}, a non-zero positive numeric value.
#'
#' @param showRaw A logical value (\code{TRUE} or \code{FALSE}) indicating whether
#' the raw data (unfitted data points) will be shown on the plot or not, respectively.
#' The default value is \code{FALSE}.
#' This argument is ignored (regarded as \code{FALSE}) in case \code{cols = -1}.
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
#' The default value is \code{3000}.}
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
#' The default value is \code{"my_baranyi"}.}
#' }
#'
#' @return A dataframe with the following columns:
#' \enumerate{
#' \item \code{colony} is the colony ID (a non-zero positive integer value)
#' or \code{-2} in case \code{cols = -2}.
#' \item \code{lag} is the \emph{\ifelse{html}{\out{lambda}}{\eqn{\lambda}}} parameter (\emph{lag time})
#' of the \emph{Baranyi and Roberts} model in \emph{hours} (a non-zero positive numeric value)
#' or \code{NA} in case the model failed to be fitted or \code{colony} has only cells which have not been divided.
#' \item \code{mumax} is the \emph{\ifelse{html}{\out{mumax}}{\eqn{\mu_{max}}}} parameter (\emph{maximum specific growth rate})
#' of the \emph{Baranyi and Roberts} model in \emph{1/hour} (a non-zero positive numeric value)
#' or \code{NA} in case the model failed to be fitted or \code{colony} has only cells which have not been divided.
#' }
#' For groups with \code{lag = NA}, no plot is generated, except for the case that \code{showRaw = TRUE}.
#' In case no cells exist, no plot is generated and \code{NULL} is returned.
#'
#' @export
#' @import igraph
#' @import data.table
#' @importFrom graphics abline axTicks legend lines mtext par plot points
#' @importFrom grDevices dev.off png
#' @importFrom stats coef lm nls
#' @importFrom utils tail

plot_baranyi <- function(LT, DT,
                         cols = -1, Ncols,
                         Nframes, frameR,
                         showRaw = FALSE,
                         save = FALSE, savePars = list(w = 3000, h = 2000, res = 300, path = getwd(), name = "my_baranyi")) {

  ################## arguments check ####################

  possible_cols <- sort(unique(V(LT)$colony))
  possible_cols <- possible_cols[possible_cols != -1]
  if (length(cols) != 1) {
    if (length(m <- setdiff(cols, possible_cols)) != 0) {
      stop(paste("Selected colony(s) ", toString(m), " do not exist\n", sep = ""))
    }
  } else {
    if (!(cols %in% c(possible_cols, -1, -2))) {
      stop(paste("Selected colony ", cols, " does not exist\n", sep = ""))
    }
  }

  #######################################

  createTable <- function(cells) {

    dt <- as.data.table(V(LT)[cells]$frame)

    # Data table with 2 columns (t: time (here in frames), LOG10N: decimal logarithm of bacterial density)
    dt <- dt[, .N, by = "V1"]
    colnames(dt) <- c("t", "LOG10N")
    dt$LOG10N <- log10(dt$LOG10N)
    dt$t <- (dt$t - 1) / frameR # t: time in hours starting from 0

    return(dt)

  }


  estimate_lag_mumax <- function(myColonies) {

    # lag time (h) : mean of life time of generation 0
    # (without cells that have not been divided --> divisionTime = Nframes)
    myLag <- V(DT)[V(DT)$colony %in% myColonies &
                     V(DT)$generation == 0 &
                     V(DT)$divisionTime != Nframes]$lifeHours

    if (length(myLag) != 0) {
      myLag <- mean(myLag)

      # mumax : slope of curve (except for lag and stationary phase)
      dt <- dt[dt$t >= myLag,]

      if (dim(dt)[1] >= 2) {
        myModel <- lm(dt$LOG10N ~ dt$t)
        myMumax <- coef(myModel)[2]
      } else {
        myMumax <- NULL
      }

    } else {
      myLag <- NULL
      myMumax <- NULL
    }

    return(list(lag = myLag, mumax = myMumax))

  }


  fitBaranyi <- function(params, title) {

    if (is.null(params$lag)) { # has only cells which have not been divided
      myModel <- NULL
    } else {
      myModel <- tryCatch(nls(formula = nlsMicrobio::baranyi, data = dt,
                              start = list(lag = params$lag, mumax = params$mumax, LOG10N0 = dt$LOG10N[1], LOG10Nmax = tail(dt$LOG10N, 1))),
                          error = function(e) {
                            # does not fit Roberts-Baranyi model
                          })
    }

    if (is.null(myModel)) {
      return(NULL)
    } else {

      myModel_vars <- all.vars(nlsMicrobio::baranyi[[3]])

      for(var in myModel_vars) {
        if (var == "t") {
          t <- dt$t
        } else {
          assign(var, unname(coef(myModel)[var]))
        }
      }

      return(list(values = eval(nlsMicrobio::baranyi[[3]]), lag = coef(myModel)["lag"], mumax = coef(myModel)["mumax"]))
    }

  }


  findYlimit <- function() {

    dt <- NULL

    for (i_colony in possible_cols) {
      colony_last_frame <- max(V(LT)[V(LT)$colony == i_colony]$frame)
      dt <- rbind(dt, as.data.table(V(LT)[V(LT)$colony == i_colony &
                                            V(LT)$frame == colony_last_frame]$colony))
    }

    # Data table with 2 columns (colony, LOG10N: decimal logarithm of bacterial density for last frame of each colony)
    dt <- dt[, .N, by = "V1"]
    colnames(dt) <- c("colony", "LOG10N")
    dt$LOG10N <- log10(dt$LOG10N)

    if ((y_upper_lim <- ceiling(max(dt$LOG10N))) < 2) {  # minimum upper ylim --> 2 (20 cells)
      y_upper_lim <- 2
    }

    return(y_upper_lim)
  }


  myFun <- function(cells, myColony, myColor, title, mainTitle, possible_cols) {

    if (save) {
      png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, "_%03d.png", sep = ""),
          width = savePars$w, height = savePars$h, res = savePars$res)
    }

    oldPar <- par()
    par(mar = c(4.5, 4.5, 4, 1))

    dt <- createTable(cells = cells)

    baranyi <- fitBaranyi(params = estimate_lag_mumax(myColonies = myColony),
                          title = title)

    if (!is.null(baranyi)) {

      yBaranyi <- baranyi$values

      if (title == "Population") {
        ylim <- c(min(yBaranyi[1], dt$LOG10N[1]), max(tail(yBaranyi, 1), tail(dt$LOG10N, 1)))
      } else {
        ylim <- c(0, findYlimit())
      }

      plot(x = dt$t, y = yBaranyi,
           main = mainTitle, xlab = "Time (h)", ylab = bquote(bold("Log"["10"] ~ "( Number of Cells )")),
           xlim = c(0, Nframes / frameR), ylim = ylim,
           cex.main = 1.5, cex.lab = 1, font = 2, font.lab = 2, las = 1,
           type = "l", lwd = 3, col = myColor)

      abline(h = axTicks(side = 2), v = axTicks(side = 1), col = "gray", lty = "dotted")

      if (showRaw == TRUE) {
        points(x = dt$t, y = dt$LOG10N, col = myColor)
        legend("topleft", cex = 1,
               legend = c("fitted Roberts and Baranyi model", "raw"),
               pch = c(NA, 1),
               lwd = c(2, NA),
               col = myColor,
               xpd = TRUE, inset = 0.01)
      } else {
        legend("topleft", cex = 1,
               legend = "fitted Roberts and Baranyi model",
               lwd = 2, col =  myColor,
               xpd = TRUE, inset = 0.01)
      }

    } else {
      if (showRaw == TRUE) {

        if (title == "Population") {
          ylim <- c(dt$LOG10N[1], tail(dt$LOG10N, 1))
        } else {
          ylim <- c(0, findYlimit())
        }

        plot(x = dt$t, y = dt$LOG10N,
             main = mainTitle, xlab = "Time (h)", ylab = bquote(bold("Log"["10"] ~ "( number of cells )")),
             xlim = c(0, Nframes / frameR), ylim = ylim,
             cex.main = 1.5, cex.lab = 1, font = 2, font.lab = 2, las = 1,
             pch = 1, col = myColor)

        abline(h = axTicks(side = 2), v = axTicks(side = 1), col = "gray", lty = "dotted")
        legend("topleft", cex = 1,
               legend = "raw",
               pch = 1, col =  myColor,
               xpd = TRUE, inset = 0.01)
      }
    }

    par(mar = oldPar$mar)

    if (save) {
      dev.off()
    }

    myData <- rbind(myData, data.frame(colony = ifelse(title == "Population", -2, myColony),
                                       lag = ifelse(is.null(baranyi), NA, baranyi$lag),
                                       mumax = ifelse(is.null(baranyi), NA, baranyi$mumax),
                                       stringsAsFactors = FALSE))

    return(myData)

  }

  #######################################

  frameR <- frameR * 60 # in frames/hour

  myData <- data.frame(colony = integer(),
                       lag = double(),
                       mumax = double(),
                       stringsAsFactors = FALSE)

  possible_cells <- get_cells(tree = LT, treeT = "LT", type = "nr")

  if (length(possible_cells) == 0) {
    return(NULL)
  }

  if (length(cols) == 1) {

    if (cols == -1) { # all colonies in one plot

      if (save) {
        png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, ".png", sep = ""),
            width = savePars$w, height = savePars$h, res = savePars$res)
      }

      oldPar <- par()
      par(mar = c(4.5, 4.5, 4, 1))

      graphColors <- getGroupColors(type = "colony", N_colors = Ncols)
      y_upper_lim <- findYlimit()

      colonies_plotted <- NULL

      for (i_colony in possible_cols) {

        dt <- createTable(cells = V(LT)[V(LT)$name %in% possible_cells &
                                          V(LT)$colony == i_colony]$name)

        baranyi <- fitBaranyi(params = estimate_lag_mumax(myColonies = i_colony),
                              title = paste("colony", i_colony))

        if (!is.null(baranyi)) {

          yBaranyi <- baranyi$values

          if (is.null(colonies_plotted)) { # first colony to plot
            plot(x = dt$t, y = yBaranyi,
                 main = "Growth Curves per Colony", xlab = "Time (h)", ylab = bquote(bold("Log"["10"] ~ "( number of cells )")),
                 xlim = c(0, Nframes / frameR), ylim = c(0, y_upper_lim),
                 cex.main = 1.5, cex.lab = 1, font = 2, font.lab = 2, las = 1,
                 type = "l", lwd = 3, col = graphColors[i_colony])
          } else {
            lines(dt$t, yBaranyi,
                  xlim = c(0, Nframes / frameR), ylim = c(0, y_upper_lim),
                  lwd = 3, col = graphColors[i_colony])
          }

          colonies_plotted <- c(colonies_plotted, i_colony)
        }

        myData <- rbind(myData, data.frame(colony = i_colony,
                                           lag = ifelse(is.null(baranyi), NA, baranyi$lag),
                                           mumax = ifelse(is.null(baranyi), NA, baranyi$mumax),
                                           stringsAsFactors = FALSE))

      }

      if (!is.null(colonies_plotted)) { # plot has been created
        abline(h = axTicks(side = 2), v = axTicks(side = 1), col = "gray", lty = "dotted")
        legend("topleft", cex = 1,
               title = expression(bold("Colony")),
               legend = colonies_plotted,
               ncol = 2, lwd = 2,
               col =  graphColors[colonies_plotted],
               xpd = TRUE, inset = 0.01)
        mtext("fitted Roberts and Baranyi models", side = 3, cex = 1)
      }

      par(mar = oldPar$mar)

      if (save) {
        dev.off()
      }

    } else if (cols == -2) { # all population

        myData <- myFun(cells = possible_cells,
                        myColony = possible_cols,
                        myColor = getGroupColors(),
                        title = "Population",
                        mainTitle = "Population Growth Curve",
                        possible_cols = NULL)

    } else { # one colony

      myData <- myFun(cells = V(LT)[V(LT)$name %in% possible_cells &
                                      V(LT)$colony == cols]$name,
                      myColony = cols,
                      myColor = getGroupColors(type = "colony", N_colors = Ncols)[cols],
                      title = paste("colony", cols),
                      mainTitle = paste("Growth Curve of colony", cols),
                      possible_cols = cols)

    }

  } else { # multiple colonies separately

    cols <- sort(cols)

    for (i_colony in cols) {

      myData <- myFun(cells = V(LT)[V(LT)$name %in% possible_cells &
                                      V(LT)$colony == i_colony]$name,
                      myColony = i_colony,
                      myColor = getGroupColors(type = "colony", N_colors = Ncols)[i_colony],
                      title = paste("colony", i_colony),
                      mainTitle = paste("Growth Curve of colony", i_colony),
                      possible_cols = cols)

    }

  }

  return(myData)

}
