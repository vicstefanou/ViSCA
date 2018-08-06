#' Compute growth curves of an attribute
#'
#' Computes the growth curve of a numeric attribute for each cell in a division tree,
#' by fitting a \emph{linear} or \emph{exponential} model to the its attribute's time-series.
#'
#' The estimated parameters as well as the RMSE of the fitted growth model
#' are added as attributes to the \code{DT}:
#' \itemize{
#' \item When \code{model = "lin"}:
#' \itemize{
#' \item \code{"<attr>_a"}, a non-zero positive numeric value in units of \code{attr} per \emph{hour}
#' \item \code{"<attr>_b"}, a positive numeric value in units of \code{attr}
#' \item \code{"<attr>_linRMSE"}
#' }
#' \item When \code{model = "exp"}:
#' \itemize{
#' \item \code{"<attr>_k"}, a non-zero positive numeric value in units of \code{attr} per \emph{hour}
#' \item \code{"<attr>_0"}, a non-zero positive numeric value in units of \code{attr}
#' \item \code{"<attr>_expRMSE"}
#' }
#' }
#' \code{NA} values are stored for cells that failed to fit the selected model
#' as well as for cells that are not included in the analysis, as returned from \code{\link{get_cells}}.
#' Information messages for cells that failed to fit the model are printed on the screen.
#'
#' @param LT The lineage tree, an object of class \code{"igraph"}.
#'
#' @param DT The corresponding division tree of the \code{LT}, an object of class \code{"igraph"}.
#'
#' @param attr The name of the attribute in the \code{LT}, a character string.
#' It can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' except for \code{"colony"}, \code{"generation"}, \code{"frame"} and \code{"age"}.
#'
#' @param model A character string naming the growth model to be fitted:
#' \itemize{
#' \item \code{"lin"} for fitting a \emph{linear} model
#' \code{y = a*t + b}
#' using the \emph{linear least squares} method provided by \code{\link[stats]{lm}}
#' \item \code{"exp"} for fitting an \emph{exponential} model
#' \code{y = y0*e^(kt)}
#' using the \emph{non-linear least squares} method provided by \code{\link[stats]{nls}}
#' }
#' where \code{y} represents the \code{attr} and \code{t} the time in \emph{hours} (starting from \code{0}).
#'
#' @param frameR Frame rate of the movie in \emph{frames} per \emph{minute}, a non-zero positive numeric value.
#'
#' @return The updated \code{DT} with the new attributes added, an object of class \code{"igraph"}.
#'
#' @export
#' @import igraph
#' @importFrom stats coef lm nls resid
#' @importFrom utils setTxtProgressBar txtProgressBar

add_attr_growth_fit_pars <- function(LT, DT, attr, model = c("lin", "exp"), frameR) {

  ############ arguments check ###########################

  if (!(model %in% c("lin", "exp"))) {
    stop("model must be \"lin\" / \"exp\"\n")
  }

  numeric_attrs <- get_attr_names(tree = LT, type = "n")

  if (!(attr %in% numeric_attrs) || attr %in% c("colony", "generation", "frame", "age")) {
    stop(paste("Wrong attr \"", attr, "\"\n", sep = ""))
  }

  ######################

  cells <- get_cells(tree = DT, treeT = "DT", type = "inc")

  if (length(cells) == 0) {
    stop("No cells in DT\n")
  }

  pb <- txtProgressBar(min = 0, max = length(cells), style = 3) ### set progress bar
  ipb <- 0

  scrM <- "" # screen message

  for (cell in cells) {

    ipb <- ipb + 1
    setTxtProgressBar(pb, ipb) ### update progress bar

    values <- as.numeric(unlist(strsplit(vertex_attr(graph = DT, name = attr, index = cell), ", ")))
    frameR <- frameR * 60 # se frames/hour
    time <- seq(0, (V(DT)[cell]$lifeFrames - 1) / frameR, 1 / frameR) # time in hours starting from 0

    if (model == "lin") {

      if (V(DT)[cell]$lifeFrames >= 2) {

        myModel <- lm(values ~ time)

        if (coef(myModel)[2] > 0 & coef(myModel)[1] >= 0) {

          DT <- set_vertex_attr(graph = DT, name = paste(attr, "a", sep = "_"), index = cell,
                                value = coef(myModel)[2])
          DT <- set_vertex_attr(graph = DT, name = paste(attr, "b", sep = "_"), index = cell,
                                value = coef(myModel)[1])
          DT <- set_vertex_attr(graph = DT, name = paste(attr, "linRMSE", sep = "_"), index = cell,
                                value = sqrt(sum(resid(myModel) ^ 2) / length(values)))

        } else if (coef(myModel)[2] <= 0) {
          newM <- paste("cell", paste("\"", cell, "\"", sep = ""), "has unchanged or decreasing", attr, "\n")
          scrM <- paste0(scrM, newM)
        } else {
          newM <- paste("cell", paste("\"", cell, "\"", sep = ""), "has negative", paste(attr, "b", sep = "_"), "\n")
          scrM <- paste0(scrM, newM)
        }
      } else {
        newM <- paste("cell", paste("\"", cell, "\"", sep = ""), "lives for", V(DT)[cell]$lifeFrames, "frame(s)\n")
        scrM <- paste0(scrM, newM)
      }
    } else { # model == "exp"

      if (V(DT)[cell]$lifeFrames >= 3) {

        df <- as.data.frame(cbind(time, values))

        myModel <- tryCatch(nls(formula = values ~ value0 * exp(k * time), data = df, start = list(k = 1, value0 = values[1])),
                            error = function(e) {
                              # does not fit an exponential model
                            })

        if (is.null(myModel)) {
          newM <- paste("cell", paste("\"", cell, "\"", sep = ""), "has", attr, "that does not fit an exponential model\n")
          scrM <- paste0(scrM, newM)
        } else {

          if (coef(myModel)[1] > 0 & coef(myModel)[2] > 0) {

            DT <- set_vertex_attr(graph = DT, name = paste(attr, "k", sep = "_"), index = cell,
                                  value = coef(myModel)[1])
            DT <- set_vertex_attr(graph = DT, name = paste(attr, "0", sep = "_"), index = cell,
                                  value = coef(myModel)[2])
            DT <- set_vertex_attr(graph = DT, name = paste(attr, "expRMSE", sep = "_"), index = cell,
                                  value = sqrt(sum(resid(myModel) ^ 2) / length(values)))

          } else if (coef(myModel)[1] <= 0) {
            newM <- paste("cell", paste("\"", cell, "\"", sep = ""), "has unchanged or decreasing", attr, "\n")
            scrM <- paste0(scrM, newM)
          } else {
            newM <- paste("cell", paste("\"", cell, "\"", sep = ""), "has zero or negative", paste(attr, "0", sep = "_"), "\n")
            scrM <- paste0(scrM, newM)
          }
        }
      } else {
        newM <- paste("cell", paste("\"", cell, "\"", sep = ""), "lives for", V(DT)[cell]$lifeFrames, "frame(s)\n")
        scrM <- paste0(scrM, newM)
      }
    }
  }

  close(pb) ### close progress bar
  cat("\n")

  cat(scrM)

  return(DT)

}
