#' Add ROC attribute
#'
#' Calculates the instantaneous Rate Of Change (ROC) of a numeric attribute of a lineage tree.
#' The ROC is calculated for each cell in the tree and can be positive or negative.
#'
#' The calculated ROC is added as an attribute to the \code{LT}:
#' \itemize{
#' \item \code{"d<attr>_norm"}, a numeric value in the range \code{[-1, 1]} in arbitrary units, when \code{norm = TRUE}
#' \item \code{"d<attr>"}, a numeric value in units of \code{attr} per \emph{hour}, when \code{norm = FALSE}
#' }
#' \code{NA} values are stored for cells that are just born
#' as well as for cells that are not included in the analysis, as returned from \code{\link{get_cells}}.
#'
#' @param LT The connected lineage tree, an object of class \code{"igraph"}.
#'
#' @param attr The name of the attribute in the \code{LT}, a character string.
#' It can be any numeric attribute, as returned from \code{\link{get_attr_names}},
#' except for \code{"colony"}, \code{"generation"}, \code{"frame"} and \code{"age"}.
#'
#' @param norm A logical value (\code{TRUE} or \code{FALSE}) indicating the type of the ROC that will be calculated.
#' When the default value \code{TRUE} is used, ROC is normalized and represents
#' the percentage of change relative to the previous frame.
#' When value \code{FALSE} is used, ROC represents the change per hour.
#'
#' @param frameR Frame rate of the movie in \emph{frames} per \emph{minute}, a non-zero positive numeric value.
#' This argument is ignored in case \code{norm = FALSE}.
#'
#' @return The updated \code{LT} with the new attribute added, an object of class \code{"igraph"}.
#'
#' @seealso \code{\link{isConnected}} for checking if a tree is connected.
#' @export
#' @import igraph

add_attr_roc <- function(LT, attr, norm = TRUE, frameR) {

  ############ arguments check ###########################

  if (!isConnected(tree = LT)) {
    stop("LT is disconnected\n")
  }

  numeric_attrs <- get_attr_names(tree = LT, type = "n")

  if (!(attr %in% numeric_attrs) || attr %in% c("colony", "generation", "frame", "age")) {
    stop(paste("Wrong attr \"", attr, "\"\n", sep = ""))
  }

  ######################

  accepted_mother_cells <- V(LT)[V(LT)$name %in% get_cells(tree = LT, treeT = "LT", type = "inc") &
                                   degree(graph = LT, v = V(LT), mode = "out") == 1]$name

  if (length(accepted_mother_cells)) {
    stop("No cells in LT\n")
  }

  cells <- V(LT)[unlist(ego(graph = LT, order = 1, nodes = accepted_mother_cells, mode = "out", mindist = 1))]$name

  dAttr_frames <- vertex_attr(graph = LT, name = attr, index = cells) -
    vertex_attr(graph = LT, name = attr, index = accepted_mother_cells)

  if (norm) {
    LT <- set_vertex_attr(graph = LT, name = paste("d", attr, "_norm", sep = ""), index = cells,
                          value = dAttr_frames / vertex_attr(graph = LT, name = attr, index = accepted_mother_cells))
  } else {
    frameR <- frameR * 60 # in frames/hour
    LT <- set_vertex_attr(graph = LT, name = paste("d", attr, sep = ""), index = cells,
                          value = dAttr_frames * frameR)
  }

  return(LT)

}
