#' Get statistic of an attribute
#'
#' Returns the statistic of a numeric attribute of a lineage or division tree.
#'
#' The statistic is calculated considering all cells that are included in the analysis,
#' as returned from \code{\link{get_cells}},
#' except for cells with \code{NA} value in \code{attr}.
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
#' @param stat A character string naming the statistic to be returned:
#' \itemize{
#' \item \code{"mean"} for the \emph{mean}
#' \item \code{"median"} for the \emph{median}
#' \item \code{"sd"} for the \emph{standard deviation}
#' \item \code{"min"} for the minimum value
#' \item \code{"max"} for the maximum value
#' }
#'
#' @return A named list with the following components:
#' \item{Ncells}{Number of cells, a positive integer value.}
#' \item{value}{The corresponding statistic of attribute \code{attr}, a numeric value,
#' or \code{NA}.}
#'
#' @export
#' @import igraph

get_attr_stats <- function(tree, treeT = c("LT", "DT"),
                           attr,
                           stat = c("mean", "median", "sd", "min", "max")) {

  ################## arguments check #######################

  if (!(treeT %in% c("LT", "DT"))) {
    stop("treeT must be \"LT\" / \"DT\"\n")
  }

  if (!(stat %in% c("mean", "median", "sd", "min", "max"))) {
    stop("stat must be \"mean\" / \"median\" / \"sd\" / \"min\" / \"max\"\n")
  }

  numeric_attrs <- get_attr_names(tree = tree, type = "n")

  if (!(attr %in% numeric_attrs) || attr %in% c("colony", "generation", "frame")) {
    stop(paste("Wrong attr \"", attr, "\"\n", sep = ""))
  }

  ###################################################

  cells <- V(tree)[V(tree)$name %in% get_cells(tree = tree, treeT = treeT, type = "inc") &
                     !is.na(vertex_attr(graph = tree, name = attr, index = V(tree)))]$name

  if (length(cells) == 0) {
    return(list(Ncells = 0, value = NA))
  } else {
    vars <- vertex_attr(graph = tree, name = attr, index = cells)
    return(list(Ncells = length(cells),
                value = eval(parse(text = paste(stat, "(vars)", sep = "")))))
  }

}
