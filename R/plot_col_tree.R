#' Plot colonies' lineage tree
#'
#' Plots the colonies' lineage tree.
#'
#' @section Prerequisites:
#' This function can be used by \emph{BaSCA} or \emph{SuperSegger} users,
#' importing the data with \code{\link{import_basca}} or \code{\link{import_ss}}, respectively.
#' \cr\cr
#' Users of \emph{Oufti} who imported the data with \code{\link{import_oufti}}
#' are \bold{excluded} from using this function, as no colony list was returned.
#' \cr\cr
#' For other users, it is necessary that a colony list was imported with \code{\link{import_json}}.
#' In other case, no colony list exists and this function cannot be used.
#'
#' Nodes represent the colony instants of the movie and are colored based on the ID of the corresponding colony.
#' The imaginary \emph{root} colony instant is colored white.
#' Gray nodes are colony instants which have arised from merged colonies.
#'
#' @param col_list A list containing all the colony instants of the movie.
#'
#' @param Ncols Number of colonies in the movie, a non-zero positive integer value.
#'
#' @param sizeV Size of vertices, a non-zero positive numeric value.
#' The default value is \code{2}.
#'
#' @param sizeE Width of edges, a non-zero positive numeric value.
#' The default value is \code{0.7}.
#'
#' @param sizeL Size of explanatory legends and title, a non-zero positive numeric value.
#' The default value is \code{1}.
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
#' The default value is \code{2000}.}
#' \item{\code{h}}{The height of the image file in \emph{pixels}, a non-zero positive integer value.
#' The default value is \code{2000}.}
#' \item{\code{res}}{The resolution of the image file in \emph{pixels} per \emph{inch} (ppi), a non-zero positive integer value.
#' The smaller this value, the larger the plot area in inches, and the smaller the text relative to the graph itself.
#' The default value is \code{250}.}
#' \item{\code{path}}{A character string naming the directory where the image file will be saved.
#' If it does not contain an absolute path, the image file will be saved relative to the current working directory \code{getwd()}.
#' The default value is the current working directory \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{/} (not \code{\\}) on Windows.}
#' \item{\code{name}}{The image file name, a character string.
#' The suffix \code{".png"} is added automatically.
#' The default value is \code{"my_col_tree"}.}
#' }
#'
#' @export
#' @import igraph
#' @importFrom graphics par text
#' @importFrom grDevices dev.off png

plot_col_tree <- function(col_list, Ncols,
                          sizeV = 2, sizeE = 0.7, sizeL = 1,
                          save = FALSE, savePars = list(w = 2000, h = 2000, res = 250, path = getwd(), name = "my_col_tree")) {

  if (save) {
    png(filename = paste(savePars$path, "/", savePars$name, "_", savePars$w, "x", savePars$h, "_", savePars$res, ".png", sep = ""),
        width = savePars$w, height = savePars$h, res = savePars$res)
  }

  edges <- NULL

  # tree root
  for (i_colony in 1:Ncols) {
    edges <- rbind(edges, c(1, i_colony + 1))
  }

  for (i_colony in 1:length(col_list)) {
    next_col_Name <- col_list[[i_colony]]$next_colName
    next_col_Id <- which(sapply(col_list, function(x) x$colName == next_col_Name))
    if (length(next_col_Id) != 0) {  # next col found (next col = "f<Nframes+1>_c0" for last frame colonies in BaSCA)
      edges <- rbind(edges, c(i_colony + 1, next_col_Id + 1))
    }
  }

  nodes <- as.matrix(unique(c(edges[, 1], edges[, 2])))
  nodes <- as.data.frame(nodes[order(nodes), ], stringsAsFactors = FALSE)

  overallTree <- graph_from_data_frame(d = as.data.frame(edges), directed = TRUE, vertices = nodes)

  branches <- decompose.graph(overallTree) # motherless colony branches do not have tracked cell instants (as created by Basca)
  i_wanted_gr <- which.max(clusters(overallTree)$csize)
  colTree <- branches[[i_wanted_gr]]
  branches[[i_wanted_gr]] <- NULL

  graphColors <- getGroupColors(type = "colony", N_colors = Ncols)

  for (i_colony in 1:Ncols) {
    V(colTree)[subcomponent(graph = colTree, v = i_colony + 1, mode = "out")]$color <- graphColors[i_colony]
  }

  V(colTree)["1"]$color <- "snow3"

  colMerged <- V(colTree)[degree(graph = colTree, v = V(colTree), mode = "in") == 2]$name
  for (col in colMerged) {
    V(colTree)[subcomponent(graph = colTree, v = col, mode = "out")]$color <- "gray30"
  }

  ## plot colonies' tree
  oldPar <- par()
  par(mar = c(0, 0, 2 * sizeL, 6), bg = "black", fg = "white")

  plot.igraph(colTree,
              vertex.shape = "circle",
              layout = layout_as_tree(colTree),
              vertex.label = NA,
              vertex.size = sizeV,
              vertex.color = V(colTree)$color,
              vertex.frame.color = NA,
              edge.width = sizeE,
              edge.arrow.size = 0)

  corners <- par("usr")
  text(x = mean(corners[1:2]), y = corners[4],
       labels = expression(bold("Colonies' Tree")),
       cex = 1.5 * sizeL)

  plotColorLegend(attr = "colony", N_colors = Ncols, size_lab = sizeL)

  par(mar = oldPar$mar, bg = oldPar$bg, fg = oldPar$fg)

  if (save) {
    dev.off()
  }

}
