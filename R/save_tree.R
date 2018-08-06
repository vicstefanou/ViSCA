#' Save a tree on disc
#'
#' Saves a lineage or division tree on disc in nodes and edges \code{.csv} files.
#'
#' @param tree The linegae or division tree to be saved, an object of class \code{"igraph"}.
#'
#' @param savePath A character string naming the directory where the \code{"nodes.csv"} and \code{"edges.csv"} files will be saved.
#' If it does not contain an absolute path, the files will be saved relative to the current working directory, \code{getwd()}.
#' The default value is the current working directory \code{getwd()}.
#' \cr\cr
#' NOTE: The components should be separated by \code{/} (not \code{\\}) on Windows.
#'
#' @param prefix A prefix that will be added to the name of the \code{"nodes.csv"} and \code{"edges.csv"} files,
#' a character string.
#' The default value is the empty string \code{""}.
#'
#' @param sep The field separator character.
#' Values on each line of the \code{"nodes.csv"} and \code{"edges.csv"} files will be separated by this character.
#' The default value is the tab separator \code{"\t"}.
#'
#' @export
#' @import igraph
#' @importFrom utils write.table

save_tree <- function(tree, savePath = getwd(), prefix = "", sep = "\t") {

  tree_nodes <- as_data_frame(x = tree, what = "vertices")
  tree_edges <- as_data_frame(x = tree, what = "edges")

  write.table(tree_nodes, file = paste(savePath, "/", prefix, "nodes.csv", sep = ""), quote = TRUE, sep = sep, na = "NA", row.names = FALSE, col.names = TRUE)
  write.table(tree_edges, file = paste(savePath, "/", prefix, "edges.csv", sep = ""), quote = TRUE, sep = sep, na = "NA", row.names = FALSE, col.names = TRUE)

  cat(paste("Graph was saved in ", savePath, "/", prefix, "nodes.csv", "\n",
            "and ", savePath, "/", prefix, "edges.csv", "\n", sep = ""))

}
