#' Select a subtree
#'
#' Selects a subtree of a lineage or division tree
#' based on multiple selection criteria combined with logical AND operator(s).
#'
#' The selection is made on the nodes of the tree.
#' Only edges whose endpoints are both included in the resulting subtree are included.
#' Therefore, the resulting subtree may be disconnected.
#'
#' @param tree The lineage or division tree from which to select a subtree, an object of class \code{"igraph"}.
#'
#' @param criteria A list containing the selection criteria.
#' Each criterion (element of the list) should be a named list with the following components:
#' \describe{
#' \item{\code{attr}}{A character string naming the attribute in the \code{tree} on which the selection will be made.
#' It can be any numeric or boolean attribute, as returned from \code{\link{get_attr_names}}.}
#' \item{\code{val}}{The value of comparison.
#' It should be a numeric value when \code{attr} is numeric and
#' a logical value (\code{TRUE} or \code{FALSE}) when \code{attr} is boolean.}
#' \item{\code{op}}{The comparison operator.
#' It should be \code{"=="}, \code{"!="}, \code{"<"}, \code{"<="}, \code{">"} or \code{">="} when \code{attr} is numeric and
#' \code{"=="} or \code{"!="} when \code{attr} is boolean.}
#' }
#'
#' @return The selected subtree, an object of class \code{"igraph"}.
#'
#' @seealso \code{\link{isConnected}} for checking if a tree is connected,
#' \code{\link{unite_trees}} for combining multiple selection criteria with logical OR operator(s).
#' @export
#' @import igraph

select_subtree <- function(tree, criteria) {

  ## subtree of tree myDT with cells of generation <= 5
  ## which are considered in the analysis
  # myCriteria <- list()
  # myCriteria[[1]] <- list(attr = "generation", val = 5, op = "<=")
  # myCriteria[[2]] <- list(attr = "isConsidered", val = TRUE, op = "==")
  # DTnew <- select_subtree(tree = myDT, criteria = myCriteria)

  ################# arguments check also in code ##############

  for (i in 1:length(criteria)) {
    d <- setdiff(c("attr", "val", "op"), names(criteria[[i]]))
    if (length(d) != 0){
      stop(paste("Missing element(s)", toString(d), "from criteria", i, "\n"))
    }
  }

  ######################################

  numeric_attrs <- get_attr_names(tree = tree, type = "n")
  boolean_attrs <- get_attr_names(tree = tree, type = "b")

  for (i in 1:length(criteria)) {

    if (criteria[[i]]$attr %in% numeric_attrs) {

      if (is.numeric(criteria[[i]]$val)) {

        if (criteria[[i]]$op %in% c("==", "!=", "<", "<=", ">", ">=")) {

          cells <- eval(parse(text = paste("V(tree)[vertex_attr(graph = tree, name = criteria[[i]]$attr, index = V(tree))",
                                           criteria[[i]]$op,
                                           "criteria[[i]]$val &
                                           !is.na(vertex_attr(graph = tree, name = criteria[[i]]$attr, index = V(tree)))]$name")))

          tree <- induced_subgraph(graph = tree, vids = cells)

        } else {
          stop(paste("Wrong op \"", criteria[[i]]$op, "\" in criteria ", i, "\n",
                     "Accepted op : \"==\" / \"!=\" / \"<\" / \"<=\" / \">\" / \">=\"\n", sep = ""))
        }

      } else {
        stop(paste("Wrong val \"", criteria[[i]]$val, "\" in criteria ", i, "\n",
                   "Accepted val : numeric\n", sep = ""))
      }

    } else if (criteria[[i]]$attr %in% boolean_attrs) {

      if (is.logical(criteria[[i]]$val)) {

        if (criteria[[i]]$op %in% c("==", "!=")) {

          cells <- eval(parse(text = paste("V(tree)[vertex_attr(graph = tree, name = criteria[[i]]$attr, index = V(tree))",
                                           criteria[[i]]$op,
                                           "criteria[[i]]$val")))

          tree <- induced_subgraph(graph = tree, vids = cells)

        } else {
          stop(paste("Wrong op \"", criteria[[i]]$op, "\" in criteria ", i, "\n",
                     "Accepted op : \"==\" / \"!=\"\n", sep = ""))
        }

      } else {
        stop(paste("Wrong val \"", criteria[[i]]$val, "\" in criteria ", i, "\n",
                   "Accepted val : TRUE / FALSE\n", sep = ""))
      }

    } else {
      stop(paste("Wrong attr in criteria ", i, "\n",
                 "attr \"", criteria[[i]]$attr, "\" is not numeric or boolean\n", sep = ""))
    }
  }

  return(tree)

}
