% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isSubtree.R
\name{isSubtree}
\alias{isSubtree}
\title{Check if a tree is a subtree of another tree}
\usage{
isSubtree(subtree, tree, treeT = c("LT", "DT"), type = c("all", "nr",
  "inc"))
}
\arguments{
\item{subtree}{The lineage or division tree which is supposed to be the subtree of \code{tree},
an object of class \code{"igraph"}.}

\item{tree}{The lineage or division tree which is supposed to be the supertree of \code{subtree},
an object of class \code{"igraph"}.}

\item{treeT}{A character string naming the type of \code{subtree} and \code{tree}:
\itemize{
\item \code{"LT"} if \code{subtree} and \code{tree} are lineage trees
\item \code{"DT"} if \code{subtree} and \code{tree} are division trees
}}

\item{type}{A character string naming the type of cells regarding which the statement will be checked:
\itemize{
\item \code{"all"} for all cells (including any existing imaginary \emph{root} cells)
\item \code{"nr"} for all non-root cells (excluding any existing imaginary \emph{root} cells)
\item \code{"inc"} for all cells included in the analysis
}}
}
\value{
A logical value (\code{TRUE} or \code{FALSE})
indicating whether the \code{subtree} is a subtree of the \code{tree} or not.
}
\description{
Checks if a lineage or division tree is a subtree of another lineage or division tree,
regarding only the nodes of the trees.
}
