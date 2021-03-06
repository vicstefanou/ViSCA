% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_attr_roc.R
\name{add_attr_roc}
\alias{add_attr_roc}
\title{Add ROC attribute}
\usage{
add_attr_roc(LT, attr, norm = TRUE, frameR)
}
\arguments{
\item{LT}{The connected lineage tree, an object of class \code{"igraph"}.}

\item{attr}{The name of the attribute in the \code{LT}, a character string.
It can be any numeric attribute, as returned from \code{\link{get_attr_names}},
except for \code{"colony"}, \code{"generation"}, \code{"frame"} and \code{"age"}.}

\item{norm}{A logical value (\code{TRUE} or \code{FALSE}) indicating the type of the ROC that will be calculated.
When the default value \code{TRUE} is used, ROC is normalized and represents
the percentage of change relative to the previous frame.
When value \code{FALSE} is used, ROC represents the change per hour.}

\item{frameR}{Frame rate of the movie in \emph{frames} per \emph{minute}, a non-zero positive numeric value.
This argument is ignored in case \code{norm = FALSE}.}
}
\value{
The updated \code{LT} with the new attribute added, an object of class \code{"igraph"}.
}
\description{
Calculates the instantaneous Rate Of Change (ROC) of a numeric attribute of a lineage tree.
The ROC is calculated for each cell in the tree and can be positive or negative.
}
\details{
The calculated ROC is added as an attribute to the \code{LT}:
\itemize{
\item \code{"d<attr>_norm"}, a numeric value in the range \code{[-1, 1]} in arbitrary units, when \code{norm = TRUE}
\item \code{"d<attr>"}, a numeric value in units of \code{attr} per \emph{hour}, when \code{norm = FALSE}
}
\code{NA} values are stored for cells that are just born
as well as for cells that are not included in the analysis, as returned from \code{\link{get_cells}}.
}
\seealso{
\code{\link{isConnected}} for checking if a tree is connected.
}
