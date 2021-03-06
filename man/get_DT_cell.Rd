% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_DT_cell.R
\name{get_DT_cell}
\alias{get_DT_cell}
\title{Get the label of the cell of a cell instant}
\usage{
get_DT_cell(LT, cell, DT)
}
\arguments{
\item{LT}{The lineage tree where the cell specified in \code{cell} belongs,
an object of class \code{"igraph"}.}

\item{cell}{The label of the cell in the \code{LT}, a character string.
It can be any non-root cell, as returned from \code{\link{get_cells}}.}

\item{DT}{The corresponding division tree of the \code{LT},
an object of class \code{"igraph"}.}
}
\value{
A named list with the following components:
\item{cell}{The label of the \code{cell} in the \code{DT}, a character string.}
\item{instant}{The instant of the cell, a non-zero positive integer value.}
In case \code{cell} is not found in the \code{DT}, \code{NULL} is returned.
}
\description{
Returns the label of a cell in a division tree,
given the label of an instant of the cell in the corresponding lineage tree.
}
\seealso{
\code{\link{get_LT_cell}} for the reverse.
}
