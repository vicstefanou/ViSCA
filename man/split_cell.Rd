% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_cell.R
\name{split_cell}
\alias{split_cell}
\title{Split a cell}
\usage{
split_cell(LT, cell, Nsplit = 2, cell_list, col_list, Ncols, pixelR,
  matFolder = getwd(), matFileName, exeFolder, mcrFolder, show = TRUE)
}
\arguments{
\item{LT}{The lineage tree where the cell specified in \code{cell} belongs, an object of class \code{"igraph"}.}

\item{cell}{The label of the cell in the \code{LT} to be splitted, a character string.
It can be any non-root cell, as returned from \code{\link{get_cells}}.}

\item{Nsplit}{Number of cells for the \code{cell} to be splitted, an integer value \code{>=2}.
The default value is \code{2}.}

\item{cell_list}{A list containing all the cell instants of the movie.}

\item{col_list}{A list containing all the colony instants of the movie.}

\item{Ncols}{Number of colonies in the movie, a non-zero positive integer value.}

\item{pixelR}{The pixel ratio in units of length, a non-zero positive numeric value.}

\item{matFolder}{A character string naming the absolute path of the directory where
the \code{.mat} file generated by \emph{BaSCA} is saved (excluding the last \code{"/"}).
The default value is the current working directory \code{getwd()}.
\cr\cr
NOTE: The components should be separated by \code{"/"} on Windows.}

\item{matFileName}{A character string naming the \code{.mat} file generated by \emph{BaSCA}
(including the suffix \code{".mat"}).
The filename is relative to the \code{matFolder}.}

\item{exeFolder}{A character string naming the absolute path of the installation folder of the MATLAB executable
(excluding the last \code{"/"}).
\cr\cr
NOTE: The components should be separated by \code{"/"} on Windows.}

\item{mcrFolder}{A character string naming the absolute path of the installation folder of the Matlab Compiler Runtime (MCR)
(excluding the last \code{"/"}).
\cr\cr
NOTE: The components should be separated by \code{"/"} on Windows.}

\item{show}{A logical value (\code{TRUE} or \code{FALSE}) indicating whether \code{\link{view_cell}}
will be called for the mother of the \code{cell} (if it exists), the \code{cell} before the split operation
and the resulting cells after the split operation.
This capability is useful in order to see the result of the function.
The default value is \code{TRUE}.}
}
\value{
A named list with the following components:
\item{LT}{The updated LT with the \code{cell} replaced, an object of class \code{"igraph"}.}
\item{cell_list}{The updated cell_list with the \code{cell} replaced and the rest resulting cells added.}
\item{branches}{A list with the single-node (connected) lineage trees
for each one of the rest resulting cells
Each branch (element of the list) is an object of class \code{"igraph"}.}
}
\description{
Splits a cell in a lineage tree.
}
\details{
There is no guarantee that \code{cell} will be indeed splitted into \code{Nsplit} cells.
After the split operation, the resulting cells are sorted by ascending distance
between their centroid and the centroid of their mother, if it exists.
Otherwise, the order of the resulting cells is random.
The original \code{cell} is replaced in both \code{LT} and \code{cell_list} by the first listed cell, which implies best matching.
For each one of the remaining resulting cells a single-node (connected) lineage tree (motherless branch) is automatically created.
}
\section{Prerequisites}{

This function can be used by \emph{BaSCA} users \bold{only},
importing the data with \code{\link{import_basca}}.
}

\seealso{
\code{\link{merge_cells}} for the reverse,
\code{\link{add_branch}} for connecting a motherless branch to a lineage tree.
}
