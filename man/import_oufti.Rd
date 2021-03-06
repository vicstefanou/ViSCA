% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_oufti.R
\name{import_oufti}
\alias{import_oufti}
\title{Import \emph{Oufti} data}
\usage{
import_oufti(file, pixelR)
}
\arguments{
\item{file}{A character string naming the \code{.csv} file generated by \emph{Oufti}
(including the suffix \code{".csv"}) from which the data is to be imported.
If it does not contain an absolute path, the file name is relative to the current working directory, \code{getwd()}.
\cr\cr
NOTE: The components should be separated by \code{"/"} on Windows.
\cr\cr
Before using this function, the user has to open the \code{.csv} file in a text editor and
manually delete or comment (by adding the character \code{"\%"} at the beginning)
all lines of the file until the line with the \code{"\% parameter values"} string.}

\item{pixelR}{The pixel ratio in units of length, a non-zero positive numeric value.}
}
\value{
A named list with the following components:
\item{col_list}{\code{NULL}}
\item{cell_list}{A list containing all the cell instants of the movie.
Each element of the list is a named list with the following components:
\itemize{
\item \code{cellName} is the name of the cell, a character string in the format \code{"c<cellId>_f<frame>"}
\item \code{frame} is the ID of the frame of the cell, a non-zero positive integer number
\item \code{colony} is the ID of the colony of the cell in the \code{frame}, equal to \code{1}
\item \code{daughterIds} is a vector of character strings containing the \code{cellName}
of the linked cell(s) in the next frame,
or \code{NULL} in case no such cells exist
\item \code{length} is the length of the cell in units of length, a non-zero positive numeric value
\item \code{area} is the area of the cell in squared units of length, a non-zero positive numeric value
\item \code{polarity} is a numeric value defining the polarity of the cell
}
}
\item{Nframes}{Number of frames in the movie, a non-zero positive integer value.
IDs of frames are in the range \code{[1, Nframes]}.}
\item{Ncols}{Number of colonies in the movie, equal to \code{1}.
IDs of colonies are in the range \code{[1, Ncols]}.}
\item{frameH}{\code{NULL}}
\item{frameW}{\code{NULL}}
}
\description{
Imports a \code{.csv} file generated by \emph{Oufti} and converts it into
a cell list containing all the cell instants of the movie.
}
\details{
\emph{Oufti} treats the whole frame as a single colony.
\cr\cr
\emph{Oufti} does not export any information about either the pixels of each cell or the pixels each colony.
\cr\cr
\emph{Oufti} users importing the data with this function are \bold{excluded} from using
\code{\link{plot_col_tree}},
\code{\link{split_cell}}, \code{\link{get_cand_merge_cells}}, \code{\link{merge_cells}},
\code{\link{get_cand_mother_cells}}, \code{\link{add_branch}},
\code{\link{create_movie}}, \code{\link{create_cell_life}} and \code{\link{view_cell}}.
}
\references{
\url{http://oufti.org/}
}
