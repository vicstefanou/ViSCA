% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ViSCA-package.R
\docType{package}
\name{ViSCA}
\alias{ViSCA}
\alias{ViSCA-package}
\title{Overview of the ViSCA package}
\description{
ViSCA is an R package for the statistical analysis and visualization of
single-cell data derived from the analysis of time-lapse cell-movies.
The package enables users to explore the spatiotemporal trends of single-cell attributes,
discover possible epigenetic effects across generations and
even identify and correct errors.
}
\details{
ViSCA was initially developed to extend the \emph{BaSCA} pipeline
with analytics, visualization and error correction capabilities.
However, most of these capabilities are general and independent.
Therefore, the package can be used for data derived from any other software,
such as \emph{SuperSegger} and \emph{Oufti}.
}
\section{Data import}{

ViSCA supports various input file formats for the analyzed cell-movies.
\cr\cr
Use \code{\link{import_basca}}, \code{\link{import_oufti}} or \code{\link{import_ss}}
to import data exported by the named software.
These functions automatically convert the input file(s) into a cell list (and colony list),
containing all the cell instants (and colony instants) of the movie, respectively.
\cr\cr
Such structures can be also directly imported from custom-made \code{.json} files,
provided that they have the appropriate format.
See \code{\link{import_json}} for details.
}

\section{FLT/FDT representation of cell-movies}{

Once the cell list of the movie is loaded,
use \code{\link{createFLT}} to transform it into a Forest of Lineage Trees (FLT) data structure.
A lineage tree (LT) node represents a cell instant (i.e. cell at a specific frame/instant of its lifespan).
A continuous segment (sequence) of LT nodes between two successive cell divisions represents the lifespan of a cell.
\cr\cr
If one reduces LT cell segments down to a single node,
he obtains the Forest of Division Trees (FDT) data structure,
capturing only each cell’s division event and summarizing its lifespan.
A division tree (DT) node represents a cell (i.e. cell at its full lifespan).
See \code{\link{createFDT}} for details.
\cr\cr
These tree data structures are objects of class \code{"igraph"} and are the core structures of the package.
}

\section{Attributes}{

Single-cell attributes are divided into two broad categories:
\itemize{
\item cell instant attributes, that may change value at each time point.
These attributes are extracted by the software, loaded into the cell list and
finally stored as node attributes in the FLT by \code{\link{createFLT}}.
\cr\cr
Some other values as also stored as node attributes in the FLT when
\code{\link{add_attr_roc}} or \code{\link{createFDT}} are called.
See the documentation of each function for more details.
\item cell life attributes characterize a cell’s whole lifespan.
These attributes are estimated by \code{\link{createFDT}} and are stored as node attributes in the FDT.
\cr\cr
Some other values as also stored as node attributes in the FDT
when \code{\link{add_attr_growth_fit_pars}} is called.
See the documentation of the function for more details.
}
}

\section{Analytics}{

ViSCA allows users to perform statistical analysis of single-cell attributes at multiple levels of organization
(whole community, sub-population, colonies, generations, subtrees of individual colonies, etc.).
Analytics capabilities are categorized into:
\itemize{
\item statistics (\code{\link{get_attr_stats}}, \code{\link{plot_hist_attr}}, \code{\link{plot_pdf_attr}},
\code{\link{plot_viobox_attr}}, \code{\link{plot_time_attr}})
\item scatterplots for correlating attributes (\code{\link{plot_dot_attr2}}, \code{\link{plot_dot_attr3}},
\code{\link{plot_dot_time_attr}}, \code{\link{plot_dot_attr2_gen2}}, \code{\link{plot_dot_attr_fam}})
\item estimation of growth curves (\code{\link{plot_baranyi}}, \code{\link{add_attr_growth_fit_pars}},
\code{\link{plot_growth_attr_fit}}, \code{\link{plot_growth_attr}}, \code{\link{plot_growth_attr_cell}})
}
}

\section{Visualization}{

ViSCA provides two different ways for visualization:
\enumerate{
\item \code{\link{plot_tree}} for visualizing a lineage or generation tree
\item \code{\link{create_movie}} for animating the segmented cells by creating videos
}
Color can be used to map a cell instant/life attribute and
capture its variability across cells, colonies, frames or generations.
\cr\cr
The user can also monitor how the life of a selected cell evolves in the movie using \code{\link{create_cell_life}}.
}

\section{Error correction}{

ViSCA allows users to correct
\itemize{
\item tracking errors with \code{\link{extract_branch}}, \code{\link{get_cand_mother_cells}} and \code{\link{add_branch}}.
\item segmantation errors with \code{\link{split_cell}}, \code{\link{get_cand_merge_cells}} and \code{\link{merge_cells}}.
This capability is offered to \emph{BaSCA} users \bold{only}.
}
}

\section{Notes}{

Some functions have prerequisites in order to be used.
See the \emph{Prerequisites} field of each function for more details.
}

\references{
A. Balomenos, P. Tsakanikas, Z. Aspridou, A. Tampakaki, K. Koutsoumanis and E. Manolakos,
\emph{“Image analysis driven single-cell analytics for systems microbiology”},
BMC Systems Biology, vol. 11, no. 1, 2017.
\cr\cr
\url{http://oufti.org/}
\cr\cr
\url{http://mtshasta.phys.washington.edu/website/SuperSegger.php}
}
