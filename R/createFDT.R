#' Create cells' FDT
#'
#' Creates the cells' Forest of Division Trees (FDT) given the corresponding Forest of Lineage Trees (FLT).
#'
#' A continuous segment (sequence) of LT nodes between two successive cell \emph{divisions} represents the lifespan of a cell.
#' A cell is \emph{divided} if
#' \itemize{
#' \item it gives birth to 2 daughter cells in the next frame,
#' \item it disappears from the field of view in the next frame or
#' \item it is not linked to any cell instant in the next frame due to tracking errors.
#' }
#' The fuction creates the Forest of Division Tree (FDT) by reducing all LT cell segments down to single nodes (DT nodes),
#' except for the imaginary \emph{root} cells.
#' \cr\cr
#' Each node of the FDT represents a cell at its full lifespan, having as character string attributes
#' the concatenation of all attributes of the \emph{collapsed} LT nodes by \code{", "},
#' with the following exceptions:
#' \itemize{
#' \item The \code{"colony"} attribute is not concatenated,
#' since the ID of the colony from which the cell emanated characterizes both cell and its instants.
#' The \code{"colony"} attribute is again a non-zero positive integer number.
#' \item Attributes \code{"age"} and/or \code{"generation"} are also not concatenated in case they exist,
#' since they are (re)evaluated by the function and updated/stored in the FLT.
#' See below for more details.
#' \item Each boolean attribute in the FLT forms also a boolean attribute in the FDT,
#' with value (\code{TRUE} or \code{FALSE}) based on the majority vote of the corresponding values
#' of the \emph{collapsed} LT nodes.
#' }
#' The concatenated \code{"name"} attribute represents the labels of the \emph{collapsed} cells (instants of the cell)
#' and is renamed to \code{"cellInstants"}.
#' The \code{"name"} attribute is again a non-zero positive integer number stored as a character string,
#' denoting the label of the cell in the FDT.
#' Value \code{"1"} corresponds to the main \emph{root} cell.
#' Values \code{"1+<i>"} correspond to the colonies' \emph{root} cells, where \code{"<i>"} is the colony ID.
#' The rest values correspond to the cells.
#' \cr\cr
#' For each numeric attribute in the FLT
#' (except for \code{"colony"} and \code{"frame"}, plus \code{"age"} and/or \code{"generation"} in case they exist),
#' the concatenation represents the cell's time-series of the attribute.
#' Given each cell's time-series of an attribute \code{"<attr>"},
#' the following numeric life attributes are estimated and stored as attributes in the corresponding FDT node:
#' \itemize{
#' \item \code{"<attr>_birth"} is the \code{"<attr>"} value of the first instant of the cell
#' \item \code{"<attr>_division"} is the \code{"<attr>"} value of the last instant of the cell
#' \item \code{"<attr>_mean"} is the \emph{mean} of \code{"<attr>"}
#' \item \code{"<attr>_sd"} is the \emph{standard deviation} of \code{"<attr>"},
#' or \code{NA} in case the cell has only one instant
#' \item \code{"<attr>_min"} is the minimum value of \code{"<attr>"}
#' \item \code{"<attr>_max"} is the maximum value of \code{"<attr>"}
#' }
#' These numeric attributes are in units of \code{"<attr>"}.
#' \cr\cr
#' The following attributes are also life attributes and are stored in each DT node:
#' \itemize{
#' \item \code{"generation"} is the ID of the generation of the cell, a positive integer value.
#' This value is also updated/stored in the \code{"generation"} attribute of the instants of the cell in the FLT,
#' since the ID of the generation characterizes both cell and its instants.
#' \item \code{"birthTime"} is the ID of the frame at which the cell is firstly spotted (born),
#' a non-zero positive integer value.
#' This value is the \code{"frame"} value of the first instant of the cell.
#' \item \code{"divisionTime"} is the ID of the frame at which the cell is lastly spotted (a frame before its \emph{division}),
#' a non-zero positive integer value.
#' This value is the \code{"frame"} value of the last instant of the cell.
#' \item \code{"lifeFrames"} is the duration of the cell life in \emph{frames}, a non-zero positive integer value.
#' This value is computed as
#' \code{\ifelse{html}{\out{lifeFrames = divisionTime - birthTime + 1}}{\deqn{lifeFrames = divisionTime - birthTime + 1}}}.
#' \item \code{"lifeHours"} is the duration of the cell life in \emph{hours}, a non-zero positive numeric value.
#' This value is computed as
#' \code{\ifelse{html}{\out{lifeHours = lifeFrames / (60 * frameR)}}{\deqn{lifeHours = \frac{lifeFrames}{60 \cdot frameR}}}}.
#' \item \code{"isConsidered"} is a logical value (\code{TRUE} or \code{FALSE})
#' indicating whether the cell will be included in the analysis or not.
#' This value is \code{FALSE} for the imaginary \emph{root} cells.
#' The value for the cells is computed based on the \code{minLife} argument.
#' It is \code{FALSE} for all leaf cells and \code{TRUE} for the rest, in case \code{minLife = 0}, or
#' \code{FALSE} for all leaf cells with \code{"lifeFrames" <= minLife} and \code{TRUE} for the rest, in case \code{minLife != 0}.
#' }
#'
#'
#' @param LTmain The main part of the overall FLT,
#' a connected lineage tree containing the imaginary \emph{root} cells
#' (object of class \code{"igraph"}).
#'
#' @param minLife Minimum life in \emph{frames} for a cell to be included in the analysis, a positive integer value.
#' The default value is \code{5}. Use value \code{0} to exclude all leaf cells from the analysis.
#'
#' @param frameR Frame rate of the movie in \emph{frames} per \emph{minute}, a non-zero positive numeric value.
#'
#' @return A named list with the following components:
#' \item{DTmain}{The corresponding main part of the overall FDT,
#' a connected division tree containing the imaginary \emph{root} cells (object of class \code{"igraph"}).}
#' \item{LTmain}{The updated LTmain with the attributes \code{"generation"} and \code{"age"} updated/added,
#' an object of class \code{"igraph"}.
#' The \code{"age"} attribute denotes the age of each cell instant in \emph{frames}.}
#' \item{Ngens}{Number of generations in the movie, a non-zero positive integer value.
#' IDs of generations are in the range \code{[0, Ngens-1]}.}
#'
#' @seealso \code{\link{isConnected}} for checking if a tree is connected,
#' \code{\link{save_tree}} for saving a tree on disc.
#' @export
#' @import igraph
#' @importFrom stats sd
#' @importFrom utils setTxtProgressBar tail txtProgressBar


createFDT <- function(LTmain, minLife = 5, frameR) {

  ################# arguments check ##########

  if (!("1" %in% V(LTmain)$name)) {
    stop("Root is missing from LTmain\n")
  }

  if (!isConnected(tree = LTmain)) {
    stop("LTmain is disconnected\n")
  }

  ################################### find cell instants to "collapse"

  roots <- V(LTmain)[!(V(LTmain)$name %in% get_cells(tree = LTmain, treeT = "LT", type = "nr"))]$name

  V(LTmain)$DTid <- c(1:gorder(LTmain))

  V(LTmain)$degree <- degree(graph = LTmain, v = V(LTmain), mode = "out")
  V(LTmain)[roots]$degree <- -1   # don't "scan" roots

  n <- V(LTmain)[V(LTmain)$degree != 1]$name
  n_d2 <- V(LTmain)[V(LTmain)$degree >= 2]$name

  cat("Creating FDT...\n")
  pb <- txtProgressBar(min = 0, max = length(n), style = 3) ### set progress bar
  ipb <- 0

  for (n_i in n) {

    ipb <- ipb + 1
    setTxtProgressBar(pb, ipb) ### update progress bar

    n_i_ancestors <- subcomponent(graph = LTmain, v = n_i, mode = "in")$name

    n_i_ancestors <- setdiff(n_i_ancestors, roots) # don't merge cells to/with roots

    if (length(n_i_ancestors[-1]) != 0) { # n_i_ancestors include n_i

      possible_mothers <- intersect(n_i_ancestors, n_d2)

      if (V(LTmain)[n_i]$degree == 2) { # possible_mothers include n_i
        if (length(possible_mothers) >= 2) { # merge to mother
          n_i_tobemerged <- match(possible_mothers[2], n_i_ancestors) - 1
        } else {  # merge to cell instant
          n_i_tobemerged <- length(n_i_ancestors)
        }
      } else { # V(LTmain)[n_i]$degree == 0 # possible_mothers don't include n_i
        if (length(possible_mothers) >= 1) { # merge to mother
          n_i_tobemerged <- match(possible_mothers[1], n_i_ancestors) - 1
        } else { # merge to cell instant
          n_i_tobemerged <- length(n_i_ancestors)
        }
      }

      V(LTmain)[V(LTmain)$name %in% n_i_ancestors[1:n_i_tobemerged]]$DTid <- V(LTmain)[n_i_ancestors[n_i_tobemerged]]$DTid

    } # else "no ancestors"
  }

  close(pb) ### close progress bar
  cat("\n")

  ########################## remove from LT

  # generation
  if ("generation" %in% vertex_attr_names(LTmain)) {
    LTmain <- delete_vertex_attr(graph = LTmain, name = "generation")
  }

  # age
  if ("age" %in% vertex_attr_names(LTmain)) {
    LTmain <- delete_vertex_attr(graph = LTmain, name = "age")
  }

  ######################################### construct DT

  #  collapse all instant attributes separated by ", "
  DTmain <- contract(graph = LTmain, mapping = V(LTmain)$DTid, vertex.attr.comb = toString)
  DTmain <- delete_vertices(graph = DTmain, v = V(DTmain)[V(DTmain)$name == ""])
  DTmain <- simplify(DTmain)

  ######################## correct attributes in DT
  cat("Adding attributes to FDT...\n")

  # cellIinstants
  V(DTmain)$cellInstants <- V(DTmain)$name
  V(DTmain)$name <- as.character(c(1:gorder(DTmain)))
  DTmain <- delete_vertex_attr(graph = DTmain, name = "DTid")

  # colony
  values <- as.numeric(unname(sapply(V(DTmain)$colony, function(x) unlist(strsplit(x, ", "))[1] )))
  DTmain <- delete_vertex_attr(graph = DTmain, name = "colony")
  V(DTmain)$colony <- values

  ######################################### extra attributes in DT

  # generation
  V(DTmain)$generation <- unname(sapply(V(DTmain)$name,
                                        function(x) length(subcomponent(graph = DTmain, v = x, mode = "in")$name))) - 2

  # lifeFrames -> duration of life in frames (frames in which the cell appears)
  V(DTmain)$lifeFrames <- unname(sapply(V(DTmain)$frame, function(x) length(unlist(strsplit(x, ", ")))))

  # lifeHours
  frameR <- frameR * 60 # in frames/hour
  V(DTmain)$lifeHours <- V(DTmain)$lifeFrames / frameR # duration of life in hours

  # birthTime (frame)
  V(DTmain)$birthTime <- as.numeric(unname(sapply(V(DTmain)$frame,
                                                  function(x) unlist(strsplit(x, ", "))[1])))

  # divisionTime (frame)
  V(DTmain)$divisionTime <- as.numeric(unname(sapply(V(DTmain)$frame,
                                                     function(x) tail(unlist(strsplit(x, ", ")), 1))))

  # isConsidered
  V(DTmain)$isConsidered <- TRUE
  V(DTmain)[!(V(DTmain)$name %in% get_cells(tree = DTmain, treeT = "DT", type = "nr"))]$isConsidered <- FALSE # root cells
  V(DTmain)$cellLastDegree <- as.numeric(unname(sapply(V(DTmain)$degree,
                                                       function(x) tail(unlist(strsplit(x, ", ")), n = 1))))
  if (minLife == 0) {
    V(DTmain)[V(DTmain)$cellLastDegree == 0]$isConsidered <- FALSE # all leaves
  } else {
    V(DTmain)[V(DTmain)$cellLastDegree == 0 &
                V(DTmain)$lifeFrames < minLife]$isConsidered <- FALSE # all leaves with "short life"
  }
  DTmain <- delete_vertex_attr(graph = DTmain, name = "degree")
  DTmain <- delete_vertex_attr(graph = DTmain, name = "cellLastDegree")

  ################################# numeric and boolean attributes of LT in DT

  numeric_attrs <- get_attr_names(tree = LTmain, type = "n")
  boolean_attrs <- get_attr_names(tree = LTmain, type = "b")

  for (attr in vertex_attr_names(DTmain)) {

    if (attr %in% numeric_attrs) { # statistics

      if (!(attr %in% c("frame", "colony"))) { # possible "generation" and "age" where deleted

        ######### replace NA in first cell of a possible ROC attribute with 0 (for future fit)

        values_all <- unname(sapply(vertex_attr(graph = DTmain, name = attr, index = V(DTmain)),
                                    function(x) unlist(strsplit(x, ", "))))

        values_all_new <- sapply(values_all, function(x) {
          if (x[1] == "NA") {
            if (length(x) == 1) {
              x <- 0
            } else { # for ROC
              x <- c(0, x[2:length(x)])
            }
          } else { # do nothing
            x <- x
          }
        })

        DTmain <- set_vertex_attr(graph = DTmain, name = attr, index = V(DTmain),
                                  value = sapply(values_all_new, function(x) toString(x)))

        ###### exclude NA in first cell of a possible ROC attribute

        values_all <- sapply(values_all, function(x) {
          if (x[1] == "NA") {
            if (length(x) == 1) {
              x <- NA # -------------> all statistics will be NA
            } else {
              x <- as.numeric(x[2:length(x)])
            }
          } else {
            x <- as.numeric(x)
          }
        })

        ####################### compute statistics

        values <- sapply(values_all, function(x) x[1])
        DTmain <- set_vertex_attr(graph = DTmain, name = paste(attr, "birth", sep = "_"), index = V(DTmain), value = values)

        values <- sapply(values_all, function(x) tail(x, 1))
        DTmain <- set_vertex_attr(graph = DTmain, name = paste(attr, "division", sep = "_"), index = V(DTmain), value = values)

        values <- sapply(values_all, function(x) mean(x))
        DTmain <- set_vertex_attr(graph = DTmain, name = paste(attr, "mean", sep = "_"), index = V(DTmain), value = values)

        values <- sapply(values_all, function(x) sd(x)) # NA if 1 value or NA
        DTmain <- set_vertex_attr(graph = DTmain, name = paste(attr, "sd", sep = "_"), index = V(DTmain), value = values)

        values <- sapply(values_all, function(x) min(x))
        DTmain <- set_vertex_attr(graph = DTmain, name = paste(attr, "min", sep = "_"), index = V(DTmain), value = values)

        values <- sapply(values_all, function(x) max(x))
        DTmain <- set_vertex_attr(graph = DTmain, name = paste(attr, "max", sep = "_"), index = V(DTmain), value = values)

      }

    } else if (attr %in% boolean_attrs) { # majority vote

      values <- as.logical(unname(sapply(vertex_attr(graph = DTmain, name = attr, index = V(DTmain)),
                                         function(x) names(which.max(table(unlist(strsplit(x, ", "))))))))
      DTmain <- delete_vertex_attr(graph = DTmain, name = attr)
      DTmain <- set_vertex_attr(graph = DTmain, name = attr, index = V(DTmain), value = values)

    }
  }

  ####################### update LT
  cat("Updating LTmain...\n")

  # generation
  for (i in unique(V(DTmain)$generation)) {
    cells <- unlist(unname(sapply(V(DTmain)[V(DTmain)$generation == i]$cellInstants,
                                  function(x) strsplit(x, split = ", "))))
    V(LTmain)[cells]$generation <- i
  }

  # age
  a <- unname(sapply(V(DTmain)$cellInstants, function(x) unlist(strsplit(x, split = ", "))))
  for (i_cell in 1:length(a)) {
    V(LTmain)[a[[i_cell]]]$age <- c(1:length(a[[i_cell]]))
  }

  LTmain <- delete_vertex_attr(graph = LTmain, name = "DTid")
  LTmain <- delete_vertex_attr(graph = LTmain, name = "degree")

  ##########################

  return(list(DTmain = DTmain, LTmain = LTmain, Ngens = max(V(DTmain)$generation) + 1)) # generations start from 0

}
