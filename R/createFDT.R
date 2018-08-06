#' Create cells' FDT
#'
#' Creates the cells' Forest of Division Trees (FDT) given the corresponding Forest of Lineage Trees (FLT).
#'
#' @param LTmain The main part of the overall FLT,
#' a connected lineage tree containing the imaginary \emph{root} cells
#' (object of class \code{"igraph"}).
#'
#' @param minLife Minimum life in \emph{frames} for a cell to be considered in the analysis, a positive integer value.
#' The default value is \code{5}. Use value \code{0} to exclude all leaf cells from the analysis.
#'
#' @param frameR Frame rate of the movie in \emph{frames} per \emph{minute}, a non-zero positive numeric value.
#'
#' @return A named list with the following components:
#' \item{DTmain}{The corresponding main part of the overall FDT,
#' a connected division tree containing the imaginary \emph{root} cells (object of class \code{"igraph"}).}
#' \item{LTmain}{The updated LTmain with the attributes \code{"generation"} and \code{"age"} added/updated,
#' an object of class \code{"igraph"}.}
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

        values <- sapply(values_all, function(x) min(x))
        DTmain <- set_vertex_attr(graph = DTmain, name = paste(attr, "min", sep = "_"), index = V(DTmain), value = values)

        values <- sapply(values_all, function(x) max(x))
        DTmain <- set_vertex_attr(graph = DTmain, name = paste(attr, "max", sep = "_"), index = V(DTmain), value = values)

        values <- sapply(values_all, function(x) mean(x))
        DTmain <- set_vertex_attr(graph = DTmain, name = paste(attr, "mean", sep = "_"), index = V(DTmain), value = values)

        values <- sapply(values_all, function(x) sd(x)) # NA if 1 value or NA
        DTmain <- set_vertex_attr(graph = DTmain, name = paste(attr, "sd", sep = "_"), index = V(DTmain), value = values)

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
