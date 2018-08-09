createAxisLab <- function(attr, unit = "", string = "") {

  attr <- unlist(strsplit(attr, "_"))

  if (length(attr) == 1) {
    attr <- paste(toupper(substr(attr, 1, 1)), substr(attr, 2, nchar(attr)), sep = "")
  } else if (length(attr) == 2) {
    if (nchar(attr[2]) == 1) {
      attr <- paste(attr[2], " ", toupper(substr(attr[1], 1, 1)), substr(attr[1], 2, nchar(attr[1])), sep="")
    } else {
      attr <- paste(toupper(substr(attr[2], 1, 1)), substr(attr[2], 2, nchar(attr[2])), " ",
                    toupper(substr(attr[1], 1, 1)), substr(attr[1], 2, nchar(attr[1])), sep = "")
    }
  } else {
    warning("attr must contain up to one \"_\"\n") ## => LT attributes must not contain "_"
    attr <- paste(attr, collapse = "_")
  }

  unit <- unlist(strsplit(unit, ","))

  if (length(unit) == 0) {
    lab <- bquote(bold(.(attr)))
  } else if (length(unit) == 1) {
    lab <- bquote(bold(.(attr) ~ ( .(unit) )))
  } else if (length(unit) == 2) {
    lab <- bquote(bold(.(attr) ~ ( .(unit[1]) ^ .(unit[2]) )))
  } else {
    warning("unit must be in format <string,number>\n")
    unit <- paste(unit, collapse = ",")
    lab <- bquote(bold(.(attr) ~ ( .(unit) )))
  }

  if (string != "") {
    lab <- bquote(bold(.(lab) ~ .(string)))
  }

  return(lab)

}


##################################################################################################

#' @importFrom grDevices rainbow

getGroupColors <- function(type = "", N_colors = NULL) {

  if (type == "colony") {
    groupColors <- rainbow(N_colors, s = 0.8)
  } else if (type == "generation") {
    groupColors <- rainbow(N_colors, s = 0.7, v = 0.8)
  } else { # population
    groupColors <- "slateblue4"
  }

  return(groupColors)

}

##################################################################################################

#' @import igraph
#' @importFrom grDevices colorRampPalette

getCellColors <- function(tree, cells, attr, N_colors = NULL) {

  if (length(cells) != 0) {

    values <- vertex_attr(graph = tree, name = attr, index = cells)

    if (attr == "colony") {

      paletteColors <- getGroupColors(type = "colony", N_colors = N_colors)
      cellColors <- paletteColors[values]

    } else if (attr == "generation") {

      paletteColors <- getGroupColors(type = "generation", N_colors = N_colors)
      cellColors <- paletteColors[values + 1]

    } else {

      palette <- colorRampPalette(c("skyblue", "slateblue4", "darkred"))

      if (is.logical(values)) { # boolean

        paletteColors <- c("darkred", "skyblue")
        cellColors <- paletteColors[as.numeric(values) + 1] # FALSE -> 1, TRUE -> 2

      } else { # numeric

        if (length(unique(values)) >= 2) {

          if (length(values[values == round(values)]) == length(values)) { # integer
            N_colors <- length(c(min(values):max(values)))
          } else {
            N_colors <- 30
          }

          paletteColors <- palette(N_colors)
          cellColors <- paletteColors[as.numeric(cut(values, breaks = N_colors))] # >= 2 unique values

        } else {
          warning("Less than 2 unique values for coloring\n")
          cellColors <- NULL
        }

      }

    }
  } else {
    warning("No cells for coloring\n")
    cellColors <- NULL
  }

  return(cellColors)

}


##################################################################################################

#' @importFrom stats na.omit

getCellColorsDot <- function(df, myPlot, tree, attr, unit = "", N_colors = NULL) {

  df <- na.omit(df) # remove na colors to find colors with the valid values
  values <- df$color
  colored <- TRUE

  if (length(values) != 0) {

    if (attr %in% c("colony", "generation") || attr %in% get_attr_names(tree = tree, type = "b")) {

      if (attr == "colony") {
        graphColors <- getGroupColors(type = attr, N_colors = N_colors)
        myColors <- graphColors[as.numeric(levels(values))]
      } else if (attr == "generation") {
        graphColors <- getGroupColors(type = attr, N_colors = N_colors)
        myColors <- graphColors[as.numeric(levels(values)) + 1]
      } else { # boolean
        myColors <- c("darkred", "skyblue") # FALSE, TRUE
      }

      myPlot <- myPlot + scale_color_manual(values = myColors)

    } else { # numeric

      if (length(unique(values)) >= 2) {

        if (length(values[values == round(values)]) == length(values)) { # integer -> density here

          myPlot <- myPlot + scale_color_gradientn(colors = c("skyblue", "slateblue4", "darkred"),
                                                   breaks = c(min(values), max(values)),
                                                   labels = c(min(values), max(values)),
                                                   na.value = "gray30")
        } else {
          myPlot <- myPlot + scale_color_gradientn(colors = c("skyblue", "slateblue4", "darkred"),
                                                   breaks = c(min(values), max(values)),
                                                   # 3 decimal digits precision
                                                   labels = c(format(round(min(values), 3), nsmall = 3),
                                                              format(round(max(values), 3), nsmall = 3)),
                                                   na.value = "gray30")
        }

        if (attr == "") { # density
          myPlot <- myPlot + labs(color = "Counts")
        } else {
          myPlot <- myPlot + labs(color = createAxisLab(attr = attr, unit = unit))
        }

      } else {
        warning("Less than 2 unique values for coloring\n")
        colored <- FALSE
      }

    }
  } else {
    warning("No cells for coloring\n")
    colored <- FALSE
  }

  return(list(myPlot = myPlot, colored = colored))

}

##################################################################################################

#' @importFrom graphics legend par

plotColorLegend <- function(attr, whichColors = NULL, N_colors, size_lab) {

  if (attr == "colony") {
    title <- expression(bold("Colony"))
    if (is.null(whichColors)) {
      legend <- c(1:N_colors)
      colors <-  getGroupColors(type = "colony", N_colors = N_colors)
    } else {
      legend <- whichColors
      colors <-  getGroupColors(type = "colony", N_colors = N_colors)[whichColors]
    }
  } else if (attr == "generation") {
    title <- expression(bold("Generation"))
    if (is.null(whichColors)) {
      legend <- c(0:(N_colors-1))
      colors <-  getGroupColors(type = "generation", N_colors = N_colors)
    } else {
      legend <- whichColors
      colors <-  getGroupColors(type = "generation", N_colors = N_colors)[whichColors + 1]
    }
  } else { # boolean attr
    title <- as.expression(bquote(bold(.(attr))))
    legend <- c("FALSE", "TRUE")
    colors <- c("darkred", "skyblue")
  }

  corners <- par("usr")

  legend(x = corners[2], y = mean(corners[3:4]),
         cex = size_lab, title = title, legend = legend, lwd = 3, col = colors,
         yjust = 0.5, bty = "n", xpd = TRUE)

}

##################################################################################################

#' @importFrom graphics mtext par plot text
#' @importFrom grDevices colorRampPalette

plotColormap <- function(values, attr, unit, size_lab, size_values, dens = FALSE) {

  if (dens == FALSE) {
    palette <- colorRampPalette(c("skyblue", "slateblue4", "darkred"))
  } else {
    palette <- colorRampPalette(c("darkred", "yellow"))
  }

  up_space <- 10
  down_space <- 13
  par(mar = c(down_space, 0, up_space, 2))

  plot(x = rep(0, 10000), y = seq(0, 10, length = 10000), pch = 15,
       axes = FALSE, xlab = "", ylab = "",
       col = palette(10000)
  )

  if (length(values[values == round(values)]) == length(values)) { # values is integer
    mtext(max(values), side = 3, cex = size_values)
    mtext(min(values), side = 1, cex = size_values)
  } else {
    mtext(format(round(max(values), 3), nsmall = 3), side = 3, cex = size_values)
    mtext(format(round(min(values), 3), nsmall = 3), side = 1, cex = size_values)
  }

  corners <- par("usr")
  text(x = corners[2], y = mean(corners[3:4]), createAxisLab(attr = attr, unit = unit), srt = 270, xpd = TRUE, cex = size_lab)

}

#################################################################################################

updateCellInBascaList <- function(x, frame, colony, colId, col_list, pixelRatio, updateDaughterIds = FALSE) {

  x$frame <- frame
  x$colony <- colony
  x$colId <- colId # numeric

  x$pixelList <- cbind(x$pixelList[, 2], x$pixelList[, 1]) # (row,col) in colImage
  x$boundaryPixelList <- cbind(x$boundaryPixelList[, 2], x$boundaryPixelList[, 1]) # (row,col) in colImage

  x$centroid <- matrix(colMeans(x$pixelList), nrow = 1, ncol = 2) #cbind(x$centroid[, 2], x$centroid[, 1]) # (row,col) in colImage

  x$distFromCentroid <- norm(x$centroid - col_list[[colId]]$colCentroid, "f") * pixelRatio

  colBoundaryPixels <- col_list[[colId]]$colBoundaryPixels
  overlapVector <- paste(x$pixelList[, 1], x$pixelList[, 2], sep = " ") %in% paste(colBoundaryPixels[, 1], colBoundaryPixels[, 2], sep = " ")
  overlapPercent <- sum(overlapVector) / dim(x$pixelList)[1]
  if (overlapPercent > 0.66) {
    x$isOnBoundary <- TRUE
  } else {
    x$isOnBoundary <- FALSE
  }

  x$cellName <- paste(x$cellName, "_f", x$frame, sep = "")

  if (updateDaughterIds) {
    if (!is.null(x$daughterIds)) {
      if (length(x$daughterIds) == 0) {
        x$daughterIds <- NULL
      } else {
        x$daughterIds <- paste(x$daughterIds, "_f", x$frame + 1, sep = "")
      }
    }
  }

  x$length <- as.numeric(x$length) * pixelRatio
  x$width <- as.numeric(x$width) * pixelRatio
  x$LW <-  x$length / x$width
  x$area <- dim(x$pixelList)[1] * pixelRatio * pixelRatio
  x$perimeter <- dim(x$boundaryPixelList)[1] * pixelRatio

  x$lengthInPixels <- NULL
  x$widthInPixels <- NULL
  x$covariance <- NULL

  # values from matrix format
  for (attr in names(x)) {
    if (is.matrix(x[[attr]])) {
      if (dim(x[[attr]])[1] == 1 & dim(x[[attr]])[2] == 1) {
        x[[attr]] <- as.vector(x[[attr]])
      }
    }
  }

  return(x)
}

#################################################################################################

#' @import igraph

updateCellInBascaLT <- function(LT, cell, cell_list, Ncols, colony) {

  cell_list_ID <- as.numeric(cell) - (Ncols + 1)

  # cellName
  V(LT)[cell]$cellName <- cell_list[[cell_list_ID]]$cellName

  # numeric attributes
  numeric_attrs <- names(cell_list[[1]])[sapply(cell_list[[1]], function(x) class(x) == "numeric" || class(x) == "integer")]
  for (attr in numeric_attrs) {
    if (attr == "colId") {
      V(LT)[cell]$colId <- as.character(cell_list[[cell_list_ID]]$colId)
    } else if (attr == "colony") {
      V(LT)[cell]$colony <- colony # correct colony
    } else {
      LT <- set_vertex_attr(graph = LT, name = attr, index = V(LT)[cell]$name, value = unlist(unname(cell_list[[cell_list_ID]][attr])))
    }
  }

  # boolean attributes
  boolean_attrs <- names(cell_list[[1]])[sapply(cell_list[[1]], function(x) is.logical(x))]
  for (attr in boolean_attrs) {
    LT <- set_vertex_attr(graph = LT, name = attr, index = V(LT)[cell]$name, value = unlist(unname(cell_list[[cell_list_ID]][attr])))
  }

  return(LT)

}

#################################################################################################

callBascaMatlabExe <- function(Nsplit, cell_list_IDs,
                          cell_list, col_list,
                          mat_folder, mat_fileName,
                          install_folder_exe, install_folder_MCR) {

  # to run the shell script, type
  # ./run_mat_correct.sh <mcr_directory> <argument_list>
  # < argument_list> :
  #   mat_folder -> without ending "/"
  #   mat_fileName -> fileName of BaSCA.mat
  #   cellIDs -> cellName of cells (without suffix "_fi", i the frame), concatanated by ","
  #   frameID -> frame
  #   colID -> colony in frame
  #   Nsplit -> 0 for merge, >=2 for split

  exe_path <- paste(install_folder_exe, "application", sep = "/")
  mcr_directory <- paste(install_folder_MCR, system(paste("cd", install_folder_MCR, "&& ls | head -1"), intern = TRUE), sep = "/")

  cellIDs <- sapply(cell_list, function(x) x$cellName)[cell_list_IDs]
  cellIDs <- paste(unname(sapply(cellIDs, function(x) unlist(strsplit(x, split = "_f"))[1])), collapse = ",")
  a <- unlist(strsplit(col_list[[cell_list[[cell_list_IDs[1]]]$colId]]$colName, split = "_c"))
  frameID <- unlist(strsplit(a[1], split = "f"))[2]
  colID <- a[2]

  system(paste(paste(exe_path, "run_mat_correct.sh", sep = "/"),
               mcr_directory, mat_folder, mat_fileName, cellIDs, frameID, colID, Nsplit))

  newCellProps <- R.matlab::readMat(paste(mat_folder, "mat_correction.mat", sep = "/"))

  file.remove(paste(mat_folder, "mat_correction.mat", sep = "/"))

  return(newCellProps)

}

#################################################################################################

#' @import ggplot2

createBlankPlot <- function() {

  blankPlot <- ggplot() +
    geom_blank(aes(1,1)) +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
    )

  return(blankPlot)

}

#################################################################################################

getAttrForTitle <- function(attr) {

  attr <- unlist(strsplit(attr, "_"))
  if (length(attr) == 1) {
    attr <- paste(toupper(substr(attr, 1, 1)), substr(attr, 2, nchar(attr)), sep = "")
  } else if (length(attr) == 2) {
    if (nchar(attr[2]) == 1) {
      attr <- paste(attr[2], " ", toupper(substr(attr[1], 1, 1)), substr(attr[1], 2, nchar(attr[1])), sep = "")
    } else {
      attr <- paste(toupper(substr(attr[2], 1, 1)), substr(attr[2], 2, nchar(attr[2])), " ",
                    toupper(substr(attr[1], 1, 1)), substr(attr[1], 2, nchar(attr[1])), sep = "")
    }
  }

  return(attr)

}

#################################################################################################

#' @import ggplot2
#' @importFrom stats lm coef

addRegressionLine <- function(myPlot, df) {

  myModel <- lm(df$var2 ~ df$var1)
  a <- unname(coef(myModel)[2])
  if (!is.na(a)) { # y = a*x + b
    b <- unname(coef(myModel)[1])
    r2 <- summary(myModel)$r.squared
    myPlot <- myPlot + geom_abline(intercept = b, slope = a, color = "gray14", size = 0.7)
    regression <- list(a = a, b = b, r2 = r2)
  } else { # y = b
    regression <- NULL # less than 2 unique data points
  }

  return(list(myPlot = myPlot, regression = regression))

}

#################################################################################################

findNbins <- function(values, Nbins) {

  if (length(unique(values)) >= 2) { # min(values) != max(values)

    if (length(values[values == round(values)]) == length(values)) { # integer
      if (Nbins >= length(c(min(values):max(values)))) {
        Nbins_old <- Nbins
        Nbins <- length(c(min(values):max(values))) - 1 # hist
        Nbins <- length(c(min(values):max(values))) # geom_histogram
        warning(paste("Less than", Nbins_old, "unique integer values\nNbins was set to", Nbins))
      }
    }

  } else {
    warning("Less than two unique values\nNbins was set to 1")
    Nbins <- 1
  }

  return(Nbins)
}

#################################################################################################

checkField <- function(myList, fieldName) {

  m <- ""
  vals <- unname(sapply(myList, function(x) x[fieldName]))

  if (is.null(unlist(vals))) {
    m <- paste("Component \"", fieldName, "\" does not exist\n", sep = "")
  }

  if (fieldName != "daughterIds") {
    if (length(which(sapply(vals, function(x) is.null(x)) == TRUE)) != 0) {
      m <- paste("Missing values in component \"", fieldName, "\"\n", sep = "")
    }
  }

  if (!is.matrix(vals[[1]]) & fieldName != "daughterIds") { # single value
    if (length(which(sapply(vals, function(x) is.na(x)) == TRUE)) != 0) {
      m <- paste("NA values in component \"", fieldName, "\"\n", sep = "")
    }
  }

  return(m)

}

#################################################################################################

checkColList <- function(col_list, frameH, frameW) {

  # colName is the name of the colony instant, a character string in the format "f<frame>_c<colony>",
  # "<frame>" and "<colony>" is the ID of the frame and colony (in the frame) of the colony instant, respectively.

  if ((m <- checkField(myList = col_list, fieldName = "colName")) != "") {
    stop(paste("colony list:", m))
  } else {
    vals <- sapply(col_list, function(x) x$colName)
    if (length(grep(pattern = "^f\\d+_c\\d+$", x = vals, value=FALSE)) != length(vals)) {
      stop("colony list: Wrong values in component \"colName\"\n")
    }
  }

  # prev_colName is a vector of character strings containing the colName
  # of the corresponding colony instant(s) in the previous frame.
  # For colony instants that do not have a corresponding colony instant in the previous frame, this is equal to "f0_c0".

  if ((m <- checkField(myList = col_list, fieldName = "prev_colName")) != "") {
    stop(paste("colony list:", m))
  } else {
    vals <- unlist(unname(sapply(col_list, function(x) x$prev_colName)))
    if (length(setdiff(vals, c(unlist(unname(sapply(col_list, function(x) x$colName))), "f0_c0"))) != 0) {
      stop("colony list: Wrong values in component \"prev_colName\"\n")
    }
  }

  # next_colName is a character string containing the colName
  # of the corresponding colony instant in the next frame.
  # For colony instants that do not have a corresponding colony instant in the next frame, this is equal to "f00_c0".

  if ((m <- checkField(myList = col_list, fieldName = "next_colName")) != "") {
    stop(paste("colony list:", m))
  } else {
    vals <- sapply(col_list, function(x) x$next_colName)
    if (length(setdiff(vals, c(unlist(unname(sapply(col_list, function(x) x$colName))), "f00_c0"))) != 0) {
      stop("colony list: Wrong values in component \"next_colName\"\n")
    }
  }

  ####################### not necessary ##################################################################
  # colImage is the mask of the box surrounding the colony instant, a matrix of 0 and 1.
  # 1s denote the pixels of cells and 0s the background pixels.
  #
  # ULcorner is an 1x2 matrix of non-zero integer values
  # denoting the upper-left pixel of the box surrounding the colony instant in global (frame) coordinates.
  # The first integer represents the row and the second the column of the pixel.

  if ((m1 <- checkField(myList = col_list, fieldName = "colImage")) == "") {
    if ((m2 <- checkField(myList = col_list, fieldName = "ULcorner")) != "") {
      stop(paste("colony list: Component \"colImage\" exists\n", m2))
    } else {

      a <- sapply(col_list, function(x) {
        x$ULcorner[1, 1] >= 1 &
          x$ULcorner[1, 2] >= 1 &
          x$ULcorner[1, 1] + nrow(x$colImage) - 1 <= frameH &
          x$ULcorner[1, 2] + ncol(x$colImage) - 1 <= frameW
          })

      if (length(which(a == FALSE)) != 0) {
        stop("colony list: Wrong values in component \"ULcorner\" or
             wrong dimensions of component \"colImage\"\n
             Out of frame boundaries")
      }

      a <- sapply(col_list, function(x) length(which((x$colImage %in% c(0, 1)) == FALSE)))
      if (length(which(a != 0)) != 0) {
        stop("colony list: Wrong values (not 0 or 1) in component \"colImage\"\n")
      }

    }
  } else if (m1 == "Component \"colImage\" does not exist\n") {
    if (checkField(myList = col_list, fieldName = "ULcorner") == "") {
      stop(paste("colony list:", m1,
                 "Component \"ULcorner\" should be omitted\n"))
    }
  } else {
    stop(paste("colony list:", m1))
  }

}

#################################################################################################

checkCellList <- function(cell_list, col_list) {

  # cellName is the name of the cell, a character string in the format \code{"<cellId>_f<frame>"}.

  if ((m <- checkField(myList = cell_list, fieldName = "cellName")) != "") {
    stop(paste("cell list:", m))
  } else {
    vals <- sapply(cell_list, function(x) x$cellName)
    if (length(grep(pattern = "_f\\d+$", x = vals, value = FALSE)) != length(vals)) {
      stop("cell list: Wrong values in component \"cellName\"\n")
    }
  }

  # frame is the ID of the frame of the cell, a non-zero positive integer number.

  if ((m <- checkField(myList = cell_list, fieldName = "frame")) != "") {
    stop(paste("cell list:", m))
  } else {
    vals <- sapply(cell_list, function(x) x$frame)
    if (class(vals) != "integer") {
      stop("cell list: Values in component \"frame\" are not integers\n")
    }
    if (length(which(vals <= 0)) != 0) {
      stop("cell list: Wrong values in component \"frame\"\n")
    }
    my_cellNames <- sapply(cell_list, function(x) x$cellName)
    if (!all.equal(vals, unname(sapply(my_cellNames, function(x) as.numeric(unlist(strsplit(x, split = "_f"))[2]))))) {
      stop("cell list: Values in component \"frame\" do not correspond to \"_f<frame>\" of component \"cellName\"\n")
    }
  }

  # colony is the ID of the colony of the cell in the frame, a non-zero positive integer number.

  if ((m <- checkField(myList = cell_list, fieldName = "colony")) != "") {
    stop(paste("cell list:", m))
  } else {
    vals <- sapply(cell_list, function(x) x$colony)
    if (class(vals) != "integer") {
      stop("cell list: Values in component \"colony\" are not integers\n")
    }
    if (length(which(vals <= 0)) != 0) {
      stop("cell list: Wrong values in component \"colony\"\n")
    }
  }

  # daughterIds is a vector of character strings containing the cellName
  # of the linked cell(s) in the next frame,
  # or NULL in case no such cells exist.

  if ((m <- checkField(myList = cell_list, fieldName = "daughterIds")) == "Component \"daughterIds\" does not exist\n") {
    stop(paste("cell list:", m))
  } else {
    # vals <- unlist(unname(sapply(cell_list, function(x) x$daughterIds)))
    # if (length(setdiff(vals, unlist(unname(sapply(cell_list, function(x) x$cellName))))) != 0) {
    #   stop("cell list: Wrong values in component \"daughterIds\"\n")
    # }
    # if (length(which(sapply(cell_list, function(x) length(x$daughterIds) > 2) == TRUE)) != 0) {
    #   stop("cell list: More than 2 elements in values of component \"daughterIds\"\n")
    # }
    cells <- unlist(unname(sapply(cell_list, function(x) x$cellName)))
    for (i_cell in 1:length(cell_list)) {
      if (!is.null(cell_list[[i_cell]]$daughterIds)) {
        if (length(setdiff(cell_list[[i_cell]]$daughterIds, cells)) != 0) {
          warning(paste("cell list: Wrong values in component \"daughterIds\" of cell", i_cell, "\n"))
        }
        if (length(cell_list[[i_cell]]$daughterIds) > 2) {
          warning(paste("cell list: More than 2 elements in values of component \"daughterIds\" of cell", i_cell, "\n"))
        }
      }
    }
  }

  # colId is a pointer to the corresponding colony instant of the cell in the col_list,
  # a non-zero positive integer value.
	# Colonies that entered the field of view at a time point and did not exist from the beginning of the movie
  # (i.e. from the first frame)
	# should not have tracked cells, until they merge (if this is the case) with another existing colony.
	# This means that no object should point to such colony instants.

  if (!is.null(col_list)) {
    if ((m <- checkField(myList = cell_list, fieldName = "colId")) != "") {
      stop(paste("colony list is imported\n
                 cell list:", m))
    } else {
      vals <- sapply(cell_list, function(x) x$colId)
      if (length(setdiff(vals, c(1:length(col_list)))) != 0) {
        stop("cell list: Values in component \"colId\" do not point to existing colony instants\n")
      }
      my_colNames <- sapply(col_list, function(x) x$colName)[vals]
      my_frames <- sapply(cell_list, function(x) x$frame)
      my_cols <- sapply(cell_list, function(x) x$colony)
      if (!all.equal(my_colNames, paste("f", my_frames, "_c", my_cols, sep = ""))) {
        stop("cell list: Values in component \"colId\" point to wrong colony instants\n")
      }
    }
  } else {
    if (checkField(myList = cell_list, fieldName = "colId") == "") {
      warning("colony list is not imported\n
              cell list: Component \"colId\" is useless and should be ommited\n")
    }
  }

  # pixelList is a Nx2 matrix of non-zero integer values
  # denoting the pixels of the cell in colony coordinates
  # (i.e. relative to the colImage of the colId^th element in the col_list).
  # Each one of the N rows indicates a pixel of the cell.
  # The first column represents the row and the second the column of the pixel.

  if (!is.null(col_list)) {
    if ((m <- checkField(myList = col_list, fieldName = "colImage")) == "") {
      if (checkField(myList = cell_list, fieldName = "pixelList") == "") {
        for (i_cell in 1:length(cell_list)) {
          if (ncol(cell_list[[i_cell]]$pixelList) != 2) {
            stop("cell list: Wrong number of columns in component \"pixelList\"\n")
          }
          if (length(which(!(cell_list[[i_cell]]$pixelList[, 1] %in% c(1:nrow(col_list[[cell_list[[i_cell]]$colId]]$colImage))))) != 0 ||
              length(which(!(cell_list[[i_cell]]$pixelList[, 2] %in% c(1:ncol(col_list[[cell_list[[i_cell]]$colId]]$colImage))))) != 0) {
            stop("cell list: Pixels in component \"pixelList\" are out of the dimensions of component \"colImage\" in colony list\n")
          }
          vals <- col_list[[cell_list[[i_cell]]$colId]]$colImage[cell_list[[i_cell]]$pixelList]
          if (length(unique(vals)) == 1) {
            if (unique(vals) != 1) {
              warning(paste("cell list: Pixels in component \"pixelList\" of cell", i_cell, "do not denote cells of the colony\n"))
            }
          } else {
            warning(paste("cell list: Pixels in component \"pixelList\" of cell", i_cell, "do not denote cells of the colony\n"))
          }
          # compute and store centroid
          cell_list[[i_cell]]$centroid <- matrix(colMeans(cell_list[[i_cell]]$pixelList), nrow = 1, ncol = 2)
        }
      }
    } else {
      if (checkField(myList = cell_list, fieldName = "pixelList") == "") {
        warning(paste("colony list: ", m,
                      "cell list: Component \"pixelList\" is useless and should be ommited\n", sep = ""))
      }
    }
  } else {
    if (checkField(myList = cell_list, fieldName = "pixelList") == "") {
      warning("colony list is not imported\n
              cell list: Component \"pixelList\" is useless and should be ommited\n")
    }
  }

  ##### rest attributes

  attributes <- unique(unlist(sapply(cell_list, function(x) names(x))))
  attributes <- attributes[!(attributes %in%
                               c("cellName", "frame", "colony", "daughterIds", "colId", "pixelList", "centroid", "boundaryPixelList"))]

  for (attr in attributes) {
    m <- checkField(myList = cell_list, fieldName = attr)
    if (!is.matrix(unname(cell_list[[1]][attr][[1]]))) { # vector or value
      vals <- unname(unlist(sapply(cell_list, function(x) x[attr])))
      if (length(vals) != length(cell_list)) { # vector
        warning(paste("cell list: Component \"", attr, "\" is a vector and should be ommited\n", sep = ""))
      } else { # value
        if (m != "") {
          stop(paste("cell list:", m))
        }
        if (is.character(vals)) {
          warning(paste("cell list: Component \"", attr, "\" is not numeric or boolean and should be ommited\n", sep = ""))
        }
      }
    } else {
      warning(paste("cell list: Component \"", attr, "\" is a matrix and should be ommited\n", sep = ""))
    }
  }

  return(cell_list = cell_list)

}
