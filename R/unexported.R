createAxisLab <- function(attr, unit = "", string = "") {

  attr <- unlist(strsplit(attr, "_"))

  if (length(attr) == 1) {
    attr <- paste(toupper(substr(attr, 1, 1)), substr(attr, 2, nchar(attr)), sep = "")
  } else if (length(attr) == 2) {
    if (nchar(attr[2]) == 1) {
      attr <- paste(attr[2], " ", toupper(substr(attr[1], 1, 1)), substr(attr[1], 2, nchar(attr[1])), sep="")
    } else {
      attr <- paste(toupper(substr(attr[2], 1, 1)), substr(attr[2], 2, nchar(attr[2])), " ", toupper(substr(attr[1], 1, 1)), substr(attr[1], 2, nchar(attr[1])), sep = "")
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
                                                   labels = c(format(round(min(values), 3), nsmall = 3), format(round(max(values), 3), nsmall = 3)), # 3 decimal digits precision
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

  legend(x = corners[2], y = mean(corners[3:4]), cex = size_lab, title = title, legend = legend, lwd = 3, col = colors, yjust = 0.5, bty = "n", xpd = TRUE)

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

updateCellInBascaLT <- function(LT, cell, cell_list, Ncols) {

  cell_list_ID <- as.numeric(cell) - (Ncols + 1)

  # cellName
  V(LT)[cell]$cellName <- cell_list[[cell_list_ID]]$cellName

  # numeric attributes
  numeric_attrs <- names(cell_list[[1]])[sapply(cell_list[[1]], function(x) class(x) == "numeric" || class(x) == "integer")]
  for (attr in numeric_attrs) {
    if (attr == "colId") {
      V(LT)[cell]$colId <- as.character(cell_list[[cell_list_ID]]$colId)
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

  system(paste(paste(exe_path, "run_mat_correct.sh", sep = "/"), mcr_directory, mat_folder, mat_fileName, cellIDs, frameID, colID, Nsplit))

  newCellProps <- R.matlab::readMat(paste(mat_folder, "mat_correction.mat", sep = "/"))

  # for (i_cell in 1:dim(newCellProps$cellProps)[3]) {
  #   newCellProps$cellProps[, , i_cell]$pixelList <- cbind(newCellProps$cellProps[, , i_cell]$pixelList[, 2], newCellProps$cellProps[, , i_cell]$pixelList[, 1])
  # }

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
      attr <- paste(toupper(substr(attr[2], 1, 1)), substr(attr[2], 2, nchar(attr[2])), " ", toupper(substr(attr[1], 1, 1)), substr(attr[1], 2, nchar(attr[1])), sep = "")
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
  vals <- unlist(unname(sapply(myList, function(x) x[fieldName])))

  if (is.null(vals)) {
    m <- paste("Component \"", fieldName, "\" does not exist\n", sep = "")
  }

  if (length(which(is.na(vals))) != 0) {
    m <- paste("NA values in component \"", fieldName, "\"\n", sep = "")
  }

  return(m)

}
