#' An example binary raster
#'
#' A raster file covering an area of 50x50 cells. The raster file contains two values: -1 and 1.
#' `system.file("raster/r_start.tif", package = "spatialising")`
#'
#' @format A raster file
#' @name r_start.tif
NULL

#' An example matrix object
#'
#' A matrix has 50 columns and 50 rows. The matrix contains two values: -1 and 1.
#'
#' @format A matrix
#' @usage data(r_start)
"r_start"

#' An example matrix object
#'
#' A matrix has 50 columns and 50 rows. The matrix contains two values: -1 and 1.
#' This object was created with `r_end = kinetic_ising(r_start, B = -0.3, J = 0.7)`
#'
#' @format A matrix
#' @usage data(r_end)
"r_end"
