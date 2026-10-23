#' Composition imbalance index
#'
#' Calculates composition imbalance index (also known as the m index) -- a sum of cell’s values
#' over the entire site divided by the number of cell in the site.
#' m has a range from -1 (site completely dominated by the -1 values) to
#' 1 (site completely dominated by the 1 values). Missing cells are ignored.
#'
#' @param x SpatRaster or matrix containing two values: -1 and 1
#'
#' @return A numeric vector
#' @seealso [spatialising::kinetic_ising()]
#' @export
#'
#' @examples
#' data(r_start, package = "spatialising")
#' composition_index(r_start)
#' ts1 = kinetic_ising(r_start, B = -0.3, J = 0.7)
#' composition_index(ts1)
#' ts2 = kinetic_ising(r_start, B = -0.3, J = 0.7, updates = 2)
#' composition_index(ts2)
#'
#' \donttest{
#'   library(terra)
#'   r1 = rast(system.file("raster/r_start.tif", package = "spatialising"))
#'   composition_index(r1)
#'   r2 = kinetic_ising(r1, B = -0.3, J = 0.7)
#'   composition_index(r2)
#' }
composition_index = function(x){
  if (inherits(x, "matrix")){
    mean(x, na.rm = TRUE)
  } else if (inherits(x, "array")) {
    apply(x, 3, mean, na.rm = TRUE)
  } else {
    terra::global(x, "mean", na.rm = TRUE)$mean
  }
}
#' Texture index
#'
#' Calculates the mean product of neighboring, non-missing cell values.
#' The texture index ranges from -1 (unlike values are always adjacent) through
#' 0 (no spatial association) to 1 (like values are always adjacent).
#'
#' @param x SpatRaster or matrix containing two values: -1 and 1
#' @param ... Arguments for [comat::get_coma()]
#'
#' @return A numeric vector
#' @seealso [spatialising::kinetic_ising()]
#' @export
#'
#' @examples
#' data(r_start, package = "spatialising")
#' texture_index(r_start)
#' ts1 = kinetic_ising(r_start, B = -0.3, J = 0.7)
#' texture_index(ts1)
#' ts2 = kinetic_ising(r_start, B = -0.3, J = 0.7, updates = 2)
#' texture_index(ts2)
#'
#' \donttest{
#'   library(terra)
#'   r1 = rast(system.file("raster/r_start.tif", package = "spatialising"))
#'   texture_index(r1)
#'   r2 = kinetic_ising(r1, B = -0.3, J = 0.7)
#'   texture_index(r2)
#' }
texture_index = function(x, ...){
  if (inherits(x, "matrix")){
    coma = comat::get_coma(x, ...)
    if (length(coma) == 1){
      return(1)
    } else {
      return((coma[1] - coma[2] - coma[3] + coma[4]) / sum(coma))
    }
  } else {
    x = terra::as.array(x)
    apply(x, 3, texture_index, ...)
  }
}
