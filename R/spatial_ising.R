#' Ising model for spatial data
#'
#' Performs simulations based on the given parameters of the Ising model
#'
#' @param x SpatRaster or matrix containing two values: -1 and 1
#' @param B External pressure
#' @param J Peer pressure - it regulates a degree of local harmonization
#' @param timesteps Specifies how many iterations are performed on the input object.
#' The output of this function has as many layers as the `timesteps` value.
#' @param updates Specifies how many iterations are performed on the input object.
#' By default it equals to the number of values in the input object.
#' @param version By default, `1`, the `x` object is converted into a matrix
#' (fast, but can be memory consuming); `version = 2` has a lower RAM impact, but
#' is much slower
#'
#' @references Ising, E., 1924. Beitrag zur theorie des ferro-und paramagnetismus. Ph.D. thesis, Grefe & Tiedemann.
#' @references Onsager, L., 1944. Crystal statistics. I. A two-dimensional model with an order-disorder transition. Physical Review 65 (3-4), 117.
#' @references Brush, S. G., 1967. History of the Lenz-Ising model. Reviews of modern physics 39 (4), 883.
#' @references Cipra, B. A., 1987. An introduction to the Ising model. The American Mathematical Monthly 94 (10), 937–959.
#'
#' @return Object of the same class as `x` with the number of layers specified by `timesteps`
#' @export
#'
#' @examples
#' data(r_start, package = "spatialising")
#' ts1 = spatial_ising(r_start, B = -0.3, J = 0.7)
#' ts10 = spatial_ising(r_start, B = -0.3, J = 0.7, timesteps = 10)
#'
#' library(terra)
#' r1 = rast(system.file("raster/r_start.tif", package = "spatialising"))
#' plot(r1)
#' r2 = spatial_ising(r1, B = -0.3, J = 0.7)
#' plot(r2)
#'
#' ri1 = spatial_ising(r1, B = -0.3, J = 0.7, timesteps = 9)
#' plot(ri1)
#'
#' ri2 = spatial_ising(r1, B = 0.3, J = 0.7, timesteps = 9)
#' plot(ri2)
#'
#' ri3 = spatial_ising(r1, B = -0.3, J = 0.4, timesteps = 9)
#' plot(ri3)
spatial_ising = function(x, B, J, timesteps = 1, updates, version = 1){
  if (timesteps > 1){
    y = vector(mode = "list", length = timesteps)
    y[[1]] = x
    for (i in seq_len(timesteps)[-1]){
      y[[i]] = spatial_ising(y[[i - 1]], B, J, timesteps = 1, updates, version)
    }
    if (version == 1){
      is_output_not_matrix = !inherits(x, "matrix")
      if (is_output_not_matrix){
        x = terra::rast(y[-1])
      } else {
        x = simplify2array(y[-1])
      }
    }
  } else if (timesteps == 1){
    if (version == 1){
      is_output_not_matrix = !inherits(x, "matrix")
      if (is_output_not_matrix){
        x_ext = terra::ext(x)
        x_crs = terra::crs(x)
        x = terra::as.matrix(x, wide = TRUE)
      }
    }
    n_rows = nrow(x); n_cols = ncol(x)
    if (missing(updates)){
      updates = n_rows * n_cols
    }
    rxs = round(stats::runif(updates, min = 1, max = n_rows))
    rys = round(stats::runif(updates, min = 1, max = n_cols))
    for (i in seq_len(updates)){
      x = single_flip2(x, B, J, rxs[i], rys[i], n_rows, n_cols)
    }
    if (version == 1){
      if (is_output_not_matrix){
        x = terra::rast(x, crs = x_crs, extent = x_ext)
      }
    }
  }
  return(x)
}
energy_diff = function(focal, neigh, B, J){
  if (focal == 1){
    if (neigh == 4){
      2 * (B + 4 * J)
    } else if (neigh == 2){
      2 * (B + 2 * J)
    } else if (neigh == 0){
      2 * (B + 0 * J)
    } else if (neigh == -2){
      2 * (B - 2 * J)
    } else if (neigh == -4){
      2 * (B - 4 * J)
    }
  } else if (focal == -1){
    if (neigh == 4){
      -2 * (B + 4 * J)
    } else if (neigh == 2){
      -2 * (B + 2 * J)
    } else if (neigh == 0){
      -2 * (B + 0 * J)
    } else if (neigh == -2){
      -2 * (B - 2 * J)
    } else if (neigh == -4){
      -2 * (B - 4 * J)
    }
  }
}
single_flip = function(input_matrix, B, J) {
  n_rows = nrow(input_matrix)
  n_cols = ncol(input_matrix)
  # choose random spin
  x = round(stats::runif(1, min = 1, max = n_rows))
  y = round(stats::runif(1, min = 1, max = n_cols))
  # neighbor sum
  nb = input_matrix[(x %% n_rows) + 1, y] + input_matrix[((x - 2) %% n_rows) + 1, y] +
    input_matrix[x, (y %% n_cols) + 1] + input_matrix[x, ((y - 2) %% n_cols) + 1]
  fo = input_matrix[x, y]
  if (energy_diff(fo, nb, B, J) <= 0){ #<= or <?
    input_matrix[x, y] = -input_matrix[x, y]
  } else if (stats::runif(1) < exp(-energy_diff(fo, nb, B, J))){
    input_matrix[x, y] = -input_matrix[x, y]
  }
  return(input_matrix)
}
single_flip2 = function(input_matrix, B, J, rx, ry, n_rows, n_cols) {
  # choose random spin
  if (missing(rx)){
    rx = round(stats::runif(1, min = 1, max = n_rows))
  }
  if (missing(ry)){
    ry = round(stats::runif(1, min = 1, max = n_cols))
  }
  # neighbor sum
  nb = input_matrix[(rx %% n_rows) + 1, ry] + input_matrix[((rx - 2) %% n_rows) + 1, ry] +
    input_matrix[rx, (ry %% n_cols) + 1] + input_matrix[rx, ((ry - 2) %% n_cols) + 1]
  fo = input_matrix[rx, ry]
  en_diff = energy_diff(fo, nb, B, J)
  if (en_diff <= 0){ #<= or <?
    input_matrix[rx, ry] = -fo
  } else if (stats::runif(1) < exp(-en_diff)){
    input_matrix[rx, ry] = -fo
  }
  return(input_matrix)
}
