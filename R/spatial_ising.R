#' Ising model for spatial data
#'
#' Performs simulations based on the given parameters of the Ising model
#'
#' @param x SpatRaster or matrix containing two values: -1 and 1
#' @param B External pressure
#' @param J Peer pressure - it regulates a degree of local harmonization
#' @param updates Specifies how many sets of iterations are performed on the input object.
#' The output of this function has as many layers as the `updates` value.
#' @param iter Specifies how many iterations are performed on the input object.
#' By default it equals to the number of values in the input object.
#' @param version By default, `1`, the `x` object is converted into a matrix
#' (fast, but can be memory consuming); `version = 2` has a lower RAM impact, but
#' is much slower
#' @param progress TRUE/FALSE
#'
#' @references Ising, E., 1924. Beitrag zur theorie des ferro-und paramagnetismus. Ph.D. thesis, Grefe & Tiedemann.
#' @references Onsager, L., 1944. Crystal statistics. I. A two-dimensional model with an order-disorder transition. Physical Review 65 (3-4), 117.
#' @references Brush, S. G., 1967. History of the Lenz-Ising model. Reviews of modern physics 39 (4), 883.
#' @references Cipra, B. A., 1987. An introduction to the Ising model. The American Mathematical Monthly 94 (10), 937–959.
#'
#' @return Object of the same class as `x` with the number of layers specified by `updates`
#' @export
#'
#' @examples
#' data(r_start, package = "spatialising")
#' ts1 = spatial_ising(r_start, B = -0.3, J = 0.7)
#' ts10 = spatial_ising(r_start, B = -0.3, J = 0.7, updates = 10)
#'
#' library(terra)
#' r1 = rast(system.file("raster/r_start.tif", package = "spatialising"))
#' plot(r1)
#' r2 = spatial_ising(r1, B = -0.3, J = 0.7)
#' plot(r2)
#'
#' # ri1 = spatial_ising(r1, B = -0.3, J = 0.7, updates = 9)
#' # plot(ri1)
#'
#' # ri2 = spatial_ising(r1, B = 0.3, J = 0.7, updates = 9)
#' # plot(ri2)
#'
#' # ri3 = spatial_ising(r1, B = -0.3, J = 0.4, updates = 9)
#' # plot(ri3)
spatial_ising = function(x, B, J, updates = 1, iter, version = 1, progress = TRUE){
  if (is.character(x)){
    is_char = TRUE
    x = terra::rast(x)
  }
  if (version == 1){
    is_output_not_matrix = !inherits(x, "matrix")
    if (is_output_not_matrix){
      x_ext = terra::ext(x)
      x_crs = terra::crs(x)
      x = terra::as.matrix(x, wide = TRUE)
    }
    x = spatial_ising_matrix(x = x, B = B, J = J, updates = updates,
                               iter = iter, progress = progress)
    if (is_output_not_matrix){
      x = terra::rast(x, crs = x_crs, extent = x_ext)
      names(x) = paste0("update", seq_len(updates))
    }

  } else if (version == 2){
    x = spatial_ising_terra(x = x, B = B, J = J, updates = updates,
                               iter = iter, progress = progress)
    names(x) = paste0("update", seq_len(updates))
  }
  # if (is_char){
  #   x = wrap(x)
  # }
  return(x)
}

spatial_ising_matrix = function(x, B, J, updates = 1, iter, progress = TRUE){
  if (updates > 1){
    y = vector(mode = "list", length = updates + 1)
    y[[1]] = x
    if (progress) pb = utils::txtProgressBar(min = 2, max = updates + 1, style = 3)
    for (i in seq_len(updates + 1)[-1]){
      y[[i]] = spatial_ising_matrix(y[[i - 1]], B, J, updates = 1, iter, progress = FALSE)
      if (progress) utils::setTxtProgressBar(pb, i)
    }
    if (progress) close(pb)
    x = simplify2array(y[-1])
  } else if (updates == 1){
    n_rows = nrow(x); n_cols = ncol(x)
    if (missing(iter)){
      iter = n_rows * n_cols
    }
    rxs = round(stats::runif(iter, min = 1, max = n_rows))
    rys = round(stats::runif(iter, min = 1, max = n_cols))
    runif_1 = stats::runif(iter)
    for (i in seq_len(iter)){
      x = single_flip2(x, B, J, rxs[i], rys[i], runif_1[i], n_rows, n_cols)
    }
  }
  return(x)
}

spatial_ising_terra = function(x, B, J, updates = 1, iter, progress = TRUE){
  if (updates > 1){
    y = vector(mode = "list", length = updates + 1)
    y[[1]] = x
    if (progress) pb = utils::txtProgressBar(min = 2, max = updates + 1, style = 3)
    for (i in seq_len(updates + 1)[-1]){
      y[[i]] = spatial_ising_terra(y[[i - 1]], B, J, updates = 1, iter, progress = FALSE)
      if (progress) utils::setTxtProgressBar(pb, i)
    }
    if (progress) close(pb)
    x = terra::rast(y[-1])
  } else if (updates == 1){
    n_rows = nrow(x); n_cols = ncol(x)
    if (missing(iter)){
      iter = n_rows * n_cols
    }
    rxs = round(stats::runif(iter, min = 1, max = n_rows))
    rys = round(stats::runif(iter, min = 1, max = n_cols))
    runif_1 = stats::runif(iter)
    for (i in seq_len(iter)){
      x = single_flip2(x, B, J, rxs[i], rys[i], runif_1[i], n_rows, n_cols)
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
# > bench::mark({is1 = spatial_ising(r1, B = -0.3, J = 0.7, updates = 250)})

energy_diff2 = function(focal, neigh, B, J) {
  if (neigh == 4) {
    2 * (B + 4 * J) * focal
  } else if (neigh == 2) {
    2 * (B + 2 * J) * focal
  } else if (neigh == 0) {
    2 * (B + 0 * J) * focal
  } else if (neigh == -2) {
    2 * (B - 2 * J) * focal
  } else if (neigh == -4) {
    2 * (B - 4 * J) * focal
  }
}
single_flip2 = function(input_matrix, B, J, rx, ry, rn, n_rows, n_cols) {
  # choose random spin
  # if (missing(rx)){
  #   rx = round(stats::runif(1, min = 1, max = n_rows))
  # }
  # if (missing(ry)){
  #   ry = round(stats::runif(1, min = 1, max = n_cols))
  # }
  # neighbor sum
  nb = input_matrix[(rx %% n_rows) + 1, ry] + input_matrix[((rx - 2) %% n_rows) + 1, ry] +
    input_matrix[rx, (ry %% n_cols) + 1] + input_matrix[rx, ((ry - 2) %% n_cols) + 1]
  fo = input_matrix[rx, ry]
  en_diff = energy_diff2(fo, nb, B, J)
  if (en_diff <= 0){ #<= or <?
    input_matrix[rx, ry] = -fo
  } else if (rn < exp(-en_diff)){
    input_matrix[rx, ry] = -fo
  }
  return(input_matrix)
}
