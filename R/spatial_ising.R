#' Title
#'
#' @param input_matrix
#' @param B
#' @param J
#' @param iter
#'
#' @return
#' @export
#'
#' @examples
spatial_ising = function(x, B, J, iter){
  if (!inherits(x, "matrix")){
    input_matrix = terra::as.matrix(x, wide = TRUE)
    x_ext = terra::ext(x)
    x_crs = terra::crs(x)
  } else {
    input_matrix = x
  }
  if (missing(iter)){
    iter = nrow(input_matrix) * ncol(input_matrix)
  }
  n_rows = nrow(input_matrix)
  n_cols = ncol(input_matrix)
  rxs = round(stats::runif(iter, min = 1, max = n_rows))
  rys = round(stats::runif(iter, min = 1, max = n_cols))
  for (i in seq_len(iter)){
    input_matrix = single_flip2(input_matrix, B, J, rxs[i], rys[i], n_rows, n_cols)
  }
  if (!inherits(x, "matrix")){
    result = terra::rast(input_matrix, crs = x_crs, extent = x_ext)
  } else {
    result = input_matrix
  }
  return(result)
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
  # n_rows = nrow(input_matrix)
  # n_cols = ncol(input_matrix)
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
    input_matrix[rx, ry] = -input_matrix[rx, ry]
  } else if (stats::runif(1) < exp(-en_diff)){
    input_matrix[rx, ry] = -input_matrix[rx, ry]
  }
  return(input_matrix)
}
