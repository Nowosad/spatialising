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
  input_matrix = terra::as.matrix(x, wide = TRUE)
  x_ext = terra::ext(x)
  x_crs = terra::crs(x)
  if (missing(iter)){
    iter = nrow(input_matrix) ^ 2
  }
  for (i in seq_len(iter)){
    input_matrix = single_flip(input_matrix, B, J)
  }
  result = terra::rast(input_matrix, crs = x_crs, extent = x_ext)
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

