#' Ising model for spatial data
#'
#' Performs simulations based on the given parameters of the Ising model
#'
#' @param x SpatRaster or matrix containing two values: -1 and 1
#' @param B External pressure (positive or negative): it tries to align cells' values with its sign
#' @param J Strength of the local autocorrelation tendency (always positive): it tries to align signs of neighboring cells
#' @param updates Specifies how many sets of iterations are performed on the input object.
#' The output of this function has as many layers as the `updates` value.
#' @param iter Specifies how many iterations are performed on the input object.
#' By default it equals to the number of values in the input object.
#' @param rule IM temporal evolution rule: either `"glauber"` (default) or `"metropolis"`
#' @param inertia Represents the modification of the algorithm aimed at suppressing the salt-and-pepper noise of the focus category present when simulating evolution of a coarse-textured pattern. With Q > 0, small patches of the focus category are not generated, thus eliminating the salt-and-pepper noise of the focus category
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
#' @importFrom Rcpp evalCpp
#'
#' @return Object of the same class as `x` with the number of layers specified by `updates`
#' @export
#'
#' @examples
#' data(r_start, package = "spatialising")
#' ts1 = kinetic_ising(r_start, B = -0.3, J = 0.7)
#' ts10 = kinetic_ising(r_start, B = -0.3, J = 0.7, updates = 10)
#'
#' \donttest{
#'   r1 = terra::rast(system.file("raster/r_start.tif", package = "spatialising"))
#'   terra::plot(r1)
#'   r2 = kinetic_ising(r1, B = -0.3, J = 0.7)
#'   terra::plot(r2)
#'
#'   library(terra)
#'   ri1 = kinetic_ising(r1, B = -0.3, J = 0.7, updates = 9)
#'   plot(ri1)
#'
#'   ri2 = kinetic_ising(r1, B = 0.3, J = 0.7, updates = 9)
#'   plot(ri2)
#'
#'   ri3 = kinetic_ising(r1, B = -0.3, J = 0.4, updates = 9)
#'   plot(ri3)
#' }
kinetic_ising = function(x, B, J, updates = 1, iter, rule = "glauber",
                         inertia = 0, version = 1, progress = FALSE){
  if (is.character(x)){
    is_char = TRUE
    x = terra::rast(x)
  }
  check_if_proper_binary(x)
  if (version == 1){
    is_output_not_matrix = !inherits(x, "matrix")
    if (is_output_not_matrix){
      x_ext = terra::ext(x)
      x_crs = terra::crs(x)
      result = terra::as.matrix(x, wide = TRUE)
    } else {
      result = x
    }
    result = kinetic_ising_matrix(x = result, B = B, J = J, updates = updates,
                               iter = iter, rule = rule, inertia = inertia, progress = progress)
    if (is_output_not_matrix){
      result = terra::rast(result, crs = x_crs, extent = x_ext)
      names(result) = paste0("update", seq_len(updates))
    }

  } else if (version == 2){
    result = x
    result = kinetic_ising_terra(x = result, B = B, J = J, updates = updates,
                               iter = iter, rule = rule, inertia = inertia, progress = progress)
    names(result) = paste0("update", seq_len(updates))
  }
  # if (is_char){
  #   result = wrap(result)
  # }
  return(result)
}

kinetic_ising_matrix = function(x, B, J, updates = 1, iter, rule, inertia, progress = TRUE){
  if (updates > 1){
    y = vector(mode = "list", length = updates + 1)
    y[[1]] = x
    if (progress) pb = utils::txtProgressBar(min = 2, max = updates + 1, style = 3)
    for (i in seq_len(updates + 1)[-1]){
      y[[i]] = kinetic_ising_matrix(y[[i - 1]], B, J, updates = 1, iter, rule, inertia, progress = FALSE)
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
    rns = stats::runif(iter)
    if (rule == "glauber"){
      x = flip_glauber_rcpp(x, B, J, rxs, rys, rns, n_rows, n_cols, inertia)
    } else if (rule == "metropolis"){
      x = flip_metropolis2_rcpp(x, B, J, rxs, rys, rns, n_rows, n_cols, inertia)
    } else {
      stop()
    }
  }
  return(x)
}


# rcpp is used instead here
# flip_glauber = function(input_matrix, B, J, rxs, rys, rns, n_rows, n_cols, inertia) {
#   for (i in seq_along(rns)){
#     rx = rxs[i]; ry = rys[i]; rn = rns[i]
#     # neighbor sum
#     nb = input_matrix[(rx %% n_rows) + 1, ry] + input_matrix[((rx - 2) %% n_rows) + 1, ry] +
#       input_matrix[rx, (ry %% n_cols) + 1] + input_matrix[rx, ((ry - 2) %% n_cols) + 1]
#     fo = input_matrix[rx, ry]
#     en_diff = energy_diff2(fo, nb, B, J, inertia)
#     P = 1 / (1 + exp(en_diff))
#     if (P > rn){
#       input_matrix[rx, ry] = -fo
#     }
#   }
#   return(input_matrix)
# }
#

# rcpp is used instead here
# flip_metropolis2 = function(input_matrix, B, J, rx, ry, rn, n_rows, n_cols, inertia) {
#   for (i in seq_along(rns)){
#     rx = rxs[i]; ry = rys[i]; rn = rns[i]
#     # neighbor sum
#     nb = input_matrix[(rx %% n_rows) + 1, ry] + input_matrix[((rx - 2) %% n_rows) + 1, ry] +
#       input_matrix[rx, (ry %% n_cols) + 1] + input_matrix[rx, ((ry - 2) %% n_cols) + 1]
#     fo = input_matrix[rx, ry]
#     en_diff = energy_diff2(fo, nb, B, J, inertia)
#     if (en_diff <= 0){ #<= or <?
#       input_matrix[rx, ry] = -fo
#     } else if (rn < exp(-en_diff)){
#       input_matrix[rx, ry] = -fo
#     }
#   }
#   return(input_matrix)
# }

kinetic_ising_terra = function(x, B, J, updates = 1, iter, rule, inertia, progress = TRUE){
  if (updates > 1){
    y = vector(mode = "list", length = updates + 1)
    y[[1]] = x
    if (progress) pb = utils::txtProgressBar(min = 2, max = updates + 1, style = 3)
    for (i in seq_len(updates + 1)[-1]){
      y[[i]] = kinetic_ising_terra(y[[i - 1]], B, J, updates = 1, iter, rule, inertia, progress = FALSE)
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
    rns = stats::runif(iter)
    for (i in seq_len(iter)){
      if (rule == "glauber"){
        x = single_flip_glauber(x, B, J, rxs[i], rys[i], rns[i], n_rows, n_cols, inertia)
      } else if (rule == "metropolis"){
        x = single_flip_metropolis2(x, B, J, rxs[i], rys[i], rns[i], n_rows, n_cols, inertia)
      } else {
        stop()
      }
    }
  }
  return(x)
}

energy_diff2 = function(focal, neigh, B, J, inertia) {
  en_diff = 2 * (B + neigh * J) * focal
  # if ((fo == -1 && nb == -4) || (fo == 1 && nb == 4)){
  if (focal == -1 && neigh == -4){
    en_diff = en_diff + inertia
  }
  return(en_diff)
}

single_flip_glauber = function(input_matrix, B, J, rx, ry, rn, n_rows, n_cols, inertia) {
  # neighbor sum
  nb = input_matrix[(rx %% n_rows) + 1, ry] +
    input_matrix[((rx - 2) %% n_rows) + 1, ry] +
    input_matrix[rx, (ry %% n_cols) + 1] +
    input_matrix[rx, ((ry - 2) %% n_cols) + 1]
  fo = input_matrix[rx, ry]
  en_diff = energy_diff2(fo, nb, B, J, inertia)
  P = 1 / (1 + exp(en_diff))
  if (P > rn){
    input_matrix[rx, ry] = -fo
  }
  return(input_matrix)
}

single_flip_metropolis2 = function(input_matrix, B, J, rx, ry, rn, n_rows, n_cols, inertia) {
  # neighbor sum
  nb = input_matrix[(rx %% n_rows) + 1, ry] +
    input_matrix[((rx - 2) %% n_rows) + 1, ry] +
    input_matrix[rx, (ry %% n_cols) + 1] +
    input_matrix[rx, ((ry - 2) %% n_cols) + 1]
  fo = input_matrix[rx, ry]
  en_diff = energy_diff2(fo, nb, B, J, inertia)
  if (en_diff <= 0){ #<= or <?
    input_matrix[rx, ry] = -fo
  } else if (rn < exp(-en_diff)){
    input_matrix[rx, ry] = -fo
  }
  return(input_matrix)
}

# slow, legacy code:
# energy_diff_metropolis = function(focal, neigh, B, J){
#   if (focal == 1){
#     if (neigh == 4){
#       2 * (B + 4 * J)
#     } else if (neigh == 2){
#       2 * (B + 2 * J)
#     } else if (neigh == 0){
#       2 * (B + 0 * J)
#     } else if (neigh == -2){
#       2 * (B - 2 * J)
#     } else if (neigh == -4){
#       2 * (B - 4 * J)
#     }
#   } else if (focal == -1){
#     if (neigh == 4){
#       -2 * (B + 4 * J)
#     } else if (neigh == 2){
#       -2 * (B + 2 * J)
#     } else if (neigh == 0){
#       -2 * (B + 0 * J)
#     } else if (neigh == -2){
#       -2 * (B - 2 * J)
#     } else if (neigh == -4){
#       -2 * (B - 4 * J)
#     }
#   }
# }

# single_flip_metropolis = function(input_matrix, B, J) {
#   n_rows = nrow(input_matrix)
#   n_cols = ncol(input_matrix)
#   # choose random spin
#   x = round(stats::runif(1, min = 1, max = n_rows))
#   y = round(stats::runif(1, min = 1, max = n_cols))
#   # neighbor sum
#   nb = input_matrix[(x %% n_rows) + 1, y] + input_matrix[((x - 2) %% n_rows) + 1, y] +
#     input_matrix[x, (y %% n_cols) + 1] + input_matrix[x, ((y - 2) %% n_cols) + 1]
#   fo = input_matrix[x, y]
#   if (energy_diff_metropolis(fo, nb, B, J) <= 0){ #<= or <?
#     input_matrix[x, y] = -input_matrix[x, y]
#   } else if (stats::runif(1) < exp(-energy_diff_metropolis(fo, nb, B, J))){
#     input_matrix[x, y] = -input_matrix[x, y]
#   }
#   return(input_matrix)
# }

check_if_proper_binary = function(x){
    if (inherits(x, "matrix")){
        sample_size = c(length(x), 10000)
        s = sample(x, size = sample_size[which.min(c(length(x), 10000))])
    } else {
        x_ncell = terra::ncell(x)
        sample_size = c(x_ncell, 10000)
        s_id = sample(seq_len(x_ncell), size = sample_size[which.min(c(x_ncell, 10000))])
        s = terra::extract(x, s_id)[[1]]
    }
    unique_s = unique(s)
    if (!all(unique_s %in% c(-1, 1))){
        stop("The input raster can only contain values of -1 and 1", call. = FALSE)
    }
}

