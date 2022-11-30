#' Ensemble of Ising models for spatial data
#'
#' Creates an ensemble of simulations based on the given parameters of the Ising model
#'
#' @param runs A number of simulations to perform
#' @param ... Arguments for [spatial_ising()]
#'
#' @return A list of objects of the same class as `x`
#' @export
#'
#' @examples
#' data(r_start, package = "spatialising")
#' l = spatial_ising_ensemble(100, r_start, B = -0.3, J = 0.7)
#'
#' # library(terra)
#' # r1 = rast(system.file("raster/r_start.tif", package = "spatialising"))
#' # plot(r1)
#' # r2 = spatial_ising_ensemble(100, r1, B = -0.3, J = 0.7)
spatial_ising_ensemble = function(runs, ...){
  l = lapply(seq_len(runs), \(i) spatial_ising(...))
  return(l)
}


#' Exemplar of an Ising model for spatial data
#'
#' Creates an ensemble of simulations based on the given parameters
#' of the Ising model and selects an exemplar
#' (a model that is closest to the average of the ensemble)
#'
#' @param runs A number of simulations to perform
#' @param ... Arguments for [spatial_ising()]
#'
#' @return Object of the same class as `x`
#' @export
#'
#' @examples
#' data(r_start, package = "spatialising")
#' l = spatial_ising_exemplar(100, r_start, B = -0.3, J = 0.7)
#'
#' # library(terra)
#' # r1 = rast(system.file("raster/r_start.tif", package = "spatialising"))
#' # plot(r1)
#' # r2 = spatial_ising_exemplar(100, r1, B = -0.3, J = 0.7)
#' # plot(r2)
spatial_ising_exemplar = function(runs, ...){
  l = spatial_ising_ensemble(runs = runs, ...)
  nlayers = ifelse(is.na(dim(l[[1]])[3]), 1, dim(l[[1]])[3])
  lc = vapply(l, magnetization, numeric(nlayers))
  la = vapply(l, texture_index, numeric(nlayers))
  if(inherits(lc, "matrix")){ #updates >1
    lc = t(lc)
    la = t(la)
    avg_metrics = c(colMeans(lc), colMeans(la))
  } else { #updates == 1
    avg_metrics = cbind(mean(lc), mean(la))
  }
  l_metrics = cbind(lc, la)
  all_metrics = rbind(avg_metrics, l_metrics)
  dist_to_avg = as.matrix(stats::dist(all_metrics))[, 1][-1]
  l_close_to_avg = l[[which.min(dist_to_avg)]]
  return(l_close_to_avg)
}
