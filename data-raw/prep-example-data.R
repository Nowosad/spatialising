RandomFieldsUtils::RFoptions(install = "no")
library(terra)
library(NLMR)
library(landscapetools)
set.seed(2022-03-23)
r1a = NLMR::nlm_fbm(50, 50, fract_dim = 1)
r2a = landscapetools::util_classify(r1a, weighting = c(0.95, 0.05))
rcl = matrix(c(1, 2, -1, 1), ncol = 2)
r_start = rast(r2a) * -1 + 3
r_start = classify(r_start, rcl)
plot(r_start)
writeRaster(r_start, "inst/raster/r_start.tif", gdal = "of=COG", overwrite = TRUE)
r_start = as.matrix(r_start, wide = TRUE)
usethis::use_data(r_start, overwrite = TRUE)

set.seed(2022-04-04)
r_end = spatial_ising(r_start, B = -0.3, J = 0.7)
usethis::use_data(r_end, overwrite = TRUE)


