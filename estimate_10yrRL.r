# Guillaume Evin
# 01/05/2026
# This script estimates 10-year return level from different RCM simulations and compares them to the return level estimated from observations.

library(extRemes)
library(ncdf4)

setwd("D:/EUROCORDEX_extremes/evaluation_precip_extremes_eurocordex/")

# list models and obs
name_models = c('ERA5_ALARO1-SFX','ERA5_CCLM6-0-1','ERA5_CNRM-ALADIN64E1',
                'ERA5_HCLIM43-ALADIN','ERA5_ICON-CLM-202407-1-1','ERA5_RACMO23E',
                'ERA5_RegCM5-0','ERA5_REMO2020-2-2','ERA5_WRF451Q')


# first loop over the models
for(model in name_models){
  print(paste0("Processing model ", model))
  list_nc = list()
  for(season in c("DJF", "MAM", "JJA", "SON")){
    fnetcdf = paste0("D:/EUROCORDEX_extremes/output_juelich_apr/", model, "1980-2020_", season, "_max_scale_S1.nc")
    list_nc[[season]] = nc_open(fnetcdf)
  }

  # get dimensions
  nx = list_nc[[season]]$dim$x$len
  ny = list_nc[[season]]$dim$y$len
  timescale = list_nc[[season]]$dim$timescale$vals


  #loop over all grid points and seasons to retrieve the variable of interest
  for(itimescale in 1:2){
    ts = timescale[itimescale]
    
    # prepare matrices
    mat_10yrRL = mat_shape = matrix(NA, nrow = ny, ncol = nx)
    
    # loop over all grid points
    for(ix in 1:nx){
      for(iy in 1:ny){
        list_maxpr = list()
        for(season in c("DJF", "MAM", "JJA", "SON")){
          # retrieve the variable of interest only for one time scale
          maxpr = ncvar_get(list_nc[[season]], "maxpr", start = c(ix, iy, itimescale, 1), count = c(1, 1, 1, -1))
          list_maxpr[[season]] = maxpr
        }

        mat_max <- matrix(unlist(list_maxpr), ncol = 4)
        annual_max <- apply(mat_max, 1, max)

        # Fit GEV using L-moments
        fit_lmom <- fevd(annual_max, type = "GEV", method = "Lmoments")

        mat_shape[iy,ix] = fit_lmom$results[3]

        # 10-year return level
        mat_10yrRL[iy,ix] = qevd(1 - 1/10, loc = fit_lmom$results[1], scale = fit_lmom$results[2], shape = fit_lmom$results[3])
      }
    }

    write.table(mat_10yrRL, paste0("D:/EUROCORDEX_extremes/10yrRL_BenPosch/matrix_10yrRL_", ts, "hr_", model, ".csv"), sep = ",", col.names = FALSE, row.names = FALSE)
    write.table(mat_shape, paste0("D:/EUROCORDEX_extremes/10yrRL_BenPosch/matrix_shape_", ts, "hr_", model, ".csv"), sep = ",", col.names = FALSE, row.names = FALSE)
  }
}