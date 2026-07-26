# Guillaume Evin
#
# 14/07/2026
# This script estimates 20-year return level from different RCM simulations
# for:
# - all seasons
# - all durations: 1h, 3h, 6h, 12h, 24h, 72h
# - all spatial scales: S1, S2, S3, S4, S5

library(extRemes)
library(ncdf4)

############################## MODELS ##########################################

# list models and obs
name_models = c() # c('ERA5_CNRM-ALADIN64E1','ERA5_ALARO1-SFX','ERA5_CCLM6-0-1',
#   'ERA5_HCLIM43-ALADIN','ERA5_ICON-CLM-202407-1-1','ERA5_RACMO23E',
#    'ERA5_RegCM5-0','ERA5_REMO2020-2-2','ERA5_WRF451Q')


# first loop over the models
for(model in name_models){
  print(paste0("Processing model ", model))

  for(season in c("DJF", "MAM", "JJA", "SON")){
    print(paste0("Processing season ", season))
    for(spatialscale in c("S1","S2","S3","S4","S5")){
      print(paste0("Processing spatial scale ", spatialscale))
      fnetcdf = paste0("D:/EUROCORDEX_extremes/DATA/SIM/",model, "1980-2020_",
                       season, "_max_scale_", spatialscale,".nc")
      nc = nc_open(fnetcdf)

      # get dimensions
      nx = nc$dim$x$len
      ny = nc$dim$y$len

      timescale = nc$dim$timescale$vals
      ntimescale = length(timescale)

      #loop over all grid points and seasons to retrieve the variable of interest
      for(itimescale in 1:ntimescale){
        ts = timescale[itimescale]

        # prepare matrices
        mat_20yrRL = mat_shape = matrix(NA, nrow = ny, ncol = nx)

        # loop over all grid points
        for(ix in 1:nx){
          for(iy in 1:ny){
            annual_max = ncvar_get(nc, "maxpr", start = c(ix, iy, itimescale, 1),
                                   count = c(1, 1, 1, -1))

            if(!any(annual_max==0)){
              # Fit GEV using L-moments
              fit_lmom <- fevd(annual_max, type = "GEV", method = "Lmoments")
              par_hat = fit_lmom$results

              mat_shape[iy,ix] = par_hat[3]

              # 10-year return level
              mat_20yrRL[iy,ix] = qevd(1 - 1/20, loc = par_hat[1], scale = par_hat[2], shape = par_hat[3])
            }
          }
        }

        write.table(mat_20yrRL,
                    paste0("D:/EUROCORDEX_extremes/20yrRL/matrix_10yrRL_",
                           ts, "hr_", spatialscale,"_",season,"_", model, ".csv"),
                    sep = ",", col.names = FALSE, row.names = FALSE,na="-9999")
        write.table(mat_shape,
                    paste0("D:/EUROCORDEX_extremes/20yrRL/matrix_shape_",
                           ts, "hr_", spatialscale,"_",season,"_", model, ".csv"),
                    sep = ",", col.names = FALSE, row.names = FALSE,na="-9999")
      }
    }
  }
}



##############################  OBS   ##########################################

# list models and obs
name_obs = c('OBS_EURADCLIM_011EUi_2013-2022',
             'OBS_COMEPHORE_011EUi_1997-2022')#,'OBS_CERRALAND_011EUi_1986-2020')
# 'OBS_RADKLIM_011EUi_2001-2022',
#'OBS_GRIPHO_011EUi_2001-2016', 'OBS_ERA5_011EUi_1980-2022',

# first loop over the models
for(obs in name_obs){
  print(paste0("Processing obs ", obs))

  for(season in c("DJF", "MAM", "JJA", "SON")){
    print(paste0("Processing season ", season))
    for(spatialscale in c("S1","S2","S3","S4","S5")){
      print(paste0("Processing spatial scale ", spatialscale))
      fnetcdf = paste0("D:/EUROCORDEX_extremes/DATA/OBS/",obs,'_',
                       season, "_max_scale_", spatialscale,".nc")
      nc = nc_open(fnetcdf)

      # get dimensions
      nx = nc$dim$x$len
      ny = nc$dim$y$len

      timescale = nc$dim$timescale$vals
      ntimescale = length(timescale)

      #loop over all grid points and seasons to retrieve the variable of interest
      for(itimescale in 1:ntimescale){
        ts = timescale[itimescale]

        # prepare matrices
        mat_20yrRL = mat_shape = matrix(NA, nrow = ny, ncol = nx)

        # loop over all grid points
        for(ix in 1:nx){
          for(iy in 1:ny){
            annual_max = ncvar_get(nc, "maxpr", start = c(ix, iy, itimescale, 1),
                                   count = c(1, 1, 1, -1))

            has_zero = any(annual_max==0)
            has_na = any(is.na(annual_max))

            if(!has_na){
              if(!has_zero){
                # Fit GEV using L-moments
                fit_lmom <- fevd(annual_max, type = "GEV", method = "Lmoments")
                par_hat = fit_lmom$results

                mat_shape[iy,ix] = par_hat[3]

                # 10-year return level
                mat_20yrRL[iy,ix] = qevd(1 - 1/20, loc = par_hat[1], scale = par_hat[2], shape = par_hat[3])
              }
            }
          }
        }

        write.table(mat_20yrRL,
                    paste0("D:/EUROCORDEX_extremes/20yrRL/matrix_10yrRL_",
                           ts, "hr_", spatialscale,"_",season,"_", obs, ".csv"),
                    sep = ",", col.names = FALSE, row.names = FALSE, na="-9999")
        write.csv2(mat_shape,
                    paste0("D:/EUROCORDEX_extremes/20yrRL/matrix_shape_",
                           ts, "hr_", spatialscale,"_",season,"_", obs, ".csv"),
                    sep = ",", col.names = FALSE, row.names = FALSE, na="-9999")
      }
    }
  }
}
