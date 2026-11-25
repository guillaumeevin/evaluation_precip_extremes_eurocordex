# Guillaume Evin
# 27/07/2026
# This script:
# - estimates annual 10-year return levels at different temporal scales from different RCM simulations and compares them to the return level
# estimated from observations (Poschlod et al., 2021, 10.5194/essd-13-983-2021).
# - estimates 20-year return levels for different seasons, time scales and spatial scales from different RCM simulations
# and compares them to the return level estimated from observations

library(extRemes)
library(ncdf4)
library(parallel) #detectCores
library(foreach) #parallelization
library(doParallel) #parallelization
library(doSNOW)
library(boot)

# list models
name_models = c("ERA5_CCLM6-0-1-URB", "ERA5_CCLM6-0-1-URB-ESG", "ERA5_REMO2020-2-2-MR2",
                "ERA5_REMO2020-2-2-TEB", "ERA5_REMO2020-2-2-iMOVE", "ERA5_REMO2020-2-2-iMOVE-LUC",
                'ERA5_WRF451Q-r2','ERA5_WRF451Q-r1','ERA5_TSMP1-140-E',
                'ERA5_ROAM-NBS','ERA5_GCOAST-AHOIB1-1', 'ERA5_ALARO1-SFX',
                'ERA5_CCLM6-0-1','ERA5_CNRM-ALADIN64E1',
                'ERA5_HCLIM43-ALADIN','ERA5_ICON-CLM-202407-1-1','ERA5_RACMO23E',
                'ERA5_RegCM5-0','ERA5_REMO2020-2-2','ERA5_WRF451Q')

# list obs
name_obs = c('OBS_KNMI_011EUi_1998-2025','OBS_CPC_011EUi_2005-2025',"OBS_INCA_011EUi_2013-2025",
  'OBS_CEH-GEAR_011EUi_1990-2016', 'OBS_KLIMGRID_011EUi_2011-2019','OBS_COMEPHORE_011EUi_1997-2025',
  'OBS_RADKLIM_011EUi_2001-2025', 'OBS_GRIPHO_011EUi_2001-2016', 
  'OBS_EMO_011EUi_1990-2022', 'OBS_ERA5_011EUi_1980-2020', 'OBS_CERRALAND_011EUi_1986-2020')

#####################################################################################################################
# STEP 1: annual 10-year return level, different temporal scales
#####################################################################################################################

#________________________________________________
# first loop over the models
#________________________________________________
if(FALSE){
  for(model in name_models){
    print(paste0("Processing model ", model))
    list_nc = list()
    for(season in c("DJF", "MAM", "JJA", "SON")){
      fnetcdf = paste0("/jarre/gevin/EUROCORDEX/DATA/SIM/", model, "1980-2020_", season, "_max_scale_S1.nc")
      list_nc[[season]] = nc_open(fnetcdf)
    }

    # get dimensions
    nx = list_nc[[season]]$dim$x$len
    ny = list_nc[[season]]$dim$y$len
    timescale = list_nc[[season]]$dim$timescale$vals
    ntimescale = length(timescale)


    #loop over all grid points and seasons to retrieve the variable of interest
    for(itimescale in 1:ntimescale){
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

          has_zero = any(annual_max==0)
          has_na = any(is.na(annual_max))

          if(!has_na){
            if(!has_zero){
              # Fit GEV using L-moments
              fit_lmom <- fevd(annual_max, type = "GEV", method = "Lmoments")

              mat_shape[iy,ix] = fit_lmom$results[3]

              # 10-year return level
              mat_10yrRL[iy,ix] = qevd(1 - 1/10, loc = fit_lmom$results[1], scale = fit_lmom$results[2], shape = fit_lmom$results[3])
            }
          }
        }
      }

      write.table(mat_10yrRL, paste0("/jarre/gevin/EUROCORDEX/RL/10YR/matrix_10yrRL_", ts, "hr_", model, ".csv"),
                  sep = ",", col.names = FALSE, row.names = FALSE)
    }

    # close netcdf
    for(fnetcdf in list_nc){
      nc_close(fnetcdf)
    }
  }
}

#________________________________________________
# second loop over the obs
#________________________________________________

if(FALSE){
  # first loop over the models
  for(obs in name_obs){
    print(paste0("Processing obs ", obs))

    list_nc = list()
    for(season in c("DJF", "MAM", "JJA", "SON")){
      fnetcdf = paste0("/jarre/gevin/EUROCORDEX/DATA/OBS/", obs,"_", season, "_max_scale_S1.nc")
      list_nc[[season]] = nc_open(fnetcdf)
    }

    # get dimensions
    nx = list_nc[[season]]$dim$x$len
    ny = list_nc[[season]]$dim$y$len
    timescale = list_nc[[season]]$dim$timescale$vals
    ntimescale = length(timescale)


    #loop over all grid points and seasons to retrieve the variable of interest
    for(itimescale in 1:ntimescale){
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

          has_zero = any(annual_max==0)
          has_na = any(is.na(annual_max))

          if(!has_na){
            if(!has_zero){
              # Fit GEV using L-moments
              fit_lmom <- fevd(annual_max, type = "GEV", method = "Lmoments")

              mat_shape[iy,ix] = fit_lmom$results[3]

              # 10-year return level
              mat_10yrRL[iy,ix] = qevd(1 - 1/10, loc = fit_lmom$results[1], scale = fit_lmom$results[2], shape = fit_lmom$results[3])
            }
          }
        }
      }

      write.table(mat_10yrRL, paste0("/jarre/gevin/EUROCORDEX/RL/10YR/matrix_10yrRL_", ts, "hr_", obs, ".csv"),
                  sep = ",", col.names = FALSE, row.names = FALSE)
    }

    # close netcdf
    for(fnetcdf in list_nc){
      nc_close(fnetcdf)
    }
  }
}

#####################################################################################################################
# STEP 2: seasonal 20-year return level, different temporal  and spatial scales
#####################################################################################################################

#Number of cores for parallelization
nbcore=20
cl <- makeSOCKcluster(nbcore)
registerDoSNOW(cl)

# number of bootstrap samples
nboot = 100

bootFunc <- function(mydata, i){
  rnddata <- mydata[i]
  tryCatch(
    expr = {
      # Fit GEV using L-moments
      fit_lmom <- fevd(rnddata, type = "GEV", method = "Lmoments")
      par_hat = fit_lmom$results
      qevd(1 - 1/20, loc = par_hat[1], scale = par_hat[2], shape = par_hat[3])
    },
    error = function(cond){
      NA
    },
    warning = function(cond){
      NA
    }
  )   
}

# first loop over the models

# list models
name_models = c('ERA5_WRF451Q-r1','ERA5_TSMP1-140-E','ERA5_ROAM-NBS',
                'ERA5_GCOAST-AHOIB1-1', 'ERA5_ALARO1-SFX', 'ERA5_CCLM6-0-1',
                'ERA5_CNRM-ALADIN64E1', 'ERA5_HCLIM43-ALADIN','ERA5_ICON-CLM-202407-1-1',
                'ERA5_RACMO23E','ERA5_RegCM5-0','ERA5_REMO2020-2-2')

if(TRUE){
  for(model in name_models){
    print(paste0("Processing model ", model))

    for(season in c("DJF", "MAM", "JJA", "SON")){
      print(paste0("Processing season ", season))
      for(spatialscale in c("S1","S2","S3","S4","S5")){
        print(paste0("Processing spatial scale ", spatialscale))
        fnetcdf = paste0("/jarre/gevin/EUROCORDEX/DATA/SIM/",model, "1980-2020_",
                         season, "_max_scale_", spatialscale,".nc")
        nc = nc_open(fnetcdf)

        ################################################
        # get dimensions
        nx = nc$dim$x$len
        ny = nc$dim$y$len
        
        timescale = nc$dim$timescale$vals
        ntimescale = length(timescale)
        
        #loop over all grid points and seasons to retrieve the variable of interest
        for(itimescale in 1:ntimescale){
          ts = timescale[itimescale]
          print(paste0("Processing temporal scale ", ts,"hr"))
        
          # loop only where there is data
          grid_annual_max = ncvar_get(nc, "maxpr", start = c(1, 1, itimescale, 1),
                                      count = c(-1, -1, 1, -1))
          index_notna = which(!is.na(grid_annual_max[,,5]),arr.ind = T)
          npixels = nrow(index_notna)
          
          pb <- txtProgressBar(min=1, max=npixels, style=3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress=progress)
          
          # loop over all grid points
          mat_boot=foreach(ipixel=1:npixels, .packages = c('extRemes','ncdf4','boot'), 
                          .options.snow=opts, .combine = 'rbind'
          ) %dopar% {
            ix = index_notna[ipixel,1]
            iy = index_notna[ipixel,2]
            vecmax = grid_annual_max[ix,iy,]
            
            vout = rep(NA,nboot+1)
            if(!all(is.na(vecmax))){
              vecmax = vecmax[!is.na(vecmax)]
              has_zero = any(vecmax==0)
              
              if(!has_zero){
                b <- boot(vecmax, bootFunc, R = nboot)
                vout = c(b$t0,b$t)
              }
            }
            vout
          }
          
          # rearrange in a matrix
          out_best = matrix(nrow=ny,ncol=nx)
          out_boot = array(dim=c(ny,nx,nboot))
          for(ipixel in 1:npixels){
              ix = index_notna[ipixel,1]
              iy = index_notna[ipixel,2]
              out_best[iy,ix] = mat_boot[ipixel,1]
              out_boot[iy,ix,] = mat_boot[ipixel,2:101]
          }
                    
          # Write to a netcdf file
          filename=paste0("/jarre/gevin/EUROCORDEX/RL/20YR/matrix_20yrRL_",
                          ts, "hr_", spatialscale,"_",season,"_", model, ".nc")

          # define dimensions
          dimlon <- ncdim_def("rlon", "rotated longitude", nc$dim$x$vals)
          dimlat <- ncdim_def("rlat", "rotated latitude", nc$dim$y$vals)
          dimboot <- ncdim_def("boot", "bootstrap sample", 1:100)

          # add best estimate
          var_best <- ncvar_def("RL20_best", "mm", list(dimlat, dimlon), longname="RL20 bestestimates", -9999)
          var_boot <- ncvar_def("RL20_boot", "mm", list(dimlat, dimlon, dimboot), longname="RL20 bootstrap estimates", -9999)
          ncnew <- nc_create(filename, list(var_best,var_boot))

          # add data
          ncvar_put(ncnew, var_best, out_best, start=c(1,1), count=c(ny,nx))
          ncvar_put(ncnew, var_boot, out_boot, start=c(1,1,1), count=c(ny,nx,nboot))

          # Don't forget to close the file
          nc_close(ncnew)
          #########################################################################################
        }
        nc_close(nc)
      }
    }
  }
}


##############################  OBS   ##########################################
name_obs =c('OBS_KNMI_011EUi_1998-2025','OBS_CPC_011EUi_2005-2025',"OBS_INCA_011EUi_2013-2025",
  'OBS_CEH-GEAR_011EUi_1990-2016', 'OBS_KLIMGRID_011EUi_2011-2019','OBS_COMEPHORE_011EUi_1997-2025',
  'OBS_RADKLIM_011EUi_2001-2025', 'OBS_GRIPHO_011EUi_2001-2016')

if(TRUE){
  # first loop over the models
  for(obs in name_obs){
    print(paste0("Processing obs ", obs))
    
    for(season in c("DJF", "MAM", "JJA", "SON")){
      print(paste0("Processing season ", season))
      for(spatialscale in c("S1","S2","S3","S4","S5")){
        print(paste0("Processing spatial scale ", spatialscale))
        fnetcdf = paste0("/jarre/gevin/EUROCORDEX/DATA/OBS/",obs,'_',
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
          
          # loop only where there is data
          grid_annual_max = ncvar_get(nc, "maxpr", start = c(1, 1, itimescale, 1),
                                      count = c(-1, -1, 1, -1))
          index_notna = which(!is.na(grid_annual_max[,,5]),arr.ind = T)
          npixels = nrow(index_notna)
          
          pb <- txtProgressBar(min=1, max=npixels, style=3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress=progress)
          
          # loop over all grid points
          mat_boot=foreach(ipixel=1:npixels, .packages = c('extRemes','ncdf4','boot'), 
                          .options.snow=opts, .combine = 'rbind'
          ) %dopar% {
            ix = index_notna[ipixel,1]
            iy = index_notna[ipixel,2]
            vecmax = grid_annual_max[ix,iy,]
            
            vout = rep(NA,nboot+1)
            if(!all(is.na(vecmax))){
              vecmax = vecmax[!is.na(vecmax)]
              has_zero = any(vecmax==0)
              
              if(!has_zero){
                b <- boot(vecmax, bootFunc, R = nboot)
                vout = c(b$t0,b$t)
              }
            }
            vout
          }
          
          # rearrange in a matrix
          out_best = matrix(nrow=ny,ncol=nx)
          out_boot = array(dim=c(ny,nx,nboot))
          for(ipixel in 1:npixels){
              ix = index_notna[ipixel,1]
              iy = index_notna[ipixel,2]
              out_best[iy,ix] = mat_boot[ipixel,1]
              out_boot[iy,ix,] = mat_boot[ipixel,2:101]
          }
                    
          # Write to a netcdf file
          filename=paste0("/jarre/gevin/EUROCORDEX/RL/20YR/matrix_20yrRL_",
                          ts, "hr_", spatialscale,"_",season,"_", obs, ".nc")

          # define dimensions
          dimlon <- ncdim_def("lon", "longitude", nc$dim$x$vals)
          dimlat <- ncdim_def("lat", "latitude", nc$dim$y$vals)
          dimboot <- ncdim_def("boot", "bootstrap sample", 1:100)

          # add best estimate
          var_best <- ncvar_def("RL20_best", "mm", list(dimlat, dimlon), longname="RL20 bestestimates", -9999)
          var_boot <- ncvar_def("RL20_boot", "mm", list(dimlat, dimlon, dimboot), longname="RL20 bootstrap estimates", -9999)
          ncnew <- nc_create(filename, list(var_best,var_boot))

          # add data
          ncvar_put(ncnew, var_best, out_best, start=c(1,1), count=c(ny,nx))
          ncvar_put(ncnew, var_boot, out_boot, start=c(1,1,1), count=c(ny,nx,nboot))

          # Don't forget to close the file
          nc_close(ncnew)
        }
        nc_close(nc)
      }
    }
  }
}

# Fermer le cluster
stopCluster(cl)