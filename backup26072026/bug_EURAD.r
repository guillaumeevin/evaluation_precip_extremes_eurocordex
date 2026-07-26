# Guillaume Evin
# 01/05/2026
# This script estimates 10-year return level from different RCM simulations and compares them to the return level estimated from observations.

library(extRemes)
library(ncdf4)

#################################################################################
# EURADCLIM
#obs = 'OBS_EURADCLIM_011EUi_2013-2022'
obs = 'OBS_RADKLIM_011EUi_2001-2022'

season = "JJA"
fnetcdf = paste0("D:/EUROCORDEX_extremes/DATA/OBS/", obs, "_", season, "_max_scale_S1.nc")
nc = nc_open(fnetcdf)

# get dimensions
nx = nc$dim$x$len
ny = nc$dim$y$len
timescale = nc$dim$timescale$vals
nyear = length(nc$dim$year$vals)

array_ratio_3to1 = array(dim=c(nx,ny,nyear))

# loop over all grid points
for(ix in 1:nx){
  for(iy in 1:ny){
      max1h = ncvar_get(nc, "maxpr", start = c(ix, iy, 1, 1), count = c(1, 1, 1, -1))
      max3h = ncvar_get(nc, "maxpr", start = c(ix, iy, 2, 1), count = c(1, 1, 1, -1))
      array_ratio_3to1[ix,iy,] = max3h/max1h
  }
}

mean(array_ratio_3to1>1,na.rm=T)
mean(array_ratio_3to1>3,na.rm=T)
arrind = which(array_ratio_3to1>3,arr.ind = T)
unique(arrind[,3])



array_ratio_6to3 = array(dim=c(nx,ny,nyear))

# loop over all grid points
for(ix in 1:nx){
  for(iy in 1:ny){
    max1h = ncvar_get(nc, "maxpr", start = c(ix, iy, 2, 1), count = c(1, 1, 1, -1))
    max3h = ncvar_get(nc, "maxpr", start = c(ix, iy, 3, 1), count = c(1, 1, 1, -1))
    array_ratio_6to3[ix,iy,] = max3h/max1h
  }
}

mean(array_ratio_6to3>1,na.rm=T)
mean(array_ratio_6to3>2,na.rm=T)
