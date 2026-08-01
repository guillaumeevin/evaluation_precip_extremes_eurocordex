import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import xesmf as xe
import numpy as np
np.float_ = np.float64
from metpy.units import units
import matplotlib.colors as mcolors
from clisops.core import subset
from metpy.plots import ctables
import matplotlib.ticker as ticker
import cmaps
import geopandas as gpd
import seaborn as sns
from shapely.geometry import Point, Polygon
import rioxarray
from math import nan
from scipy import stats
import statsmodels.api as sm

# function that estimates simple scaling parameters
#
# INPUT: maxSeas is a four-dimension array containing the seasonal maxima
# - dim 1: years
# - dim 2: temporal scale from 1h to 72h (6 agg.)
# - dim 3: longitude
# - dim 4: latitude
#
# OUTPUT: return a matrix [lon,lat] with scaling parameters between 0 and 1
def getScalingParameter(maxSeas):
    # pas de temps
    vec_pdt = [1,3,6,12,24,72]
    n_pdt = len(vec_pdt)

    # ordre moments
    vec_order = [i for i in np.arange(0.5, 4.5, 0.5)]
    n_order = len(vec_order)

    # get number of longitude and latitude points
    nLon = maxSeas.shape[2]
    nLat = maxSeas.shape[3]

    # STEP 0: rescale maxima to get the same unit for all pdt (mm/hour)
    for i_pdt in range(n_pdt):
        pdt = vec_pdt[i_pdt]
        
        # fill array
        maxSeas[:,i_pdt,:,:] = maxSeas[:,i_pdt,:,:] / pdt

    # matrix of moments for different order
    moment_mat = np.empty((nLon,nLat,n_pdt,n_order))
    scalingPar = np.empty((nLon,nLat))
    for iLon in range(nLon):
        for iLat in range(nLat):
        
            # STEP 1: compute moments for differents ordres
            for i_pdt in range(n_pdt):
                # pas de temps
                pdt = vec_pdt[i_pdt]
                
                # max pour tous les mois et tous les ans
                max_vec = maxSeas[:,i_pdt,iLon,iLat]
                
                # moments non-centre pour differents ordres: E[I^q(D,0)]
                for i_order in range(n_order):
                    order = vec_order[i_order]
                    moment_mat[iLon,iLat,i_pdt,i_order] = np.mean(max_vec**order)
                
            # STEP 2: apply linear regressions log(E_h^k) ~ a + kq*log(h) for each order
            kq = []
            for i_order in range(n_order):
                order = vec_order[i_order]
            
                # log of the moments and log of the pas de temps
                x = np.log(vec_pdt)
                y = np.log(moment_mat[iLon,iLat,:,i_order])

                # Perform regression
                reg_linear = stats.linregress(x, y)
                kq.append(reg_linear.slope)
            
            # STEP 3: apply linear regression  kq ~ k
            x = vec_order
            y = kq
            reg_linear = sm.OLS(y, x).fit()
            scalingPar[iLon,iLat] = -reg_linear.params

    return scalingPar

#______________________________________________________________________________________________________________
# estimate scaling parameters from observations
#______________________________________________________________________________________________________________
folder_obs = 'D:/EUROCORDEX_extremes/DATA/OBS/'
name_obs = ['OBS_COMEPHORE_011EUi_1997-2022','OBS_RADKLIM_011EUi_2001-2022', 'OBS_GRIPHO_011EUi_2001-2016']

if(False):
    for season in ["JJA","SON","DJF","MAM"]:
        for obs in name_obs:
        
            fObs = f'{obs}_{season}_max_scale_S1.nc'
            dsObs = xr.load_dataset(folder_obs+fObs, engine="netcdf4")

            scalingPar = getScalingParameter(dsObs.maxpr.values)
            dsObs = dsObs.cf.add_bounds(["lon","lat"])

            dsObs_cor = dsObs.drop_vars("maxpr")

            dsObs_cor["scalingPar"] = (("y", "x"), scalingPar)

            dsObs_cor.to_netcdf(f'D:/EUROCORDEX_extremes/DATA/SIMPLESCALING/simplescaling_{obs}_{season}.nc')


#______________________________________________________________________________________________________________
# estimate scaling parameters from simulations
#______________________________________________________________________________________________________________
folder_model = 'D:/EUROCORDEX_extremes/DATA/SIM/'

# reference file for the regridding
fMod = 'ERA5_WRF451Q1980-2020_SON_all_24h_perc.nc'
dsMod = xr.load_dataset(folder_model+fMod, engine="netcdf4")
rlat_424_412 = dsMod.rlat.values
rlon_424_412 = dsMod.rlon.values

# list models and obs
name_models = ['ERA5_ALARO1-SFX','ERA5_CCLM6-0-1','ERA5_CNRM-ALADIN64E1',
              'ERA5_HCLIM43-ALADIN','ERA5_ICON-CLM-202407-1-1','ERA5_RACMO23E',
              'ERA5_RegCM5-0','ERA5_REMO2020-2-2','ERA5_WRF451Q']

for season in ["JJA","SON","DJF","MAM"]:
  
  # loop over model
  for model in name_models:
    print("process data for model " + model + " and season " + season)
  
    fMod = f'{model}1980-2020_{season}_max_scale_S1.nc'
    dsMod = xr.load_dataset(folder_model+fMod, engine="netcdf4")

    if model == 'ERA5_RegCM5-0':
        dsMod = dsMod.assign_coords(rlat=rlat_424_412, rlon=rlon_424_412)

    scalingPar = getScalingParameter(dsMod.maxpr.values)

    dsMod = dsMod.cf.add_bounds(["lon","lat"])

    dsMod_cor = dsMod.drop_vars("maxpr")

    dsMod_cor["scalingPar"] = (("y", "x"), scalingPar)

    dsMod_cor.to_netcdf(f'D:/EUROCORDEX_extremes/DATA/SIMPLESCALING/simplescaling_{model}_{season}.nc')