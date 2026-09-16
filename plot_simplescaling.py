import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import os
os.environ['ESMFMKFILE'] = 'C:/Users/eving/AppData/Local/miniconda3/envs/evalprecextremes/Library/lib/esmf.mk'
import xesmf as xe
import numpy as np
np.float_ = np.float64
import matplotlib.colors as mcolors
from clisops.core import subset
from metpy.plots import ctables
import matplotlib.ticker as ticker
import cmaps
import geopandas as gpd
import seaborn as sns
import rioxarray
from math import nan
from scipy import stats
import statsmodels.api as sm
import cartopy

###############################################################################
# Start computation for the plot 
###############################################################################

# Relative biases for 10-year return levels
# list models and obs
name_models = ['ERA5_CNRM-ALADIN64E1','ERA5_ALARO1-SFX','ERA5_CCLM6-0-1',
            'ERA5_HCLIM43-ALADIN','ERA5_ICON-CLM-202407-1-1','ERA5_RACMO23E',
            'ERA5_RegCM5-0','ERA5_REMO2020-2-2','ERA5_WRF451Q']

label_models = ['ALADIN64E1','ALARO1-SFX','CCLM6-0-1',
            'HCLIM43','ICON-CLM','RACMO23E',
            'RegCM5-0','REMO2020','WRF451Q']

# Original Projection
projorig_name  = "PlateCarree"
projorig_dict  = dict()
projorig       = ccrs.__dict__[projorig_name].__call__(**projorig_dict)

# Graphic projection
projplot_name       = "LambertConformal"
projplot_dict  = dict(central_longitude=10.5, central_latitude=49.5, standard_parallels=(33,45))
projplot       = ccrs.__dict__[projplot_name].__call__(**projplot_dict)

# vector of observational products
# list models and obs
vec_obs = ["OBS_CEH-GEAR_011EUi_1990-2016", "OBS_COMEPHORE_011EUi_1997-2022",
           "OBS_CPC_011EUi_2005-2022", "OBS_GRIPHO_011EUi_2001-2016",
           "OBS_KLIMGRID_011EUi_2011-2019", "OBS_RADKLIM_011EUi_2001-2022",
           "OBS_KNMI_011EUi_1998-2022"]
zorder_obs = [7, 6, 1, 3, 4, 5, 2]  # Define the desired order of the legend entries


###############################################################################
# First plot: simple scaling estimates from obs. products
###############################################################################

folder = "D:/EUROCORDEX_extremes/DATA/SIMPLESCALING/"

# Initialiser la figure et les axes
fig, axes = plt.subplots(nrows=2, ncols=2, subplot_kw={'projection': projplot},figsize=(8,9))

vec_season = ["JJA","MAM","SON","DJF"]

discrete_cmap = plt.colormaps['terrain'].resampled(8)


for ax, season in zip(axes.flat, vec_season):
    for obs in vec_obs:
        # loadtxt is fast for simple numeric files
        fObs = f"simplescaling_{obs}_{season}.nc"
        dsObs = xr.load_dataset(folder+fObs, engine="netcdf4")

        X = dsObs.lon.values
        Y = dsObs.lat.values
        Z = dsObs.scalingPar
        
        plot = ax.pcolormesh(X, Y, Z, cmap=discrete_cmap, vmin=0.2, vmax=1,
                            transform=projorig, zorder=zorder_obs[vec_obs.index(obs)])

    # draw coastlines and borders
    ax.set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.OCEAN, color='white', zorder=10)
    ax.add_feature(cartopy.feature.LAND, color='lightgray',zorder=0,
        linewidth=0.5, edgecolor='black')
        
    # Add fine resolution borders
    # Category: 'cultural', Name: 'admin_0_boundary_lines_land', Scale: '10m'
    fine_borders = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_0_boundary_lines_land',
        scale='10m',
        facecolor='none'
    )
    ax.add_feature(fine_borders, edgecolor='black', linewidth=0.5,zorder=12)

    # Optional: Add coastlines at 10m
    ax.coastlines(resolution='10m', linewidth=0.5,zorder=12)

    # add model name
    ax.set_title(season, fontsize=14, pad=1.0)

    # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(right=0.9)
    
    # Add a colorbar axis at the bottom of the graph
    cbar_ax = fig.add_axes([0.92, 0.25, 0.03, 0.5]) #[left, bottom, width, height]
    cbar = fig.colorbar(plot, cax=cbar_ax,orientation='vertical',shrink=0.85, pad=0.05, extend = "both")
    cbar.set_label("Scaling parameter (-)", fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

# save the figure
plt.savefig(f'D:/EUROCORDEX_extremes/PAPER/simplescaling_obs.png', dpi=300, bbox_inches='tight')

###############################################################################
# Second plot: Absolute differents between simulated and observed scale parameters
###############################################################################
folder = "D:/EUROCORDEX_extremes/DATA/SIMPLESCALING/"

for season in ["JJA","MAM","SON","DJF"]:

    # Initialiser la figure et les axes
    fig, axes = plt.subplots(nrows=3, ncols=3, subplot_kw={'projection': projplot},figsize=(8,9))

    for ax, model in zip(axes.flat, name_models):
        # loadtxt is fast for simple numeric files
        fMod = f"simplescaling_{model}_{season}.nc"
        dsMod = xr.load_dataset(folder+fMod, engine="netcdf4")
        dsMod_clean = xr.Dataset(
            data_vars={"scalingPar": (("y", "x"), dsMod.scalingPar.values)},
            coords={
                "x": dsMod.x,
                "y": dsMod.y,
                "lon": (("y", "x"), dsMod.lon.values),
                "lat": (("y", "x"), dsMod.lat.values),
            }
        )
        dsMod_clean = dsMod_clean.cf.add_bounds(["lon","lat"])

        for obs in vec_obs:
            # loadtxt is fast for simple numeric files
            fObs = f"simplescaling_{obs}_{season}.nc"
            dsObs = xr.load_dataset(folder+fObs, engine="netcdf4")

            # 1. Define your 1D coordinate vectors
            lat_vec = dsObs.lat.values
            lon_vec = dsObs.lon.values

            # 2. Create the 2D coordinate matrices
            # Use indexing='ij' for matrix/geospatial indexing (Rows=Lat, Cols=Lon)
            lon, lat = np.meshgrid(lon_vec, lat_vec, indexing='xy')

            dsObs_clean = xr.Dataset(
                data_vars={"scalingPar": (("y", "x"), dsObs.scalingPar.values)},
                coords={
                    "x": dsObs.x.values,
                    "y": dsObs.y.values,
                    "lon": (("y", "x"), lon),
                    "lat": (("y", "x"), lat),
                }
            )
            dsObs_clean = dsObs_clean.cf.add_bounds(["lon","lat"])

            # regrid obs
            regridder = xe.Regridder(
                dsMod_clean, dsObs_clean, "conservative", unmapped_to_nan=True, ignore_degenerate=True
            )
            dsMod_regrid = regridder(dsMod_clean, keep_attrs=True)

            X = dsMod_regrid.lon
            Y = dsMod_regrid.lat
            ZOBS = dsObs_clean.scalingPar
            ZMOD = dsMod_regrid.scalingPar
            
            Z = ZMOD - ZOBS
            plot = ax.pcolormesh(X, Y, Z, cmap=cmaps.precip_diff_12lev, vmin=-0.26, vmax=0.26,
                                transform=projorig, zorder=zorder_obs[vec_obs.index(obs)])

        # draw coastlines and borders
        ax.set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
        ax.add_feature(cartopy.feature.OCEAN, color='white', zorder=10)
        ax.add_feature(cartopy.feature.LAND, color='lightgray',zorder=0,
            linewidth=0.5, edgecolor='black')
            
        # Add fine resolution borders
        # Category: 'cultural', Name: 'admin_0_boundary_lines_land', Scale: '10m'
        fine_borders = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_0_boundary_lines_land',
            scale='10m',
            facecolor='none'
        )
        ax.add_feature(fine_borders, edgecolor='black', linewidth=0.5, zorder=12)

        # Optional: Add coastlines at 10m
        ax.coastlines(resolution='10m', linewidth=0.5, zorder=12)

        # add model name
        ax.set_title(label_models[name_models.index(model)], fontsize=14, pad=1.0)

        # Adjust the location of the subplots on the page to make room for the colorbar
        fig.subplots_adjust(right=0.9)
        
        # Add a colorbar axis at the bottom of the graph
        cbar_ax = fig.add_axes([0.92, 0.25, 0.03, 0.5]) #[left, bottom, width, height]
        cbar = fig.colorbar(plot, cax=cbar_ax,orientation='vertical',shrink=0.85, pad=0.05, extend = "both")
        cbar.set_label("Absolute difference (-)", fontsize=12)
        cbar.ax.tick_params(labelsize=12)
        fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

    # save the figure
    plt.savefig(f'D:/EUROCORDEX_extremes/PAPER/simplescaling_{season}.png', dpi=300, bbox_inches='tight')