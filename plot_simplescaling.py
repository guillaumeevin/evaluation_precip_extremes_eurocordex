import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import os
os.environ['ESMFMKFILE'] = '/exec/gevin/.conda/envs/myxscen/lib/esmf.mk'
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
from shapely.geometry import mapping
cartopy.config['data_dir'] = '/exec/gevin/'

###############################################################################
# Start computation for the plot 
###############################################################################

# list models and labels
vec_model = ["ERA5_ALARO1-SFX","ERA5_CNRM-ALADIN64E1","ERA5_HCLIM43-ALADIN","ERA5_CCLM6-0-1",
             "ERA5_GCOAST-AHOIB1-1", "ERA5_ICON-CLM-202407-1-1",
             "ERA5_RACMO23E","ERA5_RegCM5-0","ERA5_REMO2020-2-2",
             "ERA5_ROAM-NBS","ERA5_TSMP1-140-E","ERA5_WRF451Q-r1"]
vec_label = ["ALARO1-SFX","ALADIN64E1","HCLIM43","CCLM6-0-1",
             "GC-AHOIB1", "ICON-CLM", "RACMO23E","RegCM5-0","REMO2020-2-2",
             "ROAM-NBS","TSMP1-140-E","WRF451Q-r1"]

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
vec_highres_obs = ["OBS_CEH-GEAR_011EUi_1990-2016", "OBS_COMEPHORE_011EUi_1997-2022",
           "OBS_CPC_011EUi_2005-2022", "OBS_GRIPHO_011EUi_2001-2016",
           "OBS_KLIMGRID_011EUi_2011-2019", "OBS_RADKLIM_011EUi_2001-2022",
           "OBS_KNMI_011EUi_1998-2022","OBS_INCA_011EUi_2013-2022"]
namecountry_highres_obs = ["United Kingdom","France","Switzerland","Italy","Denmark",
                   "Germany","Netherlands","Austria"]

# folder for simple scaling
folder_SS = "/jarre/gevin/EUROCORDEX/SIMPLESCALING/"

# season
vec_season = ["JJA","MAM","SON","DJF"]

# https://ec.europa.eu/eurostat/web/gisco/geodata/administrative-units/countries
folder_main = "/exec/gevin/evaluation_precip_extremes_eurocordex/"
df_border = gpd.read_file(folder_main+"DATA_FIG/CNTR_RG_20M_2024_4326.geojson")


###############################################################################
# First step: format outputs from obs. (merge high res datasets) 
# and model (reproject to obs)
###############################################################################
if(False):
    # Step 1.1: merge high-resolution datasets
    for season in vec_season:
        #------------- Crop and merge datasets
        for i_obs, obs in enumerate(vec_highres_obs):
            # read raw data and format to be clipped
            fHR = f'{folder_SS}simplescaling_{obs}_{season}.nc'
            dsHR = xr.load_dataset(fHR, engine="netcdf4")

            dsHR_clean = xr.Dataset(
            data_vars={"scalingPar": (("lat", "lon"), dsHR.scalingPar.values)},
            coords={
                "lon": dsHR.x.values,
                "lat": dsHR.y.values,
            }
            )

            dsHR_clean.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
            dsHR_clean.rio.write_crs("epsg:4326", inplace=True)

            # for the first product we retrieve he raw data
            if i_obs==0:
                scalingPar_Merge = dsHR_clean.scalingPar.values

            # border for the corresponding country
            border_country = df_border[df_border["NAME_ENGL"] == namecountry_highres_obs[i_obs]]

            # clip
            dsClip = dsHR_clean.rio.clip(border_country.geometry.apply(mapping), border_country.crs, drop=False)

            scalingPar_country = dsClip.scalingPar.values
            mask2add = ~np.isnan(scalingPar_country)
            scalingPar_Merge[mask2add] = scalingPar_country[mask2add]

        #------------------ Define a clean netcdf file
        lat_vec = dsHR.lat.values
        lon_vec = dsHR.lon.values

        # 2. Create the 2D coordinate matrices
        # Use indexing='ij' for matrix/geospatial indexing (Rows=Lat, Cols=Lon)
        lon, lat = np.meshgrid(lon_vec, lat_vec, indexing='xy')

        dsObs_clean = xr.Dataset(
            data_vars={"scalingPar": (("y", "x"), scalingPar_Merge)},
            coords={
                "x": dsHR.x.values,
                "y": dsHR.y.values,
                "lon": (("y", "x"), lon),
                "lat": (("y", "x"), lat),
            }
        )
        dsObs_clean = dsObs_clean.cf.add_bounds(["lon","lat"])
        fobs = f"simplescaling_HIGHRES_{season}.nc"
        dsObs_clean.to_netcdf(folder_SS+fobs, engine="netcdf4")

    # Step 1.2: reproject for models
    for season in vec_season:

        for imodel, model in enumerate(vec_model):
            # loadtxt is fast for simple numeric files
            fMod = f"simplescaling_{model}_{season}.nc"
            dsMod = xr.load_dataset(folder_SS+fMod, engine="netcdf4")
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

            # load high-res obs
            fobs = f"simplescaling_HIGHRES_{season}.nc"
            dsObs = xr.load_dataset(folder_SS+fobs, engine="netcdf4")

            # regrid obs
            regridder = xe.Regridder(
                dsMod_clean, dsObs, "conservative", unmapped_to_nan=True, ignore_degenerate=True
            )
            dsMod_regrid = regridder(dsMod_clean, keep_attrs=True)

            fMod_clean = f"simplescaling_{model}_{season}_reproj.nc"
            dsMod_regrid.to_netcdf(folder_SS+fMod_clean, engine="netcdf4")

###############################################################################
# First plot: simple scaling estimates from obs. products
###############################################################################
if(True):
    # Initialiser la figure et les axes
    fig, axes = plt.subplots(nrows=2, ncols=2, subplot_kw={'projection': projplot},figsize=(8,9))

    discrete_cmap = plt.colormaps['terrain'].resampled(8)

    for ax, season in zip(axes.flat, vec_season):
        # loadtxt is fast for simple numeric files
        fObs = f"simplescaling_HIGHRES_{season}.nc"
        dsObs = xr.load_dataset(folder_SS+fObs, engine="netcdf4")

        X = dsObs.lon.values
        Y = dsObs.lat.values
        Z = dsObs.scalingPar
        
        plot = ax.pcolormesh(X, Y, Z, cmap=discrete_cmap, vmin=0.2, vmax=1,
                            transform=projorig, zorder=1)

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
    plt.savefig('/exec/gevin/evaluation_precip_extremes_eurocordex/PAPER/simplescaling_obs.png', dpi=300, bbox_inches='tight')

###############################################################################
# Second plot: Absolute differents between simulated and observed scale parameters
###############################################################################
if(False):
    for season in vec_season:
        # load high-res obs
        fobs = f"simplescaling_HIGHRES_{season}.nc"
        dsObs = xr.load_dataset(folder_SS+fobs, engine="netcdf4")

        # Initialiser la figure et les axes
        fig, axes = plt.subplots(nrows=3, ncols=4, subplot_kw={'projection': projplot},figsize=(6.5,6))

        for ax, model in zip(axes.flat, vec_model):
            fMod_clean = f"simplescaling_{model}_{season}_reproj.nc"
            dsMod = xr.load_dataset(folder_SS+fMod_clean, engine="netcdf4")
        
            X = dsMod.lon
            Y = dsMod.lat
            ZOBS = dsObs.scalingPar
            ZMOD = dsMod.scalingPar
            
            Z = ZMOD - ZOBS
            plot = ax.pcolormesh(X, Y, Z, cmap=cmaps.precip_diff_12lev, vmin=-0.26, vmax=0.26,
                                transform=projorig, zorder=1)

            # draw coastlines and borders
            ax.set_extent([-5, 20, 35, 60], crs=ccrs.PlateCarree())
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
            ax.set_title(vec_label[vec_model.index(model)], fontsize=10, pad=1.0)
            
            # Add a colorbar axis at the bottom of the graph
            cbar_ax = fig.add_axes([0.15, 0.06, 0.7, 0.02]) #[left, bottom, width, height]
            cbar = fig.colorbar(plot, cax=cbar_ax,orientation='horizontal',aspect=20,shrink=0.85, pad=0, extend = "both")
            cbar.set_label("Absolute difference (-)", fontsize=10)
            cbar.ax.tick_params(labelsize=12)


        # save the figure
        plt.savefig(f'/exec/gevin/evaluation_precip_extremes_eurocordex/PAPER/simplescaling_{season}.png', dpi=300, bbox_inches='tight')