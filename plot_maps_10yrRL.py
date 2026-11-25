import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import os
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
cartopy.config['data_dir'] = '/exec/gevin/'

# Compare 10-year return levels at an annual scale from various
# reference products (Poschlo et al. 2021) to the 10-year return levels
# obtained from the regional simulations

# nc files provided by Ben Poschlod 04/05/2026

### Poschlod, B., Ludwig, R., and Sillmann, J.: Ten-year return levels of 
# sub-daily extreme precipitation over Europe, Earth Syst. Sci. Data, 13, 
# 983–1003, https://doi.org/10.5194/essd-13-983-2021, 2021

# folder with the 10-yr RL
folder = "/jarre/gevin/EUROCORDEX/RL/10YR/"

# vector of time scales
vec_time_scale = ['1','6','24']


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

for type_diff in ["abs","rel"]:
    #_________________________________________________________________
    # Start figure here
    #_________________________________________________________________

    # Initialiser la figure et les axes
    fig = plt.figure(figsize=(6.5,8.5))

    gs = fig.add_gridspec(3, 1, hspace=0.25)  # spacing between the two groups

    gs0 = gs[0].subgridspec(2, 6, hspace=0.2, wspace=0.05)  # spacing within the groups
    gs1 = gs[1].subgridspec(2, 6, hspace=0.2, wspace=0.05)
    gs2 = gs[2].subgridspec(2, 6, hspace=0.2, wspace=0.05)

    for it, time_scale in enumerate(vec_time_scale):
        if(it == 0):
            subaxes = gs0.subplots(subplot_kw={'projection': projplot})
            ax_title = fig.add_subplot(gs0[:])
        elif(it == 1):
            subaxes = gs1.subplots(subplot_kw={'projection': projplot})
            ax_title = fig.add_subplot(gs1[:])
        else:
            subaxes = gs2.subplots(subplot_kw={'projection': projplot})
            ax_title = fig.add_subplot(gs2[:])

        # Add ghost axes and titles on gs1 and gs2
        ax_title.axis('off')
        ax_title.set_title(f"{time_scale}-hour time scale", fontsize=12, pad=10.0, fontweight='bold')

        # load dsRef
        fRef = f"DS10RL_REF_{time_scale}h.nc"
        dsRef = xr.load_dataset(folder+fRef, engine="netcdf4")

        for im, model in enumerate(vec_model):
            ax = subaxes[im//6, im%6]

            lab_sim = vec_label[vec_model.index(model)]

            fmodel = f"DS10RL_{model}_withHershCorr_{time_scale}h.nc"
            dsModel = xr.load_dataset(folder+fmodel, engine="netcdf4")

            # compute relative difference
            lon = dsModel.lon.values
            lat = dsModel.lat.values
            Zref = dsRef.RL.values
            Zmodel = dsModel.RL.values
            if(type_diff=="rel"):
                Zdiff = (Zmodel - Zref)/Zref * 100
                plot = ax.pcolormesh(lon, lat, Zdiff, cmap=cmaps.precip_diff_12lev, vmin=-65, vmax=65,
                                transform=projorig, zorder=1)
            else:
                Zdiff = Zmodel - Zref
                plot = ax.pcolormesh(lon, lat, Zdiff, cmap=cmaps.precip_diff_12lev, vmin=-32.5, vmax=32.5,
                                transform=projorig, zorder=1)
            
            # draw coastlines and borders
            ax.set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
            ax.add_feature(cartopy.feature.OCEAN, color='white', zorder=0)
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
            ax.add_feature(fine_borders, edgecolor='black', linewidth=0.5)

            # Optional: Add coastlines at 10m
            ax.coastlines(resolution='10m', linewidth=0.5)

            # add model name
            ax.set_title(lab_sim, fontsize=8, pad=1.0)
        
    # Add a colorbar axis at the bottom of the graph
    cbar_ax = fig.add_axes([0.15, 0.06, 0.7, 0.02]) #[left, bottom, width, height]
    cbar = fig.colorbar(plot, cax=cbar_ax,orientation='horizontal',aspect=20,shrink=0.85, pad=0, extend = "both")
    if type_diff == "rel":
        cbar.set_label("Relative difference (%)", fontsize=10)
    else:
        cbar.set_label("Absolute difference (mm)", fontsize=10)

    cbar.ax.tick_params(labelsize=10)

    # save the figure
    plt.savefig(f'/exec/gevin/evaluation_precip_extremes_eurocordex/PAPER/maps_10yrRL_{type_diff}diff.png', dpi=300, bbox_inches='tight')
    plt.close()