import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import os
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
cartopy.config['data_dir'] = '/exec/gevin/'

from utils import unrot_lon, unrot_lat

# Compare 10-year return levels at an annual scale from various
# reference products (Poschlo et al. 2021) to the 10-year return levels
# obtained from the regional simulations

# nc files provided by Ben Poschlod 04/05/2026

### Poschlod, B., Ludwig, R., and Sillmann, J.: Ten-year return levels of 
# sub-daily extreme precipitation over Europe, Earth Syst. Sci. Data, 13, 
# 983–1003, https://doi.org/10.5194/essd-13-983-2021, 2021

# main folder
folder = "/jarre/gevin/EUROCORDEX/RL/10YR/"

# Original Projection
projorig_name  = "PlateCarree"
projorig_dict  = dict()
projorig       = ccrs.__dict__[projorig_name].__call__(**projorig_dict)

# Graphic projection
projplot_name       = "LambertConformal"
projplot_dict  = dict(central_longitude=10.5, central_latitude=49.5, standard_parallels=(33,45))
projplot       = ccrs.__dict__[projplot_name].__call__(**projplot_dict)

# limits for the plots
levels = {}
levels['1h'] = np.linspace(12, 42, 16)
levels['6h'] = np.linspace(30, 80, 11)
levels['24h'] = np.linspace(40, 210, 18)

# Add fine resolution borders
# Category: 'cultural', Name: 'admin_0_boundary_lines_land', Scale: '10m'
fine_borders = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_0_boundary_lines_land',
    scale='10m',
    facecolor='none'
)

# vector of time scales
vec_time_scale = ['1','6','24']

# type of difference
for type_diff in ["abs","rel"]:

    ############################################################################
    # Figure: 10-yr return levels for 1h, 6h, and 24h durations
    # each row is a temporal aggregation, each column is a reference product
    # - col 1: Reference 10-year return levels from Poschlod et al. (2021)
    # - col 2: High-resolution products (CPC, CEH-GEAR, COMEPHORE, GRIPHO, KLIMGRID, RADKLIM, KNMI)
    # - col 3: ERA5
    # - col 4: EMO
    # - col 5: CERRA
    ############################################################################

    fig, axes = plt.subplots(figsize=(8.5, 7), nrows = 3, ncols = 5, subplot_kw={'projection': projplot}, constrained_layout=True)

    for i_time_scale, time_scale in enumerate(vec_time_scale):
        ax = axes[i_time_scale,0]

        # retrieve reference return levels
        fobs = f"DS10RL_REF_{time_scale}h.nc"
        dsRef = xr.load_dataset(folder+fobs, engine="netcdf4")

        lon = dsRef.lon.values
        lat = dsRef.lat.values
        Zref = dsRef.RL.values
        
        # plot the reference return levels
        limits = levels[str(time_scale)+'h']
        ax.set_global()
        ax.set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
        ax.add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
        ax.add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                        linewidth=0.5, edgecolor='black')
        gls = ax.gridlines(linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
                    x_inline=False,y_inline=False,rotate_labels = False)
        gls.right_labels = False
        p0 = ax.contourf(lon, lat, Zref, cmap="YlGnBu", levels=limits, transform=projorig, extend = "both", zorder=1)
        ax.add_feature(fine_borders, edgecolor='black', linewidth=0.5, zorder=1.2)
        ax.coastlines(resolution='10m', linewidth=0.5, zorder=1.2)

        cbar = plt.colorbar(p0, ax=ax, orientation='horizontal', aspect=20, label = "(mm)")
        cbar.ax.locator_params(nbins=4)


        # add title
        if i_time_scale==0:
            ax.set_title("Reference", fontweight = "bold")

        #________________________________________________________________
        # i,1 Diff highres /ref
        #________________________________________________________________
        ax = axes[i_time_scale,1]
        ax.set_global()
        ax.set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
        ax.add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
        ax.add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                        linewidth=0.5, edgecolor='black')
        gls = ax.gridlines(linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
                    x_inline=False,y_inline=False,rotate_labels = False)
        gls.right_labels = False

        fobs = f"DS10RL_HIGHRES_withHershCorr_{time_scale}h.nc"
        dsObs = xr.load_dataset(folder+fobs, engine="netcdf4")
        Zmodel = dsObs.RL.values
        
        # compute relative difference
        if(type_diff=="rel"):
            Zdiff = (Zmodel - Zref)/Zref * 100
            p0 = ax.pcolormesh(lon, lat, Zdiff, cmap=cmaps.precip_diff_12lev, vmin=-65, vmax=65,
                            transform=projorig, zorder=1)
            label_colorbar = "(%)"
        else:
            Zdiff = Zmodel - Zref
            p0 = ax.pcolormesh(lon, lat, Zdiff, cmap=cmaps.precip_diff_12lev, vmin=-32.5, vmax=32.5,
                            transform=projorig, zorder=1)
            label_colorbar = "(mm)"
                
        ax.add_feature(fine_borders, edgecolor='black', linewidth=0.5,zorder=1.2)
        ax.coastlines(resolution='10m', linewidth=0.5,zorder=1.2)

        # add color bar at the bottom for all columns
        plt.colorbar(p0, ax=ax, extend = "both", orientation='horizontal', aspect=20, label = label_colorbar)

        # add title
        if i_time_scale==0:
            ax.set_title("High-res.", fontweight = "bold")

        #________________________________________________________________
        # i,1 Diff ERA5 /ref
        #________________________________________________________________
        ax = axes[i_time_scale,2]

        ax.set_global()
        ax.set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
        ax.add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
        ax.add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                        linewidth=0.5, edgecolor='black')
        gls = ax.gridlines(linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
                    x_inline=False,y_inline=False,rotate_labels = False)
        gls.right_labels = False

        # Relative biases for 10-year return levels
        # list models and obs
        obs = "OBS_ERA5_011EUi_1980-2020"
        fobs = f"DS10RL_{obs}_withHershCorr_{time_scale}h.nc"
        dsObs = xr.load_dataset(folder+fobs, engine="netcdf4")
        Zmodel = dsObs.RL.values
        
        # compute relative difference
        if(type_diff=="rel"):
            Zdiff = (Zmodel - Zref)/Zref * 100
            p0 = ax.pcolormesh(lon, lat, Zdiff, cmap=cmaps.precip_diff_12lev, vmin=-65, vmax=65,
                            transform=projorig, zorder=1)
            label_colorbar = "(%)"
        else:
            Zdiff = Zmodel - Zref
            p0 = ax.pcolormesh(lon, lat, Zdiff, cmap=cmaps.precip_diff_12lev, vmin=-32.5, vmax=32.5,
                            transform=projorig, zorder=1)
            label_colorbar = "(mm)"
            
        ax.add_feature(fine_borders, edgecolor='black', linewidth=0.5,zorder=1.2)
        ax.coastlines(resolution='10m', linewidth=0.5,zorder=1.2)

        # add color bar at the bottom for all columns
        plt.colorbar(p0, ax=ax, extend = "both", orientation='horizontal', aspect=20, label = label_colorbar)
        
        # add title
        if i_time_scale==0:
            ax.set_title("ERA5", fontweight = "bold")

        #________________________________________________________________
        # i,1 Diff EMO /ref
        #________________________________________________________________
        ax = axes[i_time_scale,3]


        if(i_time_scale>0):
            ax.set_global()
            ax.set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
            ax.add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
            ax.add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                            linewidth=0.5, edgecolor='black')
            gls = ax.gridlines(linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
                        x_inline=False,y_inline=False,rotate_labels = False)
            gls.right_labels = False

            # Relative biases for 10-year return levels
            # list models and obs
            obs = "OBS_EMO_011EUi_1990-2022"
            fobs = f"DS10RL_{obs}_withHershCorr_{time_scale}h.nc"
            dsObs = xr.load_dataset(folder+fobs, engine="netcdf4")
            Zmodel = dsObs.RL.values
            
            # compute relative difference
            if(type_diff=="rel"):
                Zdiff = (Zmodel - Zref)/Zref * 100
                p0 = ax.pcolormesh(lon, lat, Zdiff, cmap=cmaps.precip_diff_12lev, vmin=-65, vmax=65,
                                transform=projorig, zorder=1)
                label_colorbar = "(%)"
            else:
                Zdiff = Zmodel - Zref
                p0 = ax.pcolormesh(lon, lat, Zdiff, cmap=cmaps.precip_diff_12lev, vmin=-32.5, vmax=32.5,
                                transform=projorig, zorder=1)
                label_colorbar = "(mm)"
                
            ax.add_feature(fine_borders, edgecolor='black', linewidth=0.5,zorder=1.2)
            ax.coastlines(resolution='10m', linewidth=0.5,zorder=1.2)

            # add color bar at the bottom for all columns
            plt.colorbar(p0, ax=ax, extend = "both", orientation='horizontal', aspect=20, label=label_colorbar)
        else:
            ax.axis('off')

        # add title
        if i_time_scale==0:
            ax.set_title("EMO-1", fontweight = "bold")

        #________________________________________________________________
        # i,1 Diff CERRA /ref
        #________________________________________________________________
        ax = axes[i_time_scale,4]

        if(i_time_scale>1):
            ax.set_global()
            ax.set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
            ax.add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
            ax.add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                            linewidth=0.5, edgecolor='black')
            gls = ax.gridlines(linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
                        x_inline=False,y_inline=False,rotate_labels = False)
            gls.right_labels = False

            # Relative biases for 10-year return levels
            # list models and obs
            obs = "OBS_CERRALAND_011EUi_1986-2020"
            fobs = f"DS10RL_{obs}_withHershCorr_{time_scale}h.nc"
            dsObs = xr.load_dataset(folder+fobs, engine="netcdf4")
            Zmodel = dsObs.RL.values
            
            # compute relative difference
            if(type_diff=="rel"):
                Zdiff = (Zmodel - Zref)/Zref * 100
                p0 = ax.pcolormesh(lon, lat, Zdiff, cmap=cmaps.precip_diff_12lev, vmin=-65, vmax=65,
                                transform=projorig, zorder=1)
                label_colorbar = "(%)"
            else:
                Zdiff = Zmodel - Zref
                p0 = ax.pcolormesh(lon, lat, Zdiff, cmap=cmaps.precip_diff_12lev, vmin=-32.5, vmax=32.5,
                                transform=projorig, zorder=1)
                label_colorbar = "(mm)"
                
            ax.add_feature(fine_borders, edgecolor='black', linewidth=0.5,zorder=1.2)
            ax.coastlines(resolution='10m', linewidth=0.5,zorder=1.2)

            # add color bar at the bottom for all columns
            plt.colorbar(p0, ax=ax, extend = "both", orientation='horizontal', aspect=20, label=label_colorbar)
        else:
            ax.axis('off')

        # add title
        if i_time_scale==0:
            ax.set_title("CERRA-Land", fontweight = "bold")

        # add row title
        rows = ['1h', '6h', '24h']
        for ax, row in zip(axes[:,0], rows):
            # Place text directly to the left of the subplot box
            ax.annotate(row, xy=(0, 0.5), 
                xytext=(-10, 0), # Move text 40 points left
                xycoords='axes fraction', textcoords='offset points',
                fontweight='bold',
                ha='right', va='center', rotation=90) # Flipped 90 degrees

    fig.align_titles()  # requires matplotlib 3.9+

    # save the figure
    plt.savefig(f'/exec/gevin/evaluation_precip_extremes_eurocordex/PAPER/10yrRL_CompareObs_{type_diff}.png', dpi=300, bbox_inches='tight')
