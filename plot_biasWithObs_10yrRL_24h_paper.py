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
import cartopy.feature as cf
from google_crc32c import extend
from narwhals import col

# hershfield correction
hersh_correction = 1.14

# Original Projection
projorig_name  = "PlateCarree"
projorig_dict  = dict()
projorig       = ccrs.__dict__[projorig_name].__call__(**projorig_dict)

# Graphic projection
projplot_name       = "LambertConformal"
projplot_dict  = dict(central_longitude=10.5, central_latitude=49.5, standard_parallels=(33,45))
projplot       = ccrs.__dict__[projplot_name].__call__(**projplot_dict)

# folder with the 10-yr RL
folder = "D:/EUROCORDEX_extremes/DATA/10yrRL_BenPosch/"

# limits for the plots
levels = {}
levels['1h'] = np.linspace(12, 42, 16)
levels['3h'] = np.linspace(20, 60, 11)
levels['6h'] = np.linspace(30, 80, 11)
levels['12h'] = np.linspace(30, 110, 17)
levels['24h'] = np.linspace(40, 210, 18)

# Add fine resolution borders
# Category: 'cultural', Name: 'admin_0_boundary_lines_land', Scale: '10m'
fine_borders = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_0_boundary_lines_land',
    scale='10m',
    facecolor='none'
)

def unrot_lon(rlat, rlon, pole_lat, pole_lon):
    """
    Transform rotated longitude to regular non-rotated longitude lon(2D).
    """
    nrlat = np.shape(rlat)
    nrlon = np.shape(rlon)

    nrlat_rank = np.ndim(nrlat)
    nrlon_rank = np.ndim(nrlon)

    if(np.any(nrlat != nrlon) and (nrlat_rank != 1 or nrlon_rank != 1)):
        print("Function unrot_lon: rlat and rlon dimensions do not match")
        exit()

    if(nrlat_rank == 1 and nrlon_rank == 1):
        rlo = np.tile(rlon, (nrlat[0],1))
        rla = np.transpose([rlat]*nrlon[0])
    else:
        rla = rlat
        rlo = rlon

    rla = np.deg2rad(rla)
    rlo = np.deg2rad(rlo)

    s1 = np.sin(np.deg2rad(pole_lat))
    c1 = np.cos(np.deg2rad(pole_lat))
    s2 = np.sin(np.deg2rad(pole_lon))
    c2 = np.cos(np.deg2rad(pole_lon))

    tmp1 = s2*(-s1*np.cos(rlo)*np.cos(rla)+c1*np.sin(rla))-c2*np.sin(rlo)*np.cos(rla)
    tmp2 = c2*(-s1*np.cos(rlo)*np.cos(rla)+c1*np.sin(rla))+s2*np.sin(rlo)*np.cos(rla)

    lon = np.rad2deg(np.arctan(tmp1/tmp2))

    print('Function unrot_lon: min/max  %f / %f' % (np.min(lon[0,:]), np.max(lon[0,:])) )

    return lon

def unrot_lat(rlat, rlon, pole_lat, pole_lon):
    """
    Transform rotated latitude to regular non-rotated latitude lat(2D)
    """
    nrlat = np.shape(rlat)
    nrlon = np.shape(rlon)

    nrlat_rank = np.ndim(nrlat)
    nrlon_rank = np.ndim(nrlon)

    if(np.any(nrlat != nrlon) and (nrlat_rank != 1 or nrlon_rank != 1)):
        print("Function unrot_lat: rlat and rlon dimensions do not match")
        exit()

    if(nrlat_rank == 1 and nrlon_rank == 1):
        rlo = np.tile(rlon, (nrlat[0],1))
        rla = np.transpose([rlat]*nrlon[0])
    else:
        rla = rlat
        rlo = rlon

    rla = np.deg2rad(rla)
    rlo = np.deg2rad(rlo)

    s1 = np.sin(np.deg2rad(pole_lat))
    c1 = np.cos(np.deg2rad(pole_lat))

    lat = s1*np.sin(rla)+c1*np.cos(rla)*np.cos(rlo)
    lat = np.rad2deg(np.arcsin(lat))

    print('Function unrot_lat: min/max  %f / %f' % (np.min(lat[0,:]), np.max(lat[0,:])) )

    return lat

############################################################################
# Figure 1: 10-yr return levels for 1h duration
############################################################################

fig, axes = plt.subplots(figsize=(9, 18), nrows = 4, ncols = 2, subplot_kw={'projection': projplot})

limits = levels['24h']

#________________________________________________________________
# 0,0 Ref dataset
#________________________________________________________________

fObs = "return_levels_10a_24h.nc"
dsObs = xr.load_dataset(folder+fObs, engine="netcdf4")
var  = dsObs.return_levels_10a_24h
rlat = dsObs.rlat
rlon = dsObs.rlon
pole_lon = -162
pole_lat = 39.25

lon = unrot_lon(rlat, rlon, pole_lat, pole_lon)
lat = unrot_lat(rlat, rlon, pole_lat, pole_lon)

dsRef = xr.Dataset(
    data_vars={"RL": (("y", "x"), var.values)},
    coords={
        "x": rlon.values,
        "y": rlat.values,
        "lon": (("y", "x"), lon),
        "lat": (("y", "x"), lat),
    }
)
dsRef = dsRef.cf.add_bounds(["lon","lat"])

axes[0,0].set_global()
axes[0,0].set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
axes[0,0].add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
axes[0,0].add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                linewidth=0.5, edgecolor='black')
gls = axes[0,0].gridlines(draw_labels=True, linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
             x_inline=False,y_inline=False,rotate_labels = False)
gls.right_labels = False
p0 = axes[0,0].contourf(lon, lat, var, cmap="YlGnBu", levels=limits, transform=projorig, extend = "both", zorder=1)
axes[0,0].add_feature(fine_borders, edgecolor='black', linewidth=0.5, zorder=1.2)
axes[0,0].coastlines(resolution='10m', linewidth=0.5, zorder=1.2)

plt.colorbar(p0, ax=axes[0,0], orientation='vertical', label='(mm)') 

axes[0,0].text(
    0.05,
    0.95,
    "a)",
    transform=axes[0,0].transAxes,
    fontsize=14,
    fontweight='bold',
    va='top',
    ha='left',
)

#________________________________________________________________
# 1,0 HighRes products
#________________________________________________________________

# Relative biases for 10-year return levels
# list models and obs
vec_obs = ["OBS_CEH-GEAR_011EUi_1990-2016", "OBS_COMEPHORE_011EUi_1997-2022",
           "OBS_CPC_011EUi_2005-2022", "OBS_GRIPHO_011EUi_2001-2016",
           "OBS_KLIMGRID_011EUi_2011-2019", "OBS_RADKLIM_011EUi_2001-2022",
           "OBS_KNMI_011EUi_1998-2022"]
zorder_obs = [0, 1, 0, 1, 0, 1, 0]  # Define the desired order of the legend entries

time_scale = "24"

axes[1,0].set_global()
axes[1,0].set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
axes[1,0].add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
axes[1,0].add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                linewidth=0.5, edgecolor='black')
gls = axes[1,0].gridlines(draw_labels=True, linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
             x_inline=False,y_inline=False,rotate_labels = False)
gls.right_labels = False

for obs in vec_obs:
    fMod = f"matrix_10yrRL_{time_scale}hr_{obs}.csv"
    var = np.genfromtxt(folder+fMod, delimiter=',', missing_values='NA', filling_values=np.nan)*hersh_correction

    fMod = f'D:/EUROCORDEX_extremes/DATA/SIMPLESCALING/simplescaling_{obs}_JJA.nc'
    dsMod = xr.load_dataset(fMod, engine="netcdf4")

    # 1. Define your 1D coordinate vectors
    lat_vec = dsMod.lat.values
    lon_vec = dsMod.lon.values

    # 2. Create the 2D coordinate matrices
    # Use indexing='ij' for matrix/geospatial indexing (Rows=Lat, Cols=Lon)
    lon, lat = np.meshgrid(lon_vec, lat_vec, indexing='xy')
    p0 = axes[1,0].contourf(lon, lat, var, cmap="YlGnBu", levels=limits, transform=projorig, extend = "both", 
                            zorder=zorder_obs[vec_obs.index(obs)])
    
axes[1,0].add_feature(fine_borders, edgecolor='black', linewidth=0.5, zorder=1.2)
axes[1,0].coastlines(resolution='10m', linewidth=0.5, zorder=1.2)

plt.colorbar(p0, ax=axes[1,0], orientation='vertical', label='(mm)') 

axes[1,0].text(
    0.05,
    0.95,
    "b)",
    transform=axes[1,0].transAxes,
    fontsize=14,
    fontweight='bold',
    va='top',
    ha='left',
)


#________________________________________________________________
# 2,0 ERA5
#________________________________________________________________

# Relative biases for 10-year return levels
# list models and obs
obs = "OBS_ERA5_011EUi_1980-2022"

time_scale = "24"

axes[2,0].set_global()
axes[2,0].set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
axes[2,0].add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
axes[2,0].add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                linewidth=0.5, edgecolor='black')
gls = axes[2,0].gridlines(draw_labels=True, linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
             x_inline=False,y_inline=False,rotate_labels = False)
gls.right_labels = False

fMod = f"matrix_10yrRL_{time_scale}hr_{obs}.csv"
var = np.genfromtxt(folder+fMod, delimiter=',', missing_values='NA', filling_values=np.nan)*hersh_correction

fMod = f'D:/EUROCORDEX_extremes/DATA/OBS/{obs}_JJA_all_24h_perc.nc'
dsMod = xr.load_dataset(fMod, engine="netcdf4")

# 1. Define your 1D coordinate vectors
lat_vec = dsMod.lat.values
lon_vec = dsMod.lon.values

# 2. Create the 2D coordinate matrices
# Use indexing='ij' for matrix/geospatial indexing (Rows=Lat, Cols=Lon)
lon, lat = np.meshgrid(lon_vec, lat_vec, indexing='xy')
p0 = axes[2,0].contourf(lon, lat, var, cmap="YlGnBu", levels=limits, 
                        transform=projorig, extend = "both",zorder=1)

axes[2,0].add_feature(fine_borders, edgecolor='black', linewidth=0.5, zorder=1.2)
axes[2,0].coastlines(resolution='10m', linewidth=0.5, zorder = 1.2)

plt.colorbar(p0, ax=axes[2,0], orientation='vertical', label='(mm)')

axes[2,0].text(
    0.05,
    0.95,
    "d)",
    transform=axes[2,0].transAxes,
    fontsize=14,
    fontweight='bold',
    va='top',
    ha='left',
)

#________________________________________________________________
# 3,0 CERRA
#________________________________________________________________

# Relative biases for 10-year return levels
# list models and obs
obs = "OBS_CERRALAND_011EUi_1986-2020"

time_scale = "24"

axes[3,0].set_global()
axes[3,0].set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
axes[3,0].add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
axes[3,0].add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                linewidth=0.5, edgecolor='black')
gls = axes[3,0].gridlines(draw_labels=True, linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
             x_inline=False,y_inline=False,rotate_labels = False)
gls.right_labels = False

fMod = f"matrix_10yrRL_{time_scale}hr_{obs}.csv"
var = np.genfromtxt(folder+fMod, delimiter=',', missing_values='NA', filling_values=np.nan)*hersh_correction

fMod = f'D:/EUROCORDEX_extremes/DATA/OBS/{obs}_JJA_all_24h_perc.nc'
dsMod = xr.load_dataset(fMod, engine="netcdf4")

# 1. Define your 1D coordinate vectors
lat_vec = dsMod.lat.values
lon_vec = dsMod.lon.values

# 2. Create the 2D coordinate matrices
# Use indexing='ij' for matrix/geospatial indexing (Rows=Lat, Cols=Lon)
lon, lat = np.meshgrid(lon_vec, lat_vec, indexing='xy')
p0 = axes[3,0].contourf(lon, lat, var, cmap="YlGnBu", levels=limits, 
                        transform=projorig, extend = "both",zorder=1)

axes[3,0].add_feature(fine_borders, edgecolor='black', linewidth=0.5, zorder=1.2)
axes[3,0].coastlines(resolution='10m', linewidth=0.5, zorder = 1.2)

plt.colorbar(p0, ax=axes[3,0], orientation='vertical', label='(mm)')

axes[3,0].text(
    0.05,
    0.95,
    "f)",
    transform=axes[3,0].transAxes,
    fontsize=14,
    fontweight='bold',
    va='top',
    ha='left',
)

#________________________________________________________________
# 1,1 Diff highres /ref
#________________________________________________________________

# Relative biases for 10-year return levels
# list models and obs

axes[1,1].set_global()
axes[1,1].set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
axes[1,1].add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
axes[1,1].add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                linewidth=0.5, edgecolor='black')
gls = axes[1,1].gridlines(draw_labels=True, linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
             x_inline=False,y_inline=False,rotate_labels = False)
gls.right_labels = False

for obs in vec_obs:
    fMod = f"matrix_10yrRL_{time_scale}hr_{obs}.csv"
    var = np.genfromtxt(folder+fMod, delimiter=',', missing_values='NA', filling_values=np.nan)*hersh_correction

    fMod = f'D:/EUROCORDEX_extremes/DATA/SIMPLESCALING/simplescaling_{obs}_JJA.nc'
    dsMod = xr.load_dataset(fMod, engine="netcdf4")

    # 1. Define your 1D coordinate vectors
    lat_vec = dsMod.lat.values
    lon_vec = dsMod.lon.values

    # 2. Create the 2D coordinate matrices
    # Use indexing='ij' for matrix/geospatial indexing (Rows=Lat, Cols=Lon)
    lon, lat = np.meshgrid(lon_vec, lat_vec, indexing='xy')

    # 3; buil clean dataset
    dsObs2regrid = xr.Dataset(
        data_vars={"RL": (("y", "x"), var)},
        coords={
            "x": dsMod.lon.values,
            "y": dsMod.lat.values,
            "lon": (("y", "x"), lon),
            "lat": (("y", "x"), lat),
        }
    )
    dsObs2regrid = dsObs2regrid.cf.add_bounds(["lon","lat"])

    # 4. Regrid to the reference dataset
    # regrid obs
    regridder = xe.Regridder(
        dsObs2regrid, dsRef, "conservative", unmapped_to_nan=True, ignore_degenerate=True
    )
    dsObs = regridder(dsObs2regrid, keep_attrs=True)
    
    X = dsObs.lon
    Y = dsObs.lat
    ZREF = dsRef.RL
    ZOBS = dsObs.RL
    Z = (ZOBS - ZREF) / ZREF * 100
    p0 = axes[1,1].pcolormesh(X, Y, Z, cmap=cmaps.precip_diff_12lev, vmin=-65, vmax=65,
                            transform=projorig, zorder=zorder_obs[vec_obs.index(obs)])
    
axes[1,1].add_feature(fine_borders, edgecolor='black', linewidth=0.5,zorder=1.2)
axes[1,1].coastlines(resolution='10m', linewidth=0.5,zorder=1.2)

plt.colorbar(p0, ax=axes[1,1], orientation='vertical', label='(%)') 

axes[1,1].text(
    0.05,
    0.95,
    "c)",
    transform=axes[1,1].transAxes,
    fontsize=14,
    fontweight='bold',
    va='top',
    ha='left',
)

#________________________________________________________________
# 2,1 Diff ERA5 /ref
#________________________________________________________________

# Relative biases for 10-year return levels
# list models and obs
obs = "OBS_ERA5_011EUi_1980-2022"

axes[2,1].set_global()
axes[2,1].set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
axes[2,1].add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
axes[2,1].add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                linewidth=0.5, edgecolor='black')
gls = axes[2,1].gridlines(draw_labels=True, linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
             x_inline=False,y_inline=False,rotate_labels = False)
gls.right_labels = False

fMod = f"matrix_10yrRL_{time_scale}hr_{obs}.csv"
var = np.genfromtxt(folder+fMod, delimiter=',', missing_values='NA', filling_values=np.nan)*hersh_correction

fMod = f'D:/EUROCORDEX_extremes/DATA/OBS/{obs}_JJA_all_24h_perc.nc'
dsMod = xr.load_dataset(fMod, engine="netcdf4")

# 1. Define your 1D coordinate vectors
lat_vec = dsMod.lat.values
lon_vec = dsMod.lon.values

# 2. Create the 2D coordinate matrices
# Use indexing='ij' for matrix/geospatial indexing (Rows=Lat, Cols=Lon)
lon, lat = np.meshgrid(lon_vec, lat_vec, indexing='xy')

# 3; buil clean dataset
dsObs2regrid = xr.Dataset(
    data_vars={"RL": (("y", "x"), var)},
    coords={
        "x": dsMod.lon.values,
        "y": dsMod.lat.values,
        "lon": (("y", "x"), lon),
        "lat": (("y", "x"), lat),
    }
)
dsObs2regrid = dsObs2regrid.cf.add_bounds(["lon","lat"])

# 4. Regrid to the reference dataset
# regrid obs
regridder = xe.Regridder(
    dsObs2regrid, dsRef, "conservative", unmapped_to_nan=True, ignore_degenerate=True
)
dsObs = regridder(dsObs2regrid, keep_attrs=True)

X = dsObs.lon
Y = dsObs.lat
ZREF = dsRef.RL
ZOBS = dsObs.RL
Z = (ZOBS - ZREF) / ZREF * 100
p0 = axes[2,1].pcolormesh(X, Y, Z, cmap=cmaps.precip_diff_12lev, vmin=-65, vmax=65,
                        transform=projorig, zorder=1)
    
axes[2,1].add_feature(fine_borders, edgecolor='black', linewidth=0.5, zorder=1.2)
axes[2,1].coastlines(resolution='10m', linewidth=0.5, zorder = 1.2)

plt.colorbar(p0, ax=axes[2,1], orientation='vertical', label='(%)') 

axes[2,1].text(
    0.05,
    0.95,
    "e)",
    transform=axes[2,1].transAxes,
    fontsize=14,
    fontweight='bold',
    va='top',
    ha='left',
)

axes[0,1].set_visible(False)

#________________________________________________________________
# 3,1 Diff CERRA /ref
#________________________________________________________________

# Relative biases for 10-year return levels
# list models and obs
obs = "OBS_CERRALAND_011EUi_1986-2020"

axes[3,1].set_global()
axes[3,1].set_extent([-10, 25, 35, 65], crs=ccrs.PlateCarree())
axes[3,1].add_feature(cartopy.feature.OCEAN, color='white', zorder=1.1)
axes[3,1].add_feature(cartopy.feature.LAND, color='lightgray',zorder=-0.5,
                linewidth=0.5, edgecolor='black')
gls = axes[3,1].gridlines(draw_labels=True, linewidth=0.5, color='gray', xlocs=range(-180,180,10), ylocs=range(-90,90,10),
             x_inline=False,y_inline=False,rotate_labels = False)
gls.right_labels = False

fMod = f"matrix_10yrRL_{time_scale}hr_{obs}.csv"
var = np.genfromtxt(folder+fMod, delimiter=',', missing_values='NA', filling_values=np.nan)*hersh_correction

fMod = f'D:/EUROCORDEX_extremes/DATA/OBS/{obs}_JJA_all_24h_perc.nc'
dsMod = xr.load_dataset(fMod, engine="netcdf4")

# 1. Define your 1D coordinate vectors
lat_vec = dsMod.lat.values
lon_vec = dsMod.lon.values

# 2. Create the 2D coordinate matrices
# Use indexing='ij' for matrix/geospatial indexing (Rows=Lat, Cols=Lon)
lon, lat = np.meshgrid(lon_vec, lat_vec, indexing='xy')

# 3; buil clean dataset
dsObs2regrid = xr.Dataset(
    data_vars={"RL": (("y", "x"), var)},
    coords={
        "x": dsMod.lon.values,
        "y": dsMod.lat.values,
        "lon": (("y", "x"), lon),
        "lat": (("y", "x"), lat),
    }
)
dsObs2regrid = dsObs2regrid.cf.add_bounds(["lon","lat"])

# 4. Regrid to the reference dataset
# regrid obs
regridder = xe.Regridder(
    dsObs2regrid, dsRef, "conservative", unmapped_to_nan=True, ignore_degenerate=True
)
dsObs = regridder(dsObs2regrid, keep_attrs=True)

X = dsObs.lon
Y = dsObs.lat
ZREF = dsRef.RL
ZOBS = dsObs.RL
Z = (ZOBS - ZREF) / ZREF * 100
p0 = axes[3,1].pcolormesh(X, Y, Z, cmap=cmaps.precip_diff_12lev, vmin=-65, vmax=65,
                        transform=projorig, zorder=1)
    
axes[3,1].add_feature(fine_borders, edgecolor='black', linewidth=0.5, zorder=1.2)
axes[3,1].coastlines(resolution='10m', linewidth=0.5, zorder = 1.2)

plt.colorbar(p0, ax=axes[3,1], orientation='vertical', label='(%)') 

axes[3,1].text(
    0.05,
    0.95,
    "g)",
    transform=axes[3,1].transAxes,
    fontsize=14,
    fontweight='bold',
    va='top',
    ha='left',
)

axes[0,1].set_visible(False)


# save the figure
plt.savefig(f'D:/EUROCORDEX_extremes/PAPER/10yrRL_CompareObs_24h.png', dpi=300, bbox_inches='tight')
